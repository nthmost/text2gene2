"""
Article metadata enrichment for CitationTable.

Fetch order (fastest to slowest):
  1. Redis cache — already-enriched PMIDs
  2. pubmed.article local DB — 40M articles with structured columns
  3. NCBI efetch — fallback for any remaining misses

Results are cached in Redis per PMID with a long TTL (metadata rarely changes).
This is run automatically for web UI requests but skipped for the JSON API
(use ?enrich=true to opt in).
"""
import asyncio
import logging
import os
from xml.etree import ElementTree

import httpx

from text2gene2.cache import cache_get, cache_set
from text2gene2.models import Citation, CitationTable
from text2gene2 import rate_limit
from text2gene2.config import settings

log = logging.getLogger(__name__)

_EUTILS  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_BATCH   = 100
_TTL     = 60 * 60 * 24 * 30   # 30 days — article metadata is stable


# ── Local DB ──────────────────────────────────────────────────────────────────

_db_conn = None

def _get_db_conn():
    """Lazy singleton psycopg2 connection to the local medgen DB."""
    global _db_conn
    db_url = os.environ.get("METAPUB_DB_URL")
    if not db_url:
        return None
    try:
        import psycopg2
        if _db_conn is None or _db_conn.closed:
            _db_conn = psycopg2.connect(db_url)
            _db_conn.set_session(readonly=True, autocommit=True)
        return _db_conn
    except Exception as e:
        log.warning("enrich: local DB connect failed: %s", e)
        return None


_SELECT_LOCAL = """
SELECT pmid, title, authors, journal, year, doi, pmc_id, abstract
FROM pubmed.article
WHERE pmid = ANY(%s)
"""

def _fetch_meta_local_sync(pmids: list[int]) -> dict[int, dict]:
    """Query pubmed.article for structured metadata. Returns {pmid: meta_dict}."""
    conn = _get_db_conn()
    if conn is None or not pmids:
        return {}
    try:
        with conn.cursor() as cur:
            cur.execute(_SELECT_LOCAL, (pmids,))
            rows = cur.fetchall()
    except Exception as e:
        log.warning("enrich: local DB query error: %s", e)
        global _db_conn
        _db_conn = None
        return {}

    result = {}
    for pmid, title, authors, journal, year, doi, pmc_id, abstract in rows:
        meta: dict = {}
        if title:
            meta["title"] = title.rstrip(".")
        if authors:
            meta["authors"] = authors
        if journal:
            meta["journal"] = journal
        if year:
            meta["year"] = int(year)
        if doi:
            meta["doi"] = doi
        if pmc_id:
            meta["pmc"] = pmc_id
        if abstract:
            meta["abstract_snippet"] = abstract[:400] + ("…" if len(abstract) > 400 else "")
        result[int(pmid)] = meta
    return result


async def _fetch_meta_local(pmids: list[int]) -> dict[int, dict]:
    return await asyncio.to_thread(_fetch_meta_local_sync, pmids)


# ── NCBI efetch fallback ───────────────────────────────────────────────────────

def _api_params() -> dict:
    p: dict = {}
    if settings.ncbi_api_key:
        p["api_key"] = settings.ncbi_api_key
    return p


def _first_author(article) -> str:
    authors = article.findall(".//Author")
    if not authors:
        return ""
    a = authors[0]
    last     = (a.findtext("LastName") or "").strip()
    initials = (a.findtext("Initials") or "").strip()
    first    = (a.findtext("ForeName") or initials).strip()
    name = f"{last} {first}".strip() if first else last
    return f"{name} et al." if len(authors) > 1 else name


def _parse_article(article) -> dict:
    """Extract metadata fields from a PubmedArticle XML element."""
    meta: dict = {}

    title_el = article.find(".//ArticleTitle")
    if title_el is not None:
        meta["title"] = "".join(title_el.itertext()).strip().rstrip(".")

    meta["authors"] = _first_author(article)

    j = article.find(".//Journal")
    if j is not None:
        meta["journal"] = (
            j.findtext("ISOAbbreviation") or j.findtext("Title") or ""
        ).strip()

    for xpath in [".//PubDate/Year", ".//PubMedPubDate[@PubStatus='pubmed']/Year"]:
        yr = article.findtext(xpath)
        if yr and yr.isdigit():
            meta["year"] = int(yr)
            break

    for id_el in article.findall(".//ArticleId"):
        id_type = id_el.get("IdType", "")
        val = (id_el.text or "").strip()
        if id_type == "doi" and val:
            meta["doi"] = val
        elif id_type == "pmc" and val:
            meta["pmc"] = val

    parts = [el.text for el in article.findall(".//AbstractText") if el.text]
    if parts:
        full = " ".join(parts)
        meta["abstract_snippet"] = full[:400] + ("…" if len(full) > 400 else "")

    return meta


async def _fetch_meta_ncbi(pmids: list[int]) -> dict[int, dict]:
    """Fetch metadata from NCBI efetch for a batch of PMIDs."""
    result: dict[int, dict] = {}
    params = {
        **_api_params(),
        "db": "pubmed",
        "id": ",".join(str(p) for p in pmids),
        "rettype": "abstract",
        "retmode": "xml",
    }
    await rate_limit.clinvar.acquire()
    try:
        async with httpx.AsyncClient() as client:
            resp = await client.get(f"{_EUTILS}/efetch.fcgi", params=params, timeout=30.0)
            resp.raise_for_status()
            tree = ElementTree.fromstring(resp.text)
            for article in tree.findall(".//PubmedArticle"):
                pmid_el = article.find(".//PMID")
                if pmid_el is None or not pmid_el.text:
                    continue
                pmid = int(pmid_el.text)
                result[pmid] = _parse_article(article)
    except Exception as e:
        log.warning("enrich: NCBI efetch error for batch starting %d: %s", pmids[0], e)
    return result


# ── Main entry point ───────────────────────────────────────────────────────────

async def enrich_citations(table: CitationTable) -> CitationTable:
    """
    Populate Citation metadata fields for all citations in the table.

    Checks Redis cache first, then local pubmed.article DB, then NCBI efetch.
    """
    if not table.citations:
        return table

    meta_map: dict[int, dict] = {}
    uncached: list[int] = []

    for c in table.citations:
        key = f"meta:{c.pmid}"
        cached = await cache_get(key)
        if cached is not None:
            meta_map[c.pmid] = cached
        else:
            uncached.append(c.pmid)

    if uncached:
        # Tier 1: local DB
        local = await _fetch_meta_local(uncached)
        for pmid, meta in local.items():
            meta_map[pmid] = meta
            await cache_set(f"meta:{pmid}", meta, ttl=_TTL)

        ncbi_needed = [p for p in uncached if p not in local]
        log.debug("enrich: %d from local DB, %d need NCBI", len(local), len(ncbi_needed))

        # Tier 2: NCBI efetch for misses
        for i in range(0, len(ncbi_needed), _BATCH):
            batch = ncbi_needed[i : i + _BATCH]
            fetched = await _fetch_meta_ncbi(batch)
            for pmid, meta in fetched.items():
                meta_map[pmid] = meta
                await cache_set(f"meta:{pmid}", meta, ttl=_TTL)
            for pmid in batch:
                if pmid not in fetched:
                    await cache_set(f"meta:{pmid}", {}, ttl=_TTL)

    for citation in table.citations:
        meta = meta_map.get(citation.pmid, {})
        if meta.get("title"):
            citation.title = meta["title"]
        if meta.get("authors"):
            citation.authors = meta["authors"]
        if meta.get("journal"):
            citation.journal = meta["journal"]
        if meta.get("year"):
            citation.year = meta["year"]
        if meta.get("doi"):
            citation.doi = meta["doi"]
        citation.pmc              = meta.get("pmc")
        citation.abstract_snippet = meta.get("abstract_snippet")

    return table
