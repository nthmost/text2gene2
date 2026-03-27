"""
Abstract-level validation: check that a returned PMID actually mentions
the specific variant (not just the gene).

For each PMID we fetch the PubMed abstract via Entrez efetch and search for
the variant's key identifiers. This catches false positives like the SMAD4
case where a degenerate query matched unrelated variants.

This is run as a post-filter on the CitationTable. It is optional and off
by default (costs NCBI API calls). Enable per-request via the API or CLI.

Validation tiers (citation.validation_tier):
  "confirmed"  — variant notation found verbatim in abstract or title
  "probable"   — gene + position found (notation may differ, e.g. older paper)
  "unverified" — abstract not available or no mention found (not discarded,
                 just flagged — may still be relevant supplementary data)
"""
import logging
import re
from xml.etree import ElementTree

import httpx

from text2gene2.models import Citation, CitationTable, LVGResult
from text2gene2 import rate_limit
from text2gene2.config import settings

log = logging.getLogger(__name__)

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_BATCH  = 200   # efetch supports up to 200 IDs at once


def _api_params() -> dict:
    p: dict = {}
    if settings.ncbi_api_key:
        p["api_key"] = settings.ncbi_api_key
    return p


def _build_patterns(lvg: LVGResult) -> list[re.Pattern]:
    """Build regex patterns from the LVG to search for in abstracts."""
    patterns = []

    # Coding changes: c.1659A>C, c.1659A→C, c.1659 A>C
    for h in lvg.hgvs_c:
        change = h.split(":")[-1] if ":" in h else h
        # Allow whitespace and arrow variants (e.g. → instead of >)
        escaped = re.escape(change).replace(r"\>", r"[>→]").replace(r"\ ", r"\s*")
        patterns.append(re.compile(escaped, re.IGNORECASE))

    # Protein changes: p.R361H, R361H, Arg361His
    for h in lvg.hgvs_p:
        change = h.split(":")[-1] if ":" in h else h
        change = change.lstrip("p.(").rstrip(")")
        if change and re.search(r"[A-Za-z].*\d", change):
            patterns.append(re.compile(re.escape(change), re.IGNORECASE))

    # rsIDs
    for rsid in lvg.rsids:
        patterns.append(re.compile(re.escape(rsid), re.IGNORECASE))

    return patterns


def _position_patterns(lvg: LVGResult) -> list[re.Pattern]:
    """Fallback: look for numeric position from any HGVS form."""
    positions = set()
    for h in lvg.hgvs_c + lvg.hgvs_p:
        for m in re.finditer(r"\d+", h):
            if len(m.group()) >= 3:   # positions < 100 are too common to be signal
                positions.add(m.group())
    return [re.compile(r"\b" + p + r"\b") for p in positions]


async def _fetch_abstracts(pmids: list[int]) -> dict[int, str]:
    """Batch-fetch PubMed abstracts. Returns pmid → abstract text."""
    result: dict[int, str] = {}
    for i in range(0, len(pmids), _BATCH):
        batch = pmids[i : i + _BATCH]
        await rate_limit.clinvar.acquire()   # reuse NCBI rate limiter
        params = {
            **_api_params(),
            "db": "pubmed",
            "id": ",".join(str(p) for p in batch),
            "rettype": "abstract",
            "retmode": "xml",
        }
        try:
            async with httpx.AsyncClient() as client:
                resp = await client.get(f"{_EUTILS}/efetch.fcgi", params=params, timeout=30.0)
                resp.raise_for_status()
                tree = ElementTree.fromstring(resp.text)
                for article in tree.findall(".//PubmedArticle"):
                    pmid_el = article.find(".//PMID")
                    if pmid_el is None:
                        continue
                    pmid = int(pmid_el.text)
                    # Collect title + abstract text
                    parts = []
                    title = article.find(".//ArticleTitle")
                    if title is not None and title.text:
                        parts.append(title.text)
                    for at in article.findall(".//AbstractText"):
                        if at.text:
                            parts.append(at.text)
                    result[pmid] = " ".join(parts)
        except Exception as e:
            log.warning("efetch error for batch starting %d: %s", batch[0], e)
    return result


def _score_mention(text: str, exact_patterns: list[re.Pattern],
                   pos_patterns: list[re.Pattern]) -> str:
    """Return validation tier for a single abstract."""
    if not text:
        return "unverified"
    if any(p.search(text) for p in exact_patterns):
        return "confirmed"
    if any(p.search(text) for p in pos_patterns):
        return "probable"
    return "unverified"


async def validate_citations(table: CitationTable, max_pmids: int = 500) -> CitationTable:
    """
    Fetch abstracts and score each citation. Modifies citations in-place,
    adding a validation_tier field. Caps at max_pmids to limit API cost.

    Citations are re-sorted: confirmed > probable > unverified, then by
    confidence (source count) within each tier.
    """
    if not table.lvg or not table.citations:
        return table

    exact_pats = _build_patterns(table.lvg)
    pos_pats   = _position_patterns(table.lvg)

    if not exact_pats and not pos_pats:
        log.debug("validate: no patterns extractable from LVG, skipping")
        return table

    pmids_to_fetch = [c.pmid for c in table.citations[:max_pmids]]
    abstracts = await _fetch_abstracts(pmids_to_fetch)

    tier_order = {"confirmed": 0, "probable": 1, "unverified": 2}
    for citation in table.citations:
        text = abstracts.get(citation.pmid, "")
        citation.validation_tier = _score_mention(text, exact_pats, pos_pats)

    # Re-sort: tier first, then confidence desc
    table.citations.sort(
        key=lambda c: (tier_order.get(getattr(c, "validation_tier", "unverified"), 2),
                       -c.confidence)
    )

    confirmed = sum(1 for c in table.citations if getattr(c, "validation_tier", "") == "confirmed")
    probable  = sum(1 for c in table.citations if getattr(c, "validation_tier", "") == "probable")
    log.info("Validation: %d confirmed, %d probable, %d unverified",
             confirmed, probable, len(table.citations) - confirmed - probable)
    return table
