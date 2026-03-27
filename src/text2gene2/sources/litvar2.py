"""
LitVar2 source — NCBI's text-mined variant-to-literature database.

Two-step workflow:
  1. autocomplete: free-text variant → VarID (anchored to rsID)
  2. publications: VarID → list of PMIDs

We try every HGVS string in the LVG, plus any rsIDs found by VariantValidator.
VarID format: litvar@{rsid}##   URL-encoded: litvar%40{rsid}%23%23

API: https://www.ncbi.nlm.nih.gov/research/litvar2-api/
No API key required. No documented rate limit — we self-limit to 5 req/sec.
"""
import logging
import re
from urllib.parse import quote

import httpx

# A meaningful short protein/coding change must contain at least one digit
# (position) and at least one letter. "p." or "p.(" alone is not meaningful.
_MEANINGFUL_RE = re.compile(r"[A-Za-z].*\d|\d.*[A-Za-z]")

from text2gene2.cache import cache_get, cache_set
from text2gene2.config import settings
from text2gene2.models import LVGResult, Source, SourceResult
from text2gene2.sources.base import PMIDSource
from text2gene2 import rate_limit

log = logging.getLogger(__name__)

_BASE = "https://www.ncbi.nlm.nih.gov/research/litvar2-api"


def _varid_url(rsid: str) -> str:
    varid = f"litvar@{rsid}##"
    return f"{_BASE}/variant/get/{quote(varid, safe='')}/publications"


async def _autocomplete(query: str, client: httpx.AsyncClient) -> list[str]:
    """Return rsIDs matched by the autocomplete endpoint for a query string."""
    await rate_limit.litvar2.acquire()
    try:
        resp = await client.get(f"{_BASE}/variant/autocomplete/", params={"query": query}, timeout=15.0)
        resp.raise_for_status()
        hits = resp.json()
        return [h["rsid"] for h in hits if h.get("rsid")]
    except Exception as e:
        log.debug("LitVar2 autocomplete error for %r: %s", query, e)
        return []


async def _publications(rsid: str, client: httpx.AsyncClient) -> list[int]:
    """Return PMIDs for a given rsID via LitVar2 publications endpoint."""
    await rate_limit.litvar2.acquire()
    try:
        resp = await client.get(_varid_url(rsid), timeout=15.0)
        resp.raise_for_status()
        data = resp.json()
        return [int(p) for p in data.get("pmids", [])]
    except Exception as e:
        log.debug("LitVar2 publications error for %s: %s", rsid, e)
        return []


class LitVar2Source(PMIDSource):
    source = Source.LITVAR2

    async def query(self, lvg: LVGResult) -> SourceResult:
        cache_key = f"litvar2:{lvg.input_hgvs}"
        cached = await cache_get(cache_key)
        if cached is not None:
            return SourceResult(source=self.source, pmids=cached, cached=True)

        # Gather candidate rsIDs from two routes:
        #   a) rsIDs already found by VariantValidator
        #   b) autocomplete disambiguation for each HGVS form
        rsids: set[str] = set(lvg.rsids)

        async with httpx.AsyncClient() as client:
            # Autocomplete for input HGVS + gene-qualified forms.
            # LitVar2 works best with short, clean queries.
            queries = [lvg.input_hgvs]
            if lvg.gene_symbol:
                for p in lvg.hgvs_p[:2]:
                    short_p = p.split(":")[-1] if ":" in p else p
                    # Trim frameshift/nonsense suffix after the AA change position.
                    # But stop-codon extensions like p.(*553Cext*40) start with *,
                    # so splitting on * leaves "p.(" — meaningless. Guard against this.
                    short_p = short_p.split("*")[0].split("Ter")[0].rstrip("(").rstrip(".")
                    if _MEANINGFUL_RE.search(short_p):
                        queries.append(f"{lvg.gene_symbol} {short_p}")
                    else:
                        log.debug("LitVar2: skipping degenerate protein query %r", short_p)
                # Also try gene + coding change (useful for intronic/non-coding variants)
                for c in lvg.hgvs_c[:1]:
                    short_c = c.split(":")[-1] if ":" in c else c
                    if _MEANINGFUL_RE.search(short_c):
                        queries.append(f"{lvg.gene_symbol} {short_c}")

            for q in queries:
                found = await _autocomplete(q, client)
                rsids.update(found)

            if not rsids:
                log.debug("LitVar2: no rsIDs found for %s", lvg.input_hgvs)
                return SourceResult(source=self.source, pmids=[])

            # Fetch publications for each rsID
            all_pmids: set[int] = set()
            for rsid in rsids:
                pmids = await _publications(rsid, client)
                all_pmids.update(pmids)

        result = sorted(all_pmids)
        await cache_set(cache_key, result, ttl=settings.cache_ttl_litvar2)
        return SourceResult(source=self.source, pmids=result)
