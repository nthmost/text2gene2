"""
ClinVar source — curated variant→PMID associations via NCBI Entrez eutils.

Strategy:
  esearch ClinVar for each HGVS → get VariationIDs
  elink VariationID → PubMed → get PMIDs

With an NCBI API key: 10 req/sec.  Without: 3 req/sec.
Key is read from settings.ncbi_api_key (optional but recommended).
"""
import logging
from xml.etree import ElementTree

import httpx

from text2gene2.cache import cache_get, cache_set
from text2gene2.config import settings
from text2gene2.models import LVGResult, Source, SourceResult
from text2gene2.sources.base import PMIDSource
from text2gene2 import rate_limit

log = logging.getLogger(__name__)

_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _api_params() -> dict:
    p: dict = {"retmode": "json"}
    if settings.ncbi_api_key:
        p["api_key"] = settings.ncbi_api_key
    return p


async def _esearch_clinvar(term: str, client: httpx.AsyncClient) -> list[str]:
    """Search ClinVar for a term; returns list of ClinVar UIDs."""
    await rate_limit.clinvar.acquire()
    params = {**_api_params(), "db": "clinvar", "term": term, "retmax": "50"}
    try:
        resp = await client.get(f"{_EUTILS}/esearch.fcgi", params=params, timeout=15.0)
        resp.raise_for_status()
        data = resp.json()
        return data.get("esearchresult", {}).get("idlist", [])
    except Exception as e:
        log.debug("ClinVar esearch error for %r: %s", term, e)
        return []


async def _elink_to_pubmed(clinvar_ids: list[str], client: httpx.AsyncClient) -> list[int]:
    """Convert ClinVar UIDs → PubMed PMIDs via elink."""
    if not clinvar_ids:
        return []
    await rate_limit.clinvar.acquire()
    params = {
        **_api_params(),
        "dbfrom": "clinvar",
        "db": "pubmed",
        "id": ",".join(clinvar_ids),
        "retmode": "xml",
    }
    try:
        resp = await client.get(f"{_EUTILS}/elink.fcgi", params=params, timeout=15.0)
        resp.raise_for_status()
        tree = ElementTree.fromstring(resp.text)
        pmids = [el.text for el in tree.findall(".//Link/Id") if el.text]
        return [int(p) for p in pmids]
    except Exception as e:
        log.debug("ClinVar elink error: %s", e)
        return []


class ClinVarSource(PMIDSource):
    source = Source.CLINVAR

    async def query(self, lvg: LVGResult) -> SourceResult:
        cache_key = f"clinvar:{lvg.input_hgvs}"
        cached = await cache_get(cache_key)
        if cached is not None:
            if isinstance(cached, dict):
                prov = {int(k): v for k, v in cached.get("prov", {}).items()}
                return SourceResult(source=self.source, pmids=cached["pmids"],
                                    pmid_provenance=prov, cached=True)
            return SourceResult(source=self.source, pmids=cached, cached=True)

        all_pmids: set[int] = set()
        provenance: dict[int, str] = {}

        async with httpx.AsyncClient() as client:
            # Build search terms from the LVG.
            # ClinVar indexes by HGVS, gene symbol, and rsID.
            terms: list[str] = []

            for h in lvg.hgvs_c[:3]:      # try up to 3 coding forms
                terms.append(f'"{h}"[Variant Name]')
            for h in lvg.hgvs_g[:2]:
                terms.append(f'"{h}"[Variant Name]')
            for rsid in lvg.rsids[:3]:
                terms.append(f'"{rsid}"[RS# (All)]')

            # Fall back to the raw input if we have nothing
            if not terms:
                terms.append(f'"{lvg.input_hgvs}"[Variant Name]')

            for term in terms:
                ids = await _esearch_clinvar(term, client)
                if ids:
                    pmids = await _elink_to_pubmed(ids, client)
                    for pmid in pmids:
                        if pmid not in all_pmids:
                            # Strip ClinVar field tags for display
                            display_term = term.split('"')[1] if '"' in term else term
                            provenance[pmid] = display_term
                        all_pmids.add(pmid)

        result = sorted(all_pmids)
        prov_str = {str(k): v for k, v in provenance.items()}
        await cache_set(cache_key, {"pmids": result, "prov": prov_str},
                        ttl=settings.cache_ttl_clinvar)
        return SourceResult(source=self.source, pmids=result,
                            pmid_provenance=provenance)
