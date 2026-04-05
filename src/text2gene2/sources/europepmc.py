"""
Europe PMC source — open-access literature including preprints.

REST API: https://europepmc.org/RestfulWebService
No API key required. Covers PubMed + PMC + preprints + gray literature.
Especially useful for recent papers not yet indexed in PubMed.

Query strategy: use the expanded query builder from pipeline/expand.py which
generates gene synonym + multi-form variant OR-queries. Falls back to raw
input HGVS if expansion produces nothing.

After fetching results, does post-hoc attribution: checks which expanded
variant forms appear in each result's title/abstract to determine which
query form actually matched.
"""
import logging

import httpx

from text2gene2.cache import cache_get, cache_set
from text2gene2.config import settings
from text2gene2.models import LVGResult, Source, SourceResult
from text2gene2.sources.base import PMIDSource
from text2gene2.pipeline.expand import build_europepmc_query, expand_variant, get_gene_synonyms
from text2gene2 import rate_limit

log = logging.getLogger(__name__)

_BASE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def _attribute_hit(title: str, abstract: str, search_forms: list[str],
                   gene_synonyms: list[str]) -> str | None:
    """
    Post-hoc attribution: which search form most likely matched this result?

    Checks title+abstract for each expanded form, returns the most specific
    match. Prefers variant forms over gene-only matches.
    """
    text = f"{title} {abstract}".lower()

    # Check variant forms (most specific first — these are ordered by specificity
    # in ExpandedVariant.all_search_forms)
    for form in search_forms:
        if form.lower() in text:
            # Find which gene synonym also appears for a complete attribution
            gene_match = None
            for syn in gene_synonyms:
                if syn.lower() in text:
                    gene_match = syn
                    break
            if gene_match:
                return f"{gene_match} + {form}"
            return form

    # If no specific variant form found, check gene synonyms
    for syn in gene_synonyms:
        if syn.lower() in text:
            return f"{syn} (gene only)"

    return None


class EuropePMCSource(PMIDSource):
    source = Source.EUROPEPMC

    async def query(self, lvg: LVGResult) -> SourceResult:
        cache_key = f"europepmc:{lvg.input_hgvs}"
        cached = await cache_get(cache_key)
        if cached is not None:
            if isinstance(cached, dict):
                prov = {int(k): v for k, v in cached.get("prov", {}).items()}
                return SourceResult(source=self.source, pmids=cached["pmids"],
                                    pmid_provenance=prov,
                                    query_used=cached.get("query"), cached=True)
            return SourceResult(source=self.source, pmids=cached, cached=True)

        query = await build_europepmc_query(lvg)
        log.debug("Europe PMC query: %r", query)

        # Pre-compute expansion forms for post-hoc attribution
        exp = expand_variant(lvg)
        gene_syns = await get_gene_synonyms(lvg.gene_symbol) if lvg.gene_symbol else []

        await rate_limit.europepmc.acquire()
        params = {
            "query": query,
            "format": "json",
            "resultType": "core",
            "pageSize": 50,
            "sort": "CITED desc",
        }

        try:
            async with httpx.AsyncClient() as client:
                resp = await client.get(_BASE, params=params, timeout=20.0)
                resp.raise_for_status()
                data = resp.json()
        except Exception as e:
            log.warning("Europe PMC error: %s", e)
            return SourceResult(source=self.source, pmids=[], error=str(e))

        pmids: list[int] = []
        provenance: dict[int, str] = {}
        for hit in data.get("resultList", {}).get("result", []):
            pmid = hit.get("pmid")
            if not pmid:
                continue
            try:
                pmid_int = int(pmid)
            except ValueError:
                continue

            if pmid_int not in provenance:
                pmids.append(pmid_int)
                # Post-hoc attribution from title + abstract
                title = hit.get("title", "")
                abstract = hit.get("abstractText", "")
                match = _attribute_hit(title, abstract, exp.all_search_forms, gene_syns)
                if match:
                    provenance[pmid_int] = match

        result = sorted(set(pmids))
        prov_str = {str(k): v for k, v in provenance.items()}
        await cache_set(cache_key, {"pmids": result, "prov": prov_str, "query": query},
                        ttl=settings.cache_ttl_europepmc)
        return SourceResult(source=self.source, pmids=result, query_used=query,
                            pmid_provenance=provenance)
