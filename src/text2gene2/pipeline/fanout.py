"""
Core pipeline: given an HGVS string, fan out to all sources simultaneously.

This is the heart of text2gene2. The flow is:
  1. Expand HGVS → LVGResult  (VariantValidator, cached)
  2. asyncio.gather() across all sources in parallel
  3. Merge + deduplicate PMIDs into a CitationTable
"""
import asyncio
import logging

from text2gene2.lvg import get_lvg
from text2gene2.models import CitationTable, LVGResult, SourceResult
log = logging.getLogger(__name__)

_SOURCES = None

def _get_sources():
    global _SOURCES
    if _SOURCES is None:
        from text2gene2.sources import (
            ClinVarSource, EuropePMCSource, GoogleCSESource,
            LitVar2Source, LOVDSource, PubTatorSource,
        )
        _SOURCES = [
            LitVar2Source(),
            ClinVarSource(),
            GoogleCSESource(),
            EuropePMCSource(),
            PubTatorSource(),
            LOVDSource(),
        ]
    return _SOURCES


async def query_variant(hgvs: str) -> CitationTable:
    """
    Full pipeline: HGVS → CitationTable with PMIDs from all sources.

    All sources run concurrently via asyncio.gather.
    Individual source failures do not abort the pipeline.
    """
    # Step 1: LVG expansion (with cache)
    log.info("LVG expansion for %s", hgvs)
    lvg: LVGResult = await get_lvg(hgvs)

    # Step 2: Fan out to all sources in parallel
    sources = _get_sources()
    log.info("Querying %d sources in parallel for %s", len(sources), hgvs)
    tasks = [source.query(lvg) for source in sources]
    raw_results = await asyncio.gather(*tasks, return_exceptions=True)

    # Step 3: Normalize — convert exceptions to error SourceResults
    results: list[SourceResult] = []
    for source, outcome in zip(sources, raw_results):
        if isinstance(outcome, Exception):
            log.warning("Source %s failed: %s", source.source, outcome)
            results.append(SourceResult(source=source.source, pmids=[], error=str(outcome)))
        else:
            results.append(outcome)

    # Step 4: Merge
    table = CitationTable.from_source_results(hgvs, lvg, results)
    log.info(
        "Result for %s: %d unique PMIDs from %d sources",
        hgvs,
        len(table.citations),
        sum(1 for r in results if r.pmids),
    )
    return table
