"""
LOVD local source — curated variant-level literature references from Leiden
Open Variation Database instances.

Queries the local lovd.variant + lovd.variant_ref tables harvested from
LOVD instances worldwide. Matches the queried variant against LOVD's HGVS
notation using position-based matching (same approach as PubTator).

This is unique data — LOVD curators manually link variants to literature
references. These are not text-mined associations (like LitVar2/PubTator)
but human-curated citations.

Data source: medgen-stacks stacks/lovd/harvest.py + resolve_citations.py
"""
import asyncio
import logging
import re

from text2gene2.cache import cache_get, cache_set
from text2gene2.config import settings
from text2gene2.db import get_medgen_conn, reset_medgen_conn
from text2gene2.models import LVGResult, Source, SourceResult
from text2gene2.sources.base import PMIDSource

log = logging.getLogger(__name__)


def _extract_positions(lvg: LVGResult) -> list[str]:
    """Extract numeric positions from LVG for matching against LOVD HGVS."""
    positions = set()
    for h in lvg.hgvs_c:
        short = h.split(":")[-1] if ":" in h else h
        m = re.match(r"c\.(\d+)", short)
        if m:
            positions.add(m.group(1))
    for h in lvg.hgvs_p:
        short = h.split(":")[-1] if ":" in h else h
        m = re.search(r"[A-Z][a-z]{0,2}(\d+)", short)
        if m:
            positions.add(m.group(1))
    return list(positions)


def _extract_cdna_short(lvg: LVGResult) -> list[str]:
    """Extract short c.DNA forms for exact matching."""
    forms = []
    for h in lvg.hgvs_c:
        short = h.split(":")[-1] if ":" in h else h
        if short.startswith("c."):
            forms.append(short)
    return forms


def _query_lovd_sync(gene: str, cdna_forms: list[str],
                     positions: list[str]) -> tuple[list[int], dict[int, str]]:
    """
    Find PMIDs from lovd.variant_ref for variants matching this gene + HGVS.

    Two-tier matching:
      1. Exact c.DNA match (highest confidence)
      2. Position-based match (catches notation variants)
    """
    conn = get_medgen_conn()
    if conn is None:
        return [], {}

    pmids: list[int] = []
    provenance: dict[int, str] = {}
    seen = set()

    try:
        with conn.cursor() as cur:
            # Tier 1: Exact c.DNA match
            if cdna_forms:
                placeholders = ",".join(["%s"] * len(cdna_forms))
                cur.execute(f"""
                    SELECT DISTINCT vr.ref_id::int, v.hgvs_cdna, v.source_host
                    FROM lovd.variant_ref vr
                    JOIN lovd.variant v ON v.gene = vr.gene
                        AND v.hgvs_cdna = vr.hgvs_cdna
                        AND v.source_host = vr.source_host
                    WHERE vr.gene = %s
                      AND vr.ref_type = 'pmid'
                      AND v.hgvs_cdna IN ({placeholders})
                """, [gene] + cdna_forms)

                for pmid, hgvs, host in cur.fetchall():
                    if pmid not in seen:
                        seen.add(pmid)
                        pmids.append(pmid)
                        provenance[pmid] = f"LOVD exact: {hgvs}"

            # Tier 2: Position-based match (for notation variants)
            if positions:
                like_clauses = " OR ".join(
                    "v.hgvs_cdna LIKE %s" for _ in positions
                )
                like_params = [f"c.{pos}%" for pos in positions]

                cur.execute(f"""
                    SELECT DISTINCT vr.ref_id::int, v.hgvs_cdna, v.source_host
                    FROM lovd.variant_ref vr
                    JOIN lovd.variant v ON v.gene = vr.gene
                        AND v.hgvs_cdna = vr.hgvs_cdna
                        AND v.source_host = vr.source_host
                    WHERE vr.gene = %s
                      AND vr.ref_type = 'pmid'
                      AND ({like_clauses})
                """, [gene] + like_params)

                for pmid, hgvs, host in cur.fetchall():
                    if pmid not in seen:
                        seen.add(pmid)
                        pmids.append(pmid)
                        provenance[pmid] = f"LOVD position: {hgvs}"

    except Exception as e:
        log.warning("lovd: query error for %s: %s", gene, e)
        reset_medgen_conn()
        return [], {}

    return pmids, provenance


class LOVDSource(PMIDSource):
    source = Source.LOVD

    async def query(self, lvg: LVGResult) -> SourceResult:
        cache_key = f"lovd:{lvg.input_hgvs}"
        cached = await cache_get(cache_key)
        if cached is not None:
            if isinstance(cached, dict):
                prov = {int(k): v for k, v in cached.get("prov", {}).items()}
                return SourceResult(source=self.source, pmids=cached["pmids"],
                                    pmid_provenance=prov, cached=True)
            return SourceResult(source=self.source, pmids=cached, cached=True)

        if not lvg.gene_symbol:
            return SourceResult(source=self.source, pmids=[])

        cdna_forms = _extract_cdna_short(lvg)
        positions = _extract_positions(lvg)

        if not cdna_forms and not positions:
            return SourceResult(source=self.source, pmids=[])

        log.info("LOVD: gene=%s, cdna=%s, positions=%s", lvg.gene_symbol, cdna_forms, positions)
        pmids, provenance = await asyncio.to_thread(
            _query_lovd_sync, lvg.gene_symbol, cdna_forms, positions
        )

        prov_str = {str(k): v for k, v in provenance.items()}
        await cache_set(cache_key, {"pmids": pmids, "prov": prov_str},
                        ttl=settings.cache_ttl_litvar2)
        return SourceResult(source=self.source, pmids=pmids,
                            pmid_provenance=provenance)
