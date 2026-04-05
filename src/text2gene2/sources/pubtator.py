"""
PubTator local source — text-mined mutation-to-PMID associations from NCBI PubTator Central.

Queries the local `pubtator.mutation_mention` + `pubtator.gene_mention` tables
on the medgen Postgres DB. This is essentially LitVar2's underlying data but
queryable locally by position, not just rsID — catching papers that LitVar2's
autocomplete misses.

Strategy:
  1. Resolve gene symbol → Entrez gene_id via gene.info
  2. Extract key positions from the LVG (coding + protein)
  3. Find PMIDs where gene_mention matches the gene AND mutation_mention
     has a tmVar concept_id matching the variant position
  4. Track which tmVar concept matched for provenance

tmVar concept_id format examples:
  tmVar:c|SUB|G|1102|A    → c.1102G>A (coding substitution)
  tmVar:p|SUB|E|289|G     → p.E289G (protein substitution)
  tmVar:c|DEL|993+1_993+2 → c.993+1_993+2del
  tmVar:p|FS|K|239||21    → p.K239fsX21
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

# Cache gene symbol → entrez gene_id
_gene_id_cache: dict[str, int | None] = {}


def _get_gene_id_sync(symbol: str) -> int | None:
    """Resolve HGNC symbol to Entrez gene_id via gene.info."""
    if symbol in _gene_id_cache:
        return _gene_id_cache[symbol]

    conn = get_medgen_conn()
    if conn is None:
        return None
    try:
        with conn.cursor() as cur:
            cur.execute(
                "SELECT gene_id FROM gene.info WHERE symbol = %s AND tax_id = 9606",
                (symbol,),
            )
            row = cur.fetchone()
            gene_id = row[0] if row else None
            _gene_id_cache[symbol] = gene_id
            return gene_id
    except Exception as e:
        log.warning("pubtator: gene_id lookup failed for %s: %s", symbol, e)
        reset_medgen_conn()
        return None


def _extract_positions(lvg: LVGResult) -> list[str]:
    """
    Extract numeric positions from LVG HGVS forms for matching against
    tmVar concept_ids. Returns position strings like "1102", "289", "717".
    """
    positions = set()

    for h in lvg.hgvs_c:
        short = h.split(":")[-1] if ":" in h else h
        # Extract leading position from c.1102G>A, c.1573del, c.4710G>A, etc.
        m = re.match(r"c\.(\d+)", short)
        if m:
            positions.add(m.group(1))

    for h in lvg.hgvs_p:
        short = h.split(":")[-1] if ":" in h else h
        # Extract position from p.(Trp239Ter), p.(Gly289Glu), etc.
        m = re.search(r"[A-Z][a-z]{0,2}(\d+)", short)
        if m:
            positions.add(m.group(1))

    return list(positions)


def _query_pubtator_sync(gene_id: int, positions: list[str]) -> tuple[list[int], dict[int, str]]:
    """
    Find PMIDs where PubTator detected both the gene and a mutation at one of
    the given positions. Returns (pmids, provenance_dict).
    """
    conn = get_medgen_conn()
    if conn is None or not positions:
        return [], {}

    # Build concept_id patterns: match any position in the tmVar concept
    # e.g. position "1102" matches "tmVar:c|SUB|G|1102|A" or "tmVar:p|SUB|A|1102|T"
    # We use OR across positions, and require gene co-mention
    position_clauses = " OR ".join(
        "mm.concept_id LIKE %s" for _ in positions
    )
    # Pattern: %|{pos}|% matches position field in tmVar concept_ids
    # Also try %|{pos}% for end-of-string positions (deletions, etc.)
    params: list = [gene_id]
    like_patterns = []
    for pos in positions:
        like_patterns.append(f"%|{pos}|%")
        like_patterns.append(f"%|{pos}")
    params.extend(like_patterns)

    clause_parts = " OR ".join("mm.concept_id LIKE %s" for _ in like_patterns)

    sql = f"""
    SELECT DISTINCT mm.pmid, mm.concept_id, mm.mentions
    FROM pubtator.gene_mention gm
    JOIN pubtator.mutation_mention mm ON mm.pmid = gm.pmid
    WHERE gm.gene_id = %s
      AND mm.concept_id LIKE 'tmVar:%%'
      AND ({clause_parts})
    ORDER BY mm.pmid
    LIMIT 200
    """

    try:
        with conn.cursor() as cur:
            cur.execute(sql, params)
            rows = cur.fetchall()
    except Exception as e:
        log.warning("pubtator: query error: %s", e)
        reset_medgen_conn()
        return [], {}

    pmids = []
    provenance: dict[int, str] = {}
    seen = set()
    for pmid, concept_id, mentions in rows:
        if pmid not in seen:
            seen.add(pmid)
            pmids.append(pmid)
            # Use the mention text for provenance (more readable than concept_id)
            display = mentions.split("|")[0] if mentions else concept_id
            provenance[pmid] = display or concept_id

    return pmids, provenance


class PubTatorSource(PMIDSource):
    source = Source.PUBTATOR

    async def query(self, lvg: LVGResult) -> SourceResult:
        cache_key = f"pubtator:{lvg.input_hgvs}"
        cached = await cache_get(cache_key)
        if cached is not None:
            if isinstance(cached, dict):
                prov = {int(k): v for k, v in cached.get("prov", {}).items()}
                return SourceResult(source=self.source, pmids=cached["pmids"],
                                    pmid_provenance=prov, cached=True)
            return SourceResult(source=self.source, pmids=cached, cached=True)

        if not lvg.gene_symbol:
            return SourceResult(source=self.source, pmids=[])

        gene_id = await asyncio.to_thread(_get_gene_id_sync, lvg.gene_symbol)
        if gene_id is None:
            log.debug("pubtator: no gene_id for %s", lvg.gene_symbol)
            return SourceResult(source=self.source, pmids=[])

        positions = _extract_positions(lvg)
        if not positions:
            log.debug("pubtator: no positions extracted from LVG for %s", lvg.input_hgvs)
            return SourceResult(source=self.source, pmids=[])

        log.info("PubTator: gene_id=%d, positions=%s", gene_id, positions)
        pmids, provenance = await asyncio.to_thread(_query_pubtator_sync, gene_id, positions)

        prov_str = {str(k): v for k, v in provenance.items()}
        await cache_set(cache_key, {"pmids": pmids, "prov": prov_str},
                        ttl=settings.cache_ttl_litvar2)  # same TTL as litvar2
        return SourceResult(source=self.source, pmids=pmids,
                            pmid_provenance=provenance)
