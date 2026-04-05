"""
LVG (Lexical Variant Group) expansion via the VariantValidator public REST API.

For a given HGVS string, VariantValidator returns all equivalent representations
across transcript versions, sequence types (c., g., p., n.), and genome builds.
This is the core differentiator vs. LitVar2 — it surfaces older literature that
used deprecated transcript versions or alternate notations.

Public API: https://rest.variantvalidator.org
Rate limit:  ≤ 2 req/sec (we use 1.5 to be safe)
"""
import json
import logging
import re

import httpx

from text2gene2.cache import cache_get, cache_set
from text2gene2.config import settings
from text2gene2.models import LVGResult
from text2gene2 import rate_limit

log = logging.getLogger(__name__)

_BASE = settings.variant_validator_url.rstrip("/")
_BUILD = settings.variant_validator_genome_build


def _cache_key(hgvs: str) -> str:
    return f"lvg:vv:{hgvs}"


def _parse_response(hgvs_input: str, data: dict) -> LVGResult:
    """Parse a VariantValidator JSON response into an LVGResult."""
    lvg = LVGResult(input_hgvs=hgvs_input)

    # Flag-only or error responses
    flag = data.get("flag", "")
    if flag not in ("gene_variant", "warning"):
        lvg.warnings.append(f"VariantValidator flag: {flag}")
        return lvg

    # Each top-level key (except 'flag' and 'metadata') is a transcript result
    for key, val in data.items():
        if key in ("flag", "metadata") or not isinstance(val, dict):
            continue

        # Gene symbol
        if not lvg.gene_symbol and val.get("gene_symbol"):
            lvg.gene_symbol = val["gene_symbol"]

        # Collect warnings
        for w in val.get("validation_warnings", []):
            if w and w not in lvg.warnings:
                lvg.warnings.append(w)

        # Collect transcript HGVS
        tx = val.get("hgvs_transcript_variant", "")
        if tx and tx not in lvg.hgvs_c:
            lvg.hgvs_c.append(tx)

        # Genomic HGVS — prefer the selected assembly
        loci = val.get("primary_assembly_loci", {})
        build_key = _BUILD.lower().replace("grch", "grch")  # grch38 / grch37
        for bk in [build_key, "grch38", "grch37"]:
            genomic = loci.get(bk, {}).get("hgvs_genomic_description", "")
            if genomic and genomic not in lvg.hgvs_g:
                lvg.hgvs_g.append(genomic)
                break

        # Protein HGVS
        prot = val.get("hgvs_predicted_protein_consequence", {})
        for pkey in ("slr", "tlr"):
            p = prot.get(pkey, "")
            if p and p not in lvg.hgvs_p:
                lvg.hgvs_p.append(p)

        # rsIDs embedded in cross-references (not always present)
        db_xref = val.get("annotations", {}).get("db_xref", {})
        for v in db_xref.values():
            if isinstance(v, str) and re.match(r"^rs\d+$", v):
                if v not in lvg.rsids:
                    lvg.rsids.append(v)

    return lvg


async def get_lvg(hgvs: str, use_cache: bool = True) -> LVGResult:
    """
    Expand an HGVS string into all equivalent representations.

    Caches results for settings.cache_ttl_lvg seconds (default 7 days).
    Falls back to a minimal LVGResult on API errors so the pipeline keeps running.
    """
    key = _cache_key(hgvs)

    if use_cache:
        cached = await cache_get(key)
        if cached:
            return LVGResult.model_validate(cached)

    await rate_limit.variant_validator.acquire()

    url = f"{_BASE}/VariantValidator/variantvalidator/{_BUILD}/{hgvs}/all"
    try:
        async with httpx.AsyncClient(timeout=httpx.Timeout(10.0, connect=5.0)) as client:
            resp = await client.get(url, params={"content-type": "application/json"})
            resp.raise_for_status()
            data = resp.json()
    except httpx.HTTPStatusError as e:
        log.warning("VariantValidator HTTP error for %s: %s", hgvs, e)
        return LVGResult(input_hgvs=hgvs, warnings=[str(e)])
    except Exception as e:
        log.warning("VariantValidator error for %s: %s", hgvs, e)
        return LVGResult(input_hgvs=hgvs, warnings=[str(e)])

    lvg = _parse_response(hgvs, data)

    if use_cache:
        await cache_set(key, lvg.model_dump(), ttl=settings.cache_ttl_lvg)

    return lvg
