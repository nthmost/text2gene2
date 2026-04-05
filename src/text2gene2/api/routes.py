import asyncio
import json
import logging
from pathlib import Path

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import HTMLResponse
from starlette.responses import StreamingResponse
from fastapi.templating import Jinja2Templates

from text2gene2.lvg import get_lvg
from text2gene2.models import CitationTable, Citation, LVGResult, Source, SourceResult
from text2gene2.pipeline import query_variant
from text2gene2.pipeline.enrich import enrich_citations
from text2gene2.pipeline.expand import expand_variant, get_gene_synonyms
from text2gene2.pipeline.validate import validate_citations
from text2gene2.sources import ClinVarSource, EuropePMCSource, GoogleCSESource, LitVar2Source, PubTatorSource

log = logging.getLogger(__name__)
router = APIRouter()

_templates_dir = Path(__file__).parent.parent / "templates"
templates = Jinja2Templates(directory=str(_templates_dir))


# ── Web UI ──────────────────────────────────────────────────────────────────

# Curated demo variants: rare/orphan disease P/LP, diverse conditions
_SAMPLE_VARIANTS = [
    {"hgvs": "NM_000020.3:c.717G>A",      "gene": "ACVRL1", "disease": "Hereditary hemorrhagic telangiectasia"},
    {"hgvs": "NM_000492.4:c.1573del",     "gene": "CFTR",   "disease": "Cystic fibrosis"},
    {"hgvs": "NM_000138.5:c.4710G>A",     "gene": "FBN1",   "disease": "Marfan syndrome"},
    {"hgvs": "NM_000277.3:c.866G>A",      "gene": "PAH",    "disease": "Phenylketonuria"},
    {"hgvs": "NM_000169.3:c.160C>T",      "gene": "GLA",    "disease": "Fabry disease"},
    {"hgvs": "NM_000152.5:c.1846dup",     "gene": "GAA",    "disease": "Pompe disease"},
    {"hgvs": "NM_000053.4:c.3270C>G",     "gene": "ATP7B",  "disease": "Wilson disease"},
    {"hgvs": "NM_001110792.2:c.337C>T",   "gene": "MECP2",  "disease": "Rett syndrome"},
    {"hgvs": "NM_000520.6:c.465del",      "gene": "HEXA",   "disease": "Tay-Sachs disease"},
    {"hgvs": "NM_000350.3:c.4312C>G",     "gene": "ABCA4",  "disease": "Stargardt disease"},
    {"hgvs": "NM_000551.4:c.444del",      "gene": "VHL",    "disease": "Von Hippel-Lindau syndrome"},
    {"hgvs": "NM_001114753.3:c.715dup",   "gene": "ENG",    "disease": "Hereditary hemorrhagic telangiectasia"},
]

@router.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return templates.TemplateResponse(
        request=request, name="index.html",
        context={"sample_variants": _SAMPLE_VARIANTS},
    )


@router.get("/search", response_class=HTMLResponse)
async def search_page(request: Request, hgvs: str = ""):
    """Render the search page shell. Results stream in via SSE."""
    return templates.TemplateResponse(
        request=request, name="search.html", context={"hgvs": hgvs},
    )


def _sse(event: str, data: dict) -> str:
    """Format a Server-Sent Event."""
    return f"event: {event}\ndata: {json.dumps(data)}\n\n"


_SOURCES = [
    LitVar2Source(),
    ClinVarSource(),
    GoogleCSESource(),
    EuropePMCSource(),
    PubTatorSource(),
]


@router.get("/api/v2/stream/{hgvs_text:path}")
async def stream_search(hgvs_text: str):
    """SSE endpoint: streams pipeline results as they complete."""

    async def event_generator():
        hgvs = hgvs_text.strip()

        # Step 1: LVG expansion
        try:
            lvg = await get_lvg(hgvs)
        except Exception as e:
            yield _sse("error", {"message": f"LVG error: {e}"})
            yield _sse("done", {})
            return

        lvg_data = lvg.model_dump()
        yield _sse("lvg", lvg_data)

        # Step 2: Query expansion
        exp = expand_variant(lvg)
        gene_syns = []
        if lvg.gene_symbol:
            gene_syns = await get_gene_synonyms(lvg.gene_symbol)
        yield _sse("expansion", {
            "gene_symbol": exp.gene_symbol,
            "coding_short": exp.coding_short,
            "protein_short": exp.protein_short,
            "protein_1letter": exp.protein_1letter,
            "genomic_short": exp.genomic_short,
            "slang": exp.slang,
            "rsids": exp.rsids,
            "all_search_forms": exp.all_search_forms,
            "gene_synonyms": gene_syns,
        })

        # Step 3: Fan out to sources — stream each as it completes
        all_results: list[SourceResult] = []

        async def run_source(source):
            try:
                return await source.query(lvg)
            except Exception as e:
                log.warning("Source %s failed: %s", source.source, e)
                return SourceResult(source=source.source, pmids=[], error=str(e))

        tasks = {asyncio.create_task(run_source(s)): s for s in _SOURCES}

        for coro in asyncio.as_completed(tasks):
            result = await coro
            all_results.append(result)
            yield _sse("source", {
                "source": result.source.value,
                "pmids": result.pmids,
                "error": result.error,
                "cached": result.cached,
                "query_used": result.query_used,
                "pmid_provenance": {str(k): v for k, v in result.pmid_provenance.items()},
            })

        # Step 4: Merge and enrich
        table = CitationTable.from_source_results(hgvs, lvg, all_results)
        table = await enrich_citations(table)

        # Stream final citations
        citations_data = [c.model_dump() for c in table.citations]
        # Convert Source enums to strings for JSON
        for cd in citations_data:
            cd["sources"] = [s.value if hasattr(s, "value") else s for s in cd["sources"]]
            cd["found_by"] = {
                (k.value if hasattr(k, "value") else k): v
                for k, v in cd["found_by"].items()
            }
        yield _sse("citations", {"citations": citations_data, "total": len(table.citations)})
        yield _sse("done", {})

    return StreamingResponse(
        event_generator(),
        media_type="text/event-stream",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )


# ── JSON API ─────────────────────────────────────────────────────────────────

@router.get("/api/v2/hgvs2pmid/{hgvs_text:path}", response_model=CitationTable)
async def hgvs2pmid(hgvs_text: str, validate: bool = False) -> CitationTable:
    """
    Full pipeline: HGVS → PMIDs from all sources.

    Pass `?validate=true` to fetch abstracts and tier results as
    confirmed / probable / unverified based on variant mention in text.
    Adds ~2-5s latency for the abstract fetch.
    """
    try:
        table = await query_variant(hgvs_text)
        if validate:
            table = await validate_citations(table)
        return table
    except Exception as e:
        log.exception("Pipeline error")
        raise HTTPException(status_code=500, detail=str(e))


@router.get("/api/v2/lvg/{hgvs_text:path}", response_model=LVGResult)
async def lvg(hgvs_text: str) -> LVGResult:
    """Return the Lexical Variant Group (all equivalent HGVS forms) for a variant."""
    return await get_lvg(hgvs_text)


@router.get("/api/v2/health")
async def health():
    return {"status": "ok"}
