import logging
from pathlib import Path

from fastapi import APIRouter, HTTPException, Request
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates

from text2gene2.lvg import get_lvg
from text2gene2.models import CitationTable, LVGResult
from text2gene2.pipeline import query_variant
from text2gene2.pipeline.validate import validate_citations

log = logging.getLogger(__name__)
router = APIRouter()

_templates_dir = Path(__file__).parent.parent / "templates"
templates = Jinja2Templates(directory=str(_templates_dir))


# ── Web UI ──────────────────────────────────────────────────────────────────

@router.get("/", response_class=HTMLResponse)
async def index(request: Request):
    return templates.TemplateResponse(request=request, name="index.html")


@router.get("/search", response_class=HTMLResponse)
async def search_page(request: Request, hgvs: str = ""):
    ctx: dict = {"hgvs": hgvs, "result": None, "error": None}
    if hgvs:
        try:
            table = await query_variant(hgvs.strip())
            ctx["result"] = table
        except Exception as e:
            log.exception("Pipeline error for %s", hgvs)
            ctx["error"] = str(e)
    return templates.TemplateResponse(request=request, name="search.html", context=ctx)


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
