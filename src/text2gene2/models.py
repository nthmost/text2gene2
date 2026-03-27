"""Core data models for text2gene2."""
from enum import Enum
from typing import Any
from pydantic import BaseModel, Field


class Source(str, Enum):
    LITVAR2 = "litvar2"
    CLINVAR = "clinvar"
    GOOGLE = "google"
    EUROPEPMC = "europepmc"


class LVGResult(BaseModel):
    """Lexical Variant Group — all equivalent HGVS representations of a variant."""
    input_hgvs: str
    gene_symbol: str | None = None
    hgvs_c: list[str] = Field(default_factory=list)   # coding DNA
    hgvs_g: list[str] = Field(default_factory=list)   # genomic
    hgvs_p: list[str] = Field(default_factory=list)   # protein
    hgvs_n: list[str] = Field(default_factory=list)   # non-coding RNA
    rsids: list[str] = Field(default_factory=list)     # dbSNP rsIDs found
    warnings: list[str] = Field(default_factory=list)

    @property
    def all_hgvs(self) -> list[str]:
        """All unique HGVS strings across all sequence types."""
        seen = set()
        result = []
        for h in [self.input_hgvs] + self.hgvs_c + self.hgvs_g + self.hgvs_p + self.hgvs_n:
            if h and h not in seen:
                seen.add(h)
                result.append(h)
        return result


class SourceResult(BaseModel):
    """PMIDs returned by a single source."""
    source: Source
    pmids: list[int] = Field(default_factory=list)
    error: str | None = None
    cached: bool = False


class Citation(BaseModel):
    """A single PubMed citation with source attribution."""
    pmid: int
    sources: list[Source] = Field(default_factory=list)
    # Populated lazily by metapub or similar
    title: str | None = None
    authors: str | None = None
    journal: str | None = None
    year: int | None = None
    doi: str | None = None
    # Set by pipeline/validate.py — "confirmed" | "probable" | "unverified"
    validation_tier: str | None = None

    @property
    def confidence(self) -> int:
        """Number of independent sources that found this PMID — higher is better."""
        return len(self.sources)


class CitationTable(BaseModel):
    """Aggregated results for a variant query."""
    input_hgvs: str
    lvg: LVGResult | None = None
    by_source: dict[Source, list[int]] = Field(default_factory=dict)
    citations: list[Citation] = Field(default_factory=list)
    errors: dict[Source, str] = Field(default_factory=dict)

    @property
    def all_pmids(self) -> list[int]:
        return [c.pmid for c in self.citations]

    @classmethod
    def from_source_results(cls, hgvs: str, lvg: LVGResult | None, results: list[SourceResult]) -> "CitationTable":
        pmid_sources: dict[int, list[Source]] = {}
        errors: dict[Source, str] = {}
        by_source: dict[Source, list[int]] = {}

        for r in results:
            by_source[r.source] = r.pmids
            if r.error:
                errors[r.source] = r.error
            for pmid in r.pmids:
                pmid_sources.setdefault(pmid, []).append(r.source)

        citations = [
            Citation(pmid=pmid, sources=sources)
            for pmid, sources in sorted(pmid_sources.items(), key=lambda x: -len(x[1]))
        ]

        return cls(input_hgvs=hgvs, lvg=lvg, by_source=by_source, citations=citations, errors=errors)


class BatchJob(BaseModel):
    job_id: str
    status: str = "pending"   # pending | running | done | error
    total: int = 0
    completed: int = 0
    variants: list[str] = Field(default_factory=list)
    results: dict[str, Any] = Field(default_factory=dict)
