# Data Stacks

All local data lives in the `medgen` PostgreSQL database on loki.local (192.168.0.3:5432,
user `medgen`, password `medgen`). Stacks are loaded by
[medgen-stacks](https://github.com/text2gene/medgen-stacks) — each is independent
and can be refreshed individually.

This document describes every stack, what text2gene2 uses it for today, and what
it could be used for next.

Last bulk load: **2026-03-27** (PubMed baseline through 2026-03-29).

---

## Currently Used by text2gene2

### HGNC (`hgnc` schema)

Human Gene Nomenclature Committee — the authority on gene symbols.

| Table | Rows | Size | Source |
|-------|------|------|--------|
| `hgnc.gene` | 44,981 | 13 MB | [EBI FTP](https://www.genenames.org/download/archive/) |

**Columns:** `hgnc_id`, `symbol`, `name`, `locus_group`, `locus_type`, `status`,
`location`, `alias_symbols`, `prev_symbols`, `entrez_id`, `ensembl_gene_id`,
`omim_id`, `uniprot_ids`, `refseq_accession`

**Used by:** `pipeline/expand.py` — gene synonym lookup for query expansion.
`alias_symbols` and `prev_symbols` (pipe-delimited) provide alternate gene names
so EuropePMC finds papers using old nomenclature (e.g., ACVRL1 → ALK1, HHT2).

**Last loaded:** 2026-03-27

---

### PubMed (`pubmed` schema)

MEDLINE article metadata — titles, authors, journals, abstracts for ~40M articles.

| Table | Rows | Size | Source |
|-------|------|------|--------|
| `pubmed.article` | 39,994,988 | 230 GB | [NLM MEDLINE baseline](https://ftp.ncbi.nlm.nih.gov/pubmed/baseline/) |

**Columns:** `pmid`, `title`, `authors`, `journal`, `year`, `doi`, `pmc_id`, `abstract`

**Used by:** `pipeline/enrich.py` — local metadata lookup before falling back to
NCBI efetch. Massively reduces API calls: most PMIDs are in the local DB already.

**Last loaded:** 2026-03-29 (through pubmed26n1334.xml.gz)

---

### PubTator (`pubtator` schema)

NCBI PubTator Central — NLP-extracted gene and mutation mentions in PubMed articles.

| Table | Rows | Size | Source |
|-------|------|------|--------|
| `pubtator.gene_mention` | 83,241,805 | 11 GB | [PubTator FTP](https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/) |
| `pubtator.mutation_mention` | 7,752,131 | 1.3 GB | Same |

**gene_mention columns:** `pmid`, `gene_id`, `mentions`, `resource`

**mutation_mention columns:** `pmid`, `concept_id`, `mentions`, `resource`

`concept_id` uses the tmVar format:
- `tmVar:c|SUB|G|1102|A` → c.1102G>A (coding substitution)
- `tmVar:p|SUB|E|289|G` → p.E289G (protein substitution)
- `tmVar:c|DEL|993+1_993+2|` → c.993+1_993+2del
- `tmVar:p|FS|K|239||21` → p.K239fsX21 (frameshift)
- `rs{number}` → dbSNP rsID

**Used by:** `sources/pubtator.py` — finds PMIDs where PubTator detected both the
gene and a mutation at the same coding/protein position as the queried variant.
Catches papers that LitVar2's autocomplete misses, plus papers about nearby mutations
at the same position (relevant functional context).

**Indexes:** `gene_mention(gene_id)`, `gene_mention(pmid)`, `mutation_mention(concept_id)`,
`mutation_mention(pmid)`

**Last loaded:** 2026-03-27

---

### Entrez Gene (`gene` schema)

NCBI Entrez Gene — gene records and gene-to-PubMed associations.

| Table | Rows | Size | Source |
|-------|------|------|--------|
| `gene.info` | 193,790 (human) | 42 MB | [NCBI Gene FTP](https://ftp.ncbi.nlm.nih.gov/gene/DATA/) |
| `gene.pubmed` | 2,180,227 | 179 MB | Same |
| `gene.history` | — | 13 MB | Same |

**gene.info columns:** `gene_id`, `tax_id`, `symbol`, `description`, `type_of_gene`,
`full_name`, `other_designations`, `dbxrefs`

**gene.pubmed columns:** `gene_id`, `pmid`, `tax_id`

**Used by:** `sources/pubtator.py` — resolves gene symbol → Entrez gene_id for
PubTator queries.

**Potential use:** `gene.pubmed` could serve as a confirmation signal (issue #9) —
if a PMID is in gene2pubmed for the same gene, it adds confidence. Also useful
for disease-based discovery: filter gene2pubmed papers through abstract validation.

**Last loaded:** 2026-03-27

---

### ClinVar (`clinvar` schema)

NCBI ClinVar — curated variant classifications and literature citations.

| Table | Rows | Size | Source |
|-------|------|------|--------|
| `clinvar.variant_summary` | 4,418,157 | 2.2 GB | [ClinVar FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/) |
| `clinvar.var_citations` | 4,487,623 | 542 MB | Same |

**variant_summary columns:** `allele_id`, `variation_id`, `name`, `gene_id`,
`gene_symbol`, `hgvs_cdna`, `hgvs_protein`, `molecular_consequence`,
`clinical_significance`, `clinsig_simple`, `review_status`, `last_evaluated`,
`phenotype_ids`, `phenotype_list`, `chromosome`, `start`, `stop`, `assembly`,
`dbsnp_id`, `updated_at`

**var_citations columns:** `allele_id`, `variation_id`, `rs_id`, `citation_source`,
`citation_id`

**Currently used by:** `sources/clinvar.py` queries ClinVar via the live eutils API
(esearch → elink). The local tables are **not yet queried directly** for PMID lookup.

**Potential use (issue #8):** Query `var_citations` directly by allele_id or rs_id
to bypass the API entirely. Faster, no rate limits, works offline. Would need
`variant_summary` for HGVS → allele_id resolution.

**Note:** `hgvs_cdna` format includes the gene name: `NM_000020.3(ACVRL1):c.717G>A (p.Trp239Ter)`.
Match on gene_symbol + HGVS fragment, not exact string match.

**Last loaded:** 2026-03-27

---

### LOVD (`lovd` schema)

Leiden Open Variation Database — index of gene-specific locus databases worldwide.

| Table | Rows | Size | Source |
|-------|------|------|--------|
| `lovd.lsdb` | 84,660 | 39 MB | [grenada.lumc.nl](https://grenada.lumc.nl/LSDB_list/lsdbs) |
| `lovd.scrape_log` | 44,978 | 3.6 MB | Same |

**lsdb columns:** `gene`, `name`, `url`, `curator`, `n_variants`, `db_type`, `last_seen`

**db_type breakdown:** 54K LOVD 3.X, 28K ClinVar, 1.5K Other, 500 LOVD 2.X

**Not yet used.** This is the index for LOVD integration (issue #6). Maps genes
to their LOVD instance URLs. The actual LOVD instances contain variant-level
literature references (PMIDs + DOIs) that would need to be queried via the
LOVD REST API or HTML scraping (see `metavariant/metavariant/lovd.py` for
prior art).

**Last scraped:** 2026-03-29

---

## Available But Not Yet Used

### MedGen (`medgen` schema)

NCBI MedGen — disease/condition concepts with PubMed associations.

| Table | Rows | Size | Source |
|-------|------|------|--------|
| `medgen.concept` | 227,999 | 72 MB | [MedGen FTP](https://ftp.ncbi.nlm.nih.gov/pub/medgen/) |
| `medgen.concept_name` | — | 171 MB | Same |
| `medgen.concept_rel` | — | 117 MB | Same |
| `medgen.pubmed` | 47,822,130 | 7.6 GB | Same |

**medgen.pubmed columns:** `uid`, `pmid` (MedGen concept UID → PubMed PMID)

**Potential use (issue #10):** Disease-based PMID discovery. If we know a variant
is associated with a disease (from ClinVar's `phenotype_ids`), we can find all
PubMed articles linked to that disease concept via `medgen.pubmed`, then filter
through abstract validation. This is the "wide net, strict filter" approach —
47.8M disease→PMID links, many of which will discuss specific variants in
natural language rather than HGVS notation.

**Data chain:** variant → ClinVar phenotype → MedGen CUI → medgen.pubmed → PMIDs → validate

**Last loaded:** 2026-03-27

---

### HPO (`hpo` schema)

Human Phenotype Ontology — standardized phenotype terms linked to genes and diseases.

| Table | Rows | Size | Source |
|-------|------|------|--------|
| `hpo.term` | 19,944 | 5 MB | [HPO/JAX](https://hpo.jax.org/data/annotations) |
| `hpo.term_parent` | — | 3.3 MB | Same |
| `hpo.gene_phenotype` | 1,080,403 | 99 MB | Same |
| `hpo.disease_phenotype` | 282,127 | 50 MB | Same |

**gene_phenotype columns:** `gene_id`, `gene_symbol`, `hpo_term_id`, `hpo_term_name`,
`frequency`, `disease_id`

**Potential use:** Phenotype-aware query expansion. If a variant is in a gene
associated with specific HPO phenotypes, we can search for papers mentioning
those phenotype terms + the gene. E.g., for ACVRL1: search for "arteriovenous
malformation" + "ACVRL1" to find clinical case reports that describe the
phenotype without using HGVS notation.

Also useful for disease-based discovery (issue #10) as an alternate path:
gene → HPO phenotypes → MedGen concepts → medgen.pubmed.

**Last loaded:** 2026-03-27

---

### Orphanet (`orphanet` schema)

Orphanet — European rare disease database with gene associations.

| Table | Rows | Size | Source |
|-------|------|------|--------|
| `orphanet.disorder` | 11,456 | 1.4 MB | [Orphadata](https://www.orphadata.com/genes/) |
| `orphanet.gene_association` | 8,341 | 1.5 MB | Same |

**disorder columns:** `orpha_code`, `name`, `expert_link`

**gene_association columns:** `orpha_code`, `gene_symbol`, `gene_name`, `gene_id`,
`association_type`, `association_status`

**Potential use:** Map genes to Orphanet rare disease entries. Could enhance
disease-based discovery (issue #10) with Orphanet's curated gene-disease
associations, which are more rare-disease-focused than ClinVar's phenotype list.

**Last loaded:** 2026-03-27

---

### PMC (`pmc` schema)

PubMed Central — full-text open-access articles.

| Table | Rows | Size | Source |
|-------|------|------|--------|
| `pmc.article` | 5,103,160 | 191 GB | [PMC OA bulk download](https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/) |

**Potential use (issue #11):** Full-text validation. Current validation only
searches titles and abstracts. Many papers mention specific variants only in
methods, results tables, or supplementary data. Full-text search against
`pmc.article` would catch these.

**Note:** 191 GB — this is the single largest table. Need to verify the schema
has searchable text content (vs. just metadata) before planning full-text search.

**Last loaded:** 2026-03-29

---

### DisGeNET (`disgenet` schema)

Gene-disease and variant-disease associations from multiple curated sources.

| Table | Rows | Size | Source |
|-------|------|------|--------|
| `disgenet.gene_disease` | 0 | 40 kB | [DisGeNET](https://www.disgenet.org/) |
| `disgenet.variant_disease` | 0 | 32 kB | Same |

**Status: Empty.** The tables exist but have no data loaded. DisGeNET requires
free registration to download data files. This stack needs the data files
downloaded and loaded before it can be used.

**Potential use:** `variant_disease` maps specific variants (by rsID) to disease
associations with evidence scores. Could provide an alternate confirmation
signal and disease context for variants.

---

## External APIs (not in Postgres)

These are queried live over HTTP, not from local data:

| Source | Used by | Rate limit | Purpose |
|--------|---------|------------|---------|
| VariantValidator | `lvg/expander.py` | 1.5 req/sec | HGVS → all equivalent representations |
| LitVar2 | `sources/litvar2.py` | 5 req/sec | Text-mined variant→PMID via rsID |
| ClinVar eutils | `sources/clinvar.py` | 3 req/sec | Curated variant→PMID (esearch + elink) |
| Europe PMC | `sources/europepmc.py` | 5 req/sec | Open-access literature search |
| Google CSE | `sources/google_cse.py` | 1 req/sec | Publisher sites (not configured) |
| NCBI efetch | `pipeline/enrich.py` | 3 req/sec | Article metadata fallback |

See `DATA_SOURCES.md` for full API documentation.

---

## Refresh Schedule

All stacks are loaded manually via `medgen-stacks`. Recommended refresh cadence:

| Stack | Cadence | Why |
|-------|---------|-----|
| ClinVar | Monthly | New variant classifications and citations |
| PubMed | Monthly | New articles (baseline updates monthly) |
| PubTator | Monthly | Re-processed NLP on new articles |
| HGNC | Quarterly | Gene names change slowly |
| Gene | Quarterly | Stable |
| HPO | Quarterly | Phenotype ontology updates |
| Orphanet | Quarterly | Rare disease catalog updates |
| MedGen | Quarterly | Concept→PMID links |
| LOVD | Monthly | Scrape `grenada.lumc.nl` for new instances |
| PMC | As needed | Full-text for validation only |

To refresh a stack:
```bash
cd ~/projects/git/medgen-stacks
source .env
make <stack_name>   # e.g., make clinvar, make pubtator
```
