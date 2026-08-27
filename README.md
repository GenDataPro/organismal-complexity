```

## `python-rest/`

This directory preserves the historical Python implementation associated with the published study.

### Scripts

- `ortholog_search.py` — queries the Ensembl REST API to identify orthologous genes
- `complexity_scoring_algorithm.py` — retrieves within-species paralogues, protein-coding isoforms and Prosite Profile motifs used in complexity scoring

### Data

Historical input and intermediate files used by the Python workflow are retained under `data/`.

### Results

Historical outputs generated during the original analysis are retained under `results/`.

The historical code has been preserved as closely as possible rather than retrospectively rewritten. It therefore contains dependencies, file paths and API usage reflecting the environment in which the analysis was originally performed.

## `r-biomart/`

A later implementation redeveloped the analysis in R using `biomaRt`.

The R workflow retrieves protein-coding genes and orthologues across nine species, quantifies paralogues, protein domains and isoforms, and constructs cross-species functional-diversity comparisons.

### Primary result

`results/Complexity_Score.xlsx`

Primary gene-level output from the R/BioMart analysis, containing the component measurements and derived functional-diversity scores across the analysed species.

### Downstream analysis

`downstream-analysis/` contains subsequent data cleaning, cross-species comparison and biological interpretation.

- `Complexity_Spread 0.05.xlsx` — human-centred cross-species functional-diversity matrix, including correlation analyses using different levels of species representation
- `GO term analysis.xlsx` — Gene Ontology biological-process enrichment following candidate-gene filtering
- `STRING_gene_list.csv` — curated candidate-gene list used during downstream functional interpretation

## Repository status

This repository is being reorganised to preserve the historical analysis while making its scientific provenance clearer.

The original Python implementation is retained separately from the later R/BioMart redevelopment. Historical files are being preserved without silent modification.

A future redevelopment may modernise the workflow for reproducibility using current Ensembl interfaces, dependency management and documented analysis pipelines.
