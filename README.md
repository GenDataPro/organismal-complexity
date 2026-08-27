# Organismal Complexity

Computational analysis of the relationship between protein functional diversity and organismal complexity across metazoan species.

This repository preserves the original Python implementation associated with the 2017 PLOS ONE study, together with a later R/BioMart redevelopment and downstream analyses.

## Published study

This work originated in:

**Lopes Cardoso D, Sharpe C (2017).**  
*Relating protein functional diversity to cell type number identifies genes that determine dynamic aspects of chromatin organisation as potential contributors to organismal complexity.*  
**PLOS ONE 12(9): e0185409.**  
DOI: 10.1371/journal.pone.0185409

The paper originally linked the analysis scripts to:

https://github.com/GenDataPro/GenDataPro

They have since been reorganised and preserved in this dedicated repository.

## Published workflow

The published analysis began with **2,308 human transcription-associated genes** identified in AnimalTFDB 2.0.

Comparison with Ensembl resolved these to **2,087 distinct human genes with Ensembl gene identifiers**.

A Python script then accessed **Ensembl Release 87 through the REST API** to identify orthologues in:

- *Caenorhabditis elegans*
- *Drosophila melanogaster*
- *Ciona intestinalis*
- *Takifugu rubripes*
- *Xenopus tropicalis*
- *Gallus gallus*
- *Mus musculus*
- *Macaca mulatta*

with *Homo sapiens* as the reference species.

### Functional-diversity scoring

The original Python algorithm quantified three features for each gene:

- **Paralogues (P)** — within-species paralogues
- **Isoforms (I)** — transcripts annotated as protein coding
- **Motifs (M)** — Prosite Profile annotations collected across protein isoforms

These measures were used to calculate a functional-diversity score for orthologous genes across species and to investigate its relationship with cell-type number.

## Repository structure

    organismal-complexity/
    ├── README.md
    ├── python-rest/
    │   ├── scripts/
    │   ├── data/
    │   └── results/
    └── r-biomart/
        ├── scripts/
        ├── results/
        └── downstream-analysis/

## python-rest

This directory preserves the historical Python implementation associated with the published study.

### Scripts

- `ortholog_search.py` — queries the Ensembl REST API to identify orthologous genes
- `complexity_scoring_algorithm.py` — retrieves within-species paralogues, protein-coding isoforms and Prosite Profile motifs used in complexity scoring

### Data

Historical input and intermediate files used by the Python workflow are retained under `data/`.

### Results

Historical outputs generated during the original analysis are retained under `results/`.

The historical code has been preserved as closely as possible rather than retrospectively rewritten. It therefore contains dependencies, file paths and API usage reflecting the environment in which the analysis was originally performed.

## r-biomart

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

This repository has been reorganised to preserve the historical analysis while making its scientific provenance clearer.

The original Python implementation is retained separately from the later R/BioMart redevelopment. Historical files have been preserved without silent modification wherever possible.

Some historical scripts depend on earlier Ensembl interfaces, software versions and local file structures and are therefore not expected to run unchanged in a modern environment.

A future redevelopment may modernise the workflow for reproducibility using current Ensembl interfaces, dependency management and documented analysis pipelines.
