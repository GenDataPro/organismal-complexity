# Organismal Complexity

Computational analysis of the relationship between protein functional diversity and organismal complexity across metazoan species.

This repository preserves the original Python implementation associated with the 2017 PLOS ONE study, together with a later R/BioMart redevelopment and downstream analyses.

## Published study

This work originated in:

**Lopes Cardoso D, Sharpe C (2017).**  
*Relating protein functional diversity to cell type number identifies genes that determine dynamic aspects of chromatin organisation as potential contributors to organismal complexity.*  
**PLOS ONE 12(9): e0185409.**  
DOI: `10.1371/journal.pone.0185409`

The paper originally linked the analysis scripts to `github.com/GenDataPro/GenDataPro`.  
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

```text
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
