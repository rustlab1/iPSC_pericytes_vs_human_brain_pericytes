---
editor_options: 
  markdown: 
    wrap: 72
---

# README

## Data sources

This task is based on publicly available RNA-seq data from [Comparative
transcriptomics of iPSC-derived neural crest pericytes and primary human
brain vascular pericytes] (GEO: GSE252046, BioProject: PRJNA1057204).
The dataset includes primary pericytes and NCC-derived iPSC pericytes.
Processed count data is downloaded directly from GEO for use in the
workflow.

## How to download Data acquired from GEO:

Install SRA Toolkit and supporting tools with conda, then fetch runs
with `prefetch` and build per sample FASTQs via `fasterq dump`. The
preprocessing script streams each run and subsamples deterministically
to one million reads per sample using `seqtk sample`, compresses with
`pigz`, and avoids storing full FASTQs on disk

## Pre-processing 

filtering Filter genes with fewer than 10 counts in at least the size of
the smallest group. Log transform for visualization (log10(count+1)) and
run UMAP for dimensionality reduction.

## How the workflow works 

Workflow files are stored in workflow/.

Step 1 – Load and filter\
Purpose: Remove low-count genes and prepare data for analysis\
Tools: R, data.table\
Inputs: Count matrix from GEO

Step 2 – Visualization\
Purpose: Explore sample clustering with UMAP\
Tools: R, umap\
\

Step 3 – Differential expression analysis\
Purpose: Identify genes that differ between primary PC and NCC-PC\
Tools: R, DESeq2\
Outputs: DE results table, volcano plot, MD plot, heatmap of top DE
genes

Step 4 – GO enrichment\
Purpose: Find enriched Biological Process terms among DE genes\
Tools: R, clusterProfiler\
Inputs: Significant DE genes \
Outputs: Table of enriched terms

Notes You can skip raw FASTQ processing and start from the counts matrix
if you already have one.
