This repository contains R code for the analysis and validation of RNA-sequencing–derived 
findings using spatial transcriptomics (CODEX) data and also publicly available 
single-cell RNA-sequencing datasets.

The primary aim of the codebase is to support integrative analyses that 
compare, validate, and contextualise bulk RNA-seq results within 
spatial and single-cell frameworks, with a particular focus on cellular 
composition and transcriptional programmes.


## Overview

The analyses implemented in this repository include:

- Processing and analysis of spatial transcriptomics data
- Identification and visualisation of marker overlaps and spatial patterns
- Cell–cell interaction analysis in spatial contexts
- Integration of public single-cell RNA-sequencing datasets
- Generation of reference single-cell atlases
- Cell labelling and deconvolution-based analyses
- Comparison of pseudo-bulk and bulk RNA-seq–derived signals

The workflow is organised into sequential processing scripts, supported by reusable functions.

## Repository structure

/  
 ├── Data/                     # Input data (intentionally empty)
 ├── Outputs/                  # Generated results and figures (intentionally empty)
 ├── Process/                  # Contains the main analysis scripts
 │   ├── functions/            # Reusable helper functions
 ├── renv/                     # renv infrastructure
 └── renv.lock                 # Reproducible package environment

### Data and outputs

The `Data/` and `Outputs/` directories are intentionally left blank in this 
repository to reduce storage requirements and avoid redistributing
large or restricted datasets. Users are expected to populate these 
directories locally with the appropriate input data and generated results.

---

## Workflow organisation

Scripts in the `Process/` directory are numbered to reflect the intended 
analytical order, from initial data loading through to integrative and 
comparative analyses. Individual scripts may be run independently, but later 
steps typically assume the successful completion of earlier stages.

Reusable functions are stored separately in the `functions/` directory and 
are sourced as needed by the processing scripts.

---

## Reproducibility

This project uses **renv** to manage R package dependencies and ensure reproducibility.

To restore the computational environment:

```r
renv::restore()