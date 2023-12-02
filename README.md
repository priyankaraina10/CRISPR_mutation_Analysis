# Project 
For every (mutation, gene) pair, we wish to know whether the presence of that variant is associated to a change in the outcome of the CRISPR knock-out experiment for that gene, across all cell lines.

## Nextflow Pipeline

This pipeline performs statistical test across all (gene, mutation) pairs 

## Usage
`./nextflow run ../models_mutations.groovy --root_dir ./ --target_batch_limit 2`

## Input Files

### Gene_KOs.tsv

This tab separated value table have 11 columns and 19 rows. First column have model names and First row have gene knock-out(KO)s names 

This file has data from CRISPR single gene knock-out screens across a similar collection of cell lines, and measured for each cell line the cell count fold change depending on the knocked-out gene. These counts were subsequently quantile normalised.

### Mutations.tsv

This tab separated value table have 19 columns and 11 rows. First column have gene mutation names and First row have model names. This file has data for the presence (1) or absence (0) of a set of known cancer driver mutations across a collection of cell line models.

## Output File

A .tsv file 
