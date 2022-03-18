# Steps to generate the enhancer dataset

## 0. Download resources

- Fantom5 enhancers
- Pegasus enhancers
- FOCS enhancers (GRO-seq, Fantom5, Roadmap)
- Genehancer enhancers
- GENCODEV29lifthg19
- HGNC gene table

## 1. Preprocessing

Please check:

`enhancer_datasets_merge.ipynb`

for the detail of the preprocessing of the tables.


## 2. Post-processing before merging

`DATASETS_CREATION_FINAL.cleaned.py`

## 3. Merge, split, reannotate

First, merge all the "regulatory elements" resources into a single bed file
(suggested: split per chromosome), and use bedtools to sort positions.

You can then use the script `FINSURF-model-generation/code_train/FINSURF_train/utils/bedops_split_reannotate.sh`

Which will allow the segmentation of your bedfile into tiles that are
non-overlapping; that is if two elements from two resources are partly overlaping, the resulting bed file
will contain 3 regions: 2 regions specific to each resource, and a third
corresponding to the common element.

The reannotation allows to annotate back which are the suppporting resources
for a given region.

