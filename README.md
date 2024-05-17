# Heterogeneity-of-differentiating-hPSC-derived-corneal-limbal-stem-cells-through-scRNAseq

code files related to the manuscript: Deciphering the heterogeneity of differentiating hPSC-derived corneal limbal stem cells through single-cell RNA-sequencing.

## For regenerating all the main figures from the manuscript, follow the subsequent steps:

1. Download and pre-process all the scRNAseq and bulk RNAseq data using seq2science (0.9.2). All used sample and config files are present in the respective 'seq2science' folders. Change the path to the fastq directory and other settings where needed for your specific setup. For more information regarding running seq2science, see: https://github.com/vanheeringen-lab/seq2science

2. scRNAseq analysis (conda env bulk_rna_seq)
Follow the steps written in "Rmarkdown/scRNAseq_analysis" to perform the scRNAseq seq analysis. 

3. bulk RNAseq analysis
4. Follow the steps written in "Rmarkdown/BulkRNAseq_analysis" to perform the bulk RNAseq analysis. 

## For regenerating the supplemental figures for Pseudotime, PROGENy and SCENIC analysis from the manuscript, follow the subsequent steps:
1. Create the conda environments from the respective folders in this Git repository
```
$ conda env create -f {conda_environment_file}.yml
````

2. Run the rscript for pseudotime
```
$ Rscript ipsc_pseudotime.R
```

3. Run the ipython notebooks within activated environments for PROGENy and SCENIC analyses.
Follow the steps written in the notebooks of "PySCENIC_PROGENy".

## Individual objects can also be downloaded from Zenodo
Constructed single-cell objects can be found here: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11208468.svg)](https://doi.org/10.5281/zenodo.11208468)
