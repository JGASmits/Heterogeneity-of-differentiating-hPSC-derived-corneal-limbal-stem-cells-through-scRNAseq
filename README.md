# Heterogeneity-of-differentiating-hPSC-derived-corneal-limbal-stem-cells-through-scRNAseq

code files related to the manuscript: Deciphering the heterogeneity of differentiating hPSC-derived corneal limbal stem cells through single-cell RNA-sequencing

For regenerating all the figures from the manuscript, follow the subsequent steps:

1. Download and pre-process all the scRNAseq and bulk RNAseq data using seq2science (0.9.2). All used sample and config files are present in the respective 'seq2science' folders. Change the path to the fastq directory and other settings where needed for your specific setup. For more information regarding running seq2science, see: https://github.com/vanheeringen-lab/seq2science

2. scRNAseq analysis (conda env bulk_rna_seq)
Follow the steps written in "Rmarkdown/scRNAseq_analysis" to perform the scRNAseq seq analysis. 

3. bulk RNAseq analysis
4. Follow the steps written in "Rmarkdown/BulkRNAseq_analysis" to perform the bulk RNAseq analysis. 
