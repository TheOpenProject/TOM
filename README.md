# TOM
Translational Oncology Modelling

The R script "TOP TOM 2019-02-07.R" focuses on the following two questions:
Goal 1. Which rate laws accurately recapitulate PDX tumour growth?
Goal 2. Which rate laws accurately predict PDX tumour growth?

To run the R script, you would need to obtain the Novartis mouse clinical trial data from

https://www.nature.com/articles/nm.3954#supplementary-information

Excel files
Supplementary Table 1
Genomic profiling of PDXs and raw response and curve metrics of PCTs.

Save the file locally and set the work directory to find it:
setwd('/appropriate directory to find the following file/')
dat = read_excel("Genomic profiling of PDXs and raw response and curve metrics of PCTs.xlsx",sheet="PCT raw data")

