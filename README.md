# microarray analysis - differential expression and clustering 

microarray_a1.R: Implementation of a elementary data analysis tasks on a CEL microarray datasets

(datasets not attached)

1. loading CEL datafiles

1. boxplot of gene expression distributions within samples

![results](https://github.com/peterszabo77/microarray_data_analysis_1/blob/main/output/data_boxplot.png)

1. filtering: remove genes with low expression

1. linear model fit: gene expressions ~ sidepop (treatment group)

1. write differential expressions table into 'modelfit_output_de.txt'

1. save heatmap of most differencially expressed genes with clustering

![results](https://github.com/peterszabo77/microarray_data_analysis_1/blob/main/output/heatmap_diffexp.png)

### Software / libraries

- R: Bioconductor

### References

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25191

Introduction to Bioinformatics with R: A Practical Guide for Biologists by Edward R. Curry



