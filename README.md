# microarray analysis

Implementation of a elementary data analysis tasks on a CEL microarray datasets
(datasets not included)

## loading CEL datafiles

### boxplot of gene expression distributions within samples
![results](https://github.com/peterszabo77/microarray_data_analysis/blob/master/images/results.png)

## linear model fit: gene expressions ~ sidepop (treatment group)

### write differential expressions table into 'modelfit_output_de.txt'

### save heatmap of most differencially expressed genes with clustering



## filtering: remove genes with low expression

filtering: remove genes with low expression
The script implements a randomization statistical test. Consider a fully connected graph with directed edges.
We can consider the connections among three nodes (A,B,C) as hierarchical, if an A->B and a B->C
connection implies an A->C connection (if A is dominant over B, and B is dominant over C, then
A is dominant over C). Using this concept we can calculate a hierarchy index (H) for any directed graph.
The script calculates the H index of an input graph (given by a connection matrix) and returns a p value, which is an estimate of the probability that the H index of a completely random graph is at least as large
as the H index of the input (one-sided significance of hierarchy).

## Software / libraries
- Python: numpy, matplotlib, graphviz

## Results (examples for graphs with different H-indices)

Different shades of blue indicate the overall dominance of nodes (proportion of outgoing edges).

![results](https://github.com/peterszabo77/graph-property-statistics/blob/master/images/results.png)
