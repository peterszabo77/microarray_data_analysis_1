library(affy) # to load CEL format
library(limma)
library(annotate) # to obtain gene names
library(hgu133plus2.db) # to obtain gene names
library(gplots)

outputdir = 'output'

########################################################################
# read CELL dataset
# output: data_boxplot.png
#
print('loading dataset')

normdata = justRMA() # read nomalized data
array_expr = exprs(normdata) # array of expression values, shape(N_genes, N_samples)

sample_names = sampleNames(normdata)
probeset_names = featureNames(normdata)

print('expession array shape')
print(dim(array_expr))
print('sample names:')
print(sample_names)

png(file.path(outputdir, "data_boxplot.png"))
boxplot(array_expr, main="gene expression distributions within samples")
dev.off()

#########################################################################
# filtering: remove genes with low expression
#
print('filtering: remove genes with low expression')
medianexprs = apply(array_expr, MARGIN=1, median)
lowexprs = which(medianexprs<5)
array_expr = array_expr[-lowexprs,]
print(paste(length(lowexprs), 'genes removed'))
print(paste(dim(array_expr)[1], 'genes remained'))

########################################################################
# linear model fit: gene expressions ~ sidepop
# saved: modelfit_output.txt
#
print('performing linear model fitting: expressions ~ sidepop membership')

design = cbind(c(1,1,1,1,1,1), sidepop=c(0,0,0,1,1,1)) # based on side pop membership
colnames(design) = c('Intercept', 'Sidepop')
modelfit = lmFit(array_expr, design)
modelfit = eBayes(modelfit)
#modelfit_output = topTable(modelfit, coef=c(2, 3), adjust="fdr", number=nrow(array_expr), p.value=0.05)
modelfit_output = topTable(modelfit, coef=2, adjust="fdr")
modelfit_output$ID = rownames(modelfit_output)
modelfit_output$Symbol = unlist(mget(rownames(modelfit_output),hgu133plus2SYMBOL))

write.table(modelfit_output, file=file.path(outputdir, "modelfit_output_de.txt"), sep="\t", quote=FALSE, row.names=FALSE)

########################################################################
# save heatmap of most differencially expressed genes with clustering
# saved: heatmap.pdf
#
print('save heatmap of most differencially expressed genes with clustering')

diffexp_gene_IDs = as.character(modelfit_output$ID)
de_row_idxs = which(rownames(array_expr) %in% diffexp_gene_IDs)
gene_names = unlist(mget(diffexp_gene_IDs,hgu133plus2SYMBOL))

png(file.path(outputdir, "heatmap_diffexp.png"))
heatmap(array_expr[de_row_idxs,], scale="row", labCol=c("nSP1","nSP2","nSP3","SP1","SP2","SP3"), labRow=gene_names,
main="heatmap of DE genes")
dev.off()
