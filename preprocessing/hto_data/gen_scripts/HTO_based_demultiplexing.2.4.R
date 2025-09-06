############    -------   HTO-based Deconvolution    ------    ############
cat('############    -------   HTO-based Deconvolution    ------    ############\n')

# Still need to add the doublet enrichment.


# Original author: Vicente Fajardo
# Version: 2
# Subversion: 4
# Updates.
# ---> Version updates:
#   Some manual inspection of results gotten with previous versions of this program have shown that the classification usually works great for positive cells, while negative cells are usually not well classified since they seem to convey info to appropriately being classified into one positive class. Therefore, this version adds one important step: the declassification of negative cells into one class as long as they accomplish this next criteria:
#     The ratio between the best hashtag count versus the second best hashtag count is > n (check parameter definition).
# ---> Subversion updates:
#   .1 - Filtering based on GenEx data called cells.
#   .2 - Set a threshold previous to classification when considering a raw file.
#   .3 - Threshold in subversion .2 is now applied considering the UMI value of the top hashtag only.
#      - Get to decide if the original MULTI-seq classification should be used or the one used by the laternative method.
#      - Doublet enrichment: When the original method classes are to be used, this allows for the FC threshold to be applied over singlets, then calling them doublets if the value smaller than the one seen.
#      - We've realized MULTI-seq method is prone to some obvious errors. Therefore, whenever the method assigns to a cell a class whose UMI count is smaller than a different class' one, this program corrects for that and turns it into a Negative.

### -------------------------- Description -------------------------- ###
# Hashtag data (HTO-data)-based deconvolution of samples pooled alltogether. We take advantage of the MULTI-seq algorithm wrapped in the allmighty Seurat package.

cat('\n\n')
### -------------------------- Dependencies ------------------------- ###
cat('### -------------------------- Dependencies ------------------------- ###\n')
library(optparse)
library(Seurat)
library(ggplot2)
library(Rtsne)
library(parallel)
library(stringr)
source('/home/vfajardo/scripts/functions/R_handy_functions.R')
cat('Dependencies loaded...\n')

cat('\n\n')
cat('### --------------------------- Arguments --------------------------- ###\n')
### --------------------------- Arguments --------------------------- ###
# Declaring arguments to parse from command line ----------------------->
option.list <- list(
  make_option(opt_str="--ReportsPath", type="character", default=NULL, dest="reports.path", help="Absolute path to directory to save DE genes table."),
  make_option(opt_str="--HTOData", type="character", default=NULL, dest="hto.data.file", help="Absolute path to HTO data file. Provide raw matrix if it is suppoused to be filtered according to gene expression library called cells."),
  make_option(opt_str="--PrjName", type="character", default='Prj', dest="prj.name", help="Project name."),
  make_option(opt_str="--HTONames", type="character", default=NULL, dest="hto.names", help="Character vector-style. When HTO and CITE-seq data were gotten together, antobodies not specified in this list will be filtered out."),
  make_option(opt_str="--FCThold", type="numeric", default=3, dest="fc.thold", help="UMI count Fold Change Threshold to be used for Negative desclassification. I.e., given a cell initially classified as 'Negative', its hashtags counts will be evaluated for the first two with the highest counts, then getting their ratio for a fold change. If this value is larger than the treshold (this parameter), then the cell will be classified as the hashtag with the highest UMI count."),
  make_option(opt_str="--GenExData", type="character", default=NULL, dest="genex.data.file", help="Character. If it desired to filter a HTO raw matrix based on the barcodes called as cells by the gene expression library, then provide the absolute path to the cellranger directory output by count. If no filtering should be done, leave as NULL."),
  make_option(opt_str="--UMIThold", type="integer", default=NULL, dest="umi.thold", help="Integer. Works whenever raw data must be filtered (assummed whenever a gene expression data file is provided). Then, as part of a prefilteting step, threshold for the number of UMIs for a barcode to be considered as cell. Leave as NULL if no threshold should be applied for filtering."),
  make_option(opt_str="--DoReclass", type="character", default=NULL, dest="do.reclass", help="Character, value indicating what kind of reclassification should be done if any. Valid values, none or simply NULL for no reclassification, 'negative' for negative declassification (negative values get a class according to an alternative method) or 'doublet' for doublet enrichment method (singlets classified as 'doublets' whenever they aren't over the ratio threshold value).")
)
# Getting arguments from command line and setting their values to their respective variable names.
opt.parser = OptionParser(option_list=option.list);
opt = parse_args(opt.parser);
# Moving options to their own variables
reports.path <- opt$reports.path
hto.data.file <- opt$hto.data.file
prj.name <- opt$prj.name
hto.names <- opt$hto.names
fc.thold <- opt$fc.thold
genex.data.file <- opt$genex.data.file
umi.thold <- opt$umi.thold
if(is.null(umi.thold)) umi.thold <- 100
do.reclass <- opt$do.reclass
# Reassign appropriately if no value was provided for thsi flag.
if(is.null(do.reclass)) do.reclass <- 'none'
# Just if an expression needs to be evaluated to be considered as an R object, use...
if(!is.null(hto.names)) hto.names <- eval(expr=parse(text=hto.names))
# Then, show arguments.
cat('Reports path:', reports.path, '\n')
cat('Absolute path to HTO data:', hto.data.file, '\n')
cat('Project name:', prj.name, '\n')
cat('HTO names:', paste0(hto.names, collapse=','), '\n')
cat('UMI threshold:', umi.thold, '\n')

### --------------------------- Functions --------------------------- ###
# 1 -------------------------------------------------------------------->
# Name: Get two top HTOs per cell.
# Description:
# Given a cell idx within the dimension of the matrix, this function will attempt to return the amount of UMIs for the top two UMI count hashtags, as well as their corresponding classes.
# Arguments ------------------------->
# this.cell - cell index within the hto data matrix dimensions (specifically, columns).
# Function:

get.top.htos <- function(this.cell){
  this.cell <- hto.data[, this.cell]
  this.cell <- sort(x=this.cell, decreasing=TRUE)
  to.output <- data.frame(top.umis=this.cell[1], second.umis=this.cell[2], top.class=names(this.cell)[1], second.class=names(this.cell)[2], stringsAsFactors=FALSE)
  return(to.output)
}

cat('\n\n')
cat('### --------------------------- Load data --------------------------- ###\n')
### --------------------------- Load data --------------------------- ###
hto.data <- read.10X.data(file=hto.data.file)
cat('HTO data has been loaded.\n')

cat('\n\n')
cat('### ---------------------- Data preprocessing ----------------------- ###\n')
### ---------------------- Data preprocessing ----------------------- ###

preprocess.path <- paste0(reports.path, '/preprocessing_reports')
create.dir(dir.path=preprocess.path, path.desc='Preprocessing')

# ---> According to HTO abs.
# If necessary, remove non-HTO CITE-seq data.
if(!is.null(hto.names)){
  # Make sure, we've got data for at least two hashtags.
  if(!all(hto.names %in% rownames(hto.data))) cat('=============================== Warning ===============================\n\nNot all hashtag names input are part of the HTO data. Be cautios.\n\n=======================================================================\n')
  hto.names <- hto.names[hto.names %in% rownames(hto.data)]
  if(length(hto.names)<2) stop('================================ Error ================================\n\nAt least two hashtags must be provided.\n\n=======================================================================\n')
  # Then, clean data.
  hto.data <- hto.data[hto.names, ]
}else{
  hto.names <- rownames(hto.data)
}

# ---> According to GenEx data
if(!is.null(genex.data.file)){
  cat('\n\n# ---> Filtering according to GenEx data.\n')
  # ---> Load GenEx data.
  if(!dir.exists(genex.data.file)) stop(paste0('It was required for the HTO data to be filtered according to GenEx cells, but no appropriate file was provided, instead: ', genex.data.file))
  barcodes.file <- list.files(path=genex.data.file, pattern='barcodes.tsv', full.names=TRUE)
  if(length(barcodes.file)!=1) stop(paste0('Even though a valid directory (', genex.data.file, ') for GenEx data was provided, no barcodes.tsv<.gz> file was found.'))
  if(grepl(x=barcodes.file, pattern='.gz$')){ tmp.cmd <- paste0('zcat ', barcodes.file); genex.cells <- system(command=tmp.cmd, intern=TRUE) } else genex.cells <- read.csv(file=barcodes.file, header=FALSE, stringsAsFactors=FALSE)[, 1]
  # Standardize everything (no suffix).
  genex.cells <- str_replace(string=genex.cells, pattern='-\\d+$', replacement='')
  colnames(hto.data) <- str_replace(string=colnames(hto.data), pattern='-\\d+$', replacement='')
  # ---> Reports for GenEx data.
  total.cells <- length(genex.cells)
  no.cells.in.hto <- sum(genex.cells %in% colnames(hto.data))
  prop.cells.in.hto <- no.cells.in.hto/length(genex.cells)
  tmp.output <- data.frame(total.cells=total.cells, cells.in.hto=no.cells.in.hto, prop.cells.in.hto=prop.cells.in.hto)
  # ---> Output.
  tmp.file.name <- paste0(preprocess.path, '/GenExBasedFiltering.csv')
  write.csv(file=tmp.file.name, x=tmp.output, row.names=FALSE, quote=FALSE)
  # ---> Filter according to GenEx data.
  hto.data <- hto.data[, colnames(hto.data) %in% genex.cells]
  cat('Filtering applied!\n')
}

# ---> According to UMI counts.
if(!is.null(umi.thold)){
  cat('\n\n# ---> Filtering according to UMI counts.\n')
  cat(paste0('UMI count threshold to use: ', umi.thold, '\n'))
  # For final filtering reports.
  total.cells <- ncol(hto.data)
  # ---> Filter poor cells.
  # To speed next calculations up, using a faster operation, we get rid of cells we know apriori are a total waste because won't pass the threshold for the top class.
  hto.data <- hto.data[, Matrix::colSums(hto.data)>umi.thold]
  # ---> Top classes and UMI counts.
  # Calculate top UMIs and classes values across cells.
  top.htos.per.cell <- mclapply(X=1:ncol(hto.data), FUN=get.top.htos)
  top.htos.per.cell <- Reduce(x=top.htos.per.cell, f=rbind)
  rownames(top.htos.per.cell) <- colnames(hto.data)
  # ---> Filtering reports.
  no.over.thold <- sum(top.htos.per.cell$top.umis>umi.thold)
  prop.over.thold <- no.over.thold/total.cells
  tmp.output <- data.frame(total.cells=total.cells, no.cells.over.thold=no.over.thold, prop.cells.over.thold=prop.over.thold)
  tmp.file.name <- paste0(preprocess.path, '/UMIBasedTholdFiltering.csv')
  write.csv(file=tmp.file.name, x=tmp.output, row.names=FALSE, quote=FALSE)
  # ---> Filter
  hto.data <- hto.data[, top.htos.per.cell$top.umis>umi.thold]
  top.htos.per.cell <- top.htos.per.cell[top.htos.per.cell$top.umis>umi.thold, ]
  cat('Filtering applied!\n')
}
cat('\n\nHTO table has been filtered if it was required.\n')

cat('\n\n')
cat('### ------------------------- Main program -------------------------- ###\n')
### ------------------------- Main program -------------------------- ###

cat('\n\n')
cat('### ----------------------- HTO seurat object ----------------------- ###\n')
### ----------------------- HTO seurat object ----------------------- ###
seurat.obj <- CreateSeuratObject(counts=hto.data, project=prj.name, assay='HTO', min.cells=3)
# Fix HTO names if necessary
hto.names <- str_replace_all(string=hto.names, pattern='_', replacement='-')
seurat.obj <- NormalizeData(object=seurat.obj, assay="HTO", normalization.method="CLR")
cat('HTO-assay created as a seurat object and processed accordingly.\n')
# We simply get rid of the original data file since we won't use it anymore.
rm(hto.data)

cat('\n\n')
cat('### ------------------------- Deconvolution ------------------------- ###\n')
### ------------------------- Deconvolution ------------------------- ###
# MULTI-seq deconvolution.
seurat.obj <- MULTIseqDemux(object=seurat.obj, assay="HTO", autoThresh=TRUE, maxiter=10)
cat('Deconvolution applied based on MULTIseq.\n')

cat('\n\n')
cat('### ------------------------ Reclassification ----------------------- ###\n')
### ------------------------ Reclassification ----------------------- ###
# ---> Raw data.
# Get UMI counts and classification to work with.
hto.data.and.class <- cbind(as.data.frame(t(as.matrix(seurat.obj@assays$HTO@counts))), seurat.obj$MULTI_ID)
colnames(hto.data.and.class) <- c(colnames(hto.data.and.class)[1:(ncol(hto.data.and.class)-1)], 'MULTI-seq_classification')
hto.data.and.class[, 'MULTI-seq_classification'] <- as.character(hto.data.and.class[, 'MULTI-seq_classification'])

# ---> Alternative method for classification.
# Get ratios and best class for 'Negative' declassification process..
alt.class.info <- mclapply(X=rownames(hto.data.and.class), FUN=function(tmp.barcode){
  barcode.counts <- hto.data.and.class[tmp.barcode, ]
  # For general purposes.
  hto.counts <- as.integer(barcode.counts[, 1:(length(barcode.counts[1, ])-1)])
  names(hto.counts) <- colnames(barcode.counts[, 1:(length(barcode.counts[1, ])-1)])
  hto.counts <- sort(x=hto.counts, decreasing=TRUE)
  normal.class.ratio <- hto.counts[1]/hto.counts[2]
  alt.class <- names(hto.counts)[1]
  # Taking into account the class identified by MULTI-seq.
  barcode.class <- barcode.counts[1, 'MULTI-seq_classification']
  class.counts <- hto.counts[names(hto.counts)==barcode.class]
  hto.counts <- c(class.counts, hto.counts[names(hto.counts)!=barcode.class])
  alt.class.ratio <- hto.counts[1]/hto.counts[2]
  # Format and output.
  to.return <- data.frame(normal.class.ratio=normal.class.ratio, alt.class=alt.class, alt.class.ratio=alt.class.ratio, stringsAsFactors=FALSE)
  return(to.return)
})
alt.class.info <- Reduce(x=alt.class.info, f=rbind)
# Add doublet or infered positive class according to FC threshold.
alt.class.info[alt.class.info$normal.class.ratio<fc.thold, 'alt.class'] <- 'Doublet'
# Add info.
hto.data.and.class <- cbind(hto.data.and.class, alt.class.info)

# Check concordance between normal class ratio and the alternative class ratio.
are.all.same <- all(hto.data.and.class$normal.class.ratio == hto.data.and.class$alt.class.ratio)
count.same <- sum(hto.data.and.class$normal.class.ratio == hto.data.and.class$alt.class.ratio)
prop.same <- count.same / nrow(hto.data.and.class)
tmp.output <- data.frame(all.concordant=are.all.same, concordance.freq=count.same, concordance.prop=prop.same, stringsAsFactors=FALSE)
tmp.file.name <- paste0(preprocess.path, '/HashtagsRatios.csv')
write.csv(file=tmp.file.name, x=tmp.output, row.names=FALSE, quote=FALSE)

# ---> Reclassification.
switch(EXPR=do.reclass,
  'none'=cat('No reclassification performed.\n'),
  'negative'={
    # ---> Negative declassification.
    cat('# ---> Negative declassification.\n')
    tmp.alt.class <- hto.data.and.class[Cells(seurat.obj)[seurat.obj@meta.data[, 'MULTI_ID']=='Negative'], 'alt.class']
    seurat.obj@meta.data[seurat.obj@meta.data[, 'MULTI_ID']=='Negative', 'MULTI_ID'] <- tmp.alt.class
  },
  'doublet'={
    # ---> Doublet enrichment.
    cat('# ---> Doublet enrichment.\n')
    seurat.obj@meta.data[hto.data.and.class$alt.class.ratio<fc.thold & hto.data.and.class[, 'MULTI-seq_classification']!='Negative', 'MULTI_ID'] <- 'Doublet'
  },
  stop('No appropriate reclassification value provided.\n')
)

# ---> MULTI-seq error correction.
# singlet.idxs <- seurat.obj@meta.data[, 'MULTI_ID']!='Negative' & seurat.obj@meta.data[, 'MULTI_ID']!='Doublet'
# method.errs <- seurat.obj@meta.data[, 'MULTI_ID'] != top.htos.per.cell$top.class & singlet.idxs
# seurat.obj@meta.data[method.errs, 'MULTI_ID'] <- 'Negative'

cat('\n\n')
cat('### ---------------------------- Raw data --------------------------- ###\n')
### ---------------------------- Raw data --------------------------- ###
raw.data.path <- paste0(reports.path, '/raw_data')
if(!dir.exists(raw.data.path)) dir.create(raw.data.path)

# ---> Raw data.
tmp.file.name <- paste0(raw.data.path, '/RawHTOCountsWithMULTISeqClassification.csv')
write.csv(x=hto.data.and.class, file=tmp.file.name)
cat('Raw data output to a csv file.\n')

cat('\n\n')
cat('### ------------------------ General reports ------------------------ ###\n')
### ------------------------ General reports ------------------------ ###
gen.reports.path <- paste0(reports.path, '/general_reports')
dir.create(gen.reports.path)

# ---> Classification summary
cat('# Classification summary\n')
tmp.table <- as.data.frame(table(seurat.obj@meta.data[, 'MULTI_ID']))
colnames(tmp.table) <- c('Classification', 'Frequency')
tmp.file.name <- paste0(gen.reports.path, '/GeneralClassificationSummary.csv')
write.csv(x=tmp.table, file=tmp.file.name)

# ---> Classification per cell.
cat('# Classification per cell.\n')
tmp.file.name <- paste0(gen.reports.path, '/ClassificationPerCell.csv')
write.csv(x=seurat.obj[['MULTI_ID']], file=tmp.file.name, quote=FALSE)

# ---> Classification methods agreement.
tmp.output <- sum(hto.data.and.class[, 'MULTI-seq_classification']==hto.data.and.class$alt.class)/nrow(hto.data.and.class)
tmp.file.name <- paste0(gen.reports.path, '/MethodsAgreement.csv')
write.csv(x=tmp.output, file=tmp.file.name, quote=FALSE)

cat('\n\n')
cat('### --------------------- HTO counts across IDs --------------------- ###')
### --------------------- HTO counts across IDs --------------------- ###
counts.per.ids.path <- paste0(reports.path, '/HTOs_across_IDs')
if(!dir.exists(counts.per.ids.path)) dir.create(counts.per.ids.path)

for(tmp.hto in hto.names){
  # Directory for each tag.
  tmp.counts.per.ids.path <- paste0(counts.per.ids.path, '/', tmp.hto)
  dir.create(tmp.counts.per.ids.path)
  # Density plot across IDs considering other tags.
  tmp.df <- data.frame(counts=seurat.obj@assays$HTO@data[tmp.hto, ], ID=seurat.obj@meta.data[, 'MULTI_ID'])
  tmp.ggplot <- ggplot(data=tmp.df, aes(x=counts, fill=ID)) + geom_density(alpha=0.6)
  tmp.file.name <- paste0(tmp.counts.per.ids.path, '/DensityPlotForHTO', tmp.hto, 'CountsAcrossIDsAndTags.pdf')
  pdf(file=tmp.file.name)
  print(tmp.ggplot)
  dev.off()
  # Density plot across IDs not considering other tags.
  tmp.df <- data.frame(counts=seurat.obj@assays$HTO@data[tmp.hto, ], ID=ifelse(test=seurat.obj@meta.data[, 'MULTI_ID']==tmp.hto, yes=tmp.hto, no=paste0('non-', tmp.hto)))
  tmp.ggplot <- ggplot(data=tmp.df, aes(x=counts, fill=ID)) + geom_density(alpha=0.6)
  tmp.file.name <- paste0(tmp.counts.per.ids.path, '/DensityPlotForHTO', tmp.hto, 'CountsAcrossIDs.pdf')
  pdf(file=tmp.file.name)
  print(tmp.ggplot)
  dev.off()
}
cat('Plots output.\n')

cat('\n\n')
cat('### ------------------ Doublet and negatives rate ------------------- ###\n')
### ------------------ Doublet and negatives rate ------------------- ###
tags.path <- paste0(reports.path, '/negative_and_doublet_tags')
dir.create(tags.path)

# Is doublet?
seurat.obj$is.doublet <- ifelse(test=seurat.obj$MULTI_ID=='Doublet', yes='doublet', no='s/n')
# Is negative?
seurat.obj$is.negative <- ifelse(test=seurat.obj$MULTI_ID=='Negative', yes='negative', no='classified')

# * Negative tag analysis.
# Histogram.
tmp.ggplot <- ggplot(data=seurat.obj@meta.data, aes(x=nCount_HTO, fill=is.negative)) + geom_histogram(binwidth=100, position='dodge') + xlab('HTO counts')
tmp.file.name <- paste0(tags.path, '/NegativeTagHTOCountsHistogram.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot)
dev.off()
# Density plot.
tmp.ggplot <- ggplot(data=seurat.obj@meta.data, aes(x=nCount_HTO, fill=is.negative)) + geom_density(alpha=0.8) + xlab('HTO counts')
tmp.file.name <- paste0(tags.path, '/NegativeTagHTOCountsDensity.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot)
dev.off()

# * Doublet tag analysis.
# Histogram.
tmp.ggplot <- ggplot(data=seurat.obj@meta.data, aes(x=nCount_HTO, fill=is.doublet)) + geom_histogram(binwidth=100, position='dodge') + xlab('HTO counts')
tmp.file.name <- paste0(tags.path, '/DoubletTagHTOCountsHistogram.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot)
dev.off()
# Density plot.
tmp.ggplot <- ggplot(data=seurat.obj@meta.data, aes(x=nCount_HTO, fill=is.doublet)) + geom_density(alpha=0.8) + xlab('HTO counts')
tmp.file.name <- paste0(tags.path, '/DoubletTagHTOCountsDensity.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot)
dev.off()

cat('\n\n')
cat('### ------------------- Dimensionality reduction -------------------- ###\n')
### ------------------- Dimensionality reduction -------------------- ###
dim.reduction.path <- paste0(reports.path, '/dim_reduction')
if(!dir.exists(dim.reduction.path)) dir.create(dim.reduction.path)

# We'll go through this part only if we have HTO data for more than 2 Hashtags.
if(nrow(seurat.obj@assays$HTO@data) > 2){
  # Scaling RNA data, we only scale the variable features here for efficiency
  seurat.obj <- ScaleData(seurat.obj, assay='HTO')

  # ---> PCA
  seurat.obj <- RunPCA(object=seurat.obj, assay='HTO', features=rownames(seurat.obj))
  tmp.file.name <- paste0(dim.reduction.path, '/PCAByHashtagClassification.pdf')
  pdf(file=tmp.file.name)
  print(PCAPlot(object=seurat.obj))
  dev.off()
  # ---> tSNE
  seurat.obj <- RunTSNE(object=seurat.obj, assay='HTO', features=rownames(seurat.obj), check_duplicates = FALSE)
  tmp.file.name <- paste0(dim.reduction.path, '/TSNEByHashtagClassification.pdf')
  pdf(file=tmp.file.name)
  print(TSNEPlot(object=seurat.obj))
  dev.off()
  # ---> UMAP
  seurat.obj <- RunUMAP(object=seurat.obj, assay='HTO', features=rownames(seurat.obj))
  tmp.file.name <- paste0(dim.reduction.path, '/UMAPByHashtagClassification.pdf')
  pdf(file=tmp.file.name)
  print(UMAPPlot(object=seurat.obj))
  dev.off()
  cat('Dimensionality reduction-related plots have been output.\n')
}

cat('\n\n')
cat('### -------------------------- Seurat object ------------------------ ###\n')
### -------------------------- Seurat object ------------------------ ###

seurat.objs.path <- paste0(reports.path, '/seurat_objs')
dir.create(seurat.objs.path)
# ---> Seurat objects.
tmp.file.name <- paste0(seurat.objs.path, '/HTOAssaySeuratObj.RDS')
saveRDS(object=seurat.obj, file=tmp.file.name)
cat('Seurat object for the HTO assay has been output.\n')

# ---> Sessio info.
sessionInfo()

cat('\n\nProgram finished!\n')
