############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    -----------   Tables stack    -----------    ############
###########    ------------ BOOST project 1  -------------    ###########
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############

# ---> About the script.
# Version: 0
# Subversion: 1
# Script to get all bioinformatics-based tables for the BOOST project.


############    -----------------------------------------    ############
### -------------------------- Description -------------------------- ###
############    -----------------------------------------    ############

# Script to get all bioinformatics-based tabkles for the paper from the BOOST 1 project.


############    -----------------------------------------    ############
### --------------------------- Libraries --------------------------- ###
############    -----------------------------------------    ############
library(Seurat)
library(data.table)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(fgsea)
library(immunarch)
library(Hmisc)
library(corrplot)
library(jsonlite)


############    -----------------------------------------    ############
### --------------------------- Functions --------------------------- ###
############    -----------------------------------------    ############

coll.version <- 0.1
gen.data.path <- '/mnt/BioAdHoc/Groups/vd-sette/BOOST/paper_developments/BOOST_1/paper_items'
coll.file <- paste0(gen.data.path, '/jobs_scripts/functions_collection.', coll.version, '.R')
if(file.exists(coll.file)) source(coll.file) else stop('Faulty function collection file.\n')
total.workers <- try(expr=as.integer(system(command='echo $PBS_NUM_PPN', intern=TRUE)))
if(is.na(total.workers)) total.workers <- try(expr=as.integer(system(command='echo $SLURM_JOB_CPUS_PER_NODE', intern=TRUE)))
if(!is.numeric(total.workers)) stop('Problem while defining number of processors/CPUs available per node.\n')


############    -----------------------------------------    ############
### ----------------------- General Arguments ----------------------- ###
############    -----------------------------------------    ############

# ---> General definitions.
# @ Paper.
main.prj <- 'BOOST'
this.prj <- 'BOOST_1'
this.table <- 'table_stack_v1'
# @ Seed.
set.seed(seed=1)
# @ Dataset labels.
obj.extended.names <- c(
  b1.cd4='CD4_Full',
  b1.cd8='CD8_Full'
)
gen.cell.types <- c(
    b1.cd4='CD4',
    b1.cd8='CD8'
)
# @ Cluster labels for each dataset.
clust.labs <- c(
    b1.cd4='RNA_snn_res.0.4',
    b1.cd8='RNA_snn_res.0.2'
)
clusts.lab <- clust.labs
# @ For tag analysis.
# min.no.cells.per.group <- 50 # @ Internal parameter for tag specific analysis (part of figure ?).
# colors.table <- NULL
# ---> Cluster identities
clusts.defs <- list(
  `b1.cd4`=c(
    `0`='TREG',
    `1`='TH17',
    `2`='Poly',
    `3`='TFH',
    `4`='TCM',
    `5`='TH2',
    `6`='TFR',
    `7`='Unknown'
  ),
  `b1.cd8`=c(
    `0`='GZMKhi',
    `1`='exTEFF',
    `2`='MAIT',
    `3`='TEFF',
    `4`='LTBhi',
    `5`='TCM',
    `6`='IFNR'
  )
)
# ---> Color defintions.
# @ Color per cell type.
cell.type.cols <- c(
  CD4='#00BFFF',
  CD8='#EE82EE'
)
# @ Color per cluster per dataset.
# See definitions at: https://docs.google.com/spreadsheets/d/1_lNI4M7P-LBuwfQWEgag4pHaycMZMNcmzGv_k_e7x5U/edit#gid=0
clusts.cols <- list(
  `b1.cd4`=c(
    `0`='#319272', # TREG
    `1`='#0000CD', # TH17
    `2`='#D3B2CB', # Poly
    `3`='#A6B0E3', # TFH
    `4`='#00BFFF', # TCM/Naive
    `5`='#FFB4B4', # TH2
    `6`='#32CD32', # TFR
    `7`='#676765' # Unknown
  ),
  `b1.cd8`=c(
    `0`='#4DA6FF', # GZMKhi
    `1`='#CC99FF', # exTEFF
    `2`='#676765', # MAIT
    `3`='#F2E6FF', # TEFF
    `4`='#BFBFBF', # LTBhi
    `5`='#EE82EE', # Naive/TCM
    `6`='#FFD700' # IFNR
  )
)
# @ Dose cohorts' colors
cohort.cols <- c(
    `dose.2`='#0096ED',
    `dose.3`='#894C00',
    `dose.4`='#009000'
)
# @ DEA colors
deg.class.cols <- c(
    `Negative`='#0036FE',
    `Positive`='#F60D00',
    `Other`='#D6D7D7'
)
# @ Other terms' colors.
signatures.col.scale <- c('#ffffff', '#ffffe0', '#ffffad', '#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
# @ Lane tag info.
lane.tags <- c(
  `CD4_Full`='cell.type.tag;facs.sorting.batch.tag',
  `CD8_Full`='cell.type.tag;facs.sorting.batch.tag'
)
# @ Others.
subset.n <- 50000
# ---> Path definitions.
# General data and reports paths.
prj.gen.path <- paste0('/mnt/BioAdHoc/Groups/vd-sette/', main.prj, '/paper_developments/', this.prj)
gen.data.path <- paste0(prj.gen.path, '/paper_items')
final.tables.path <- paste0(prj.gen.path, '/final_tables')
this.table.path <- paste0(final.tables.path, '/', this.table)
if(!dir.exists(this.table.path)) dir.create(this.table.path)
reports.path <- this.table.path
# ---> File definitions.
# Seurat objects.
objs.files.list <- paste0(gen.data.path, '/', 'SeuratObj_', obj.extended.names, '.RDS')
names(objs.files.list) <- names(obj.extended.names)
# TCR data.
tcr.meta.files <- c(
    `b1.cd4`='/mnt/BioAdHoc/Groups/vd-sette/BOOST/paper_developments/BOOST_1/exploratory_analyses/preliminary_process/Round-4_CD4_2023-10-19/CellBasedTCRData.csv',
    `b1.cd4`='/mnt/BioAdHoc/Groups/vd-sette/BOOST/paper_developments/BOOST_1/exploratory_analyses/preliminary_process/Round-4_CD8_2023-10-19/CellBasedTCRData.csv'
)
# Aggr tables.
aggr.files <- paste0(gen.data.path, '/', 'AggrTable_', obj.extended.names, '.csv')
names(aggr.files) <- names(obj.extended.names)
# Gene sets' files (for FGSE analyses).
module.feats.dir <- paste0(gen.data.path, '/module_signatures')
# # Cluster markers.
marker.files.list <- paste0(gen.data.path, '/ClusterMarkers_', obj.extended.names, '.csv')
names(marker.files.list) <- obj.extended.names
# Specific files.
#     Donor-specific cell counts before and after QC filtering.
# qc.inf.1.file <- paste0(gen.data.path, '/QCInf_CellCount_Donor.csv')
#     Annotated metadata for datasets before QC filtering.
# no.qc.meta.file <- paste0(gen.data.path, '/AnnotatedGExAggrMetaData.RDS')
#     Donors' metadata.
# donor.meta.file <- paste0(gen.data.path, '/DonorMetadata.csv')
#     Brief genome reference details.
gen.ref.file <- paste0(gen.data.path, '/DLCP_geneAnnot_hg38.complete.txt')
# ---> Check directories and files.
if(!all(dir.exists(gen.data.path), dir.exists(this.fig.path), dir.exists(module.feats.dir))) stop(paste0('Following paths must be already defined:\n', gen.data.path, '\n', this.fig.path, '\n', module.feats.dir, '\n'))
if(!dir.exists(reports.path)) dir.create(reports.path)
essential.files <- c(
  objs.files.list#,
#   no.qc.meta.file, msigdb.file
)
essential.files <- essential.files[!unlist(lapply(X=essential.files, FUN=file.exists))]
if(length(essential.files) > 0) stop(paste0('Next essential files are not appropriately defined -make sure you\'ve got their names right-:\n', paste0(essential.files, collapse='\n'), '\n'))
# ---> Present arguments (to do)


############    -----------------------------------------    ############
### ------------------------- Data Loading -------------------------- ###
############    -----------------------------------------    ############

# ---> Seurat obejcts.
# Main seurat objects.
srt.objs.list <- lapply(X=objs.files.list, FUN=readRDS)
names(srt.objs.list) <- names(objs.files.list)

# ---> TCR data.
b1.tcr.meta <- lapply(X=tcr.meta.files, FUN=fread)

# # ---> Cluster markers.
# markers.list <- lapply(X=marker.files.list, FUN=fread)
markers.list <- lapply(X=marker.files.list, FUN=function(x){
  if(file.exists(x)) fread(file=x) else NA
})

# ---> Specific files.
#     Metadata of unfiltered seurat objects per cell type.
# unfilt.meta <- readRDS(file=unfilt.meta.file)
#     General donors' metadata and demographics files.
# donor.meta <- fread(file=donor.meta.file)
# donor.demo <- fread(file=donor.demo.file)
#     Brief genome reference details.
gen.ref <- fread(file=gen.ref.file)
#     File listing the IDs for sequencing batches as in original lab recored.
# seq.team.trans <- fread(file=seq.team.trans.file)
#     GSEA results.
# gsea.res <- fread(file=gsea.res.file)

# ---> Gene sets' files (for FGSE analyses).
module.feats.files <- list.files(path=module.feats.dir, recursive=FALSE, include.dirs=FALSE, pattern='signature.csv$', full.names=TRUE)
modules.feats <- lapply(X=module.feats.files, FUN=function(tmp.file){
  tmp.feats <- read.csv(file=tmp.file, stringsAsFactors=FALSE)
  if(!'feature' %in% colnames(tmp.feats)) stop(paste0('File ', tmp.file, ' not appropriately defined. Column listing the gene features should be named \'feature\'.\n\n'))
  tmp.feats <- tmp.feats$feature
  # Replace underscores with points in the gene names just as Seurat does when one creates a seurat object.
  tmp.check <- any(str_detect(string=tmp.feats, pattern='_'))
  if(tmp.check){
    tmp.warning <- paste0('Feature names cannot have underscores (\'_\'), replacing with dashes (\'-\'). This for file: ', basename(tmp.file))
    warning(tmp.warning)
    tmp.feats <- str_replace_all(string=tmp.feats, pattern='_', replacement='-')
  }
  # Return.
  return(tmp.feats)
})
names(modules.feats) <- list.files(path=module.feats.dir, recursive=FALSE, include.dirs=FALSE, pattern='signature.csv$', full.names=FALSE)
names(modules.feats) <- str_replace(string=names(modules.feats), pattern='_signature.csv$', replacement='')
names(modules.feats) <- str_replace_all(string=names(modules.feats), pattern='_', replacement='.')
modules.names <- names(modules.feats)
names(modules.names) <- str_replace_all(string=modules.names, pattern='\\.', replacement='_')

# ---> Define aggr tables. ----- TO REMOVE -----
# Individual aggr table file paths
# aggr.table.path <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/sequencing_data/08-07-2023/aggr/data'
# aggr.data.path <- '/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/sequencing_data/08-07-2023/aggr'
# # Define cell and tissue types.
# aggr.pffx.lab <- 'R24_Cancer_Batches-1-to-20_'
# cell.types <- c('CD4', 'CD8', 'Myeloid', 'NKAndB')
# tissue.types <- c('Normal-lung', 'Tumor-lung')
# aggr.files <- paste0(aggr.table.path, '/', aggr.pffx.lab, cell.types)
# aggr.files <- paste0(rep(x=aggr.files, each=2), '_', tissue.types)
# aggr.files <- paste0(aggr.files, '_aggr_table_annotated.1.0.csv')
# if(!all(file.exists(aggr.files))) stop('Unexpected error.\n')
# names(aggr.files) <- paste0(aggr.pffx.lab, rep(x=cell.types, each=2), '_', tissue.types)


############    -----------------------------------------    ############
### ---------------------- Data preprocessing ----------------------- ###
############    -----------------------------------------    ############

# ---> Ex vivo lung CD4 T cells
tmp.lab <- 'b1.cd4'
# Clusters to remove: 9 & 10 (too small).
clusts.to.rm <- c('8', '9', '10')
tmp.data <- Cells(srt.objs.list[[tmp.lab]])[
    !srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]] %in% clusts.to.rm
]
srt.objs.list[[tmp.lab]] <- subset(x=srt.objs.list[[tmp.lab]], cells=tmp.data)
# Set populations as factors.
srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]] <- factor(
    x=as.character(srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]]),
    levels=names(clusts.defs[[tmp.lab]])
)

# ---> Ex vivo lung CD8 T cells
tmp.lab <- 'b1.cd8'
# Clusters to remove: 7, 8, 9 & 10 (too small).
clusts.to.rm <- c('7', '8', '9', '10')
tmp.data <- Cells(srt.objs.list[[tmp.lab]])[
    !srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]] %in% clusts.to.rm
]
srt.objs.list[[tmp.lab]] <- subset(x=srt.objs.list[[tmp.lab]], cells=tmp.data)
# Remove cells beyond 7.5 on UMAP 1 scale.
tmp.data <- Cells(srt.objs.list[[tmp.lab]])[!srt.objs.list[[tmp.lab]]@reductions$umap@cell.embeddings[, 'UMAP_1'] > 7.5]
srt.objs.list[[tmp.lab]] <- subset(x=srt.objs.list[[tmp.lab]], cells=tmp.data)
# Set populations as factors.
srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]] <- factor(
    x=as.character(srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]]),
    levels=names(clusts.defs[[tmp.lab]])
)

# ---> TCR data.
# Relevant sample info to be extracted directly from seurat object.
b1.tcr.meta <- rbindlist(l=b1.tcr.meta, use.names=TRUE)
b1.tcr.meta <- b1.tcr.meta[, .(
    cell.type.tag, barcode,
    clonotype.tag,
    cdr3b.nt.seq, cdr3a.nt.seq,
    cdr3b.aa.seq, cdr3a.aa.seq,
    trb.v, trb.j, tra.v, tra.j
)]
# intersect(x=colnames(meta.data), y=colnames(b1.tcr.meta))
b1.tcr.meta <- merge(
    x=meta.data, y=b1.tcr.meta,
    by=c('cell.type.tag', 'barcode'),
    all.x=FALSE, all.y=TRUE, sort=FALSE
)

# ---> Specific files.
#     Metadata of unfiltered seurat objects per cell type.
# data.sets <- c(
#   'CD4'='R24_Cancer_Batches-1-to-20_CD4_Normal-lung',
#   'CD8'='R24_Cancer_Batches-1-to-20_CD8_Normal-lung',
#   'Myeloid'='R24_Cancer_Batches-1-to-20_Myeloid_Normal-lung',
#   'NKAndB'='R24_Cancer_Batches-1-to-20_NKAndB_Normal-lung'
# )
# unfilt.meta <- unfilt.meta[data.sets]
# names(unfilt.meta) <- names(data.sets)
# #     General donors' metadata and demographics.
# donor.meta <- as.data.table(separate(data=donor.meta, col='donor.tag', into=c('hashtag.tag', 'chrom.batch.tag'), sep='-', convert=TRUE))
# # donor.demo[, donor.id.tag:=str_replace(string=ID, pattern='DLCP[0]*', replacement='')]
#     Brief genome reference details.
gen.ref[gene_type=='protein_coding', gene_type:='Protein-coding']
# #     File listing the IDs for sequencing batches as in original lab recored.
# seq.team.trans[, seq.team.record.id:=
#   str_replace(string=seq.team.record.id, pattern='NV0', replacement='')
# ]
# tmp.data <- seq.team.trans[, seq.team.record.id]; names(tmp.data) <- seq.team.trans[, seq.batch.id]
# seq.team.trans <- tmp.data


############    -----------------------------------------    ############
### ------------------ Summary stats for features ------------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
# None.

# ---------------------------------------------------------------------->
# Name: Retrieve expression data.
# Description:
# This function will provide a dgcMatrix for a subset of cells -taken from the main environment seurat object- where the subset is taken according to a tag of interest and a specific value of such tag (i.e., the group of cells belonging to such tag value). Plus, data can be output as LogNormalized data or CPMs.
# If no more than 10 cells, the function will return NULL value.
# Arguments ------------------------->
# tag.of.interest - Character, specific tag of interest.
# group - Character, group of interest (must be a valid value for the tag of interest).
# seurat.obj and tmp.slot variable values are taken from the main environment.
# Function:
# get.exp.data <- function(seurat.obj, tag.of.interest, group, data.slot='counts', norm.data='cpms'){
#   # Get expression data and subset keeping only group cells.
#   exp.data <- GetAssayData(object=seurat.obj, slot=data.slot, assay='RNA')
#   group.cells <- Cells(seurat.obj)[seurat.obj@meta.data[, tag.of.interest]==group & !is.na(seurat.obj@meta.data[, tag.of.interest])]
#   if(length(group.cells)<10) return(NULL)
#   exp.data <- exp.data[, group.cells]
#   # Apply normalization if needed.
#   if(norm.data=='cpms'){
#     scale.factor <- 1000000
#     size.factors <- scale.factor / Matrix::colSums(exp.data)
#     exp.data <- Matrix::t(size.factors * Matrix::t(exp.data))
#   }
#   return(exp.data)
# }
# ---------------------------------------------------------------------->

# ---------------------------------------------------------------------->
# Name: Retrieve expression stats summary for a given tag.
# Description:
# Provided a seurat object, this function will provide a table with expression stats summary for a given tag defined in the object's metadata.
# Arguments ------------------------->
# tag.of.interest - Character, specific tag of interest.
# seurat.obj and tmp.slot variable values are taken from the main environment.
# Function:
# get.exp.data.for.tag.vals <- function(seurat.obj, tag.of.interest, data.slot='counts', norm.data='cpms', report.pos.mean=TRUE){
#   # Apply calculations.
#   # @ Tag values as groups.
#   tag.values <- unique(as.character(seurat.obj@meta.data[, tag.of.interest]))
#   tag.values <- tag.values[!is.na(tag.values)]
#   # Remove NA values.
#   tag.values <- gtools::mixedsort(tag.values[!is.na(tag.values)])
#   # @ Retrieve stats per group.
#   stats.per.values <- lapply(X=tag.values, FUN=function(tmp.val){
#     cat(paste0(tmp.val, '\n'))
#     # @ Group cells' expression values.
#     exp.data <- get.exp.data(seurat.obj=seurat.obj, tag.of.interest=tag.of.interest, group=tmp.val, data.slot=data.slot, norm.data=norm.data)
#     if(is.null(exp.data)) return(NULL)
#     # @ Mean
#     feat.means <- Matrix::rowMeans(exp.data)
#     # @ Proportion of expressing cells.
#     feat.props <- Matrix::rowSums(exp.data>0)/ncol(exp.data)
#     # @ Mean of expressing cells.
#     if(report.pos.mean){
#       exp.data[exp.data==0] <- NA
#       feat.pos.means <- Matrix::rowMeans(exp.data, na.rm=TRUE) # NA result means there's no data with value different than NA (where NA=0 in our case). Then, we need to change such values.
#       feat.pos.means[is.na(feat.pos.means)] <- 0
#     }else{
#       feat.pos.means <- NA
#     }
#     # @ Pack data and deliver.
#     feat.stats <- data.table(feature=rownames(seurat.obj), mean=feat.means, mean.pos=feat.pos.means, prop=feat.props)
#     return(feat.stats)
#   })
#   # Get a tag label and name accordingly.
#   tag.lab <- ifelse(test=grepl(x=tag.of.interest, pattern='RNA_snn_res'), yes='resolution', no=str_replace(string=tag.of.interest, pattern='\\.tag', replacement=''))
#   names(stats.per.values) <- paste((tag.lab), tag.values, sep='.')
#   # Remove any NULL values (smallest groups with cell no < 10)
#   stats.per.values <- stats.per.values[!sapply(X=stats.per.values, FUN=is.null)]
#   # Condense values into a single object.
#   stats.per.values <- rbindlist(l=stats.per.values, use.names=TRUE, idcol='tag.value')
#   return(stats.per.values)
# }
# ---------------------------------------------------------------------->

# ---------------------------------------------------------------------->
# Name: Retrieve expression stats summary for a given tag across seurat objects.
# Description:
# This function will provide a table with expression stats summary for a given tag defined in the object's metadata. Process will be done across seurat objects.
# Arguments ------------------------->
# tag.of.interest - Character, specific tag of interest.
# Function:
# get.exp.summ <- function(tag.of.interest, data.slot='counts', norm.data='cpms', report.pos.mean=TRUE){
#   stats.per.values <- lapply(X=names(srt.objs.list), FUN=function(tmp.pop){
#     tmp.data <- get.exp.data.for.tag.vals(seurat.obj=srt.objs.list[[tmp.pop]], tag.of.interest=tag.of.interest, data.slot=data.slot, norm.data=norm.data, report.pos.mean=report.pos.mean)
#     return(tmp.data)
#   })
#   names(stats.per.values) <- names(srt.objs.list)
#   # Standardize format.
#   tmp.data <- lapply(X=names(srt.objs.list), FUN=function(tmp.pop){
#     # Spread data stat by stat.
#     pop.data <- stats.per.values[[tmp.pop]]
#     tmp.pttn <- paste0('resolution.|', str_replace(string=tag.of.interest, pattern='\\.tag$', replacement='.'))
#     pop.data[, tag.value:=str_replace(string=tag.value, pattern=tmp.pttn, replacement='')]
#     # For mean
#     tmp.data.1 <- pop.data[, .(mean=tag.value, feature, stat=mean)]
#     tmp.data.1 <- as.data.table(spread(data=tmp.data.1, key=mean, value=stat, sep='.'))
#     # For mean of positive cells.
#     if(report.pos.mean){
#       tmp.data.2 <- pop.data[, .(pos_mean=tag.value, feature, stat=mean.pos)]
#       tmp.data.2 <- as.data.table(spread(data=tmp.data.2, key=pos_mean, value=stat, sep='.'))
#     }
#     # For precentage of positive cells.
#     tmp.data.3 <- pop.data[, .(prop=tag.value, feature, stat=prop)]
#     tmp.data.3 <- as.data.table(spread(data=tmp.data.3, key=prop, value=stat, sep='.'))
#     # Sanity check or merge data.
#     tmp.check <- tmp.data.1[, .N] == tmp.data.3[, .N]
#     if(report.pos.mean) tmp.check <- tmp.check & (tmp.data.2[, .N] == tmp.data.3[, .N])
#     if(!tmp.check) stop(paste0('Unexpected error for pop.: ', tmp.pop))
#     tmp.data <- if(report.pos.mean) merge(x=tmp.data.1, y=tmp.data.2, by='feature') else tmp.data.1
#     tmp.data <- merge(x=tmp.data, y=tmp.data.3, by='feature')
#     # Sanity check.
#     tmp.lab <- str_replace(string=tmp.pop, pattern='tumor.|hlty.', replacement='')
#     colnames(tmp.data) <- paste0(tmp.lab, '.', colnames(tmp.data))
#     colnames(tmp.data)[colnames(tmp.data)==paste0(tmp.lab, '.feature')] <- 'gene'
#     return(tmp.data)
#   })
#   # Merge calculations across cell types.
#   tmp.data <- Reduce(x=tmp.data, f=function(x, y){
#     merge(x=x, y=y, by='gene', all=TRUE)
#   })
#   return(tmp.data)
# }
# # ---------------------------------------------------------------------->

# # ---> General clusters tag.
# for(tmp.obj in names(srt.objs.list)){
#   srt.objs.list[[tmp.obj]]@meta.data[, 'cluster.tag'] <- srt.objs.list[[tmp.obj]]@meta.data[, clusts.lab[tmp.obj]]
# }

# # ---> Summary for clusters.
# clust.tag.exp.summ <- get.exp.summ(tag.of.interest='cluster.tag', data.slot='data', norm.data='non-cpms', report.pos.mean=FALSE)

# # ---> Summary for sex tag
# sex.tag.exp.summ <- get.exp.summ(tag.of.interest='sex.tag', data.slot='data', norm.data='non-cpms', report.pos.mean=FALSE)

# ---> Summary for sex-cluster tag
# sex.clust.tag.exp.summ <- lapply(X=names(srt.objs.list), FUN=function(tmp.pop){
#   srt.objs.list[[tmp.pop]]@meta.data[, 'sex.cluster.tag'] <- paste(srt.objs.list[[tmp.pop]]@meta.data[, 'sex.tag'], srt.objs.list[[tmp.pop]]@meta.data[, clusts.lab[tmp.pop]], sep='.')
#   tmp.rows <- str_detect(string=srt.objs.list[[tmp.pop]]@meta.data[, 'sex.cluster.tag'], pattern='^NA\\.')
#   srt.objs.list[[tmp.pop]]@meta.data[tmp.rows, 'sex.cluster.tag'] <- NA
#   tmp.data <- get.exp.data.for.tag.vals(seurat.obj=srt.objs.list[[tmp.pop]], tag.of.interest='sex.cluster.tag', data.slot='data', norm.data='non-cpms', report.pos.mean=FALSE)
#   return(tmp.data)
# })
# names(sex.clust.tag.exp.summ) <- str_replace(string=names(srt.objs.list), pattern='Healthy-|Tumor-', replacement='')
# sex.clust.tag.exp.summ <- rbindlist(l=sex.clust.tag.exp.summ, use.names=TRUE, idcol='cell.type')
# sex.clust.tag.exp.summ[, mean.pos:=NULL]
# sex.clust.tag.exp.summ[, cluster:=str_extract(string=tag.value, pattern='\\d+$')]
# sex.clust.tag.exp.summ[, sex:=str_extract(string=tag.value, pattern='Male|Female')]
# sex.clust.tag.exp.summ[, tag.value:=NULL]
# tmp.data.1 <- sex.clust.tag.exp.summ[, .(cell.type, cluster, feature, sex, mean)]
# tmp.data.1 <- spread(data=tmp.data.1, key=sex, value=mean)
# tmp.data.2 <- sex.clust.tag.exp.summ[, .(cell.type, cluster, feature, sex, prop)]
# tmp.data.2 <- spread(data=tmp.data.2, key=sex, value=prop)
# sex.clust.tag.exp.summ <- merge(x=tmp.data.1, y=tmp.data.2, by=c('cell.type', 'cluster', 'feature'), suffixes=c('.mean', '.prop'))
# # Further include cohort stats.
# tmp.data <- lapply(X=names(srt.objs.list), FUN=function(tmp.pop){
#   tmp.data <- get.exp.data.for.tag.vals(seurat.obj=srt.objs.list[[tmp.pop]], tag.of.interest='cluster.tag', data.slot='data', norm.data='non-cpms', report.pos.mean=FALSE)
#   return(tmp.data)
# })
# names(tmp.data) <- str_replace(string=names(srt.objs.list), pattern='Healthy-|Tumor-', replacement='')
# tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='cell.type')
# tmp.data[, mean.pos:=NULL]
# tmp.data[, cluster:=str_extract(string=tag.value, pattern='\\d+$')]
# colnames(tmp.data)[colnames(tmp.data) %in% c('mean', 'prop')] <- paste('cohort', colnames(tmp.data)[colnames(tmp.data) %in% c('mean', 'prop')], sep='.')
# sex.clust.tag.exp.summ <- merge(x=sex.clust.tag.exp.summ, y=tmp.data, by=c('cell.type', 'feature', 'cluster'))
# sex.clust.tag.exp.summ[, gene.id:=feature]; sex.clust.tag.exp.summ[, feature:=NULL]
# sex.clust.tag.exp.summ[, cluster:=as.character(cluster)]

# ---> Summary for sex-cluster tag for the general cohort
# Further include cohort stats.
# for(tmp.obj in names(srt.objs.list)){
#   srt.objs.list[[tmp.obj]]@meta.data[, 'cluster.tag'] <- srt.objs.list[[tmp.obj]]@meta.data[, clusts.lab[tmp.obj]]
# }
# tmp.data <- lapply(X=names(srt.objs.list), FUN=function(tmp.pop){
#   tmp.data <- get.exp.data.for.tag.vals(seurat.obj=srt.objs.list[[tmp.pop]], tag.of.interest='cluster.tag', data.slot='data', norm.data='non-cpms', report.pos.mean=FALSE)
#   return(tmp.data)
# })
# names(tmp.data) <- str_replace(string=names(srt.objs.list), pattern='Healthy-|Tumor-', replacement='')
# tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='cell.type')
# tmp.data[, mean.pos:=NULL]
# tmp.data[, cluster:=str_extract(string=tag.value, pattern='\\d+$')]
# colnames(tmp.data)[colnames(tmp.data) %in% c('mean', 'prop')] <- paste('cohort', colnames(tmp.data)[colnames(tmp.data) %in% c('mean', 'prop')], sep='.')
# sex.clust.tag.exp.summ <- tmp.data
# sex.clust.tag.exp.summ[, tag.value:=NULL]
# sex.clust.tag.exp.summ[, gene.id:=feature]; sex.clust.tag.exp.summ[, feature:=NULL]


############    -----------------------------------------    ############
### ----------------------- On libraries' QCs ----------------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
table.qc.path <- paste0(reports.path, '/table_on_lib_qcs')
if(!dir.exists(table.qc.path)) dir.create(table.qc.path)

# ---> General definitions.
# Seurat objects of interest.
table.qc.objs <- obj.extended.names

# ---> Aggregation QC summaries.
# @ Load summaries.
aggr.qc.files <- paste0(gen.data.path, '/AggrSummary_', table.qc.objs, '.json')
aggr.qcs <- lapply(X=aggr.qc.files, FUN=function(tmp.file) fromJSON(txt=tmp.file))
names(aggr.qcs) <- table.qc.objs
# @ Get relevant data in a tidy format per aggregation of interest.
aggr.qcs <- lapply(X=names(aggr.qcs), FUN=process.aggr.json, lib.name.pttn='^BOOST')
names(aggr.qcs) <- table.qc.objs
# @ Get general aggregation metrics only.
tmp.data <- lapply(X=aggr.qcs, FUN=function(x) return(x[['gen.qc']]))
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='aggr')
tmp.data <- as.data.table(spread(data=tmp.data, key=group, value=value))
# @ Row order.
tmp.order <- c(`CD4_Full`=1, `CD8_Full`=2)
tmp.data[, to.order:=tmp.order[aggr]]
setorderv(x=tmp.data, cols='to.order', order=1)
# @ Column order.
cols.order <- c('aggr', 'pre.total', 'post.total', 'pre.mean', 'post.mean')
tmp.data <- tmp.data[, ..cols.order]
colnames(tmp.data) <- c(
  'Aggregation ID',
  'Pre-Normalization Total Number of Reads',
  'Post-Normalization Total Number of Reads',
  'Pre-Normalization Mean Reads per Cell',
  'Post-Normalization Mean Reads per Cell'
)
# @ Output general aggregation metrics only.
tmp.file.name <- paste0(table.qc.path, '/QCSumm_Aggr.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)
# @ Merge per-library metrics.
ind.aggr.qcs <- lapply(X=aggr.qcs, FUN=function(x) return(x[['libs.qc']]))
ind.aggr.qcs <- rbindlist(l=ind.aggr.qcs, use.names=TRUE, idcol='aggr.id')

# ---> Per-library QC summaries.
# @ General QC information for GEx and HTO data.
sample.ids <- lapply(X=table.qc.objs, FUN=list.samples.info)
sample.ids <- unique(rbindlist(l=sample.ids, use.names=TRUE, fill=TRUE))

# @ GEx summries.
# Define metric summary files.
uniq.sample.ids <- unique(sample.ids[, sample.id])
sample.qc.files <- unique(paste0(uniq.sample.ids, '/metrics_summary.csv'))
to.check <- all(file.exists(sample.qc.files))
if(!to.check) stop('Unexpected error. Faulty GEx metrics summary files.\n')
# Load data as a single list keepin track of ID.
tmp.data <- lapply(X=sample.qc.files, FUN=fread)
names(tmp.data) <- uniq.sample.ids
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='sample.id')
# Merge metadata with metrics.
tmp.data[, library_id:=basename(str_replace(string=sample.id, pattern='/outs', replacement=''))
]
# Include per-library aggr information.
tmp.check <- tmp.data[, all(library_id %in% ind.aggr.qcs[, sample.id])]; if(!tmp.check) stop('No perfect intersection between aggr-related tables.\n')
tmp.data <- merge(x=tmp.data, y=ind.aggr.qcs, by.x='library_id', by.y='sample.id', all=TRUE)
tmp.data <- merge(x=tmp.data, y=sample.ids, by.x='library_id', by.y='library_id', all=TRUE)
# Redefine sample ID.
tmp.data[, sample.id:=library_id]
# Final formatting.
# * Col order
cols.order <- c('Estimated Number of Cells',	'Mean Reads per Cell',	'Median Genes per Cell',	'Number of Reads',	'Valid Barcodes',	'Fraction Reads in Cells',	'Total Genes Detected',	'Median UMI Counts per Cell', 'Reads Mapped to Genome',	paste0('Reads Mapped Confidently to ', c('Genome',	'Intergenic Regions', 'Intronic Regions',	'Exonic Regions',	'Transcriptome')), 'Reads Mapped Antisense to Gene', 'Sequencing Saturation',	'Q30 Bases in Barcode',	'Q30 Bases in RNA Read',	'Q30 Bases in RNA Read 2',	'Q30 Bases in UMI', colnames(ind.aggr.qcs)[2:length(ind.aggr.qcs)])
# all(cols.order %in% colnames(tmp.data))
cols.order <- c(setdiff(x=colnames(tmp.data), y=cols.order), cols.order)
cols.to.rm <- c('sample.id.x', 'sample.id.y', 'sample.id', 'molecule_h5', 'hto.tag', 'lane.tag', 'vdj.path') # Columns to excluded
cols.order <- setdiff(x=cols.order, y=cols.to.rm)
tmp.data <- tmp.data[, ..cols.order]
# * Row order.
tmp.order <- c(`CD4_Full`=1, `CD8_Full`=2)
tmp.data[, to.order:=tmp.order[aggr.id]]
setorderv(x=tmp.data, cols='to.order', order=1)
tmp.data[, to.order:=NULL]
# # Harmonize sequencing batches ### NOTE: THIS IS A STEP QUITE SPECIFIC FOR THIS PROJECT.
# tmp.data$seq.batch.tag <- unlist(lapply(X=tmp.data$seq.batch.tag, FUN=function(x){
#   x <- as.character(str_split(string=x, pattern=';', simplify=TRUE))
#   x <- paste0(seq.team.trans[x], collapse=';')
#   return(x)
# }))
# Get a copy of this file to aid in the development of another table. ### NOTE: THIS IS A STEP QUITE SPECIFIC FOR THIS PROJECT.
bkup.1 <- copy(tmp.data)
# Output table.
tmp.file.name <- paste0(table.qc.path, '/QCSumm_GEx.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)

# @ HTO summary
# Define metric summary files.
uniq.sample.ids <- unique(sample.ids[, .(sample.id, hto.tag)])
for.files <- uniq.sample.ids[, hto.tag]
for.files <- str_replace(string=for.files, pattern='^.+deconvolution/HTO_based/', replacement='')
tmp.dates <- str_extract(string=for.files, pattern='^\\d{2}-\\d{2}-\\d{4}')
for.files <- str_replace(string=for.files, pattern='^\\d{2}-\\d{2}-\\d{4}', replacement=paste0(tmp.dates, '/count'))
for.files <- str_extract(string=for.files, pattern='^([^/]+/){4}')
sample.qc.files <- paste0('/mnt/BioAdHoc/Groups/vd-vijay/vfajardo/sequencing_data/', for.files, '/outs/metrics_summary.csv')
to.check <- all(file.exists(sample.qc.files))
sample.qc.files[!file.exists(sample.qc.files)] <- NA
# Load data as a single list keepin track of ID.
tmp.data <- lapply(X=sample.qc.files, FUN=function(x){
  if(is.na(x)) return(NA) else return(fread(x))
})
names(tmp.data) <- uniq.sample.ids[, sample.id]
tmp.data <- tmp.data[!is.na(tmp.data)]
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='sample.id')
# Merge metadata with metrics.
tmp.data[, library_id:=basename(str_replace(string=sample.id, pattern='/outs', replacement=''))]
# Include per-library aggr information.
tmp.check <- tmp.data[, all(library_id %in% ind.aggr.qcs[, sample.id])]; if(!tmp.check) stop('No perfect intersection between aggr-related tables.\n')
tmp.data <- merge(x=tmp.data, y=ind.aggr.qcs, by.x='library_id', by.y='sample.id', all=TRUE)
tmp.data <- merge(x=tmp.data, y=sample.ids, by.x='library_id', by.y='library_id', all=TRUE)
# Redefine sample ID.
tmp.data[, sample.id:=library_id]
# Final formatting.
# * Col order
colnames(tmp.data) <- str_replace(string=colnames(tmp.data), pattern='^Antibody:\\s', replacement='')
cols.order <- c('Mean Reads per Cell',	'Valid Barcodes', 'Fraction Antibody Reads',	'Fraction Antibody Reads Usable',	'Antibody Reads Usable per Cell',	'Fraction Unrecognized Antibody',	'Antibody Reads in Cells',	'Median UMIs per Cell (summed over all recognized antibody barcodes)', 'Sequencing Saturation',	'Q30 Bases in Barcode',	'Q30 Bases in Antibody Read',	'Q30 Bases in UMI')
cols.order <- c(setdiff(x=colnames(tmp.data), y=cols.order), cols.order)
cols.to.rm <- c('sample.id.x', 'sample.id.y', 'sample.id', 'molecule_h5', 'hto.tag', 'lane.tag', 'vdj.path') # Columns to excluded
cols.order <- setdiff(x=cols.order, y=cols.to.rm)
tmp.data <- tmp.data[, ..cols.order]
# * Row order.
tmp.order <- c(`Full_CD4`=1, `Full_CD8`=2)
tmp.data[, to.order:=tmp.order[aggr.id]]
setorderv(x=tmp.data, cols='to.order', order=1)
tmp.data[, to.order:=NULL]
# Final touches to the library IDs
tmp.data[, library_id:=str_replace(string=library_id, pattern='_Gex|_GEX|_GEx', replacement='_CITE')] # Suffix
# Output table.
tmp.file.name <- paste0(table.qc.path, '/QCSumm_HTO.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)

# # @ TCR summary
# # Define metric summary files.
# uniq.sample.ids <- unique(sample.ids[, .(sample.id, vdj.path)])
# sample.qc.files <- paste0(uniq.sample.ids[, vdj.path], '/metrics_summary.csv')
# to.check <- all(file.exists(sample.qc.files))
# # Load data as a single list keepin track of ID.
# tmp.data <- lapply(X=sample.qc.files, FUN=fread)
# names(tmp.data) <- uniq.sample.ids[, sample.id]
# tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='sample.id')
# # Merge metadata with metrics.
# tmp.data[, sample.id:=
#   str_replace(
#     string=basename(str_replace(string=sample.id, pattern='/outs', replacement='')),
#     pattern='_Gex|_GEX', replacement=''
#   )
# ]
# tmp.data <- merge(y=tmp.data, x=ind.aggr.qcs, by.y='sample.id', by.x='sample.id', all=TRUE)
# tmp.data <- merge(x=tmp.data, y=sample.ids, by.x='sample.id', by.y='library_id', all=TRUE)
# # Redefine sample ID.
# tmp.data[, tmp.tag:=peptide.pool.tag]
# tmp.data[tmp.tag=='NA', tmp.tag:='UNSTIM']
# tmp.data[, sample.id:=paste(tract.tag, cell.type.tag, tmp.tag, 'TCR', sep='_')]
# # Final formatting.
# colnames(tmp.data) <- str_replace(string=colnames(tmp.data), pattern='^Antibody:\\s', replacement='')
# cols.order <- c('Estimated Number of Cells',	'Mean Read Pairs per Cell',	'Number of Cells With Productive V-J Spanning Pair',	'Number of Read Pairs',	'Valid Barcodes',	'Fraction Reads in Cells',	'Median TRA UMIs per Cell',	'Median TRB UMIs per Cell', 'Reads Mapped to Any V(D)J Gene',	'Reads Mapped to TRA',	'Reads Mapped to TRB',	'Mean Used Read Pairs per Cell', 'Q30 Bases in Barcode',	'Q30 Bases in RNA Read 1',	'Q30 Bases in RNA Read 2',	'Q30 Bases in Sample Index',	'Q30 Bases in UMI', 'Cells With Productive V-J Spanning Pair',	'Cells With Productive V-J Spanning (TRA, TRB) Pair',	'Paired Clonotype Diversity',	'Cells With TRA Contig',	'Cells With TRB Contig',	'Cells With CDR3-annotated TRA Contig',	'Cells With CDR3-annotated TRB Contig',	'Cells With V-J Spanning TRA Contig',	'Cells With V-J Spanning TRB Contig',	'Cells With Productive TRA Contig',	'Cells With Productive TRB Contig')
# # all(cols.order %in% colnames(tmp.data))
# cols.order <- c(setdiff(x=colnames(tmp.data), y=cols.order), cols.order)
# tmp.data <- tmp.data[, ..cols.order]
# setorderv(x=tmp.data, cols=c('cell.type.tag'))
# # Output table.
# tmp.file.name <- paste0(table.2.path, '/QCSumm_TCR.csv')
# fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)

# ---> Tidy table showing in which libraries each donor was sequenced.
# Combine all data.
# tmp.data.1 <- unique(bkup.1[, .(library_id=str_replace(string=library_id, pattern='_GEx', replacement=''), chrom.batch.tag, seq.batch.tag, cell.type.tag, donors.sort.batch.tag)])
# tmp.data.2 <- lapply(X=names(srt.objs.list), FUN=function(cell.type){
#     tmp.data <- as.data.table(srt.objs.list[[cell.type]]@meta.data)
#     tmp.data <- unique(
#         tmp.data[
#             !is.na(full.donor.id.tag),
#             .(full.donor.id.tag, donors.sort.batch.tag)
#         ]
#     )
#     return(tmp.data)
# })
# tmp.data.2 <- rbindlist(l=tmp.data.2, use.names=TRUE)
# tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='donors.sort.batch.tag')
# tmp.data.1 <- donor.demo[, .(donor.id.tag=as.character(donor.id.tag), Age, Sex=Gender)]
# tmp.data <- merge(x=tmp.data.1, y=tmp.data, by.x='donor.id.tag', by.y='full.donor.id.tag')
# # Spread cell info across columns.
# tmp.data[, cell.type.tag:=str_replace(string=donors.sort.batch.tag, pattern='\\.[A-Z]', replacement='')]
# tmp.data[, donors.sort.batch.tag:=NULL]
# tmp.data <- unique(tmp.data)
# tmp.data.1 <- tmp.data[,
#     .(
#         seq.batch.tag=paste(unique(
#             sort(unlist(str_split(string=seq.batch.tag, pattern=';')))
#         ), collapse=';')
#     ),
#     by=donor.id.tag
# ]
# tmp.data[, seq.batch.tag:=NULL]; tmp.data <- unique(tmp.data)
# tmp.data <- spread(data=tmp.data, key=cell.type.tag, value=library_id, fill='-')
# tmp.data <- merge(x=tmp.data, y=tmp.data.1, by='donor.id.tag')
# setorderv(x=tmp.data, col=c('chrom.batch.tag', 'CD4.5A'))
# # Output table.
# tmp.file.name <- paste0(table.qc.path, '/DonorDistributionAcrossLibs.csv')
# fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)

# # ---> Tidy table showing hashtag information per library ID.
# # Combine all data.
# tmp.data.1 <- unique(bkup.1[, .(library_id=str_replace(string=library_id, pattern='_GEx', replacement=''), chrom.batch.tag, cell.type.tag, donors.sort.batch.tag)])
# tmp.data.2 <- lapply(X=names(srt.objs.list), FUN=function(cell.type){
#     tmp.data <- as.data.table(srt.objs.list[[cell.type]]@meta.data)
#     tmp.data <- unique(
#         tmp.data[
#             !is.na(full.donor.id.tag),
#             .(full.donor.id.tag, donors.sort.batch.tag)
#         ]
#     )
#     return(tmp.data)
# })
# tmp.data.2 <- rbindlist(l=tmp.data.2, use.names=TRUE)
# tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='donors.sort.batch.tag')
# tmp.data.1 <- donor.demo[, .(donor.id.tag=as.character(donor.id.tag), Age, Sex=Gender)]
# tmp.data <- unique(merge(x=tmp.data.1, y=tmp.data, by.x='donor.id.tag', by.y='full.donor.id.tag'))
# # Combine with hashtag information.
# tmp.data.2 <- donor.meta[, .(donor.id.tag=as.character(donor.id.tag), chrom.batch.tag, hashtag.tag)]
# tmp.data <- merge(x=tmp.data, y=tmp.data.2, by=c('donor.id.tag', 'chrom.batch.tag'))
# # Fix hashtag info.
# tmp.data[, hashtag.tag:=str_replace(string=hashtag.tag, pattern='Hashtag', replacement='Hashtag ')]
# # Add TotalSeqC info.
# tmp.data.2 <- paste0('C0', seq(from=251, to=260, by=1))
# names(tmp.data.2) <- 1:length(tmp.data.2)
# tmp.data[, totalseq.c.id:=tmp.data.2[str_replace(string=hashtag.tag, pattern='Hashtag\\s[0]*', replacement='')]]
# # Set row order.
# # * Row order.
# tmp.order <- c(`CD8`=1, `CD4`=2, `NKB`=3, `MYE`=4)
# tmp.data[, to.order:=tmp.order[cell.type.tag]]
# setorderv(x=tmp.data, cols=c('to.order', 'library_id', 'hashtag.tag'), order=1)
# tmp.data[, to.order:=NULL]
# # Fix cell type info.
# tmp.data.2 <- c(`CD4`='CD4+ T cells', `CD8`='CD8+ T cells', `NKB`='B and NK cells', `MYE`='Myeloid cells')
# tmp.data[, cell.type.tag:=tmp.data.2[cell.type.tag]]
# # * Col order
# cols.order <- c(
#   `Library`='library_id',
#   `FACS sorting run ID`='chrom.batch.tag',
#   `Predefined immune cell population`='cell.type.tag',
#   `TotalSeq-C antibody ID`='totalseq.c.id',
#   `Hashtag ID`='hashtag.tag',
#   `Donor ID`='donor.id.tag',
#   `Sex`='Sex',
#   `Age`='Age'
# )
# cols.order %in% colnames(tmp.data)
# tmp.data <- tmp.data[, ..cols.order]
# colnames(tmp.data) <- names(cols.order)
# # Output table.
# tmp.file.name <- paste0(table.qc.path, '/HashtagInfoAcrossLibs.csv')
# fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)



############    -----------------------------------------    ############
### ------------------- On donor-specific reports ------------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
table.donors.path <- paste0(reports.path, '/table_on_donor_reports')
if(!dir.exists(table.donors.path)) dir.create(table.donors.path)

# ---> General definitions.
# Seurat objects of interest.
table.donors.objs <- obj.extended.names

# ---> Cell counts per condition per donor.
# @ Obtain each table.
# sep.conv.tables <- lapply(X=names(table.donors.objs), FUN=function(set.lab) get.sep.conv.table(table.obj=table.donors.objs[set.lab], set.lab=set.lab))
# names(sep.conv.tables) <- names(table.donors.objs)
# # @ Output each table.
# lapply(X=names(sep.conv.tables), FUN=function(table.obj){
#   # Retrieve data.
#   tmp.data <- sep.conv.tables[[table.obj]]
#   # Clusters order.
#   cluster.cols <- grep(x=colnames(tmp.data), pattern='^\\d+$', value=TRUE)
#   cluster.cols <- mixedsort(cluster.cols)
#   other.cols <- setdiff(x=colnames(tmp.data), y=cluster.cols)
#   cols.order <- c(other.cols, cluster.cols)
#   tmp.data <- tmp.data[, ..cols.order]
#   # Column names.
#   tmp.names <- c(
#     other.cols,
#     paste0('Cell count in cluster ', cluster.cols)
#   )
#   colnames(tmp.data) <- tmp.names
#   # Sort
#   tmp.data[str_count(string=Donor)==2, Donor:=paste0('00', Donor)]
#   setorderv(x=tmp.data, cols=c('Donor'), order=c(1))
#   # Further info. ### NOTE: THIS IS A STEP QUITE SPECIFIC FOR THIS PROJECT.
#   tmp.data.2 <- as.data.table(rbind(
#     colMeans(tmp.data[Donor!='Unassigned', 2:ncol(tmp.data)]),
#     apply(X=tmp.data[Donor!='Unassigned', 2:ncol(tmp.data)], MARGIN=2, FUN=median),
#     apply(X=tmp.data[Donor!='Unassigned', 2:ncol(tmp.data)], MARGIN=2, FUN=min),
#     apply(X=tmp.data[Donor!='Unassigned', 2:ncol(tmp.data)], MARGIN=2, FUN=max)
#   ))
#   tmp.data.2 <- cbind(c('Average', 'Median', 'Minimum', 'Maximum'), tmp.data.2)
#   tmp.data <- rbind(tmp.data, tmp.data.2, use.names=FALSE)
#   # Output
#   tmp.file.name <- paste0(table.donors.path, '/CellCountsPerCond_Set-', obj.extended.names[[table.obj]], '.csv')
#   fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=FALSE)
# })


############    -----------------------------------------    ############
### --------------------- On clusters' markers ---------------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
table.markers.path <- paste0(reports.path, '/table_on_cluster_markers')
if(!dir.exists(table.markers.path)) dir.create(table.markers.path)

# ---> General definitions.
# Seurat objects of interest.
table.markers.objs <- obj.extended.names

# ---> Clusters' markers.
# Output definition as a function.
output.markers <- function(obj.of.int, output.path, gene.id='ensembl'){
  # Retrieve and pre-process table.
  tmp.data <- markers.list[[obj.of.int]]
  if(all(is.na(tmp.data))) return(NULL)
  tmp.data <- tmp.data[p.adj<=0.05 & avg.lfc>=0.25]
  # Merge to include further gene info.
  tmp.data.2 <- gen.ref[, .(gene.id=Geneid, gene.name=gene_name, ens.biotype=gene_type, gen.location=paste0(Chr, ':', Start, '-', End))]
  # Sanity check.
  if(gene.id=='ensembl'){
    tmp.check <- tmp.data[!gene.id %in% tmp.data.2[, gene.id], .N]
  }else{
    tmp.check <- tmp.data[!gene.name %in% tmp.data.2[, gene.name], .N]
  }
  if(tmp.check > 0) stop('Some DEGs were not found in reference.\n')
  if(gene.id=='ensembl'){
    tmp.data <- merge(x=tmp.data, y=tmp.data.2, by=c('gene.id', 'gene.name'), all.x=TRUE, all.y=FALSE)
  }else{
    tmp.data[, gene.id:=NULL]
    tmp.data <- merge(x=tmp.data, y=tmp.data.2, by=c('gene.name'), all.x=TRUE, all.y=FALSE)
  }
  # Formatting.
  cols.order <- c(
    `Cluster`='element',
    `Ensembl ID`='gene.id', `Gene ID`='gene.name', `Ensembl biotype`='ens.biotype', `Genomic location of transcript`='gen.location',
    `Log2 fold change`='avg.lfc', `Adj. P value`='p.adj'
  )
  tmp.names <- c(
    names(cols.order),
    grep(x=colnames(tmp.data), pattern='mean.resolution.|mean.t.pop', value=TRUE),
    grep(x=colnames(tmp.data), pattern='prop.resolution.|prop.t.pop', value=TRUE)
  )
  cols.order <- c(
    cols.order,
    grep(x=colnames(tmp.data), pattern='mean.resolution.|mean.t.pop', value=TRUE),
    grep(x=colnames(tmp.data), pattern='prop.resolution.|prop.t.pop', value=TRUE)
  )
  tmp.data <- tmp.data[, ..cols.order]
  setorderv(x=tmp.data, cols=c('element', 'p.adj', 'avg.lfc'), order=c(1, 1, -1))
  colnames(tmp.data) <- tmp.names
  # Output.
  tmp.file.name <- paste0(output.path, '/ClusterMarkers_Set-', obj.of.int, '.csv')
  fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=FALSE)
  return(NA)
}
# Output each table.
for(table.obj in table.markers.objs){
  output.markers(obj.of.int=table.obj, output.path=table.markers.path, gene.id='name')
}


############    -----------------------------------------    ############
### ---------------------- On gene signatures ----------------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
# table.gsea.path <- paste0(reports.path, '/table_on_gene_signatures')
# if(!dir.exists(table.gsea.path)) dir.create(table.gsea.path)

# # ---> MSigDB table.
# # Output somewhere else and saved to the general data path. Loaded and defined here in variable "msigdb"

# # ---> GSEA results
# # * Preprocesing.
# gsea.res[, cell.type:=str_replace(string=obj.extended.names[data.set], pattern='Healthy-', replacement='')]
# gsea.res[, cluster.def:=str_replace(string=cluster.def, pattern='^\\w+\\.\\d+\\.', replacement='')]
# tmp.data.1 <- gsea.res
# tmp.data.2 <- unique(msigdb[, .(pathway=name, pathway.id=gs_id)])
# tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='pathway')
# # * Row order.
# tmp.order <- c(`Myeloid`=1, `NK`=2, `B`=3, `CD4`=4, `CD8`=5)
# tmp.data[, to.order:=tmp.order[cell.type]]
# setorderv(x=tmp.data, cols=c('to.order', 'cluster', 'pathway'), order=c(1, 1, 1))
# tmp.data[, to.order:=NULL]
# # Fix cell type info.
# tmp.data.2 <- c(`CD4`='CD4+ T cells', `CD8`='CD8+ T cells', `B`='B cells', `NK`='NK cells', `Myeloid`='Myeloid cells')
# tmp.data[, cell.type:=tmp.data.2[cell.type]]
# # * Col order
# cols.order <- c(
#   `Cell type`='cell.type',
#   `Cluster`='cluster',
#   `Signature ID`='pathway.id',
#   `Signature name`='pathway',
#   `Enrichment score`='ES',
#   `Normalized enrichment score`='NES',
#   `P-value`='pval',
#   `Adjusted P-value`='padj'
# )
# tmp.data <- tmp.data[, ..cols.order]
# colnames(tmp.data) <- names(cols.order)
# # Output table.
# tmp.file.name <- paste0(table.gsea.path, '/GSEAResults.csv')
# fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)
