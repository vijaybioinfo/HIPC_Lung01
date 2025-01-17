############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    -----   Preliminary figure stack    -----    ############
###########    ------------ BOOST project 1  -------------    ###########
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############

# ---> About the script.
# Version: 0
# Subversion: 1
# Best module: R/4.2.2


############    -----------------------------------------    ############
### -------------------------- Description -------------------------- ###
############    -----------------------------------------    ############

# Script to get all bioinformatics-based panel figures for any of the final main and supplementary figures for the paper from the BOOST 1 project.



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
# library(VennDiagram)
# library(gtools)
# library(ggnewscale)
# library(pheatmap)
# library(mefa)
# library(UpSetR)
# library(magrittr)
# library(PerformanceAnalytics)
# Only required for diversity metric calculation. Available only in v4.2
library(immunarch)
# Only required for correlation analyses. Available only in v4.2
library(Hmisc)
library(corrplot)


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
this.figure <- 'figure_stack_v1'
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
    `0`='#00FF7F', # TREG
    `1`='#FF7F00', # TH17 D
    `2`='#E41A1C', # Poly D
    `3`='#984EA3', # TFH D
    `4`='#377EB8', # TCM/Naive D
    `5`='#00FFFF', # TH2 D
    `6`='#FFDD33', # TFR
    `7`='#4B0082' # Unknown D
  ),
  `b1.cd8`=c(
    `0`='#FF1493', # GZMKhi D
    `1`='#8B4513', # exTEFF D
    `2`='#666666', # MAIT D
    `3`='#D37A48', # TEFF D
    `4`='#F781BF', # LTBhi D
    `5`='#4DAF4A', # Naive/TCM D
    `6`='#FF00FF' # IFNR D
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
alt.col.scale <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
# @ Others.
subset.n <- 50000
# ---> Path definitions.
# General data and reports paths.
prj.gen.path <- paste0('/mnt/BioAdHoc/Groups/vd-sette/', main.prj, '/paper_developments/', this.prj)
gen.data.path <- paste0(prj.gen.path, '/paper_items')
final.figs.path <- paste0(prj.gen.path, '/final_figures')
this.fig.path <- paste0(final.figs.path, '/', this.figure)
if(!dir.exists(this.fig.path)) dir.create(this.fig.path)
reports.path <- paste0(this.fig.path, '/', this.figure, '_panels')
# ---> File definitions.
# Seurat objects.
objs.files.list <- paste0(gen.data.path, '/', 'SeuratObj_', obj.extended.names, '.RDS')
names(objs.files.list) <- names(obj.extended.names)
# TCR data.
tcr.meta.files <- c(
    `b1.cd4`='/mnt/BioAdHoc/Groups/vd-sette/BOOST/paper_developments/BOOST_1/exploratory_analyses/preliminary_process/Round-4_CD4_2023-10-19/CellBasedTCRData.csv',
    `b1.cd4`='/mnt/BioAdHoc/Groups/vd-sette/BOOST/paper_developments/BOOST_1/exploratory_analyses/preliminary_process/Round-4_CD8_2023-10-19/CellBasedTCRData.csv'
)
# Gene sets' files (for module scoring & FGSE analyses).
module.feats.dir <- paste0(gen.data.path, '/module_signatures')
# Specific files.
#     Donor-specific cell counts before and after QC filtering.
# qc.inf.1.file <- paste0(gen.data.path, '/QCInf_CellCount_Donor.csv')
#     Annotated metadata for datasets before QC filtering.
# no.qc.meta.file <- paste0(gen.data.path, '/AnnotatedGExAggrMetaData.RDS')
#     Donors' metadata.
# donor.meta.file <- paste0(gen.data.path, '/DonorMetadata.csv')
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
### --------------------------- Functions --------------------------- ###
############    -----------------------------------------    ############


############    -----------------------------------------    ############
### ------------------------- Data Loading -------------------------- ###
############    -----------------------------------------    ############

# ---> Seurat obejcts.
# Main seurat objects.
srt.objs.list <- lapply(X=objs.files.list, FUN=readRDS)
names(srt.objs.list) <- names(objs.files.list)

# ---> TCR data.
b1.tcr.meta <- lapply(X=tcr.meta.files, FUN=fread)

# ---> Specific files.
#     Donor-specific cell counts before and after QC filtering.
# qc.inf.1 <- fread(file=qc.inf.1.file)

# ---> Differential analysis results
# @ Two-group comparisons
dea.res.files <- list.files(path=gen.data.path, pattern='DEA_Comp', full.names=TRUE)
dea.res.files <- dea.res.files[str_detect(string=dea.res.files, pattern='-vs-')]
names(dea.res.files) <- str_replace(
    string=basename(dea.res.files),
    pattern='\\.csv', replacement=''
)
dea.res <- lapply(X=dea.res.files, FUN=fread)
tmp.cols <- c('gene.id', 'p.val', 'p.adj', 'lfc')
dea.res <- lapply(X=dea.res, FUN=function(x) return(x[, ..tmp.cols]))
# dea.res <- rbindlist(l=dea.res, use.names=TRUE, idcol='test.id')

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


############    --- --------------------------------------    ############
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
# Set populations as factors.
srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]] <- factor(
    x=as.character(srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]]),
    levels=names(clusts.defs[[tmp.lab]])
)

# ---> Retrieve general GEx metadata.
meta.data <- lapply(X=names(srt.objs.list), FUN=function(data.set){
    seurat.obj <- srt.objs.list[[data.set]]
    tmp.data <- as.data.table(seurat.obj@meta.data)
    tmp.data$barcode <- row.names(seurat.obj@meta.data)
    tmp.tag <- clust.labs[data.set]
    tmp.data[, clusters.tag:=as.character(get(tmp.tag))]
    return(tmp.data)
})
names(meta.data) <- names(srt.objs.list)
meta.data <- rbindlist(l=meta.data, use.names=TRUE, idcol='data.set')
meta.data[, gex.cell.type.tag:=toupper(str_extract(string=data.set, pattern='cd[48]'))]
# Set factors.
tmp.lvls <- c('CD4', 'CD8')
tmp.vals <- factor(x=meta.data[, gex.cell.type.tag], levels=tmp.lvls)
set(x=meta.data, j='gex.cell.type.tag', value=tmp.vals)
meta.data$donor.id.tag <- factor(
    x=as.character(meta.data$donor.id.tag),
    levels=as.character(sort(unique(meta.data$donor.id.tag)))
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
#     Annotated metadata for datasets before QC filtering.
# names(no.qc.meta) <- str_replace(string=names(no.qc.meta), pattern='R24_Cancer_Batches-1-to-20_', replacement='')
# no.qc.meta <- no.qc.meta[str_detect(string=names(no.qc.meta), pattern='Normal-lung')]
# names(no.qc.meta) <- str_replace(string=names(no.qc.meta), pattern='_Normal-lung', replacement='')
# no.qc.meta <- rbindlist(l=no.qc.meta, use.names=TRUE, idcol='data.set')

# ---> Obtain cell subset of each dataset that's equal across cell types.
subset.list <- lapply(X=srt.objs.list, FUN=function(seurat.obj){
    sample(x=Cells(seurat.obj), size=subset.n)
})

# # ---> Number of cells per immune population.
# tmp.data <- lapply(X=srt.objs.list, FUN=function(seurat.obj){
#     tmp.data <- as.data.table(seurat.obj@meta.data)
#     tmp.data <- data.table(
#         total.in.obj=tmp.data[, .N],
#         total.ex.vivo=tmp.data[culture.cond.tag=='Ex vivo', .N],
#         total.ex.vivo.with.tcr=tmp.data[culture.cond.tag=='Ex vivo' & !is.na(clonotype.tag), .N]
#     )
#     return(tmp.data)
# })
# tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='data.set')
# tmp.file.name <- paste0(reports.path, '/CellCounts.csv')
# fwrite(x=tmp.data, file=tmp.file.name)


############    -----------------------------------------    ############
### --------------------------- Figure on --------------------------- ###
### --------------------- Intro to general data --------------------- ###
############    -----------------------------------------    ############

### ------------------------- Text details -------------------------- ###

# @ Total number of good-quality single-cell transcriptomes.
meta.data[, .N]

# @ Total number of good-quality single-cell transcriptomes per cell type.
meta.data[, .N, by=gex.cell.type.tag]

# @ Total number of samples
meta.data[!is.na(exp.sample.id.tag), uniqueN(exp.sample.id.tag), by=gex.cell.type.tag][, sum(V1)]

# @ Number of samples per cohort
tmp.data.1 <- meta.data[
    !is.na(dose.2.time.point.tag),
    .(
        cohort='dose.2',
        sample.count=uniqueN(exp.sample.id.tag))
    ,
    by=(gex.cell.type.tag)
]
tmp.data.2 <- meta.data[
    !is.na(dose.3.time.point.tag),
    .(
        cohort='dose.3',
        sample.count=uniqueN(exp.sample.id.tag))
    ,
    by=(gex.cell.type.tag)
]
tmp.data.3 <- meta.data[
    !is.na(dose.4.time.point.tag),
    .(
        cohort='dose.4',
        sample.count=uniqueN(exp.sample.id.tag))
    ,
    by=(gex.cell.type.tag)
]
tmp.data <- list(
    tmp.data.1,
    tmp.data.2,
    tmp.data.3
)
tmp.data <- rbindlist(l=tmp.data, use.names=T)
tmp.data[, sum(sample.count)]


# @ Median number of cells per sample for each T cell lineage.
meta.data[,
    .(cell.count=.N),
    by=.(data.set, ext.donor.id.tag)
][, median(cell.count), by=data.set]



############    -----------------------------------------    ############
### --------------------------- Figure on --------------------------- ###
### -------------------- Intro & TCR repertoires -------------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
fig.trc.reps.path <- paste0(reports.path, '/figure_on_tcr_reps')
if(!dir.exists(fig.trc.reps.path)) dir.create(fig.trc.reps.path)
# ---> Comments: Main dataset for this section, Lower-CD8 and Lower-CD4


### ------------------------- Text details -------------------------- ###


### -------------------------- Main Figure -------------------------- ###

# ---> Donor-specific repertoire diversity per major cell type compartment (CD4 and CD8).
# Retrieve diversity metrics.
# tmp.metrics <- c('gini', 'inv.simp', 'gini.simp', 'chao1')
# tmp.metrics <- c('gini', 'inv.simp')
tmp.metrics <- c('gini')
div.metrics <- lapply(
    X=tmp.metrics,
    FUN=function(x){
        cat(x, '\n')
        tmp.data <- summ.div.metrics(
            tcr.meta=b1.tcr.meta,
            donor.meta=donor.meta, div.metric=x
        )
        return(tmp.data)
    }
)
names(div.metrics) <- tmp.metrics
div.metrics <- rbindlist(l=div.metrics, use.names=TRUE, idcol='metric')
# Add extended timepoint
div.metrics[, 
    gen.time.point.tag:=paste(
        str_extract(string=dose.cohort.tag, pattern='\\d$'), 
        time.point.tag,
        sep=' - '
    )
]
tmp.lvls <- c(0, 1, 3, 6, 9, 12)
tmp.lvls <- paste(
    rep(x=2:4, each=length(tmp.lvls)), tmp.lvls, sep=' - '
)
div.metrics$gen.time.point.tag <- factor(
    x=div.metrics$gen.time.point.tag,
    levels=tmp.lvls
)
# Plot.
for(tmp.metric in tmp.metrics){
    for(cell.type in div.metrics[, unique(cell.type.tag)]){
        tmp.data <- div.metrics[
            metric==tmp.metric & cell.type.tag==cell.type
        ]
        metric.lab <- paste0(tmp.metric, '.All')
        tmp.data[, metric.lab:=get(metric.lab)]
        tmp.data[, metric.lab:=1-metric.lab]
        # Plot
        tmp.ggplot <- ggplot(data=tmp.data, aes(x=gen.time.point.tag, y=metric.lab)) +
            geom_boxplot(alpha=0, width=0.7, linewidth=3, outlier.shape=NA, color='black') +
            geom_jitter(aes(fill=dose.cohort.tag), shape=21, size=6, stroke=2.5, color='black', width=0, height=0.01) +
            scale_y_continuous(breaks=scales::pretty_breaks(n=4)) +
            # scale_color_manual(values=cell.type.cols) +
            labs(
                x='Time point',
                y=paste0('Diversity [1-', tmp.metric, ']'),
                fill='Dose cohort'
            )
        tmp.lab <- paste0('/Metric-', tmp.metric, '_GenTimePoint_Donor_', cell.type)
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.trc.reps.path, file.name=tmp.lab, type='pdf',
            blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=16, height=10
        )
    }
}
# Output table w/ metric values.
tmp.file.name <- paste0(fig.trc.reps.path, '/DivMetric_ExtDetails.csv')
fwrite(file=tmp.file.name, x=div.metrics, na=NA, quote=TRUE)


### ------------------------- Supp. figure -------------------------- ###

# ---> QC features distribution after QC filtering.
# qc.feats <- c(
#     `UMICount`='nCount_RNA',
#     `GeneCount`='nFeature_RNA',
#     `MitUMIFreq`='percent.mt'
# )
# for(qc.feat in names(qc.feats)){
#     tmp.feat <- qc.feats[[qc.feat]]
#     meta.data[, tmp.val:=get(tmp.feat)+1]
#     tmp.ggplot <- ggplot(data=meta.data, aes(x=facs.sorting.batch.tag, y=tmp.val, fill=facs.sorting.batch.tag)) +
#         geom_violin(alpha=0.8, trim=FALSE, adjust=0.8) +
#         geom_boxplot(width=0.05, alpha=0.7, outlier.shape=NA) +
#         scale_y_log10() +
#         labs(x='FACS batch', y=paste0(qc.feat, ' (log10)')) +
#         theme_bw() + theme(legend.position='none')
#     tmp.lab <- paste0('/QCMetric-', qc.feat, '_FACSBatch')
#     publish.plot(
#         tmp.ggplot=tmp.ggplot, output.path=fig.trc.reps.path, file.name=tmp.lab, type='pdf',
#         blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=14, height=3.5
#     )
#     meta.data[, tmp.val:=NULL]
# }

# # ---> Donor-specific cell counts before after filtering.
# tmp.data <- as.data.table(gather(data=qc.inf.1, key='qc.stage', value='cell.count', -`pre.donor.id.tag`))
# tmp.vals <- c(
#     `pre.filter.count`='Pre-filter',
#     `post.filter.count`='Post-filter'
# )
# tmp.data <- tmp.data[qc.stage %in% names(tmp.vals)]
# tmp.data[, qc.stage:=factor(x=tmp.vals[qc.stage], levels=tmp.vals)]
# tmp.ggplot <- ggplot(data=tmp.data, aes(x=qc.stage, y=cell.count)) +
#     geom_boxplot(color='black', width=0.7, linewidth=9, outlier.shape=NA) +
#     geom_jitter(size=10, color='black', width=0.3) +
#     scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
#     labs(x='QC stage', y='Number of cells')
# tmp.lab <- paste0(
#     '/CellCount_QCStage_Donor'
# )
# publish.plot(
#     tmp.ggplot=tmp.ggplot, output.path=fig.trc.reps.path, file.name=tmp.lab, type='pdf',
#     blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=7, height=14
# )


############    -----------------------------------------    ############
### --------------------------- Figure on --------------------------- ###
### -------------- CD4 T cell transcriptomic phenotype -------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
fig.cd4.phen.path <- paste0(reports.path, '/figure_on_cd4_phen')
if(!dir.exists(fig.cd4.phen.path)) dir.create(fig.cd4.phen.path)
# ---> Comments: Main dataset for this section, CD4 T cell set.
main.obj.name <- 'b1.cd4'


### ------------------------- Text details -------------------------- ###

# @ Fraction accounted for by representative subsets
# Cluster 0
tmp.vals <- c('0')
meta.data[
    data.set=='b1.cd4' & clusters.tag %in% tmp.vals
    ,
.N]/meta.data[data.set=='b1.cd4', .N]
# Clusters 0, 1, 2 and 3.
tmp.vals <- c('0', '1', '2', '3')
meta.data[
    data.set=='b1.cd4' & clusters.tag %in% tmp.vals
    ,
.N]/meta.data[data.set=='b1.cd4', .N]


### -------------------------- Main Figure -------------------------- ###

# ---> UMAP plots depicting clusters.
# @ Tiff versions
# Whole dataset.
tmp.ggplot <- get.umap.gg(obj.name=main.obj.name, attempted.format='tiff', cell.subset=NULL)
tmp.lab <- paste0(obj.extended.names[main.obj.name], '_PopsOnUMAP_Whole')
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.cd4.phen.path, file.name=tmp.lab, type='tiff',
    blank.comp=blank.complement.1, comp.comp=theme_classic(), do.legend=FALSE
)
# Subset
tmp.ggplot <- get.umap.gg(obj.name=main.obj.name, attempted.format='tiff', cell.subset=subset.list[[main.obj.name]])
tmp.lab <- paste0(obj.extended.names[main.obj.name], '_PopsOnUMAP_Subset')
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.cd4.phen.path, file.name=tmp.lab, type='tiff',
    blank.comp=blank.complement.1, comp.comp=theme_classic(), do.legend=FALSE
)
# @ PDF versions
# Whole dataset.
tmp.ggplot <- get.umap.gg(obj.name=main.obj.name, attempted.format='pdf', cell.subset=NULL)
tmp.lab <- paste0(obj.extended.names[main.obj.name], '_PopsOnUMAP_Whole')
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.cd4.phen.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.1, comp.comp=theme_classic(), do.legend=FALSE
)
# Subset
tmp.ggplot <- get.umap.gg(obj.name=main.obj.name, attempted.format='pdf', cell.subset=subset.list[[main.obj.name]])
tmp.lab <- paste0(obj.extended.names[main.obj.name], '_PopsOnUMAP_Subset')
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.cd4.phen.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.1, comp.comp=theme_classic(), do.legend=FALSE
)

# ---> Proportion of cells per cluster.
# tmp.ggplot <- get.cell.fracs.gg(
#     seurat.obj=srt.objs.list[[main.obj.name]],
#     x.var=NULL, fill.var=clust.labs[main.obj.name], fill.cols=clusts.cols[[main.obj.name]]
# )
# tmp.lab <- paste0(obj.extended.names[main.obj.name], '_CellFracs_AllCells_Clusters')
# publish.plot(
#     tmp.ggplot=tmp.ggplot, output.path=fig.cd4.phen.path, file.name=tmp.lab, type='pdf',
#     blank.comp=blank.complement.3.1, do.legend=FALSE,
#     width=2, height=21
# )

# ---> Tag-specific analysis
# @ Dose-specific cohort
tmp.data <- get.tag.analysis(
  seurat.obj=srt.objs.list[[main.obj.name]],
  tmp.tag='dose.cohort.tag', clusters.tag=clust.labs[[main.obj.name]]
)
tmp.data$index.tag <- factor(x=as.character(tmp.data$index.tag), levels=rev(c('dose.2', 'dose.3', 'dose.4')))
# Plot
tmp.ggplot <- ggplot(data=tmp.data, aes(x=cluster.tag, y=scl.abs.freq, fill=index.tag)) +
    geom_bar(stat='identity', position='fill', width=0.6, color='black', linewidth=1.5) +
    scale_y_continuous(expand=expansion(add=c(0, 0)), breaks=scales::pretty_breaks(n=3)) +
    scale_fill_manual(values=cohort.cols) +
    labs(x='Cluster', y='Cell fraction', fill='Dose cohort')
# Output.
tmp.width <- tmp.data[, uniqueN(cluster.tag)] * 0.9 + 0.3
tmp.height <- 7
tmp.lab <- paste0(obj.extended.names[main.obj.name], '_NormFracts-DoseCohort')
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.cd4.phen.path, file.name=tmp.lab, type='pdf', blank.comp=blank.complement.3, do.legend=FALSE,
    width=tmp.width, height=tmp.height
)

# ---> Dot plots.
# Define genes of interest.
these.markers <- c(
    'FOXP3', 'TIGIT', 'CTLA4', 'IKZF2', # TREG
    'IL1R2', 'IL1R1', #TFR
    'TNF', 'CSF2', 'IFNG', 'IL2', # TH1
    'IL21', 'POU2AF1', 'BTLA', 'CD200', 'CXCL13', 'PDCD1', # TFH
    'RORA', 'CCR6', 'CTSH', 'IL4I1', # TH17
    'IL5', 'IL13', 'GATA3', # TH2
    'CCR7', 'TCF7', 'IL7R' # TCM
)
# Define clusters' order.
these.pops <- c('6', '0', '2', '7', '3', '1', '5', '4')
# Plot
tmp.ggplot <- dot.plot(
    seurat.obj=srt.objs.list[[main.obj.name]], features=these.markers, slot='data', do.norm=FALSE, ensembl=FALSE,
    groups.tag=clust.labs[main.obj.name], groups.order=these.pops, groups.of.int=NULL, filter.tag=NULL, groups.to.filter=NULL, keep=FALSE, na.rm=TRUE, feature.thold=NULL,
    this.color.scale=signatures.col.scale, col.min=NULL, col.max=NULL,
    scale.by='radius', dot.scale=12, size.min=NA, size.max=NA,
    file.name=NULL
)
tmp.ggplot <- tmp.ggplot + theme(legend.position='bottom')
# Output.
tmp.width <- (length(these.pops) * 0.5) + 0.3
tmp.height <- (length(these.markers) * 0.4) + 0.3
tmp.lab <- paste0(obj.extended.names[main.obj.name], '_Markers_Opt-A')
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.cd4.phen.path, file.name=tmp.lab, type='pdf', blank.comp=blank.complement.2, do.legend=TRUE,
    width=tmp.width, height=tmp.height,
    c.width=tmp.width+2, c.height=tmp.height,
    legend.height=2.6, legend.width=0.5
)

# # ---> Expression of specific genes
genes.of.int <- c(
    'FOXP3', 'CCR6'
)
scale.tholds <- c(
    FOXP3=NA,
    CCR6=NA
)
get.gene.exp.plots(
    seurat.obj=srt.objs.list[[main.obj.name]], 
    genes.of.int=genes.of.int, scale.tholds=scale.tholds,
    col.scale=alt.col.scale,
    output.path=fig.cd4.phen.path, dot.size=1
)

# ---> Module scoring calculation (if necessary).
mods.of.int <- c(
    # TH1='th1.signature1.arlehamn.score',
    TH17='th17.signature2.arlehamndan.score', TFH='tfh.signature.locci.score',
    TFR='tfr.signature.guo.score',
    # TREG='treg.signature.schmiedle.score'
    
)
if(!all(mods.of.int %in% colnames(srt.objs.list[[main.obj.name]]@meta.data))){
    mod.in.path <- paste0(module.feats.dir)
    mod.out.path <- paste0(module.feats.dir, '/module_signatures_', obj.extended.names[main.obj.name])
    if(!dir.exists(mod.out.path)) dir.create(mod.out.path)
    srt.objs.list[[main.obj.name]] <- get.module.scores(
        seurat.obj=srt.objs.list[[main.obj.name]], module.feats.dir=mod.in.path, reports.path=mod.out.path,
        signal.thold=0, is.ensembl=FALSE
    )
}

# # ---> Module scoring as UMAP heatmaps.
mods.of.int <- c(
    # TH1='th1.signature1.arlehamn.score',
    TH17='th17.signature2.arlehamndan.score', TFH='tfh.signature.locci.score',
    TFR='tfr.signature.guo.score'
    # TREG='treg.signature.schmiedle.score'
    
)
scale.tholds <- list(
    # TH1=NA,
    TH17=NA, TFH=NA, TFR=NA
    # TREG=NA
)
get.mod.score.plots(
    seurat.obj=srt.objs.list[[main.obj.name]], 
    mods.of.int=mods.of.int, scale.tholds=scale.tholds,
    output.path=fig.cd4.phen.path, dot.size=0.3
)

# ---> FGSEA
# @ TREG population.
srt.objs.list[[main.obj.name]]@meta.data[, 'treg.tag'] <- ifelse(
  test=str_detect(
    string=clusts.defs[[main.obj.name]][
        as.character(srt.objs.list[[main.obj.name]]@meta.data[, clust.labs[main.obj.name]])
    ],
    pattern='^TREG'
  ),
  yes='TREG',
  no='Rest'
)
do.fgsea(
    metrics.src=srt.objs.list[[main.obj.name]], modules.feats=modules.feats['treg.signature.schmiedle'],  metric='signal.to.noise', cell.subset=subset.list[[main.obj.name]],
    tag.of.int='treg.tag',
    output.path=fig.cd4.phen.path, vals.to.depict='TREG', files.preffix=paste0(obj.extended.names[main.obj.name])
)
# @ TFR population.
srt.objs.list[[main.obj.name]]@meta.data[, 'tfr.tag'] <- ifelse(
  test=str_detect(
    string=clusts.defs[[main.obj.name]][
        as.character(srt.objs.list[[main.obj.name]]@meta.data[, clust.labs[main.obj.name]])
    ],
    pattern='^TFR'
  ),
  yes='TFR',
  no='Rest'
)
do.fgsea(
    metrics.src=srt.objs.list[[main.obj.name]], modules.feats=modules.feats['tfr.signature.guo'],  metric='signal.to.noise', cell.subset=subset.list[[main.obj.name]],
    tag.of.int='tfr.tag',
    output.path=fig.cd4.phen.path, vals.to.depict='TFR', files.preffix=paste0(obj.extended.names[main.obj.name])
)
# # @ TFH population.
srt.objs.list[[main.obj.name]]@meta.data[, 'tfh.tag'] <- ifelse(
  test=str_detect(
    string=clusts.defs[[main.obj.name]][
        as.character(srt.objs.list[[main.obj.name]]@meta.data[, clust.labs[main.obj.name]])
    ],
    pattern='^TFH'
  ),
  yes='TFH',
  no='Rest'
)
do.fgsea(
    metrics.src=srt.objs.list[[main.obj.name]], modules.feats=modules.feats['tfh.signature.locci'],  metric='signal.to.noise', cell.subset=subset.list[[main.obj.name]],
    tag.of.int='tfh.tag',
    output.path=fig.cd4.phen.path, vals.to.depict='TFH', files.preffix=paste0(obj.extended.names[main.obj.name])
)
# @ TH17 population.
srt.objs.list[[main.obj.name]]@meta.data[, 'th17.tag'] <- ifelse(
  test=str_detect(
    string=clusts.defs[[main.obj.name]][
        as.character(srt.objs.list[[main.obj.name]]@meta.data[, clust.labs[main.obj.name]])
    ],
    pattern='^TH17'
  ),
  yes='TH17',
  no='Rest'
)
do.fgsea(
    metrics.src=srt.objs.list[[main.obj.name]], modules.feats=modules.feats['th17.signature2.arlehamndan'],  metric='signal.to.noise', cell.subset=subset.list[[main.obj.name]],
    tag.of.int='th17.tag',
    output.path=fig.cd4.phen.path, vals.to.depict='TH17', files.preffix=paste0(obj.extended.names[main.obj.name])
)

# ---> Fraction of cells per cluster per donor.
# tmp.ggplot <- get.cell.fracs.gg(
#     seurat.obj=srt.objs.list[[main.obj.name]],
#     x.var='donor.id.tag', fill.var=clust.labs[main.obj.name], fill.cols=clusts.cols[[main.obj.name]],
#     bar.width=0.6, line.width=1.3
# )
# tmp.lab <- paste0(obj.extended.names[main.obj.name], '_CellFracs_DonorID_Clusters')
# publish.plot(
#     tmp.ggplot=tmp.ggplot, output.path=fig.cd4.phen.path, file.name=tmp.lab, type='pdf',
#     blank.comp=blank.complement.3.1, do.legend=FALSE,
#     width=21, height=7
# )


############    -----------------------------------------    ############
### --------------------------- Figure on --------------------------- ###
### -------------- CD8 T cell transcriptomic phenotype -------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
fig.cd8.phen.path <- paste0(reports.path, '/figure_on_cd8_phen')
if(!dir.exists(fig.cd8.phen.path)) dir.create(fig.cd8.phen.path)
# ---> Comments: Main dataset for this section, CD8 T cell set.
main.obj.name <- 'b1.cd8'


### ------------------------- Text details -------------------------- ###

# @ Fraction accounted for by representative subsets
# Cluster 0
tmp.vals <- c('0')
meta.data[
    data.set=='b1.cd4' & clusters.tag %in% tmp.vals
    ,
.N]/meta.data[data.set=='b1.cd8', .N]


### -------------------------- Main Figure -------------------------- ###

# ---> UMAP plots depicting clusters.
# @ Retrieve temporary object.
tmp.obj <- srt.objs.list[[main.obj.name]]
# Remove cells beyond 7.5 on UMAP 1 scale.
tmp.data <- Cells(tmp.obj)[!tmp.obj@reductions$umap@cell.embeddings[, 'UMAP_1'] > 7.5]
tmp.obj <- subset(x=tmp.obj, cells=tmp.data)
# ---> Obtain cell subset of each dataset that's equal across cell types.
tmp.subset <- sample(x=Cells(tmp.obj), size=subset.n)

# @ Tiff versions
# Whole dataset.
tmp.ggplot <- get.umap.gg(
    obj.name=main.obj.name, seurat.obj=tmp.obj,
    attempted.format='tiff', cell.subset=NULL
)
tmp.lab <- paste0(obj.extended.names[main.obj.name], '_PopsOnUMAP_Whole')
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.cd8.phen.path, file.name=tmp.lab, type='tiff',
    blank.comp=blank.complement.1, comp.comp=theme_classic(), do.legend=FALSE
)
# Subset
tmp.ggplot <- get.umap.gg(
    obj.name=main.obj.name, seurat.obj=tmp.obj,
    attempted.format='tiff', cell.subset=tmp.subset
)
tmp.lab <- paste0(obj.extended.names[main.obj.name], '_PopsOnUMAP_Subset')
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.cd8.phen.path, file.name=tmp.lab, type='tiff',
    blank.comp=blank.complement.1, comp.comp=theme_classic(), do.legend=FALSE
)
# @ PDF versions
# Whole dataset.
tmp.ggplot <- get.umap.gg(
    obj.name=main.obj.name, seurat.obj=tmp.obj,
    attempted.format='pdf', cell.subset=NULL
)
tmp.lab <- paste0(obj.extended.names[main.obj.name], '_PopsOnUMAP_Whole')
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.cd8.phen.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.1, comp.comp=theme_classic(), do.legend=FALSE
)
# Subset
tmp.ggplot <- get.umap.gg(
    obj.name=main.obj.name, seurat.obj=tmp.obj,
    attempted.format='pdf', cell.subset=tmp.subset
)
tmp.lab <- paste0(obj.extended.names[main.obj.name], '_PopsOnUMAP_Subset')
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.cd8.phen.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.1, comp.comp=theme_classic(), do.legend=FALSE
)

# ---> Proportion of cells per cluster.
# tmp.ggplot <- get.cell.fracs.gg(
#     seurat.obj=srt.objs.list[[main.obj.name]],
#     x.var=NULL, fill.var=clust.labs[main.obj.name], fill.cols=clusts.cols[[main.obj.name]]
# )
# tmp.lab <- paste0(obj.extended.names[main.obj.name], '_CellFracs_AllCells_Clusters')
# publish.plot(
#     tmp.ggplot=tmp.ggplot, output.path=fig.cd8.phen.path, file.name=tmp.lab, type='pdf',
#     blank.comp=blank.complement.3.1, do.legend=FALSE,
#     width=2, height=21
# )

# ---> Tag-specific analysis
# @ Dose-specific cohort
tmp.data <- get.tag.analysis(
  seurat.obj=srt.objs.list[[main.obj.name]],
  tmp.tag='dose.cohort.tag', clusters.tag=clust.labs[[main.obj.name]]
)
tmp.data$index.tag <- factor(x=as.character(tmp.data$index.tag), levels=rev(c('dose.2', 'dose.3', 'dose.4')))
# Plot
tmp.ggplot <- ggplot(data=tmp.data, aes(x=cluster.tag, y=scl.abs.freq, fill=index.tag)) +
    geom_bar(stat='identity', position='fill', width=0.6, color='black', linewidth=1.5) +
    scale_y_continuous(expand=expansion(add=c(0, 0)), breaks=scales::pretty_breaks(n=3)) +
    scale_fill_manual(values=cohort.cols) +
    labs(x='Cluster', y='Cell fraction', fill='Dose cohort')
# Output.
tmp.width <- tmp.data[, uniqueN(cluster.tag)] * 0.9 + 0.3
tmp.height <- 7
tmp.lab <- paste0(obj.extended.names[main.obj.name], '_NormFracts-DoseCohort')
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.cd8.phen.path, file.name=tmp.lab, type='pdf', blank.comp=blank.complement.3, do.legend=FALSE,
    width=tmp.width, height=tmp.height
)

# ---> Dot plots.
# Define genes of interest.
these.markers <- c(
    'ISG15', 'MX1', 'OAS1', # IFNR
    'KLRB1', 'TRAV1-2', 'TRBV20-1', # MAIT
    # 'S1PR4', 'KLF2', 'GZMK', 'FGFBP2', 'KLRC2', # GZMK-expressing
    'GZMK', 'FGFBP2', 'KLRC2', # GZMK-expressing
    'LTB', 'LGALS1', 'LTA',  # LTBhi
    # 'GNLY', 'GZMB', 'PRF1', 'IFNG', 'TNF', # Effector
    'GNLY', 'GZMB', 'PRF1', 'IFNG', 'TNF', # Effector
    'TOX', 'PDCD1', 'HAVCR2', 'LAG3', # Exhausted
    'CCR7', 'TCF7', 'IL7R' # TCM
)
# Define clusters' order.
these.pops <- c('6', '2', '0', '4', '3', '1', '5')
# # Plot
tmp.ggplot <- dot.plot(
    seurat.obj=srt.objs.list[[main.obj.name]], features=these.markers, slot='data', do.norm=FALSE, ensembl=FALSE,
    groups.tag=clust.labs[main.obj.name], groups.order=these.pops, groups.of.int=NULL, filter.tag=NULL, groups.to.filter=NULL, keep=FALSE, na.rm=TRUE, feature.thold=NULL,
    this.color.scale=signatures.col.scale, col.min=NULL, col.max=NULL,
    scale.by='radius', dot.scale=12, size.min=NA, size.max=NA,
    file.name=NULL
)
tmp.ggplot <- tmp.ggplot + theme(legend.position='bottom')
# Output.
tmp.width <- (length(these.pops) * 0.5) + 0.3
tmp.height <- (length(these.markers) * 0.4) + 0.3
tmp.lab <- paste0(obj.extended.names[main.obj.name], '_Markers_Opt-A')
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.cd8.phen.path, file.name=tmp.lab, type='pdf', blank.comp=blank.complement.2, do.legend=TRUE,
    width=tmp.width, height=tmp.height,
    c.width=tmp.width+2, c.height=tmp.height,
    legend.height=2.6, legend.width=0.5
)

# ---> Module scoring calculation (if necessary).
mods.of.int <- c(Cytotoxicity='tcell.cytotoxicity.guo.score', Exhaustion='exhaustion.consensus.score')
if(!all(mods.of.int %in% colnames(srt.objs.list[[main.obj.name]]@meta.data))){
    mod.in.path <- paste0(module.feats.dir)
    mod.out.path <- paste0(module.feats.dir, '/module_signatures_', obj.extended.names[main.obj.name])
    if(!dir.exists(mod.out.path)) dir.create(mod.out.path)
    srt.objs.list[[main.obj.name]] <- get.module.scores(
        seurat.obj=srt.objs.list[[main.obj.name]], module.feats.dir=mod.in.path, reports.path=mod.out.path,
        signal.thold=0, is.ensembl=FALSE
    )
}

# ---> Module scoring as UMAP heatmaps.
mods.of.int <- c(Cytotoxicity='tcell.cytotoxicity.guo.score', Exhaustion='exhaustion.consensus.score')
scale.tholds <- list(
    Cytotoxicity=NA,
    Exhaustion=0.7
)
get.mod.score.plots(
    seurat.obj=srt.objs.list[[main.obj.name]], 
    mods.of.int=mods.of.int, scale.tholds=scale.tholds,
    output.path=fig.cd8.phen.path, dot.size=0.3
)

# ---> FGSEA
# @ Effector CD8 T cell populations.
srt.objs.list[[main.obj.name]]@meta.data[, 'effector.tag'] <- ifelse(
  test=str_detect(
    string=clusts.defs[[main.obj.name]][
        as.character(srt.objs.list[[main.obj.name]]@meta.data[, clust.labs[main.obj.name]])
    ],
    pattern='TEFF'
  ),
  yes='Effector',
  no='Rest'
)
do.fgsea(
    metrics.src=srt.objs.list[[main.obj.name]], modules.feats=modules.feats['tcell.cytotoxicity.guo'],  metric='signal.to.noise', cell.subset=subset.list[[main.obj.name]],
    tag.of.int='effector.tag',
    output.path=fig.cd8.phen.path, vals.to.depict='Effector', files.preffix=paste0(obj.extended.names[main.obj.name])
)
# @ T cell exhaustion consensus signature
srt.objs.list[[main.obj.name]]@meta.data[, 'exhaust.tag'] <- ifelse(
  test=str_detect(
    string=clusts.defs[[main.obj.name]][
        as.character(srt.objs.list[[main.obj.name]]@meta.data[, clust.labs[main.obj.name]])
    ],
    pattern='^ex'
  ),
  yes='Exhausted',
  no='Rest'
)
do.fgsea(
    metrics.src=srt.objs.list[[main.obj.name]], modules.feats=modules.feats['exhaustion.consensus'],  metric='signal.to.noise', cell.subset=subset.list[[main.obj.name]],
    tag.of.int='exhaust.tag',
    output.path=fig.cd8.phen.path, vals.to.depict='Exhausted', files.preffix=paste0(obj.extended.names[main.obj.name])
)
# @ Exhaustion signature, comparison between early and late timepoints for all cells.
tmp.cells <- Cells(srt.objs.list[[main.obj.name]])[!is.na(srt.objs.list[[main.obj.name]]@meta.data[, 'evl.sample.tag'])]
do.fgsea(
    metrics.src=srt.objs.list[[main.obj.name]], modules.feats=modules.feats['exhaustion.consensus'],  metric='signal.to.noise',
    cell.subset=tmp.cells,
    tag.of.int='evl.sample.tag',
    output.path=fig.cd8.phen.path, vals.to.depict='Late', files.preffix=paste0(obj.extended.names[main.obj.name])
)
# @ Exhaustion signature, comparison between early and late timepoints for all cells.
tmp.cells <- Cells(srt.objs.list[[main.obj.name]])[
    !is.na(srt.objs.list[[main.obj.name]]@meta.data[, 'evl.sample.tag']) &
    srt.objs.list[[main.obj.name]]@meta.data[, clust.labs[[main.obj.name]]]=='1'
]
do.fgsea(
    metrics.src=srt.objs.list[[main.obj.name]], modules.feats=modules.feats['exhaustion.consensus'],  metric='signal.to.noise',
    cell.subset=tmp.cells,
    tag.of.int='evl.sample.tag',
    output.path=fig.cd8.phen.path, vals.to.depict='Late',
    files.preffix=paste0(obj.extended.names[main.obj.name], '_C-1_')
)

# ---> Fraction of cells per cluster per donor.
# tmp.ggplot <- get.cell.fracs.gg(
#     seurat.obj=srt.objs.list[[main.obj.name]],
#     x.var='donor.id.tag', fill.var=clust.labs[main.obj.name], fill.cols=clusts.cols[[main.obj.name]],
#     bar.width=0.6, line.width=1.3
# )
# tmp.lab <- paste0(obj.extended.names[main.obj.name], '_CellFracs_DonorID_Clusters')
# publish.plot(
#     tmp.ggplot=tmp.ggplot, output.path=fig.cd8.phen.path, file.name=tmp.lab, type='pdf',
#     blank.comp=blank.complement.3.1, do.legend=FALSE,
#     width=21, height=7
# )

############    -----------------------------------------    ############
### --------------------------- Figure on --------------------------- ###
### --------------------- Differential analysis --------------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
fig.dea.path <- paste0(reports.path, '/figure_on_dea')
if(!dir.exists(fig.dea.path)) dir.create(fig.dea.path)
# ---> Comments: Main dataset for this section, none. All of them were used.
main.obj.name <- NULL


### -------------------------- Main Figure -------------------------- ###

# @ Custom gene labels
custom.gene.labs <- list(
    `CD4_Full_DEA_Comp-0-vs-1_Subset-C-3-DC-2`=c(
        'CXCL13', 'IL2RA', 'CXCR5',
        'IL2', 'CD200')
)
# @ Custom P-value bounds.
custom.p.bounds <- c(
    `CD4_Full_DEA_Comp-0-vs-1_Subset-C-3-DC-2`=200,
    `CD4_Full_DEA_Comp-CD4REpos-vs-CD4REneg_Subset-C-1`=70,
    `CD8_Full_DEA_Comp-CD4REpos-vs-CD4REneg_Subset-C-0-4-6`=70
)

# @ Process per test
for(test.id in names(dea.res)){
    tmp.data <- copy(dea.res[[test.id]])
    # DEG classes
    tmp.data[, gen.deg.class:=
        p.adj<=0.05 & (
            lfc<=(-0.25) |
            lfc>=0.25
        )
    ]
    tmp.data[, gen.deg.class:=ifelse(test=gen.deg.class, yes='DEG', no='Other')]
    tmp.data[, deg.class:=gen.deg.class]
    tmp.data[gen.deg.class=='DEG', deg.class:=ifelse(
        test=lfc<=(-0.25),
        yes='Negative', no='Positive'
    )]
    # LFC bound
    lfc.bound <- tmp.data[abs(lfc)>1, floor(quantile(x=abs(lfc), probs=0.99))]
    if(is.na(lfc.bound)) lfc.bound <- tmp.data[, round(x=max(abs(lfc)), digits=2)]
    tmp.data[lfc<=(-lfc.bound), lfc:=-lfc.bound]
    tmp.data[lfc>=(lfc.bound), lfc:=lfc.bound]
    # P-value bound
    tmp.data[, log.p:=-log10(p.adj)]
    tmp.data[is.infinite(log.p), log.p:=tmp.data[!is.infinite(log.p), max(log.p)]]
    if(test.id %in% names(custom.p.bounds)){
        p.bound <- custom.p.bounds[test.id]
    }else{
        p.bound <- tmp.data[, floor(quantile(x=log.p, probs=0.99))]
    }
    tmp.data[log.p>=p.bound, log.p:=p.bound]
    # DEG annotations
    setorderv(x=tmp.data, cols=c('p.adj', 'lfc'), order=c(1, -1))
    gen.genes <- tmp.data[1:16, gene.id]
    if(test.id %in% names(custom.gene.labs)){
        custom.genes <- unique(c(gen.genes, custom.gene.labs[[test.id]]))
    }else{
        custom.genes <- gen.genes
    }
    custom.genes <- data.table(gene.id=custom.genes)
    custom.genes <- merge(x=custom.genes, y=tmp.data, by='gene.id')
    # Plot.
    y.bound <- floor(log10(p.bound));
    if(y.bound==0) y.bound <- 0.5
    size.vals <- c(`DEG`=3, `Other`=0.5)
    tmp.ggplot <- ggplot(
            data=tmp.data,
            aes(x=lfc, y=log.p, fill=deg.class, size=gen.deg.class)
        ) +
        geom_point(shape=21, color='black', stroke=0.5) +
        geom_hline(yintercept=-log10(0.05), linewidth=1, linetype='dashed', color='#D6D7D7') +
        geom_vline(xintercept=-0.25, linewidth=1, linetype='dashed', color='#D6D7D7') +
        geom_vline(xintercept=0.25, linewidth=1, linetype='dashed', color='#D6D7D7') +
        scale_x_continuous(limits=c(-lfc.bound, lfc.bound)) +
        scale_y_continuous(
            expand=expansion(add=c(
                0, y.bound
            ))
        ) +
        scale_size_manual(values=size.vals) +
        scale_fill_manual(values=deg.class.cols) +
        labs(x='LFC', y='-log10(adj. P-value)') +
        theme(legend.position='none')
    tmp.lab <- test.id
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.dea.path, file.name=tmp.lab, type='pdf', blank.comp=blank.complement.3, do.legend=FALSE,
        repel.data=custom.genes, repel.label='gene.id'
    ) 
    # DEG class summary
    tmp.data <- tmp.data[, .(gene.count=.N), by=deg.class]
    tmp.file.name <- paste0(fig.dea.path, '/', tmp.lab, '_DEGSumm.csv')
    fwrite(file=tmp.file.name, x=tmp.data)
}