############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    -----   Preliminary figure stack    -----    ############
###########    ---------- HIPC Lung-1 project  -----------    ###########
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############

# ---> About the script.
# Version: 0
# Subversion: 1
# Best module: R/4.2.2


############    -----------------------------------------    ############
### -------------------------- Description -------------------------- ###
############    -----------------------------------------    ############

# Script to get all bioinformatics-based panel figures for any of the final main and supplementary figures of the HIPC Lung-1 project's paper.


############    -----------------------------------------    ############
### --------------------------- Libraries --------------------------- ###
############    -----------------------------------------    ############
library(Seurat)
library(data.table)
library(tidyr)
library(stringr)
library(fgsea)
library(gtools)
# Plotting tools.
library(grid)
library(ggplot2)
library(ggbreak)
library(ggpubr)
library(pheatmap)
library(VennDiagram)
library(circlize)
# library(ggnewscale)
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

coll.version <- 0.3
gen.data.path <- '/path/to/paper_items'
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
main.prj <- 'HIPC'
this.prj <- 'HIPC-Lung'
this.figure <- 'figure_stack_v6'
# @ Seed.
set.seed(seed=1)
# @ Dataset labels.
obj.extended.names <- c(
    c40.lg.cd4.d40='C40-CD4',
    c40.lg.cd8.d40='C40-CD8',
    c40.lg.cd8.trm='C40-CD8-TRM',
    c40.lg.cd8.tem='C40-CD8-TEM',
    c40.ln.cd4.d08='C40-LN-CD4',
    c40.ln.cd8.d08='C40-LN-CD8'
)
main.extended.names <- c(
    c40.lg.cd4.d40='C40-CD4',
    c40.lg.cd8.d40='C40-CD8'
)
gen.cell.types <- c(
    c40.lg.cd4.d40='CD4',
    c40.lg.cd8.d40='CD8',
    c40.lg.cd8.trm='CD8',
    c40.lg.cd8.tem='CD8',
    c40.ln.cd4.d08='CD4',
    c40.ln.cd8.d08='CD8'
)
gen.subsets <- list(
    `trm.tag`='TRM',
    `tcm.tag`='TCM',
    `tem.tag`=c('TEM', 'CD4-CTL')
)
# @ Cluster labels for each dataset.
clust.labs <- c(
    c40.lg.cd4.d40='RNA_snn_res.0.2',
    c40.lg.cd8.d40='RNA_snn_res.0.2',
    c40.lg.cd8.trm='RNA_snn_res.0.2',
    c40.lg.cd8.tem='RNA_snn_res.0.2',
    c40.ln.cd4.d08='RNA_snn_res.0.3',
    c40.ln.cd8.d08='RNA_snn_res.0.3'
)
# @ Aggregation set correspondence.
aggr.set.corr <- c(
    c40.lg.cd4.d40='c40.lg', # C40 stands for "Cohort of 40 donors"
    c40.lg.cd8.d40='c40.lg',
    c40.lg.cd8.trm='c40.lg',
    c40.lg.cd8.tem='c40.lg',
    c40.ln.cd4.d08='c40.ln',
    c40.ln.cd8.d08='c40.ln'
)
aggr.set.labs <- c(
    `c40.lg`='C40-LG',
    `c40.ln`='C40-LN'
)
aggr.set.main <- list(
    `c40.lg`=names(main.extended.names),
    `c40.ln`=c('c40.ln.cd4.d08', 'c40.ln.cd8.d08')
)
aggr.simp.labs <- c(
    `c40.lg`='LG',
    `c40.ln`='LN'
)
# @ Ag groups.
ag.groups <- list(
    `Broad-viruses`=c('CMV', 'EBV'),
    `Respiratory-pathogens`=c('IAV', 'SARS-CoV-2', 'MPV', 'PIV', 'RSV', 'B-pertussis-Vax', 'B-pertussis-Rest', 'Aspergillus'),
    `Respiratory-viruses`=c('IAV', 'SARS-CoV-2', 'MPV', 'PIV', 'RSV')
)
# @ TCR UCMs
tcr.ucms <- c(
    'GIANA', 'TCRdist3', 'iSMART',
    'clusTCR', 'GLIPH2'
)
# @ For tag analysis.
# min.no.cells.per.group <- 50 # @ Internal parameter for tag specific analysis (part of figure ?).
# colors.table <- NULL
# ---> Cluster identities
clusts.defs <- list(
  `c40.lg.cd4.d40`=c(
    `0`='TRM',
    `1`='TCM',
    `2`='CD4-CTL',
    `3`='TREG',    
    `4`='TFH',
    `5`='Prolif'
  ),
  `c40.lg.cd8.d40`=c(
    `0`='TRM',
    `1`='TEM',
    `2`='GZMKhi',
    `3`='TCM',
    `4`='NKG2Cpos',
    `5`='MAIT',
    `6`='Prolif'
  ),
  `c40.lg.cd8.trm`=c(
    `0`='BATFhi',
    `1`='IL7Rhi'
  ),
  `c40.lg.cd8.tem`=c(
    `0`='GZMKpos',
    `1`='ZNF683pos',
    `2`='NKG2Cpos',
    `3`='TEM 3',
    `4`='TEM 4'
  ),
  `c40.ln.cd4.d08`=c(
    `0`='TN',
    `1`='TI', # Intermediate
    `2`='IFNR',
    `3`='TREG',
    `4`='TCM',
    `5`='TFH', # (to contain 7)
    `6`='TRM'
  ),
  `c40.ln.cd8.d08`=c(
    `0`='GZMKhi',
    `1`='IFNR', # (to contain 6)
    `2`='TI',
    `3`='TRM',
    `4`='TCM',
    `5`='TN',
    # `6`='IFNR-2', 
    `7`='MAIT'
    # `8`='TEM' # (excluded)
  )
)
# ---> General aesthetic definitions.
# @ General distinction between original and bounded values.
tmp.shapes <- c('ori'=1, `bounded`=16)
# @ General dot color.
gen.dot.col <- '#929292'
# @ General dot link color and width.
gen.link.col <- '#BFBFBF'
gen.link.width <- 0.7
# ---> Specific color definitions.
# @ Color per cell type.
cell.type.cols <- c(
    CD4='#00BFFF',
    CD8='#EE82EE'
)
# @ Color per MHC restriction.
mhc.rest.cols <- c(
    `Class-II`='#00BFFF',
    `Class-I`='#EE82EE'
)
# @ General phenotypes.
gen.subset.cols <- c(
    `TRM`='#7469E0',
    `TCM`='#779CC4',
    `TEM`='#D86F75',
    `IFNR`='#FFD700',
    `IL2pos`='#754600'
)
# @ Color per cluster per dataset.
# See definitions at: https://docs.google.com/spreadsheets/d/1_lNI4M7P-LBuwfQWEgag4pHaycMZMNcmzGv_k_e7x5U/edit#gid=0
clusts.cols <- list(
  `c40.lg.cd4.d40`=c(
    `0`='#5D72E0', # TRM
    `1`='#46B2AF', # TCM
    `2`='#D4824B', # CD4-CTL
    `3`='#5FAF6A', # TREG
    `4`='#9679D1', # TFH
    `5`='#E8E8E8' # Prolif
  ),
  `c40.lg.cd8.d40`=c(
    `0`='#8B5FE0', # TRM
    `1`='#DC5C9F', # TEM
    `2`='#4E9FD1', # GZMKhi
    `3`='#A985D8', # TCM
    `4`='#63B76C', # NKG2Cpos
    `5`='#6D6D6A', # MAIT
    `6`='#E8E8E8' # Prolif
  ),
  `c40.lg.cd8.trm`=c(
    `0`='#9C80D3', # BATFhi
    `1`='#F9E6A6' # IL7Rhi
  ),
  `c40.lg.cd8.tem`=c(
    `0`='#FFD1DC', # GZMKpos
    `1`='#BF0DD9', # ZNF683pos
    `2`='#9C6376', # NKG2Cpos
    `3`='#676765', # TEM 3
    `4`='#9A9A98' # TEM 4
  ),
  `c40.ln.cd4.d08`=c(
    `0`='#A4D3C4', # TN
    `1`='#D29BA7', # Intermediate
    `2`='#F5C48D', # THIFNR
    `3`='#5FAF6A', # TREG
    `4`='#46B2AF', # TCM
    `5`='#9679D1', # TFH (to contain 7)
    `6`='#5D72E0' # TRM
  ),
  `c40.ln.cd8.d08`=c(
    `0`='#4E9FD1', # GZMKhi
    `1`='#F5C48D', # IFNR (to contain 6)
    `2`='#88C78E', # Intermediate
    `3`='#8B5FE0', # TRM
    `4`='#A985D8', # TCM
    `5`='#E3D679', # TN
    `7`='#6D6D6A' # MAIT
    # `8`='#DC5C9F', # TEM (removed)
  )
)
# @ Anatomical site colors.
an.site.cols <- c(
    `LG`='#A7C7E7',
    `LN`='#F8C3A6'
)
# @ Peptide pool colors.
pp.cols <- c(
    `CMV`='#FF888F',
    `EBV`='#FFA828',
    `SARS-CoV-2`='#DCB9FB',
    `IAV`='#5CB8FE',
    `PIV`='#FEEF74',
    `RSV`='#81E9DC',
    `MPV`='#77DD77',
    `B-pertussis-Vax`='#C3A9A2',
    `B-pertussis-Rest`='#DBAD9C',
    `Aspergillus`='#C23B22',
    `Alternaria`='#319272',
    `Cross-reactive`='#CCCC00',
    `HAV`='#e69500',
    `HBV`='#cc8400',
    `HCV`='#b37400',
    `YFV`='#C1779B',
    `HIV`='#720000',
    `Homo sapiens`='#779BC1',
    `Multiple`='#808080'
)
# @ Specificity assignment approach colors
app.cols <- c(
    `Integrative`='#F9BB80',
    `Experimental`='#CBCBCB'
)
# @ TCR UCMs.
ucm.cols <- viridisLite::viridis(n=length(tcr.ucms), option='plasma')
names(ucm.cols) <- tcr.ucms
# @ Related to subset bias analysis of clonotypes.
bias.1.cols <- c(
    `Plastic`='#E07B91',
    `Singleton`='#A7C7E7',
    `Biased to same cluster`='#A9D5C7',
    `Biased to diff. cluster`='#F2E2B3'
)
# @ Other terms' colors.
signatures.col.scale <- c('#ffffff', '#ffffe0', '#ffffad', '#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
# @ Others.
rel.spcs <- c(
    'CMV', 'EBV',
    'SARS-CoV-2', 'IAV', 'PIV', 'RSV', 'MPV',
    'B-pertussis-Vax', 'B-pertussis-Rest',
    'Aspergillus'
)
subset.n <- 50000
# Limit of detection
limit.of.detect <- 10
# Highly expanded clonoype, absolute freq. threshold.
hec.freq.thold <- 10
# ---> Path definitions.
# General data and reports paths.
prj.gen.path <- paste0('/path/to/user/', main.prj, '/paper_developments/', this.prj)
gen.data.path <- paste0(prj.gen.path, '/paper_items')
final.figs.path <- paste0(prj.gen.path, '/final_figures')
this.fig.path <- paste0(final.figs.path, '/', this.figure)
if(!dir.exists(this.fig.path)) dir.create(this.fig.path)
reports.path <- paste0(this.fig.path, '/', this.figure, '_panels')
# TCR UCM systematic comparison results.
ucm.comp.path <- '/path/to/user/R24/paper_developments/DICE-VirusSpec/exploratory_analyses/clust_by_sim_best_opts/basic_pps_v6/benchmarking/option_benchmarking/EpRef-2'
# DEA results.
dea.res.paths <- paste0(gen.data.path, '/dea_results_', main.extended.names)
names(dea.res.paths) <- names(main.extended.names)
# ---> File definitions.
# HIPC cohort seurat objects.
objs.files.list <- paste0(gen.data.path, '/', 'SeuratObj_', obj.extended.names, '.RDS')
names(objs.files.list) <- names(obj.extended.names)
# Ag specific TCR references
#   Class I-restricted TCR reference.
ref.1.meta.file <- paste0(gen.data.path, '/AgSpcTCR_Ref-1.csv')
#   Class II-restricted TCR reference.
# DICE-TCR cohort TCR data.
ref.2.meta.file <- paste0(gen.data.path, '/AgSpcTCR_Ref-2.csv')
# Proliferation assay data.
prol.assay.meta.file <- paste0(gen.data.path, '/ProlifAssay.csv')
# Prediction results.
#   Class I-restricted TCRs.
pred.ref.1.files <- c(
    `c40.lg`=paste0(gen.data.path, '/IntersectedPredictionDetails_Ref-1.csv'),
    `c40.ln`=paste0(gen.data.path, '/IntersectedPredictionDetails_LN_Ref-1.csv')
)
#   Class II-restricted TCRs.
pred.ref.2.files <- c(
    `c40.lg`=paste0(gen.data.path, '/IntersectedPredictionDetails_Ref-2.csv'),
    `c40.ln`=paste0(gen.data.path, '/IntersectedPredictionDetails_LN_Ref-2.csv')
)
# Barcode-based prediction results (includes selected background cells).
#   Class I-restricted TCRs.
pred.bc.ref.1.file <- c(paste0(gen.data.path, '/BarcodePredictions_Ref-1.csv'))
#   Class II-restricted TCRs.
pred.bc.ref.2.file <- paste0(gen.data.path, '/BarcodePredictions_Ref-2.csv')
# TCR UCM systematic comparison results.
ucm.comp.files <- paste0(ucm.comp.path, '/', tcr.ucms, '/BestOptResults.RDS')
names(ucm.comp.files) <- tcr.ucms
# TCR validation results.
val.res.file <- paste0(gen.data.path, '/HIPC_L40_TCRValidation_2025-09-01.csv')
# Gene sets' files (for module scoring & FGSE analyses).
module.feats.dir <- paste0(gen.data.path, '/module_signatures')
msigdb.file <- paste0(gen.data.path, '/MSigDB_HALLMARK_2023-03-31.csv')
# Specific files.
#     Donor nomenclature correspondence.
donor.nom.corr.file <- paste0(gen.data.path, '/NumenclatureCorrespondence.csv')
#     Donor-specific cell counts before and after QC filtering.
qc.inf.1.file <- paste0(gen.data.path, '/QCInf_CellCount_Donor.csv')
#     Previously reported TCR motifs
pub.motif.file <- paste0(gen.data.path, '/PublicMotifInfo.0.1.csv')
#     Serology data for CMV and EBV
ser.ce.data.file <- paste0(gen.data.path, '/SerologyData_CMV-EBV.csv')
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
ref.1.meta <- fread(file=ref.1.meta.file)
ref.2.meta <- fread(file=ref.2.meta.file)
prol.assay.meta <- fread(file=prol.assay.meta.file)
pred.ref.1.l <- lapply(X=pred.ref.1.files, FUN=fread)
pred.ref.2.l <- lapply(X=pred.ref.2.files, FUN=fread)
pred.bc.ref.1 <- fread(file=pred.bc.ref.1.file)
pred.bc.ref.2 <- fread(file=pred.bc.ref.2.file)

# ---> TCR UCM systematic comparison results.
ucm.comp.res <- lapply(X=ucm.comp.files, FUN=readRDS)

# --> TCR validation results.
val.res <- fread(file=val.res.file)

# ---> Specific files.
#     Donor nomenclature correspondence.
donor.nom.corr <- fread(file=donor.nom.corr.file)
donor.nom.corr <- unique(donor.nom.corr[, .(uk.donor.id, donor.id.tag=final.donor.id.tag)])
donor.nom.corr <- donor.nom.corr[!is.na(uk.donor.id)]
#     Donor-specific cell counts before and after QC filtering.
qc.inf.1 <- fread(file=qc.inf.1.file)
#     Previously reported TCR motifs
pub.motif <- fread(file=pub.motif.file, na=c(NA, 'NA', ''))
#     Serology data for CMV and EBV
ser.ce.data <- fread(file=ser.ce.data.file)

# ---> Gene sets' files (for FGSE analyses).
module.feats.files <- list.files(path=paste0(module.feats.dir, '/base_signatures'), recursive=FALSE, include.dirs=FALSE, pattern='signature.csv$', full.names=TRUE)
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
names(modules.feats) <- list.files(path=paste0(module.feats.dir, '/base_signatures'), recursive=FALSE, include.dirs=FALSE, pattern='signature.csv$', full.names=FALSE)
names(modules.feats) <- str_replace(string=names(modules.feats), pattern='_signature.csv$', replacement='')
names(modules.feats) <- str_replace_all(string=names(modules.feats), pattern='_', replacement='.')
modules.names <- names(modules.feats)
names(modules.names) <- str_replace_all(string=modules.names, pattern='\\.', replacement='_')

# ---> Gene sets from MSigDB.
msigdb.data <- fread(file=msigdb.file)

# ---> DEA results.

# @ Index specificity vs background strategy.
# Define summarized results of index-vs-bg strategy.
ivb.res.paths <- paste0(dea.res.paths, '/idx-vs-bg_comps')
names(ivb.res.paths) <- names(dea.res.paths)
ivb.res <- lapply(X=ivb.res.paths, FUN=function(x){
    # Define output files.
    x <- paste0(list.dirs(path=x), '/Summary_AllGenes_AllMetrics.csv')
    x <- x[
        file.exists(x) &
        str_detect(string=x, pattern='index-ag-vs-bg_clusters-\\d+')
    ]
    if(length(x)<1) stop('Unexpected error!\n')
    # Define T cell subsets
    tmp.vals <- str_extract(string=x, pattern='index-ag-vs-bg_clusters-\\d+')
    names(x) <- str_extract(string=tmp.vals, pattern='\\d+$')
    # Load results.
    x <- lapply(X=x, FUN=fread)
    x <- rbindlist(l=x, use.names=TRUE, idcol='cluster')
    return(x)
})

# @ Subpopulation comparisons among pathogens.
# Define and load individual results.
apc.res.paths <- paste0(dea.res.paths, '/among-pats_comps')
names(apc.res.paths) <- names(dea.res.paths)
apc.res <- lapply(X=apc.res.paths, FUN=function(x){
    x <- list.dirs(path=x, full.names=TRUE, recursive=FALSE)
    x <- x[str_detect(string=basename(path=x), pattern='label-pop-\\d+')]
    if(length(x)<1){
        stop('Failed to define summarized results directory.\n')
    }
    x <- paste0(x, '/DEA/_multiple_comparisons/MergedResults_DEGsInfo_MultipleComparisons.csv')
    names(x) <- str_replace(
        string=str_extract(string=x, pattern='pop-\\d+-[^_]+'),
        pattern='pop-', replacement=''
    )
    tmp.check <- all(file.exists(x))
    if(!tmp.check) stop('Unexpected error.\n')
    x <- lapply(X=x, FUN=fread)
    x <- rbindlist(l=x, use.names=TRUE, idcol='pop.tag', fill=TRUE)
    x[, cluster.tag:=str_extract(string=pop.tag, pattern='^\\d+')]
    tmp.cols <- c('cluster.tag', setdiff(x=colnames(x), y='cluster.tag'))
    x <- x[, ..tmp.cols]
    return(x)
})
names(apc.res) <- names(apc.res.paths)


############    -----------------------------------------    ############
### ---------------------- Data preprocessing ----------------------- ###
############    -----------------------------------------    ############

# ---> Seurat objects, basic processing
# Combination instructions.
comb.insts <- list(
    `c40.ln.cd4.d08`=list(
        `5`=c('5', '7')
    ),
    `c40.ln.cd8.d08`=list(
        `1`=c('1', '6')
    )
)
# Process per seurat object
tmp.subset <- c(gen.subsets, list(`gzmk.tag`=c('GZMKhi')))
for(tmp.lab in names(srt.objs.list)){
    # Combine clusters (if necessary)
    if(tmp.lab %in% names(comb.insts)){
        for(tmp.comb in comb.insts[[tmp.lab]]){
            tmp.new.val <- tmp.comb[1]
            tmp.data <- as.character(srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]])
            tmp.data[tmp.data %in% tmp.comb] <- tmp.new.val
            srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]] <- tmp.data
        }
    }
    # Keep only relevant clusters
    clusts.to.keep <- names(clusts.cols[[tmp.lab]])
    tmp.data <- Cells(srt.objs.list[[tmp.lab]])[
        as.character(srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]]) %in% clusts.to.keep
    ]
    srt.objs.list[[tmp.lab]] <- subset(x=srt.objs.list[[tmp.lab]], cells=tmp.data)
    # Set populations as factors.
    srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]] <- factor(
        x=as.character(srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]]),
        levels=names(clusts.defs[[tmp.lab]])
    )
    # Add population tag for common populations.
    for(tmp.tag in names(tmp.subset)){
        tmp.clusts <- names(clusts.defs[[tmp.lab]])[
            clusts.defs[[tmp.lab]] %in% tmp.subset[[tmp.tag]]
        ]
        srt.objs.list[[tmp.lab]]@meta.data[, tmp.tag] <- factor(
            x=rep(x='Rest', times=length(Cells(srt.objs.list[[tmp.lab]]))),
            levels=c(tmp.subset[[tmp.tag]][1], 'Rest')
        )
        tmp.idxs <- srt.objs.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]] %in% tmp.clusts
        srt.objs.list[[tmp.lab]]@meta.data[tmp.idxs, tmp.tag] <- tmp.subset[[tmp.tag]][1]
    }
}
# @ Custom changes.
# Wrong donor ID label in LN datasets.
tmp.vals <- c('c40.ln.cd4.d08', 'c40.ln.cd8.d08')
for(tmp.lab in tmp.vals){
    to.check <- !is.na(srt.objs.list[[tmp.lab]]@meta.data[, 'donor.id.tag']) & srt.objs.list[[tmp.lab]]@meta.data[, 'donor.id.tag']=='L024'
    srt.objs.list[[tmp.lab]]@meta.data[to.check, 'donor.id.tag'] <- 'L027'
}

# ---> Retrieve general GEx metadata.
meta.data.l <- lapply(X=aggr.set.main, FUN=function(main.sets){
    meta.data <- lapply(X=main.sets, FUN=function(data.set){
        seurat.obj <- srt.objs.list[[data.set]]
        tmp.data <- as.data.table(seurat.obj@meta.data)
        tmp.data$barcode <- row.names(seurat.obj@meta.data)
        tmp.tag <- clust.labs[data.set]
        tmp.data[, clusters.tag:=as.character(get(tmp.tag))]
        return(tmp.data)
    })
    names(meta.data) <- main.sets
    meta.data <- rbindlist(l=meta.data, use.names=TRUE, idcol='data.set')
    meta.data[, gex.lin.tag:=toupper(str_extract(string=data.set, pattern='cd[48]'))]
    # Set factors.
    tmp.lvls <- c('CD4', 'CD8')
    tmp.vals <- factor(x=meta.data[, as.character(gex.lin.tag)], levels=tmp.lvls)
    set(x=meta.data, j='gex.lin.tag', value=tmp.vals)
    return(meta.data)
})

# --> TCR validation results.
# Keep data for positive results only.
val.res <- val.res[`Used in HIPC1`=='Used']
# Keep only significant columns only w/ proper names
tmp.cols <- c(
    `Dataset`='val.data.set',
    `Clonotype validation ID`='val.clonotype.tag',
    `Cell Type`='gex.lin.tag',
    `AA CDR3B`='cdr3b.aa.seq',
    `AA CDR3A`='cdr3a.aa.seq',
    `TRBV gene`='trb.v',
    `TRBJ gene`='trb.j',
    `TRAV gene`='tra.v',
    `TRAJ gene`='tra.j',
    `Putative specificity (if known)`='val.specificity.tag'
)
tmp.vals <- names(tmp.cols)
val.res <- val.res[, ..tmp.vals]
colnames(val.res) <- tmp.cols[colnames(val.res)]
# Custom changes.
val.res[, unique(val.data.set)]
val.res[val.data.set=='HIPC L40 / Ref.-matched TCRs', val.data.set:='HIPC-040']
tmp.vals <- c(
    `SARS`='SARS-CoV-2',
    `FLU`='IAV',
    `RSPV`='PIV',
    `MPV`='MPV',
    `EBV`='EBV',
    `RSV`='RSV',
    `B. pertussis (vax)`='B-pertussis-Vax',
    `Aspergillus`='Aspergillus',
    `CMV/Alternaria`='CMV',
    `CMV`='CMV',
    `SARS-CoV-2`='SARS-CoV-2',
    `IAV`='IAV'
)
tmp.check <- val.res[, all(val.specificity.tag %in% names(tmp.vals))]
if(!tmp.check) stop('Unexpected error!\n')
val.res[, val.specificity.tag:=tmp.vals[val.specificity.tag]]
tmp.check <- val.res[, all(val.specificity.tag%in%names(pp.cols))]
if(!tmp.check) stop('Unexpected error!\n')
# Separate multi-chain clonotypes into separate rows.
val.res <- as.data.table(separate_rows(data=val.res, cdr3b.aa.seq, trb.v, trb.j, sep=';'))
val.res <- as.data.table(separate_rows(data=val.res, cdr3a.aa.seq, tra.v, tra.j, sep=';'))

# ---> TCR data
# @ TCR references. 
#   Class I-restricted TCR reference.
# Set factors.
tmp.lvls <- c(
  'CMV',
  'EBV',
  'SARS-CoV-2',
  'IAV',
  'HAV',
  'HBV',
  'HCV',
  'Homo sapiens',
  'Cross-reactive'
)
ref.1.meta <- ref.1.meta[specificity.class.tag %in% tmp.lvls]
tmp.vals <- factor(x=ref.1.meta[, specificity.class.tag], levels=tmp.lvls)
set(x=ref.1.meta, j='specificity.class.tag', value=tmp.vals)
#   Class II-restricted TCR reference.
# NOTE: THIS IS A TEMPORARY FIX. THE HIV DATA SHOULD BE EXCLUDED UPSTREAM IN THE PROCESS, RIGHT BEFORE ASSIGNING PREFERENTIAL REACTIVITIES (IDEALLY, EVEN DURING THE AGGREGATION PROCESS).
ref.2.meta <- ref.2.meta[!(!is.na(specificity.class.tag) & specificity.class.tag=='HIV')]
# Set factors.
#   Peptide pool
tmp.lvls <- names(pp.cols)[names(pp.cols) %in% ref.2.meta[, unique(specificity.class.tag)]]
tmp.vals <- factor(x=ref.2.meta[, specificity.class.tag], levels=tmp.lvls)
set(x=ref.2.meta, j='specificity.class.tag', value=tmp.vals)
#   Clonotype status.
tmp.lvls <- c('Singleton', 'Ambiguous', 'Unambiguous')
tmp.vals <- factor(x=ref.2.meta[, gen.clonotype.status], levels=tmp.lvls)
set(x=ref.2.meta, j='gen.clonotype.status', value=tmp.vals)

# @ Proliferation assay data.
# Keep only peptide pools of interest.
# NOTE: CUSTOM CHANGE. MAY MODIFY IN FUTURE VERSIONS OF THE SCRIPT/FIGURES.
custom.pps <- c(
    'CMV', 'EBV',
    'IAV', 'SARS-CoV-2', 'MPV', 'PIV', 'RSV',
    'B-pertussis-Vax', 'B-pertussis-Rest',
    'Aspergillus'
)
prol.assay.meta <- prol.assay.meta[
    !(!is.na(specificity.class.tag) &
    specificity.class.tag=='RV') &
    specificity.class.tag %in% custom.pps
]
# Set factors.
#   Peptide pool
tmp.lvls <- names(pp.cols)[names(pp.cols) %in% prol.assay.meta[, unique(specificity.class.tag)]]
tmp.vals <- factor(x=prol.assay.meta[, specificity.class.tag], levels=tmp.lvls)
set(x=prol.assay.meta, j='specificity.class.tag', value=tmp.vals)
#   Clonotype status.
tmp.lvls <- c('Singleton', 'Ambiguous', 'Unambiguous')
tmp.vals <- factor(x=prol.assay.meta[, gen.clonotype.status], levels=tmp.lvls)
set(x=prol.assay.meta, j='gen.clonotype.status', value=tmp.vals)

# @ Prediction results.
# Custom changes: Wrong donor ID label in LN datasets.
pred.ref.1.l[['c40.ln']][donor.id.tag=='L024', donor.id.tag:='L027']
pred.ref.2.l[['c40.ln']][donor.id.tag=='L024', donor.id.tag:='L027']
# General function to set factors.
set.pred.res.factors <- function(pred.ref){
    #   Donor ID
    tmp.lvls <- c(paste0('L00', 1:9), paste0('L0', 10:40))
    tmp.vals <- factor(x=pred.ref[, donor.id.tag], levels=tmp.lvls)
    set(x=pred.ref, j='donor.id.tag', value=tmp.vals)
    # NOTE: CUSTOM CHANGE. MAY MODIFY IN FUTURE VERSIONS OF THE SCRIPT/FIGURES.
    # @ Keep only predictions for custom list of antigens.
    custom.pps <- c(
        'CMV', 'EBV',
        'IAV', 'SARS-CoV-2', 'MPV', 'PIV', 'RSV',
        'B-pertussis-Vax', 'B-pertussis-Rest',
        'Aspergillus'
    )
    pred.ref[
        !is.na(consensus.pred) & !consensus.pred %in% custom.pps,
        `:=`(
            consensus.pred=NA,
            consensus.approach.class=NA
        )
    ]
    #   Consensus approach
    tmp.lvls <- c('Match', 'Experimental', 'Clustering')
    tmp.vals <- factor(x=pred.ref[, consensus.approach.class], levels=tmp.lvls)
    set(x=pred.ref, j='consensus.approach.class', value=tmp.vals)
    #   Consensus prediction.
    tmp.lvls <- names(pp.cols)[names(pp.cols) %in% pred.ref[, unique(consensus.pred)]]
    tmp.vals <- factor(x=pred.ref[, consensus.pred], levels=tmp.lvls)
    set(x=pred.ref, j='consensus.pred', value=tmp.vals)
    #   Experimental prediction.
    tmp.lvls <- names(pp.cols)[names(pp.cols) %in% pred.ref[, unique(exp.pred)]]
    tmp.vals <- factor(x=pred.ref[, exp.pred], levels=tmp.lvls)
    set(x=pred.ref, j='exp.pred', value=tmp.vals)
    #   Experimental match chain class
    tmp.lvls <- c('BC', 'Beta', 'Alpha')
    tmp.vals <- factor(x=pred.ref[, exp.chain.class], levels=tmp.lvls)
    set(x=pred.ref, j='exp.chain.class', value=tmp.vals)
    #   Reference match chain class
    tmp.lvls <- c('BC', 'IA', 'IB', 'Beta', 'Alpha')
    tmp.vals <- factor(x=pred.ref[, match.chain.class], levels=tmp.lvls)
    set(x=pred.ref, j='match.chain.class', value=tmp.vals)
    #   Reference match match class
    tmp.lvls <- c('Perfect', 'Imperfect', 'Unambiguous', 'Ambiguous')
    tmp.vals <- factor(x=pred.ref[, match.match.class], levels=tmp.lvls)
    set(x=pred.ref, j='match.match.class', value=tmp.vals)
    #   Return
    return(pred.ref)
}
#   Class I-restricted TCRs.
# Set factors.
pred.ref.1.l <- lapply(X=pred.ref.1.l, FUN=set.pred.res.factors)
#   Class II-restricted TCRs.
# Set factors.
pred.ref.2.l <- lapply(X=pred.ref.2.l, FUN=set.pred.res.factors)

# ---> TCR data addition to general metadata, HIPC cohort.
# @ Add basic prediction results. Function to be used for either cell type.
add.pred.ann <- function(qry.data, ref.data, ref.subset){
    #       For reference 1.
    tmp.data <- ref.data[
        !is.na(consensus.pred),
        .(clonotype.tag=qry.clone.id, donor.id.tag, cdr3b.aa.seq, cdr3a.aa.seq, trb.v, trb.j, tra.v, tra.j, ag.spc.tag=consensus.pred)
    ]
    merge.cols <- c('clonotype.tag', 'cdr3b.aa.seq', 'cdr3a.aa.seq', 'trb.v', 'trb.j', 'tra.v', 'tra.j')
    #           Sanity checks.
    tmp.data.1 <- unique(tmp.data[clonotype.tag %in% qry.data[, clonotype.tag], ..merge.cols])
    tmp.data.2 <- unique(qry.data[clonotype.tag %in% tmp.data.1[, clonotype.tag], ..merge.cols])
    to.check.1 <- tmp.data.2[, .N]
    to.check.2 <- merge(x=tmp.data.1, y=tmp.data.2, by=merge.cols)
    tmp.check <- to.check.1==to.check.2[, .N]
    if(!tmp.check) stop('Unexpected error while adding prediction annotations to general HIPC TCR data (1).\n')
    #           Proceed to merme.
    merge.cols <- c('gex.lin.tag', 'donor.id.tag', merge.cols)
    tmp.data[, gex.lin.tag:=ref.subset]
    tmp.data <- merge(
        x=qry.data, y=tmp.data,
        by=merge.cols,
        all.x=TRUE, all.y=FALSE,
        sort=FALSE
    )
    tmp.check <- tmp.data[, .N]==qry.data[, .N]
    if(!tmp.check) stop('Unexpected error while adding prediction annotations to general HIPC TCR data (2).\n')
    return(tmp.data)
}
# @ Add consensus predictions to GEx metadata and their corresponding seurat objects.
for(data.set in names(aggr.set.main)){
    # @ Add validation IDs (if necessary).
    if(data.set=='c40.lg'){
        merge.cols <- c(
            'gex.lin.tag',
            'cdr3b.aa.seq', 'trb.v', 'trb.j',
            'cdr3a.aa.seq', 'tra.v', 'tra.j'
        )
        meta.data.l[[data.set]] <- merge(
            x=meta.data.l[[data.set]], y=val.res,
            by=merge.cols,
            all.x=TRUE, all.y=FALSE,
            sort=FALSE
        )
        to.check <- val.res[!val.clonotype.tag %in% meta.data.l[[data.set]][, val.clonotype.tag], 'HIPC-040' %in% val.data.set]
        if(to.check) stop('Unexpected error.\n')
    }
    # @ Add basic prediction results from ref. 1
    meta.data.l[[data.set]] <- add.pred.ann(qry.data=meta.data.l[[data.set]], ref.data=pred.ref.1.l[[data.set]], ref.subset='CD8')
    # @ Add basic prediction results from ref. 2
    meta.data.l[[data.set]] <- add.pred.ann(qry.data=meta.data.l[[data.set]], ref.data=pred.ref.2.l[[data.set]], ref.subset='CD4')
    # @ Combine pred. results from both sources.
    tmp.check <- meta.data.l[[data.set]][!is.na(ag.spc.tag.x) & !is.na(ag.spc.tag.y), .N==0]
    if(!tmp.check) stop('Unexpected error while adding prediction annotations to general HIPC TCR data (merged).\n')
    meta.data.l[[data.set]][!is.na(ag.spc.tag.y), ag.spc.tag.x:=ag.spc.tag.y]
    meta.data.l[[data.set]][, `:=`(ag.spc.tag=ag.spc.tag.x, ag.spc.tag.y=NULL, ag.spc.tag.x=NULL)]
    tmp.check <- meta.data.l[[data.set]][!is.na(ag.spc.tag), all(ag.spc.tag %in% rel.spcs)]
    if(!tmp.check) stop('Unexpected error. All final specificities in metadata should be included  as relevant.\n')
    # @ Add broad pathogen type.
    tmp.groups <- c(
        `CMV`='Herpesvirus',
        `EBV`='Herpesvirus',
        `IAV`='Respiratory pat.',
        `SARS-CoV-2`='Respiratory pat.',
        `MPV`='Respiratory pat.',
        `PIV`='Respiratory pat.',
        `RSV`='Respiratory pat.',
        `B-pertussis-Vax`='Respiratory pat.',
        `B-pertussis-Rest`='Respiratory pat.',
        `Aspergillus`='Respiratory pat.'
    )
    meta.data.l[[data.set]][, broad.pat.type:=tmp.groups[ag.spc.tag]]
    # @ Add broad virus type.
    tmp.groups <- c(
        `CMV`='Herpesvirus',
        `EBV`='Herpesvirus',
        `IAV`='Respiratory virus',
        `SARS-CoV-2`='Respiratory virus',
        `MPV`='Respiratory virus',
        `PIV`='Respiratory virus',
        `RSV`='Respiratory virus'
    )
    meta.data.l[[data.set]][, broad.vir.type:=tmp.groups[ag.spc.tag]]
    # Set factors.
    tmp.vals <- factor(x=meta.data.l[[data.set]][, ag.spc.tag], levels=rel.spcs)
    set(x=meta.data.l[[data.set]], j='ag.spc.tag', value=tmp.vals)
    # Final column order (for my own sanity)
    tmp.cols <- c('data.set', 'gex.lin.tag', 'barcode', 'donor.id.tag')
    tmp.cols <- c(tmp.cols, setdiff(x=colnames(meta.data.l[[data.set]]), y=tmp.cols))
    meta.data.l[[data.set]] <- meta.data.l[[data.set]][, ..tmp.cols]
    # Add specificity annotations to seurat objects.
    for(t.subset in names(obj.extended.names)){
        if(aggr.set.corr[[t.subset]]!=data.set) next
        tmp.data <- as.data.frame(meta.data.l[[data.set]][, .(barcode, ag.spc.tag)])
        row.names(tmp.data) <- tmp.data$barcode; tmp.data$barcode <- NULL
        tmp.check <- all(Cells(srt.objs.list[[t.subset]]) %in% row.names(tmp.data))
        if(!tmp.check) stop(paste0('Could not find full list of barcodes from object ', t.subset, ' in general metadata file.\n'))
        srt.objs.list[[t.subset]]@meta.data[, 'ag.spc.tag'] <- tmp.data[
            Cells(srt.objs.list[[t.subset]]),
            'ag.spc.tag'
        ]
    }
}

# ---> Formal subset-bias determination
tmp.reports.path <- paste0(reports.path, '/subset_bias_det')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
# Define clusters to remove.
clusts.to.rm <- c(
    `c40.lg.cd4.d40`=c('5'),
    `c40.lg.cd8.d40`=c('6')
)
meta.data.l <- lapply(X=names(aggr.set.labs), FUN=function(data.set){
    tmp.vals <- meta.data.l[[data.set]][, unique(data.set)]
    tmp.data <- lapply(X=tmp.vals, FUN=function(t.subset){
        this.reports.path <- paste0(tmp.reports.path, '/', obj.extended.names[t.subset])
        if(!dir.exists(this.reports.path)) dir.create(this.reports.path)
        tmp.data.1 <- meta.data.l[[data.set]][gex.lin.tag==gen.cell.types[t.subset]]
        if(t.subset %in% names(clusts.to.rm)){
            bias.det.inp <- tmp.data.1[!clusters.tag %in% clusts.to.rm[[t.subset]]]
        }else{
            bias.det.inp <- tmp.data.1
        }
        tmp.data.2 <- det.group.biases(
            tcr.data=bias.det.inp,
            group.tag=clust.labs[t.subset], group.cols=clusts.cols[[t.subset]],
            group.fc.thold=2, fad.1.thold=1,
            this.reports.path=this.reports.path
        )
        tmp.cols <- setdiff(x=colnames(tmp.data.2), y=colnames(tmp.data.1))
        tmp.cols <- c('barcode', tmp.cols)
        tmp.data.2 <- tmp.data.2[, ..tmp.cols]
        tmp.data.1 <- merge(
            x=tmp.data.1, y=tmp.data.2,
            by='barcode', all=TRUE
        )
        return(tmp.data.1)
    })
    tmp.data <- rbindlist(l=tmp.data, use.names=TRUE)
    merge.cols <- intersect(x=colnames(meta.data.l[[data.set]]), y=colnames(tmp.data))
    # setdiff(x=colnames(meta.data.l[[data.set]]), y=merge.cols)
    # setdiff(x=colnames(tmp.data), y=merge.cols)
    to.check <- merge(
        x=meta.data.l[[data.set]], y=tmp.data,
        by=merge.cols, all=TRUE
    )
    tmp.check <- to.check[, .N]==meta.data.l[[data.set]][, .N]
    if(!tmp.check) stop('Unexpected error.\n')
    return(tmp.data)
})
names(meta.data.l) <- names(aggr.set.labs)

# ---> Obtain cell subset of each dataset that's equal across cell types.
subset.list <- lapply(X=srt.objs.list, FUN=function(seurat.obj){
    if(subset.n>=length(Cells(seurat.obj))){
        return(Cells(seurat.obj))
    }else{
        sample(x=Cells(seurat.obj), size=subset.n)
    }
})

# ---> Specific files.
#     Previously reported TCR motifs
# Filter out motifs w/out a valid PMID.
pub.motif <- pub.motif[pmid %like% '^\\d+$']
#     Serology data for CMV and EBV
# No meaningful comparisons can be performed for EBV since all donors are EBV+.
ser.cmv.data <- ser.ce.data[, .(
    uk.donor.id=`Participant.Participant ID`,
    cmv.status=`CMV interpretation`,
    cmv.igg=`CMV IgG AU/mL`
)]
tmp.vals <- c(
    `Evidence primary infection at some time`='Positive',
    `No evidence of infection`='Negative'
)
ser.cmv.data[, cmv.status:=tmp.vals[cmv.status]]
ser.cmv.data[cmv.igg=='not detected', cmv.igg:=NA]
ser.cmv.data[, cmv.igg:=as.numeric(cmv.igg)]
ser.cmv.data <- merge(x=donor.nom.corr, y=ser.cmv.data, by='uk.donor.id', all.x=FALSE, all.y=TRUE)
ser.cmv.data[, uk.donor.id:=NULL]
tmp.lvls <- c('Positive', 'Negative')
tmp.vals <- factor(x=ser.cmv.data$cmv.status, levels=tmp.lvls)
set(x=ser.cmv.data, j='cmv.status', value=tmp.vals)
#     Annotated metadata for datasets before QC filtering.
# names(no.qc.meta) <- str_replace(string=names(no.qc.meta), pattern='R24_Cancer_Batches-1-to-20_', replacement='')
# no.qc.meta <- no.qc.meta[str_detect(string=names(no.qc.meta), pattern='Normal-lung')]
# names(no.qc.meta) <- str_replace(string=names(no.qc.meta), pattern='_Normal-lung', replacement='')
# no.qc.meta <- rbindlist(l=no.qc.meta, use.names=TRUE, idcol='data.set')

# ---> Number of cells per immune population.
tmp.data <- lapply(X=names(main.extended.names), FUN=function(data.set){
    seurat.obj <- srt.objs.list[[data.set]]
    tmp.data <- as.data.table(seurat.obj@meta.data)
    tmp.data <- data.table(
        total.in.obj=tmp.data[, .N],
        total.ex.vivo=tmp.data[culture.cond.tag=='Ex vivo', .N],
        total.ex.vivo.with.tcr=tmp.data[culture.cond.tag=='Ex vivo' & !is.na(clonotype.tag), .N]
    )
    return(tmp.data)
})
names(tmp.data) <- names(main.extended.names)
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='data.set')
tmp.file.name <- paste0(reports.path, '/CellCounts.csv')
fwrite(x=tmp.data, file=tmp.file.name)


### -------------------------- Preflights --------------------------- ###
### ---------------- For anatomical site assessments ---------------- ###

# ---> Anatomical site donor overlap &, accordingly, custom donor order.
site.donor.ovlp <- intersect(
    meta.data.l[[1]][, unique(donor.id.tag)],
    meta.data.l[[2]][, unique(donor.id.tag)]
)
tmp.data <- meta.data.l[['c40.lg']][
    !is.na(donor.id.tag) & !is.na(clonotype.tag),
    .(cell.count=.N),
    by=donor.id.tag
]
setorderv(x=tmp.data, cols='cell.count', order=-1)
custom.donor.ord.1 <- c(
    tmp.data[donor.id.tag %in% site.donor.ovlp, donor.id.tag],
    tmp.data[!donor.id.tag %in% site.donor.ovlp, donor.id.tag]
)
# Manually set alternative custom donor orders. Based on stats set based on TRM enrichment applied on section titled: "Antigen-specific T cells: Repertoire comparisons between lineages".
#       CD4 order (based on mean).
custom.donor.ord.2 <- c(
    '15', '19', '01', '37', '25', '22', '27', '10', '14', '35',
    '39', '12', '13', '38', '24', '26', '40', '09', '02', '23',
    '04', '32', '16', '28', '30', '03', '05', '31', '34', '06',
    '36', '17', '20', '08', '11', '18', '07', '29', '21', '33'
)
custom.donor.ord.2 <- paste0('L0', custom.donor.ord.2)
#       CD8 order (based on max).
custom.donor.ord.3 <- c(
    '26', '35', '12', '10', '03', '36', '15', '14', '19', '32',
    '23', '02', '13', '24', '33', '40', '30', '17', '28', '22',
    '01', '08', '25', '27', '38', '34', '20', '05', '39', '31',
    '37', '21', '07', '06', '04', '16', '11', '29', '18', '09'
)
custom.donor.ord.3 <- paste0('L0', custom.donor.ord.3)

# ---> For tissue clonality overlap
tmp.data <- get.comp.clone.freqs(
  gen.tcr.data=meta.data.l,
  clone.group.tags=NULL,
  donor.subset=site.donor.ovlp
)
gen.tcr.data.1 <- copy(tmp.data[[1]])
gen.tcr.data.2 <- copy(tmp.data[[2]])
rm(tmp.data)

# ---> Donor metadata collection.
# Retrieve donor-specific metadata.
tmp.cols <- colnames(meta.data.l[['c40.lg']])[str_detect(string=colnames(meta.data.l[['c40.lg']]), pattern='^donor\\.')]
tmp.cols <- setdiff(x=tmp.cols, y='donor.tag')
donor.meta <- unique(meta.data.l[['c40.lg']][, ..tmp.cols])
# Select relevant variables to show.
tmp.cols <- c(
    'donor.id.tag'='donor.id.tag',
    'age'='donor.age.tag',
    'bmi'='donor.bmi.tag',
    'gender'='donor.gender.tag',
    'smoking.status'='donor.smoking.status.tag',
    'alcohol.units'='donor.alcohol.units.a.week.tag',
    'covid.shot.no'='donor.known.shot.no.tag',
    'covid.shot.types'='donor.vaccine.types.tag',
    'covid.shot.span'='donor.vacc.cov.sample.shortest.time.int.tag',
    'flu.shot.span'='donor.vacc.flu.sample.time.int.tag'
)
donor.meta <- donor.meta[, ..tmp.cols]
colnames(donor.meta) <- names(tmp.cols)

# ---> Donor metadata to be consistently displayed.
# Categorize variables.
donor.meta.hmap <- copy(donor.meta)
# @ Set donor order
tmp.vals <- factor(x=as.character(donor.meta.hmap$donor.id.tag), levels=rev(custom.donor.ord.1))
set(x=donor.meta.hmap, j='donor.id.tag', value=tmp.vals)
# @ Gender
donor.meta.lvls <- list(gender=c('Female', 'Male'))
# @ Age
# 4 categories. Chosen according to quartiles, then adjusted for ease in communication (according to 5-y bins). <66, 66-70, 71-75, >75
# donor.meta.hmap[, quantile(x=age, probs=c(0.25, 0.5, 0.75), na.rm=TRUE)]
donor.meta.hmap[,
    age:=ifelse(
        test=age<66,
        yes='<66',
        no=ifelse(
            test=age<71,
            yes='66-70',
            no=ifelse(
                test=age<76,
                yes='71-75',
                no='>75'
            )
        )
    )
]
donor.meta.lvls[['age']] <-  c('<66', '66-70', '71-75', '>75')
# @ BMI
donor.meta.hmap[,
    bmi:=ifelse(
        test=bmi<18.5,
        yes='<18.5',
        no=ifelse(
            test=bmi<=25,
            yes='18.5-25',
            no=ifelse(
                test=bmi<=30,
                yes='25-30',
                no='>30'
            )
        )
    )
]
donor.meta.lvls[['bmi']] <- c('<18.5', '18.5-25', '25-30', '>30')
# @ Smoking status
donor.meta.lvls[['smoking.status']] <- c('Never Smoked', 'Ex smoker', 'Current smoker')
# @ Alcohol units
donor.meta.hmap[,
    alcohol.units:=ifelse(
        test=alcohol.units==0,
        yes='Abstainer',
        no=ifelse(
            test=alcohol.units<=7,
            yes='Low consumption',
            no=ifelse(
                test=alcohol.units<=15,
                yes='Moderate consumption',
                no='High consumption'
            )
        )
    )
]
donor.meta.lvls[['alcohol.units']] <- c('Abstainer', 'Low consumption', 'Moderate consumption', 'High consumption')
# @ Covid vaccine types.
donor.meta.hmap[covid.shot.types=='', covid.shot.types:=NA]
donor.meta.hmap[str_detect(string=covid.shot.types, pattern=';'), covid.shot.types:='Hybrid']
donor.meta.lvls[['covid.shot.types']] <- c('Hybrid', 'Pfizer-BioNTech', 'AstraZeneca')
# @ Total number of Covid shots.
# donor.meta.hmap[, .N, by=covid.shot.no]
donor.meta.hmap[,
    covid.shot.no:=ifelse(
        test=covid.shot.no==0,
        yes='None',
        no=ifelse(
            test=covid.shot.no>2,
            no='2',
            yes='>2'
        )
    )
]
donor.meta.lvls[['covid.shot.no']] <- c('None', '2', '>2')
# @ Time span since last Covid shot.
# donor.meta.hmap[, summary(covid.shot.span)]
donor.meta.hmap[,
    covid.shot.span:=ifelse(
        test=covid.shot.span<30,
        yes='<1m',
        no=ifelse(
            test=covid.shot.span>90,
            no='>1m & <3m',
            yes='>3m'
        )
    )
]
donor.meta.lvls[['covid.shot.span']] <- c('<1m', '>1m & <3m', '>3m')
# @ Time span since last flu shot.
# donor.meta.hmap[, summary(covid.shot.span)]
donor.meta.hmap[,
    flu.shot.span:=ifelse(
        test=flu.shot.span<30,
        yes='<1m',
        no=ifelse(
            test=flu.shot.span>90,
            no='>1m & <3m',
            yes='>3m'
        )
    )
]
donor.meta.lvls[['flu.shot.span']] <- c('<1m', '>1m & <3m', '>3m')
# @ Tidy format.
donor.meta.hmap <- as.data.table(gather(data=donor.meta.hmap, key='donor.var', value='var.val', -`donor.id.tag`))
donor.meta.hmap[!is.na(var.val), fill.val:=paste(donor.var, var.val, sep='.')]
# Set levels.
tmp.lvls <- lapply(X=names(donor.meta.lvls), FUN=function(x) return(paste(x, donor.meta.lvls[[x]], sep='.')))
tmp.lvls <- unlist(tmp.lvls)
donor.meta.hmap[!is.na(fill.val), all(fill.val%in%tmp.lvls)]
donor.meta.hmap$fill.val <- factor(x=donor.meta.hmap$fill.val, levels=tmp.lvls)
donor.meta.hmap$donor.var <- factor(x=donor.meta.hmap$donor.var, levels=c(
    'flu.shot.span', 'covid.shot.span', 'covid.shot.no', 'covid.shot.types', 'alcohol.units', 'smoking.status', 'bmi', 'gender', 'age'
))
# @ Define colors.
donor.meta.cols <- c(
    'gender.Female'='#BE82BF',
    'gender.Male'='#8FD5C9',
    'age.<66'='#93F6FF',
    'age.66-70'='#91D2F3',
    'age.71-75'='#99C1FF',
    'age.>75'='#919EF3',
    'bmi.<18.5'='#8FD5C9',
    'bmi.18.5-25'='#C2E7DF',
    'bmi.25-30'='#E1C8E3',
    'bmi.>30'='#BE82BF',
    'smoking.status.Never Smoked'='#8FD5C9',
    'smoking.status.Ex smoker'='#FDFDFF',
    'smoking.status.Current smoker'='#BE82BF',
    'alcohol.units.Abstainer'='#8FD5C9',
    'alcohol.units.Low consumption'='#C2E7DF',
    'alcohol.units.Moderate consumption'='#E1C8E3',
    'alcohol.units.High consumption'='#BE82BF',
    'covid.shot.types.Hybrid'='#A89D77',
    'covid.shot.types.Pfizer-BioNTech'='#645E47',
    'covid.shot.types.AstraZeneca'='#F5E4AE',
    'covid.shot.no.None'='#F5BAE9',
    'covid.shot.no.2'='#DBA7D1',
    'covid.shot.no.>2'='#A880A0',
    'covid.shot.span.<1m'='#52A858',
    'covid.shot.span.>1m & < 3m'='#66D16E',
    'covid.shot.span.>3m'='#79F581',
    'flu.shot.span.<1m'='#52A858',
    'flu.shot.span.>1m & < 3m'='#66D16E',
    'flu.shot.span.>3m'='#79F581'
)
# Sanity check
to.check <- donor.meta.hmap[, levels(fill.val)]
to.check[!to.check %in% names(donor.meta.cols)]


############    -----------------------------------------    ############
### --------------------------- Figure on --------------------------- ###
### -------------------- Intro & TCR repertoires -------------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
fig.tcr.reps.path <- paste0(reports.path, '/figure_on_tcr_reps_intro')
if(!dir.exists(fig.tcr.reps.path)) dir.create(fig.tcr.reps.path)
# ---> Comments: Main dataset for this section, Lower-CD8 and Lower-CD4


### ------------------------- Text details -------------------------- ###

# ---> Total number of all T cells and T cells w/ available TCR data for the median donor.
for(data.set in names(aggr.set.main)){
    cat(data.set, '\n')
    tmp.data <- meta.data.l[[data.set]]
    tmp.capt <- tmp.data[,
        .(cell.count=.N),
        by=donor.id.tag
    ][, median(cell.count)]
    cat('Median GEx cell count across donors: ', tmp.capt, '\n')
    tmp.capt <- tmp.data[
        !is.na(clonotype.tag),
        .(cell.count=.N),
        by=donor.id.tag
    ][, paste0(median(cell.count), ', ', min(cell.count), '-', max(cell.count))]
    cat('Median and range of GEx & TCR cell count across donors: ', tmp.capt, '\n')
    tmp.capt <- tmp.data[
        !is.na(clonotype.tag),
        .(cell.count=.N),
        by=.(gex.lin.tag, donor.id.tag)
    ][, paste0(median(cell.count), ', ', min(cell.count), '-', max(cell.count))]
    cat('Median and range of GEx & TCR cell count (lineage-specific, separately for CD4s and CD8s) across donors: ', tmp.capt, '\n')
    tmp.capt <- tmp.data[
        !is.na(clonotype.tag),
        uniqueN(barcode)
    ]
    cat('Total GEx & TCR cell count across donors (confirmation 1):', tmp.capt, '\n')
    tmp.capt <- tmp.data[
        !is.na(clonotype.tag),
        .N,
    ]
    cat('Total GEx & TCR cell count across donors (confirmation 2):', tmp.capt, '\n')
    cat('\n\n')
}

# ---> Percentage of cells accounted for by top clonally expanded clonotypes in each compartment.
meta.data <- meta.data.l[['c40.lg']]
tmp.data.1 <- meta.data[
    !is.na(clonotype.tag),
    .(abs.freq=.N),
    by=.(gex.lin.tag, donor.id.tag, clonotype.tag)
]
setorderv(x=tmp.data.1, cols='abs.freq', order=-1)
tmp.data.1 <- tmp.data.1[,
    .(top.abs.freq=.SD[1:10, sum(abs.freq)]),
    by=.(gex.lin.tag, donor.id.tag)
]
tmp.data.1 <- tmp.data.1[!is.na(top.abs.freq)]
tmp.data.2 <- meta.data[
    !is.na(clonotype.tag),
    .(total.freq=.N),
    by=.(gex.lin.tag, donor.id.tag)
]
tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by=c('gex.lin.tag', 'donor.id.tag'))
tmp.data[, top.rel.freq:=top.abs.freq/total.freq]
tmp.data[,
    .(median(top.rel.freq)),
    by=gex.lin.tag
]


### -------------------------- Main Figure -------------------------- ###

# ---> Donor metadata.
# @ Version w/ full set of selected variables.
tmp.ggplot <- ggplot(data=donor.meta.hmap, aes(x=donor.var, y=donor.id.tag, fill=fill.val)) +
    geom_tile(color='black', linewidth=0.8, linetype=1, alpha=0.8) +
    coord_fixed() +
    scale_fill_manual(values=donor.meta.cols) +
    labs(x='', y='', fill='Value') +
    theme(axis.text.x=element_text(angle=45))
    # theme_void()
tmp.lab <- paste0('DonorMeta_Hmap_Std-1')
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.tcr.reps.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.1, do.legend=FALSE, do.rotate=TRUE, height=10, width=1.6
)
tmp.lab <- paste0('DonorMeta_Hmap_Std-2')
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.tcr.reps.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.1, do.legend=FALSE, do.rotate=TRUE, height=10, width=10
)
# @ Version w/ custom info.
tmp.vals <- c('gender', 'age', 'bmi', 'smoking.status', 'alcohol.units')
tmp.ggplot <- ggplot(data=donor.meta.hmap[donor.var %in% tmp.vals], aes(x=donor.var, y=donor.id.tag, fill=fill.val)) +
    geom_tile(color='black', linewidth=1, linetype=1, alpha=0.8) +
    coord_fixed() +
    scale_fill_manual(values=donor.meta.cols) +
    labs(x='', y='', fill='Value') +
    theme(axis.text.x=element_text(angle=45))
tmp.lab <- paste0('DonorMeta_Hmap_Std-3')
publish.plot(
    tmp.ggplot=tmp.ggplot+theme_void(), output.path=fig.tcr.reps.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.1, do.legend=FALSE, do.rotate=FALSE, height=10, width=1.3
)
tmp.lab <- paste0('DonorMeta_Hmap_Std-4')
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.tcr.reps.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.1, do.legend=FALSE, do.rotate=TRUE, height=10, width=7
)

# ----> Number of T cells per donor.
for(data.set in names(aggr.set.main)){
    tmp.data <- meta.data.l[[data.set]][
        !is.na(clonotype.tag) & !is.na(donor.id.tag),
        .(cell.count=.N),
        by=.(donor.id.tag, gex.lin.tag)
    ]
    to.test <- spread(data=tmp.data, key=gex.lin.tag, value=cell.count)
    tmp.test <- wilcox.test(x=to.test[, 'CD4'], y=to.test[, 'CD8'], paired=TRUE)
    tmp.caption <- paste0('P-value is ', tmp.test$p.value)
    tmp.ggplot <- ggplot(data=tmp.data, aes(x=gex.lin.tag, y=cell.count)) +
        geom_line(aes(group=donor.id.tag), color=gen.link.col, linewidth=gen.link.width) +
        geom_jitter(shape=1, stroke=5, size=10, color=gen.dot.col, width=0, height=0) +
        geom_boxplot(aes(color=gex.lin.tag, fill=gex.lin.tag), alpha=0.4, width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
        scale_color_manual(values=cell.type.cols) +
        scale_fill_manual(values=cell.type.cols) +
        labs(x='Cell type', y='Number of cells', color='', caption=tmp.caption)
    tmp.lab <- paste0(
        '/', aggr.set.labs[data.set], '_CellCount_CellType_Donor'
    )
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.tcr.reps.path, file.name=tmp.lab, type='pdf',
        blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=7, height=14
    )
}

# ----> Fraction of clonally expanded cells per donor.
for(data.set in names(aggr.set.main)){
    tmp.data <- meta.data.l[[data.set]][
        !is.na(clonotype.tag) & !is.na(donor.id.tag),
        .(clone.size=.N),
        by=.(clonotype.tag, donor.id.tag, gex.lin.tag)
    ]
    tmp.data <- merge(x=meta.data.l[[data.set]], y=tmp.data, by=c('clonotype.tag', 'donor.id.tag', 'gex.lin.tag'))
    tmp.data <- tmp.data[,
        .(exp.cell.frac=.SD[clone.size>1, .N]/.SD[, .N]),
        by=.(donor.id.tag, gex.lin.tag)
    ]
    to.test <- spread(data=tmp.data, key=gex.lin.tag, value=exp.cell.frac)
    tmp.test <- wilcox.test(x=to.test[, 'CD4'], y=to.test[, 'CD8'], paired=TRUE)
    tmp.caption <- paste0('P-value is ', tmp.test$p.value)
    tmp.ggplot <- ggplot(data=tmp.data, aes(x=gex.lin.tag, y=exp.cell.frac)) +
        geom_line(aes(group=donor.id.tag), color=gen.link.col, linewidth=gen.link.width) +
        geom_jitter(shape=1, stroke=5, size=10, color=gen.dot.col, width=0, height=0) +
        geom_boxplot(aes(color=gex.lin.tag, fill=gex.lin.tag), alpha=0.4, width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=4)) +
        scale_color_manual(values=cell.type.cols) +
        scale_fill_manual(values=cell.type.cols) +
        labs(x='Cell type', y='Fraction of expanded cells', color='', caption=tmp.caption)
    tmp.lab <- paste0(
        '/', aggr.set.labs[data.set], '_ExpandedCellFract_CellType_Donor'
    )
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.tcr.reps.path, file.name=tmp.lab, type='pdf',
        blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=7, height=14
    )
}

# ---> Cell-to-clonotype ratio, comparison between anatomical sites
# Upper bound.
tmp.thold <- 7
for(cell.lin in names(cell.type.cols)){
    # Ratio per anatomic site.
    tmp.data <- lapply(X=names(meta.data.l), FUN=function(aggr.lab){
        tmp.data <- meta.data.l[[aggr.lab]]
        tmp.data <- tmp.data[
            !is.na(donor.id.tag) &
            !is.na(clonotype.tag) &
            data.set %in% aggr.set.main[[aggr.lab]] &
            gex.lin.tag==cell.lin
        ]
        tmp.data <- tmp.data[,
            .(ratio=.N/uniqueN(clonotype.tag)),
            by=.(donor.id.tag)
        ]
        return(tmp.data)
    })
    names(tmp.data) <- aggr.simp.labs[names(meta.data.l)]
    tmp.data <- merge(
        x=tmp.data[['LG']], y=tmp.data[['LN']],
        by='donor.id.tag',
        suffixes=c('.LG', '.LN')
    )
    # Statistical test. Paired Wilcoxon signed-rank test.
    tmp.test <- wilcox.test(
        x=tmp.data[, ratio.LG],
        y=tmp.data[, ratio.LN],
        alternative='two.sided', paired=TRUE, exact=TRUE
    )
    tmp.caption <- paste0('Wilcoxon signed-rank test nominal P-value is: ', tmp.test$p.value)
    # Final format
    tmp.data <- as.data.table(gather(data=tmp.data, key='site', value='ratio', -`donor.id.tag`))
    tmp.data[, site:=str_extract(string=site, pattern='LN|LG$')]
    tmp.limits <- c(1, 7)
    # Set upper bound.
    tmp.data[, bounded.ratio:=ratio]
    tmp.data[, point.status:='ori']
    tmp.data[bounded.ratio>tmp.thold, `:=`(bounded.ratio=tmp.thold, point.status='bounded')]
    # Plot
    tmp.ggplot <- ggplot(data=tmp.data, aes(x=site)) +
        geom_line(aes(y=bounded.ratio, group=donor.id.tag), color=gen.link.col, linewidth=gen.link.width) +
        geom_jitter(aes(y=bounded.ratio, shape=point.status), stroke=5, size=10, color=gen.dot.col, width=0, height=0) +
        geom_boxplot(aes(y=ratio, color=site, fill=site), alpha=0.4, width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
        coord_cartesian(ylim=c(1, tmp.thold)) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=3)) + # expand=expansion(add=c(0, 0))
        scale_shape_manual(values=tmp.shapes) +
        scale_color_manual(values=an.site.cols) +
        scale_fill_manual(values=an.site.cols) +
        labs(x='Cell type', y='Metric', color='', caption=tmp.caption)
    tmp.lab <- paste0('/SiteComp_Metric-Ratio_Lin-', cell.lin)
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.tcr.reps.path, file.name=tmp.lab, type='pdf',
        blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=7, height=14
    )
}

# ---> TCR repertoire profile per donor and cell type. Plot that depicts clonal expansion info.
for(data.set in names(aggr.set.main)){
    # Set order of clonotypes according to clone size.
    tmp.data.1 <- meta.data.l[[data.set]][
        !is.na(clonotype.tag) & !is.na(donor.id.tag),
        .(clone.size=.N),
        by=.(
            donor.id.tag, clonotype.tag, gex.lin.tag
        )
    ]
    setorderv(x=tmp.data.1, cols=c('clone.size'), order=c(-1))
    tmp.data.1[, order.tag:=paste(donor.id.tag, clonotype.tag, gex.lin.tag, sep=';')]
    tmp.data.1$order.tag <- factor(x=tmp.data.1$order.tag, levels=tmp.data.1$order.tag)
    # Set expansion categorical variable.
    tmp.data.1[, expansion.tag:=ifelse(test=clone.size>1, yes='Expanded', no='Unexpanded')]
    # @ Set donor order
    tmp.lvls <- rev(custom.donor.ord.1[custom.donor.ord.1 %in% tmp.data.1$donor.id.tag])
    tmp.vals <- factor(x=as.character(tmp.data.1$donor.id.tag), levels=tmp.lvls)
    set(x=tmp.data.1, j='donor.id.tag', value=tmp.vals)
    # Get number of unique clonotypes per donor and cell type.
    tmp.data.2 <- meta.data.l[[data.set]][
        !is.na(clonotype.tag) & !is.na(donor.id.tag),
        .(clone.count=uniqueN(clonotype.tag)),
        by=.(donor.id.tag, gex.lin.tag)
    ]
    tmp.data.2[clone.count<100, clone.count:=NA]
    # @ Set donor order
    tmp.lvls <- rev(custom.donor.ord.1[custom.donor.ord.1 %in% tmp.data.2$donor.id.tag])
    tmp.vals <- factor(x=as.character(tmp.data.2$donor.id.tag), levels=tmp.lvls)
    set(x=tmp.data.2, j='donor.id.tag', value=tmp.vals)
    # Define file width according to donor amount.
    tmp.width <- tmp.data.1[, uniqueN(donor.id.tag)]*14.62/40 + 0.38
    # Plot.
    for(cell.type in gen.cell.types){
        tmp.cols <- c('Expanded'=cell.type.cols[[cell.type]], 'Unexpanded'='#d3d3d3')
        to.plot.1 <- tmp.data.1[gex.lin.tag==cell.type]
        to.plot.2 <- tmp.data.2[gex.lin.tag==cell.type]
        # Plot.
        tmp.ggplot <- ggplot(data=to.plot.1, aes(x=donor.id.tag, y=clone.size)) +
            geom_bar(aes(group=order.tag, color=expansion.tag), stat='identity', position='stack', width=0.8, linewidth=1, alpha=0) +
            geom_point(data=to.plot.2, aes(y=clone.count), size=7, shape=21, fill='white', color='black') +
            scale_y_continuous(expand=expansion(add=c(20, 30)), breaks=scales::pretty_breaks(n=3)) +
            scale_color_manual(values=tmp.cols) +
            labs(x='Donor ID', y='Cell count', color='')
        tmp.lab <- paste0('/', aggr.set.labs[data.set], '_StackedSizedClones_Donor_Expansion_CellType-', cell.type)
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.tcr.reps.path, file.name=tmp.lab, type='pdf',
            blank.comp=blank.complement.7, do.legend=FALSE, do.rotate=TRUE, width=tmp.width, height=5
        )
    }
}

# ---> Cummulative relative frequency according to frequency rank.
for(data.set in names(aggr.set.main)){
    # ---> Retrieve clonotype's relative frequencies.
    tmp.data.1 <- meta.data.l[[data.set]][
        !is.na(clonotype.tag) & !is.na(donor.id.tag),
        .(freq.abs=.N),
        by=.(
            gex.lin.tag, donor.id.tag, clonotype.tag
        )
    ]
    tmp.data.2 <- tmp.data.1[,
        .(freq.total=sum(freq.abs)),
        by=.(gex.lin.tag, donor.id.tag)
    ]
    tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by=c('gex.lin.tag', 'donor.id.tag'))
    tmp.data <- tmp.data[, freq.rel:=freq.abs/freq.total]
    setorderv(x=tmp.data, col='freq.rel', order=-1)
    # ---> For each sample, determine cummulative relative frequency & scaled rank.
    cell.lins <- tmp.data[, unique(gex.lin.tag)]
    donor.ids <- tmp.data[, unique(donor.id.tag)]
    tmp.data <- lapply(X=cell.lins, FUN=function(cell.lin){
        tmp.data <- lapply(X=donor.ids, FUN=function(donor.id){
            tmp.data <- tmp.data[gex.lin.tag==cell.lin & donor.id.tag==donor.id]
            tmp.data[, rank.raw:=1:.N]
            tmp.data[, rank.scl:=rank.raw*100/.N]
            tmp.data$cumm.freq.rel <- Reduce(x=tmp.data$freq.rel, f=sum, accumulate=TRUE)
            return(tmp.data)
        })
        tmp.data <- rbindlist(l=tmp.data, use.names=TRUE)
        return(tmp.data)
    })
    tmp.data <- rbindlist(l=tmp.data, use.names=TRUE)
    # ---> Find best fit curve w/ AIC
    # Polynomial degree options
    # degree.opts <- 1:10
    # tmp.lins <- unique(gen.cell.types)
    # tmp.evals <- lapply(X=tmp.lins, FUN=function(tmp.lin){
    #     tmp.evals <- lapply(X=degree.opts, FUN=function(d){
    #         tmp.evals <- lapply(X=donor.ids, FUN=function(tmp.donor.id){
    #             to.model <- tmp.data[
    #                 gex.lin.tag==tmp.lin &
    #                 donor.id.tag==tmp.donor.id,
    #                 .(x=rank.scl, y=cumm.freq.rel)
    #             ]
    #             tmp.model <- lm(formula=y ~ poly(x=x, degree=d), data=to.model)
    #             tmp.aic <- AIC(tmp.model)
    #             return(tmp.aic)
    #         })
    #         tmp.evals <- data.table(
    #             donor.id.tag=donor.ids, aic=unlist(tmp.evals)
    #         )
    #         return(tmp.evals)
    #     })
    #     names(tmp.evals) <- degree.opts
    #     tmp.evals <- rbindlist(l=tmp.evals, use.names=TRUE, idcol='degree')
    #     return(tmp.evals)
    # })
    # names(tmp.evals) <- tmp.lins
    # tmp.evals <- rbindlist(l=tmp.evals, use.names=TRUE, idcol='cell.lin')
    # tmp.evals[, .(median(aic)), by=.(cell.lin, degree)][, min(V1), by=cell.lin]
    # Best fit a 10-degree polynomial, so no.
    # ---> Determine closes point to each percent value.
    # For plotting purposes. Otherwise, too many dots to be easily handled.
    tmp.vals <- seq(from=0.01, to=1, length.out=100)
    tmp.data <- lapply(X=tmp.vals, FUN=function(tmp.pct){
        tmp.data.1 <- tmp.data[,
            .(rank.raw=.SD[
                abs(cumm.freq.rel-tmp.pct)==min(abs(cumm.freq.rel-tmp.pct)),
                rank.raw
            ]),
            by=.(gex.lin.tag, donor.id.tag)
        ]
        tmp.data.1 <- merge(x=tmp.data.1, y=tmp.data, by=c('gex.lin.tag', 'donor.id.tag', 'rank.raw'))
        return(tmp.data.1)
    })
    names(tmp.data) <- tmp.vals*100
    tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='rank.pct')
    # ---> Plot.
    for(tmp.lin in unique(gen.cell.types)){
        to.plot <- tmp.data[gex.lin.tag==tmp.lin]
        # @ Stats. Median cummulative frequency close to 10%.
        tmp.caption <- to.plot[,
            .(cumm.freq.rel=.SD[
                which.min(abs(rank.scl-10)),
                cumm.freq.rel
            ]),
            by=donor.id.tag
        ][, median(cumm.freq.rel)*100]
        tmp.caption <- paste0('Median cummulative rel. freq. at 10% of scaled rank: ', round(x=tmp.caption, digits=4))
        # Plot.
        tmp.ggplot <- ggplot(data=to.plot, aes(x=rank.scl, y=cumm.freq.rel)) +
            geom_line(aes(group=donor.id.tag), linewidth=0.8, color=cell.type.cols[tmp.lin]) +
            # geom_point(size=1, shape=1, color=cell.type.cols[tmp.lin]) +
            geom_smooth(
                method='loess', formula=y~x, se=TRUE, level=0.95,
                color='black', linewidth=4
            ) +
            geom_vline(xintercept=10, linewidth=2, color='#000080', linetype='dashed') +
            scale_y_continuous(expand=expansion(add=c(0, 0)), breaks=scales::pretty_breaks(n=3)) +
            scale_x_continuous(expand=expansion(add=c(0, 0)), breaks=scales::pretty_breaks(n=3)) +
            labs(x='Scaled rank', y='Cummulative rel. freq.', color='', caption=tmp.caption)
        tmp.lab <- paste0('/', aggr.set.labs[data.set], '_CummRelFreq_ScaledRank_CellType-', tmp.lin)
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.tcr.reps.path, file.name=tmp.lab, type='pdf',
            blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=14, height=14
        )
    }
}

# ---> Retrieve diversity metrics.
# tmp.metrics <- c('gini', 'inv.simp', 'gini.simp', 'chao1')
# tmp.metrics <- c('gini', 'inv.simp')
tmp.metrics <- c('gini')
div.metrics <- lapply(X=names(aggr.set.main), FUN=function(data.set){
    div.metrics <- lapply(X=tmp.metrics, FUN=function(x){
        cat(x, '\n')
        tmp.data <- summ.div.metrics(
            tcr.meta=meta.data.l[[data.set]],
            donor.meta=donor.meta, div.metric=x,
            main.sets=aggr.set.main[[data.set]]
        )
        return(tmp.data)
    })
    names(div.metrics) <- tmp.metrics
    div.metrics <- rbindlist(l=div.metrics, use.names=TRUE, idcol='metric')
    return(div.metrics)
})
names(div.metrics) <- names(aggr.set.main)

# ---> Donor-specific repertoire diversity per T-cell lineage (CD4 and CD8).
# Process per metric.
for(data.set in names(aggr.set.main)){
    for(tmp.metric in tmp.metrics){
        # ---> Retrieve general data
        tmp.data <- div.metrics[[data.set]][metric==tmp.metric]
        tmp.data <- tmp.data[
            CD8.cell>100 & CD4.cell>100,
            .(donor.id.tag, CD8.All, CD4.All)
        ]
        # ---> Association between lineages
        tmp.data[, tmp.var:='All']
        tmp.ggplot <- ggplot(data=tmp.data, aes(x=CD4.All, y=CD8.All)) +
            geom_point(shape=1, stroke=5, size=12) +
            geom_smooth(method='lm', formula=y~x, se=FALSE, fullrange=TRUE, color='black', linewidth=5) +
            scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
            scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
            labs(x='Diversity metric in CD4 lineage', y='Diversity metric in CD8 lineage', color='Phenotype')
        tmp.lab <- paste0('/', aggr.set.labs[data.set], '_Metric-', tmp.metric, '_Lineage_Association')
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.tcr.reps.path, file.name=tmp.lab, type='pdf',
            stat.cor=TRUE, cor.group='tmp.var',
            blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=14, height=14
        )
        tmp.data[, tmp.var:=NULL]
        # ---> Comparison between lineages.
        # Statistical test. Paired Wilcoxon signed-rank test.
        tmp.test <- wilcox.test(
            x=tmp.data[, CD8.All],
            y=tmp.data[, CD4.All],
            alternative='two.sided', paired=TRUE, exact=TRUE

        )
        tmp.caption <- paste0('Wilcoxon signed-rank test nominal P-value is: ', tmp.test$p.value)
        # Plot
        tmp.data <- as.data.table(gather(data=tmp.data, key='cell.type', value='metric', -`donor.id.tag`))
        tmp.data[, cell.type:=str_extract(string=cell.type, pattern='CD[84]')]
        # tmp.data$cell.type <- factor(x=tmp.data$cell.type, levels=gen.cell.types)
        tmp.limits <- if(tmp.metric=='gini') c(0, 1) else NULL
        tmp.ggplot <- ggplot(data=tmp.data, aes(x=cell.type, y=metric)) +
            geom_line(aes(group=donor.id.tag), color=gen.link.col, linewidth=gen.link.width) +
            geom_jitter(shape=1, stroke=5, size=10, color=gen.dot.col, width=0, height=0) +
            geom_boxplot(aes(color=cell.type, fill=cell.type), alpha=0.4, width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
            scale_y_continuous(limits=tmp.limits, expand=expansion(add=c(0, 0)), breaks=scales::pretty_breaks(n=3)) +
            scale_color_manual(values=cell.type.cols) +
            scale_fill_manual(values=cell.type.cols) +
            labs(x='Cell type', y='Metric', color='', caption=tmp.caption)
        tmp.lab <- paste0('/', aggr.set.labs[data.set], '_Metric-', tmp.metric, '_Lineage_Comp')
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.tcr.reps.path, file.name=tmp.lab, type='pdf',
            blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=7, height=14
        )
    }
    # Output table w/ metric values.
    tmp.file.name <- paste0(fig.tcr.reps.path, '/', aggr.set.labs[data.set], '_DivMetricPerPopAndDonor.csv')
    fwrite(file=tmp.file.name, x=div.metrics[[data.set]], na=NA)
}

# ---> Diversity comparison between anatomical sites
cell.thold <- 100
for(tmp.metric in tmp.metrics){
    for(cell.lin in names(cell.type.cols)){
        tmp.data <- lapply(X=div.metrics, FUN=function(tmp.data){
            tmp.cols <- c('donor.id.tag', paste0(cell.lin, c('.All', '.cell')))
            tmp.data <- tmp.data[metric==tmp.metric]
            tmp.data <- tmp.data[, ..tmp.cols]
            colnames(tmp.data) <- c('donor.id.tag', 'metric', 'cell.count')
            tmp.data <- tmp.data[cell.count>=cell.thold]
            return(tmp.data[, .(donor.id.tag, metric)])
        })
        tmp.data <- merge(
            x=tmp.data[['c40.lg']], y=tmp.data[['c40.ln']],
            by='donor.id.tag',
            suffixes=c('.LG', '.LN')
        )
        # Statistical test. Paired Wilcoxon signed-rank test.
        tmp.test <- wilcox.test(
            x=tmp.data[, metric.LG],
            y=tmp.data[, metric.LN],
            alternative='two.sided', paired=TRUE, exact=TRUE

        )
        tmp.caption <- paste0('Wilcoxon signed-rank test nominal P-value is: ', tmp.test$p.value)
        # Plot
        tmp.data <- as.data.table(gather(data=tmp.data, key='site', value='metric', -`donor.id.tag`))
        tmp.data[, site:=str_extract(string=site, pattern='LN|LG$')]
        # tmp.data$cell.type <- factor(x=tmp.data$cell.type, levels=gen.cell.types)
        tmp.limits <- if(tmp.metric=='gini') c(0, 1) else NULL
        tmp.ggplot <- ggplot(data=tmp.data, aes(x=site, y=metric)) +
            geom_line(aes(group=donor.id.tag), color=gen.link.col, linewidth=gen.link.width) +
            geom_jitter(shape=1, stroke=5, size=10, color=gen.dot.col, width=0, height=0) +
            geom_boxplot(aes(color=site, fill=site), alpha=0.4, width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
            scale_y_continuous(limits=tmp.limits, expand=expansion(add=c(0, 0)), breaks=scales::pretty_breaks(n=3)) +
            scale_color_manual(values=an.site.cols) +
            scale_fill_manual(values=an.site.cols) +
            labs(x='Cell type', y='Metric', color='', caption=tmp.caption)
        tmp.lab <- paste0('/SiteComp_Metric-', tmp.metric, '_Lin-', cell.lin)
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.tcr.reps.path, file.name=tmp.lab, type='pdf',
            blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=7, height=14
        )
    }
}

# ---> Repertoire overlap between anatomical sites, separate strategies.
# @ Clonotype-based overlap between anatomical sites.
for(tmp.lin in names(cell.type.cols)){
    to.plot <- lapply(X=gen.tcr.data.1, FUN=function(x){
        x <- x[gex.lin.tag==tmp.lin]
        x <- as.data.table(separate_rows(data=x, comb.ids, sep='/'))
        x <- x[!(site.ovlp=='Overlapped' & comb.ids%like%'^LG.|^LN')]
        x <- x[, paste(donor.id.tag, comb.ids, sep=';')]
        # x <- x[, comb.ids]
        # x <- unlist(str_split(string=x, pattern='\\/'))
        return(x)
    })
    tmp.vals <- c(`c40.lg`='LG', `c40.ln`='LN')
    names(to.plot) <- tmp.vals[names(to.plot)]
    # Complete v
    tmp.file.name <- paste0(fig.tcr.reps.path, '/SimBetSites_Ovlp-Clone', '_Lin-', tmp.lin, '.C.tiff')
    VennDiagram::venn.diagram(
        x=to.plot,
        force.unique=TRUE,
        hyper.test=TRUE, total.population=uniqueN(unlist(to.plot)),
        filename=tmp.file.name, lwd=3, fill=an.site.cols, col='black',
        disable.logging=TRUE
    )
    # Blank v
    tmp.file.name <- paste0(fig.tcr.reps.path, '/SimBetSites_Ovlp-Clone', '_Lin-', tmp.lin, '.B.tiff')
    VennDiagram::venn.diagram(
        x=to.plot,
        force.unique=TRUE,
        hyper.test=FALSE, total.population=uniqueN(unlist(to.plot)),
        filename=tmp.file.name, lwd=3, fill=an.site.cols, col='black',
        cat.cex=FALSE, cex=FALSE,
        disable.logging=TRUE
    )
}
# @ Cell-based overlap between anatomical sites. Plot that depicts clonal expansion info.
for(tmp.site in names(gen.tcr.data.1)){
    tmp.data.1 <- gen.tcr.data.1[[tmp.site]]
    setorderv(x=tmp.data.1, cols=c('site.ovlp', 'freq.abs'), order=c(1, -1))
    tmp.data.1[, order.tag:=paste(donor.id.tag, clonotype.tag, gex.lin.tag, sep=';')]
    tmp.data.1$order.tag <- factor(x=tmp.data.1$order.tag, levels=tmp.data.1$order.tag)
    # Set expansion categorical variable.
    tmp.data.1[, expansion.tag:=ifelse(test=freq.abs>0, yes='Expanded', no='Unexpanded')]
    tmp.data.1[expansion.tag=='Expanded', expansion.tag:=site.ovlp]
    # @ Set donor order
    tmp.lvls <- rev(custom.donor.ord.1[custom.donor.ord.1 %in% tmp.data.1$donor.id.tag])
    tmp.vals <- factor(x=as.character(tmp.data.1$donor.id.tag), levels=tmp.lvls)
    set(x=tmp.data.1, j='donor.id.tag', value=tmp.vals)
    # Get number of unique clonotypes per donor, cell type and overlap status.
    tmp.data.2 <- tmp.data.1[,
        .(clone.count=uniqueN(clonotype.tag)),
        by=.(site.ovlp, donor.id.tag, gex.lin.tag)
    ]
    tmp.data.2[clone.count<100, clone.count:=NA]
    # @ Set donor order
    tmp.lvls <- rev(custom.donor.ord.1[custom.donor.ord.1 %in% tmp.data.2$donor.id.tag])
    tmp.vals <- factor(x=as.character(tmp.data.2$donor.id.tag), levels=tmp.lvls)
    set(x=tmp.data.2, j='donor.id.tag', value=tmp.vals)
    # Define file width according to donor amount.
    tmp.width <- tmp.data.1[, uniqueN(donor.id.tag)]*14.62/40 + 0.38
    # Plot.
    for(cell.type in unique(gen.cell.types)){
        anti.site <- names(an.site.cols)[names(an.site.cols)!=aggr.simp.labs[tmp.site]]
        tmp.cols <- c(
            `Overlapped`=unname(an.site.cols[anti.site]),
            # `Unique`=unname(cell.type.cols[cell.type]),
            `Unique`=unname(an.site.cols[aggr.simp.labs[tmp.site]]),
            'Unexpanded'='#D3D3D3'
        )
        to.plot.1 <- tmp.data.1[gex.lin.tag==cell.type]
        to.plot.2 <- tmp.data.2[gex.lin.tag==cell.type]
        # Plot.
        tmp.ggplot <- ggplot(data=to.plot.1, aes(x=donor.id.tag, y=freq.abs)) +
            geom_bar(aes(group=order.tag, color=expansion.tag), stat='identity', position='stack', width=0.8, linewidth=1, alpha=0) +
            # geom_point(data=to.plot.2, aes(y=clone.count, fill=site.ovlp), size=7, stroke=1.5, shape=21, color='black') +
            scale_y_continuous(expand=expansion(add=c(20, 30)), breaks=scales::pretty_breaks(n=3)) +
            scale_color_manual(values=tmp.cols) +
            # scale_fill_manual(values=tmp.cols) +
            labs(x='Donor ID', y='Cell count', color='')
        tmp.lab <- paste0('SimBetSites_Ovlp-Cell', '_Lin-', cell.type, '_Site-', aggr.set.labs[tmp.site], '_StackedSizedClones')
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.tcr.reps.path, file.name=tmp.lab, type='pdf',
            blank.comp=blank.complement.7, do.legend=FALSE, do.rotate=TRUE, width=tmp.width, height=5
        )
    }
}

# ---> Repertoire comparison between anatomical sites, combined (VJ's) strategy, general function.
plot.size.aware.clone.sharing <- function(
    gen.tcr.data=NULL,
    clone.group.tags=NULL,
    groups.to.filter=NULL,
    y.log=FALSE,
    size.tholds=c(2, 5, 10, 15),
    this.reports.path, file.prfx
){
    if(is.null(gen.tcr.data) & is.null(clone.group.tags)){
        this.tcr.data <- gen.tcr.data.1
    }else{
        if(!is.null(gen.tcr.data) & !is.null(clone.group.tags)){
            # Apply recalculation of clonotype frequencies and overlap on the basis of extra columns provided by user.
            tmp.data <- get.comp.clone.freqs(
                gen.tcr.data=gen.tcr.data,
                clone.group.tags=clone.group.tags,
                donor.subset=site.donor.ovlp
            )
            this.tcr.data <- tmp.data[[1]]
            # Filter if requested.
            if(!is.null(groups.to.filter)){
                # Check if we have proper input.
                tmp.check <- is.list(groups.to.filter) & length(groups.to.filter)==length(clone.group.tags)
                if(!tmp.check) stop('Faulty filtering rules. Must a list with length matching the clonotype group vector\'s.\n')
                for(tmp.idx in 1:length(clone.group.tags)){
                    tmp.col <- clone.group.tags[tmp.idx]
                    tmp.vals <- groups.to.filter[[tmp.idx]]
                    this.tcr.data <- lapply(X=this.tcr.data, FUN=function(tmp.data){
                        tmp.data <- tmp.data[
                            !is.na(tmp.col) & get(tmp.col) %in% tmp.vals
                        ]
                        return(tmp.data)
                    })
                }
            }
        }else{
            stop('Faulty plug-ins for gen.tcr.data and clone.group.tags. When either is different than NULL, the other one is also expected to be different than NULL.\n')
        }
    }
    for(tmp.site in names(this.tcr.data)){
        # Plot.
        for(tmp.lin in unique(gen.cell.types)){
            anti.site <- names(an.site.cols)[names(an.site.cols)!=aggr.simp.labs[tmp.site]]
            anti.site <- names(aggr.simp.labs)[aggr.simp.labs==anti.site]
            # Collect data on reference site.
            tmp.data.1 <- this.tcr.data[[tmp.site]][
                gex.lin.tag==tmp.lin,
                .(clonotype.tag, donor.id.tag, freq.abs, freq.rel, comb.ids, site.ovlp)
            ]
            tmp.data.1 <- as.data.table(separate_rows(data=tmp.data.1, comb.ids, sep='/'))
            new.cols <- if(tmp.site=='c40.lg') c('ref.id', 'qry.id') else c('qry.id', 'ref.id')
            fill.opt <- if(tmp.site=='c40.lg') 'right' else 'left'
            tmp.data.1 <- as.data.table(separate(data=tmp.data.1, col='comb.ids', into=new.cols, sep=';', fill=fill.opt)) 
            tmp.data.1 <- tmp.data.1[!(site.ovlp=='Overlapped' & is.na(qry.id))]
            # tmp.data.1[, .N, by=.(clonotype.tag, donor.id.tag)][N>1] # Some reference clonotypes are expected to come up multiple times and they represent matches to multi-chain clonotypes in the inquiry dataset that had a sibling single-chain clonotype.
            tmp.data.1[, ref.id:=str_replace(string=ref.id, pattern=paste0(aggr.simp.labs[tmp.site], '.'), replacement='')]
            tmp.check <- tmp.data.1[, all(clonotype.tag==ref.id)]
            if(!tmp.check) stop('Unexpected error.\n')
            tmp.data.1[, ref.id:=NULL]
            # Collect data on inquiry site and merge accordingly.
            tmp.data.2 <- this.tcr.data[[anti.site]][
                gex.lin.tag==tmp.lin,
                .(clonotype.tag, donor.id.tag, freq.abs, freq.rel)
            ]
            tmp.data <- merge(
                x=tmp.data.1, y=tmp.data.2,
                by.x=c('donor.id.tag', 'qry.id'),
                by.y=c('donor.id.tag', 'clonotype.tag'),
                all.x=TRUE, all.y=FALSE,
                suffixes=c('.ref', '.qry')
            )
            # Force unique entries per clonotypes. Max among match options.
            tmp.warn <- tmp.data[, .N, by=.(clonotype.tag, donor.id.tag)][, .SD[N>1, .N]/.N]
            if(tmp.warn>0.1) warning(paste0('Total of ', tmp.warn, ' reference clonotypes had multiple match options. Max frequencies selected among them.'))
            tmp.data.1 <- tmp.data[, .(entry.no=.N), by=.(clonotype.tag, donor.id.tag)]
            tmp.data <- merge(
                x=tmp.data, y=tmp.data.1,
                by=c('clonotype.tag', 'donor.id.tag')
            )
            tmp.data.1 <- tmp.data[entry.no==1]
            tmp.data.2 <- tmp.data[entry.no>1]
            setorderv(x=tmp.data.2, cols='freq.abs.qry', order=-1) # To overcome top expansion ranks.
            tmp.data.2 <- tmp.data.2[,
                .SD[1, .(qry.id, freq.abs.ref, freq.rel.ref, site.ovlp, freq.abs.qry, freq.rel.qry, entry.no)],
                by=.(clonotype.tag, donor.id.tag)
            ]
            tmp.data <- rbind(tmp.data.1, tmp.data.2)
            # Final sanity checks and format.
            tmp.data[, entry.no:=NULL]
            tmp.check <- tmp.data[site.ovlp=='Overlapped', !any(is.na(freq.rel.qry))]
            if(!tmp.check) stop('Unexpected error.\n')
            tmp.check <- tmp.data[, .N, by=.(clonotype.tag, donor.id.tag)][N>1, .N==0]
            if(!tmp.check) stop('Unexpected error.\n')
            setorderv(x=tmp.data, cols=c('site.ovlp', 'freq.rel.ref'), order=c(-1, -1))
            for(tmp.thold in size.tholds){
                # Filter according to absolute frequency threshold and sort accordingly.
                to.plot <- tmp.data[
                    freq.abs.ref>=tmp.thold,
                    .(
                        rank=1:.N,
                        freq.rel.ref, freq.rel.qry, site.ovlp
                    )
                ]
                v.line.x <- to.plot[site.ovlp=='Overlapped', min(rank)]
                to.plot <- as.data.table(gather(data=to.plot, key=type, value=freq.rel, -`rank`, -`site.ovlp`))
                to.plot <- to.plot[!is.na(freq.rel)]
                to.plot[, site:=aggr.simp.labs[tmp.site]]
                to.plot[type %like% 'qry$', site:=aggr.simp.labs[anti.site]]
                to.plot[, type:=NULL]
                to.plot[site.ovlp=='Unique', site:='REF']
                # Set upper bound of y scale (99 percentile).
                tmp.bound <- ceiling(to.plot[, quantile(x=freq.rel, probs=0.99)]*100)/100
                to.plot[freq.rel>tmp.bound, freq.rel:=tmp.bound]
                # Set line order.
                tmp.lvls <- rev(c(aggr.simp.labs[tmp.site], aggr.simp.labs[anti.site], 'REF'))
                tmp.vals <- factor(x=to.plot[, as.character(site)], levels=tmp.lvls)
                set(x=to.plot, j='site', value=tmp.vals)
                # Set x axis breaks.
                # tmp.breaks.1 <- ceiling(to.plot[site=='REF', max(rank)/2])
                tmp.breaks.1 <- ceiling(seq(from=1, to=to.plot[site=='REF', max(rank)], length.out=4))
                tmp.breaks.1 <- tmp.breaks.1[1:3]
                tmp.breaks.2 <- to.plot[site!='REF', range(rank)]
                tmp.breaks.2 <- c(
                    tmp.breaks.2[1],
                    tmp.breaks.2[1] + ceiling((tmp.breaks.2[2] - tmp.breaks.2[1])/2),
                    tmp.breaks.2[2]
                )
                tmp.breaks <- c(1, tmp.breaks.1, tmp.breaks.2)
                # Plot 1.
                tmp.cols <- c(
                    `REF`=unname(an.site.cols[aggr.simp.labs[tmp.site]]),
                    an.site.cols
                )
                tmp.ggplot <- if(y.log) ggplot(data=to.plot, aes(x=rank, y=log10((freq.rel+0.01)*100))) else ggplot(data=to.plot, aes(x=rank, y=freq.rel))
                tmp.ggplot <- tmp.ggplot +
                    geom_line(aes(color=site), linewidth=3.2) +
                    # geom_vline(xintercept=v.line.x, linewidth=2, color='#000080', linetype='dashed') +
                    scale_x_continuous(expand=expansion(add=c(0, 0)), breaks=tmp.breaks) +
                    scale_y_continuous(expand=expansion(add=c(0, 0)), limits=c(0, NA), breaks=scales::pretty_breaks(n=3)) +
                    scale_color_manual(values=tmp.cols) +
                    labs(x='Clonotype rank', y='Relative freq. (% of T cells in site)', color='Type')
                tmp.lab <- paste0(file.prfx, '_Site-', aggr.simp.labs[tmp.site], '_Lin-', tmp.lin, '_Thold-', tmp.thold, '_App-1')
                publish.plot(
                    tmp.ggplot=tmp.ggplot, output.path=this.reports.path, file.name=tmp.lab, type='pdf',
                    blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=28, height=7
                )
                # Plot 2.
                # Set group order.
                tmp.lvls <- c(aggr.simp.labs[tmp.site], aggr.simp.labs[anti.site], 'REF')
                tmp.vals <- factor(x=to.plot[, as.character(site)], levels=tmp.lvls)
                set(x=to.plot, j='site', value=tmp.vals)
                tmp.cols <- c(
                    # `REF`=unname(an.site.cols[aggr.simp.labs[tmp.site]]),
                    `REF`='#D6A3E3',
                    '#FFFFFF'
                )
                names(tmp.cols)[2] <- aggr.simp.labs[tmp.site] 
                tmp.vals <- c('REF', aggr.simp.labs[tmp.site])
                tmp.ggplot <- ggplot(data=to.plot[site%in%tmp.vals], aes(x='', fill=site)) +
                    geom_bar(position='fill', width=0.6, color='black', linewidth=1) +
                    scale_y_continuous(expand=expansion(add=c(0, 0)), breaks=scales::pretty_breaks(n=3)) +
                    scale_fill_manual(values=tmp.cols) +
                    labs(x='', y='Clonotype %', fill='Type')
                tmp.lab <- paste0(file.prfx, '_Site-', aggr.simp.labs[tmp.site], '_Lin-', tmp.lin, '_Thold-', tmp.thold, '_App-2')
                publish.plot(
                    tmp.ggplot=tmp.ggplot, output.path=this.reports.path, file.name=tmp.lab, type='pdf',
                    blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=2, height=28
                )
            }
        }
    }
}

# ---> Repertoire comparison between anatomical sites, combined (VJ's) strategy, process.
tmp.reports.path <- paste0(fig.tcr.reps.path, '/site_comp_combined')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
# Run process.
plot.size.aware.clone.sharing(
    gen.tcr.data=NULL,
    clone.group.tags=NULL,
    groups.to.filter=NULL,
    y.log=FALSE,
    size.tholds=c(2, 5, 10, 15),
    this.reports.path=tmp.reports.path, file.prfx='FreqCompBetSites'
)
plot.size.aware.clone.sharing(
    gen.tcr.data=NULL,
    clone.group.tags=NULL,
    groups.to.filter=NULL,
    y.log=TRUE,
    size.tholds=c(2, 5, 10, 15),
    this.reports.path=tmp.reports.path, file.prfx='FreqCompBetSites_Log-y'
)

# ---> Repertoire overlap between anatomical sites, integrative strategy (circos plot).
tmp.reports.path <- paste0(fig.tcr.reps.path, '/site_comp_circos')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
# ---> Process per T-cell lineage and donor.
tmp.donors <- sort(gen.tcr.data.1[[1]][, unique(donor.id.tag)])
tmp.lins <- sort(gen.tcr.data.1[[1]][, unique(gex.lin.tag)])
for(tmp.donor in tmp.donors){
    for(tmp.lin in tmp.lins){
        # ---> Segment-specific info
        site.orders <- c(
            `c40.lg`=-1,
            `c40.ln`=1
        )
        tmp.data.1 <- lapply(X=names(gen.tcr.data.1), FUN=function(tmp.site){
            tmp.data <- gen.tcr.data.1[[tmp.site]]
            tmp.data.1 <- tmp.data[
                gex.lin.tag==tmp.lin &
                donor.id.tag==tmp.donor,
                .(clonotype.tag, freq.rel, freq.abs, site.ovlp)
            ]
            setorderv(x=tmp.data.1, cols='freq.rel', order=site.orders[tmp.site])
            tmp.data.2 <- Reduce(x=tmp.data.1[, freq.rel], f=sum, accumulate=TRUE)
            tmp.data.2 <- c(0, tmp.data.2[1:(length(tmp.data.2)-1)])
            tmp.data.1$cumm.rel.freq <- tmp.data.2
            return(tmp.data.1)
        })
        names(tmp.data.1) <- aggr.simp.labs[names(gen.tcr.data.1)]
        tmp.data.1 <- rbindlist(l=tmp.data.1, use.names=TRUE, idcol='site')
        # ---> Link-specific info.
        tmp.data.2 <- lapply(X=gen.tcr.data.2, FUN=function(tmp.data){
            tmp.data <- tmp.data[
                site.ovlp=='Overlapped' &
                gex.lin.tag==tmp.lin &
                donor.id.tag==tmp.donor, 
                comb.ids
            ]
            tmp.data <- str_split(string=tmp.data, pattern='\\/')
            tmp.data <- unlist(tmp.data)
            tmp.data <- tmp.data[str_detect(string=tmp.data, pattern=';')]
            return(tmp.data)
        })
        tmp.check <- all(tmp.data.2[[1]] %in% tmp.data.2[[2]]) & all(tmp.data.2[[2]] %in% tmp.data.2[[1]])
        if(!tmp.check) stop('Unexpected error!\n')
        tmp.data.2 <- tmp.data.2[[1]]
        tmp.data.2 <- data.table(link=tmp.data.2)
        tmp.data.2 <- as.data.table(separate(data=tmp.data.2, col=link, into=c('lg.id', 'ln.id'), sep=';'))

        # ---> Circos plot, option 4
        n.opts <- c(50, 100, 250)
        for(n.opt in n.opts){
            # ---> Segment-specific info
            site.orders <- c(
                `c40.lg`=-1,
                `c40.ln`=1
            )
            to.plot.1 <- lapply(X=names(gen.tcr.data.1), FUN=function(tmp.site){
                tmp.data <- gen.tcr.data.1[[tmp.site]]
                tmp.data.1 <- tmp.data[
                    gex.lin.tag==tmp.lin &
                    donor.id.tag==tmp.donor,
                    .(clonotype.tag, freq.rel, freq.abs, site.ovlp)
                ]
                setorderv(x=tmp.data.1, cols='freq.abs', order=-1)
                tmp.data.1 <- tmp.data.1[1:n.opt, ]
                tmp.data.1[, freq.total:=sum(freq.abs)]
                tmp.data.1[, freq.rel:=freq.abs/freq.total]
                tmp.data.1$cumm.rel.freq <- c(0, Reduce(x=tmp.data.1$freq.rel, f=sum, accumulate=TRUE)[1:tmp.data.1[, .N-1]])
                setorderv(x=tmp.data.1, cols='freq.rel', order=site.orders[tmp.site])
                tmp.data.2 <- Reduce(x=tmp.data.1[, freq.rel], f=sum, accumulate=TRUE)
                tmp.data.2 <- c(0, tmp.data.2[1:(length(tmp.data.2)-1)])
                tmp.data.1$cumm.rel.freq <- tmp.data.2
                return(tmp.data.1)
            })
            names(to.plot.1) <- aggr.simp.labs[names(gen.tcr.data.1)]
            to.plot.1 <- rbindlist(l=to.plot.1, use.names=TRUE, idcol='site')
            # ---> Link-specific info.
            to.plot.2 <- tmp.data.2[
                lg.id %in% to.plot.1[site=='LG', clonotype.tag] &
                ln.id %in% to.plot.1[site=='LN', clonotype.tag]
            ]
            # ---> Output
            tmp.file.name <- paste0(
                tmp.reports.path, '/',
                'SimBetSites_Circos-4_Size-', n.opt, '_Lin-', tmp.lin, '_Donor-', tmp.donor, '.B.pdf'
            )
            pdf(file=tmp.file.name)
            # Initialize the Circos plot
            circos.par('track.height'=0.3)
            circos.initialize(
                sectors=to.plot.1[, unique(site)],
                xlim=c(0, 1)
            )
            # Add sector tracks
            circos.track(ylim=c(0, 1), panel.fun=function(x, y){
                sector <- CELL_META$sector.index
                xlim <- CELL_META$xlim
                ylim <- CELL_META$ylim
                circos.rect(
                    xlim[1], 0, xlim[2], 1,
                    col=an.site.cols[sector]
                )
                # Lines showing clonotype relative frequency.
                circos.lines(
                    x=to.plot.1[site==sector, cumm.rel.freq],
                    y=rep(x=1, times=to.plot.1[site==sector, .N]),
                    type='h', col='black'
                )
                # circos.rect(
                #     xleft=0, ybottom=0, xright=1, ytop=0.5,
                #     col=adjustcolor("#FDFD96", alpha.f=0)
                # )
            }, track.height=0.22, bg.border=NA)
            # Create links between sectors (representing shared TCRs)
            for(idx in 1:nrow(to.plot.2)){
                tmp.check.1 <- to.plot.1[
                    site=='LG' &
                    clonotype.tag==to.plot.2[idx, lg.id],
                    freq.abs>1
                ]
                tmp.check.2 <- to.plot.1[
                    site=='LN' &
                    clonotype.tag==to.plot.2[idx, ln.id],
                    freq.abs>1
                ]
                tmp.check <- tmp.check.1 | tmp.check.2
                if(tmp.check){
                    to.plot.1[
                        site=='LG' &
                        clonotype.tag==to.plot.2[idx, lg.id],
                        .(start=cumm.rel.freq, end=cumm.rel.freq+freq.rel)
                    ]
                    coords.1 <- unlist(to.plot.1[
                        site=='LG' &
                        clonotype.tag==to.plot.2[idx, lg.id],
                        .(start=cumm.rel.freq, end=cumm.rel.freq+freq.rel)
                    ])
                    coords.2 <- unlist(to.plot.1[
                        site=='LN' &
                        clonotype.tag==to.plot.2[idx, ln.id],
                        .(start=cumm.rel.freq, end=cumm.rel.freq+freq.rel)
                    ])
                    circos.link(
                        sector.index1='LG', point1=coords.1,
                        sector.index2='LN', point2=coords.2,
                        col=adjustcolor("#B6E2D3", alpha.f=0.5),
                        border=adjustcolor("#A8E6CF", alpha.f=1)
                    )
                }
            }
            # Clear the plot
            dev.off()
            circos.clear()
        }
    }
}

# ---> Associations between donor-specific continuous variables and diversity metrics
for(data.set in names(aggr.set.main)){
    # @ Temporary directory
    tmp.reports.path <- paste0(fig.tcr.reps.path, '/', aggr.set.labs[data.set], '_div-met_donor-var_corrs')
    if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
    # Define donor-specific continuous variables
    tmp.vars <- c('age', 'bmi', 'alcohol.units')
    tmp.cols <- c('donor.id.tag', tmp.vars)
    tmp.meta <- donor.meta[, ..tmp.cols]
    # Process per metric.
    for(tmp.metric in tmp.metrics){
        # ---> Retrieve general data
        tmp.data <- div.metrics[[data.set]][metric==tmp.metric]
        tmp.data <- tmp.data[
            CD8.cell>100 & CD4.cell>100,
            .(donor.id.tag, CD8.All, CD4.All)
        ]
        tmp.data <- as.data.table(gather(data=tmp.data, key='cell.type', value='metric', -`donor.id.tag`))
        tmp.data[, cell.type:=str_extract(string=cell.type, pattern='CD4|CD8')]
        tmp.data <- merge(x=tmp.data, y=tmp.meta, by='donor.id.tag')
        # Process per cell type.
        tmp.types <- c('CD4', 'CD8')
        for(tmp.type in tmp.types){
            to.plot <- tmp.data[cell.type==tmp.type, ]
            # Process per donor-specific variable.
            for(tmp.var in tmp.vars){
                to.plot[, tmp.var:='All']
                tmp.lab <- str_to_title(string=str_replace_all(string=tmp.var, pattern='_', replacement=' '))
                tmp.ggplot <- ggplot(data=to.plot, aes_string(x='metric', y=tmp.var)) +
                    geom_point(shape=1, stroke=5, size=10) +
                    geom_smooth(method='lm', formula=y~x, se=FALSE, linewidth=4) +
                    scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
                    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
                    labs(
                        x='Diversity metric',
                        y=tmp.lab
                    )
                tmp.lab <- paste0('/Metric-', tmp.metric, '_Lin-', tmp.type, '_Var-', tmp.lab)
                publish.plot(
                    tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
                    stat.cor=TRUE, cor.group='tmp.var',
                    blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=14, height=14
                )
            }
        }
    }
}


### ------------------------- Supp. figure -------------------------- ###

# ---> QC features distribution after QC filtering.
qc.feats <- c(
    `UMICount`='nCount_RNA',
    `GeneCount`='nFeature_RNA',
    `MitUMIFreq`='percent.mt'
)
for(data.set in names(aggr.set.main)){
    for(qc.feat in names(qc.feats)){
        tmp.feat <- qc.feats[[qc.feat]]
        meta.data.l[[data.set]][, tmp.val:=get(tmp.feat)+1]
        tmp.ggplot <- ggplot(data=meta.data.l[[data.set]], aes(x=facs.sorting.batch.tag, y=tmp.val, fill=facs.sorting.batch.tag)) +
            geom_violin(alpha=0.8, trim=FALSE, adjust=0.8) +
            geom_boxplot(width=0.05, alpha=0.7, outlier.shape=NA) +
            scale_y_log10() +
            labs(x='FACS batch', y=paste0(qc.feat, ' (log10)')) +
            theme_bw() + theme(legend.position='none')
        tmp.lab <- paste0('/', aggr.set.labs[data.set], '_QCMetric-', qc.feat, '_FACSBatch')
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.tcr.reps.path, file.name=tmp.lab, type='pdf',
            blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=28, height=7
        )
        meta.data.l[[data.set]][, tmp.val:=NULL]
    }
}

# ---> Donor-specific cell counts before after filtering.
tmp.data <- as.data.table(gather(data=qc.inf.1, key='qc.stage', value='cell.count', -`pre.donor.id.tag`))
tmp.vals <- c(
    `pre.filter.count`='Pre-filter',
    `post.filter.count`='Post-filter'
)
tmp.data <- tmp.data[qc.stage %in% names(tmp.vals)]
tmp.data[, qc.stage:=factor(x=tmp.vals[qc.stage], levels=tmp.vals)]
tmp.ggplot <- ggplot(data=tmp.data, aes(x=qc.stage, y=cell.count)) +
    geom_boxplot(color='black', width=0.7, linewidth=4, outlier.shape=NA) +
    geom_jitter(size=10, color='black', width=0.3) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
    labs(x='QC stage', y='Number of cells')
tmp.lab <- paste0(
    '/CellCount_QCStage_Donor'
)
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.tcr.reps.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=7, height=14
)

# ---> Donor repertoire similarity between anatomical sites.
# Retrieve similarity indexes. 
tmp.data <- lapply(X=names(cell.type.cols), FUN=function(cell.lin){
    tmp.data <- lapply(X=site.donor.ovlp, FUN=function(donor.id.1){
        tmp.data <- lapply(X=site.donor.ovlp, FUN=function(donor.id.2){
            # Merge data from different sites
            merge.cols <- setdiff(x=tcr.cols, y='clonotype.tag')
            tmp.data <- merge(
                x=gen.tcr.data.1[['c40.lg']][gex.lin.tag==cell.lin & donor.id.tag==donor.id.1],
                y=gen.tcr.data.1[['c40.ln']][gex.lin.tag==cell.lin & donor.id.tag==donor.id.2],
                by=merge.cols, all=TRUE,
                suffixes=c('.LG', '.LN')
            )
            # Cosine similarity
            tmp.vals.1 <- tmp.data[, freq.rel.LG]
            tmp.vals.1[is.na(tmp.vals.1)] <- 0
            tmp.vals.2 <- tmp.data[, freq.rel.LN]
            tmp.vals.2[is.na(tmp.vals.2)] <- 0
            cosine.idx <- lsa::cosine(x=tmp.vals.1, y=tmp.vals.2)[1, 1]
            if(length(cosine.idx)>1){
                warning('Issue found.\n')
                cosine.idx <- cosine.idx[1]
            }
            tmp.data <- data.table(
                # raw.jaccard.idx=get.jaccard(x=sample.1.seqs.1, y=sample.2.seqs.1),
                cosine.idx=cosine.idx
            )
            return(tmp.data)
        })
        names(tmp.data) <- site.donor.ovlp
        tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='donor.id.2')
        return(tmp.data)
    })
    names(tmp.data) <- site.donor.ovlp
    tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='donor.id.1')
    return(tmp.data)
})
names(tmp.data) <- names(cell.type.cols)
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='cell.lin')
# Viz of similarity indexes.
idx.types <- c(
  # 'UnweightJaccard'='raw.jaccard.idx',
  # 'WeightedJaccard'='nor.jaccard.idx',
  'Cosine'='cosine.idx'
)
for(idx.type in names(idx.types)){
    for(tmp.lin in names(cell.type.cols)){
        # Matrix to plot.
        tmp.cols <- c('donor.id.1', 'donor.id.2', idx.types[idx.type])
        to.plot <- tmp.data[
            cell.lin==tmp.lin,
            ..tmp.cols
        ]
        to.plot <- spread(data=to.plot, key=donor.id.2, value=idx.types[idx.type])
        tmp.row.names <- to.plot$donor.id.1
        to.plot$donor.id.1 <- NULL
        to.plot <- as.matrix(to.plot)
        row.names(to.plot) <- tmp.row.names
        to.plot <- to.plot[site.donor.ovlp, site.donor.ovlp]
        # Set metadata.
        col.meta <- data.frame(
            row.names=colnames(to.plot),
            `Donor`=colnames(to.plot)
        )
        # Define colors for metadata tracks.
        tmp.cols <- viridisLite::viridis(n=length(site.donor.ovlp), option='plasma')
        names(tmp.cols) <- site.donor.ovlp
        ann.colors <- list(
            `Donor`=tmp.cols
        )
        # Set color scale and breaks for heatmap.
        col.breaks <- seq(from=min(range(to.plot, na.rm=TRUE)), to=max(range(to.plot, na.rm=TRUE)), length.out=100)
        mid.point <- which.min(abs(col.breaks - 0))
        hmap.col.scale.1 <- colorRampPalette(c('#ffffff', '#fffff6', '#ffffed', '#ffffb2'))(mid.point)
        hmap.col.scale.2 <- colorRampPalette(c('#ffff77', '#ffff00', '#ffffb2', '#ffff77', '#ffd280', '#ffd280', '#ffbc40', '#ffa500', '#ff5300', '#ff0000', '#8b0000', '#3f0000'))(100-(mid.point+1))
        hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
        # Plot.
        tmp.file.name <- paste0(fig.tcr.reps.path, '/SimBetSites_', idx.type, '_Lin-', tmp.lin, '.C.pdf')
        pheatmap(
            mat=to.plot, scale='none',
            color=hmap.col.scale, breaks=col.breaks, border_color='black',
            cluster_rows=FALSE, cluster_cols=FALSE,
            annotation_col=col.meta, annotation_row=col.meta, annotation_colors=ann.colors,
            show_colnames=FALSE, show_rownames=TRUE,
            legend=TRUE, annotation_legend=TRUE, annotation_names_row=TRUE, annotation_names_col=TRUE,
            filename=tmp.file.name, heigh=5, width=5.5
        )
        # Blank version.
        tmp.file.name <- paste0(fig.tcr.reps.path, '/SimBetSites_', idx.type, '_Lin-', tmp.lin, '.B.pdf')
        pheatmap(
            mat=to.plot, scale='none',
            color=hmap.col.scale, breaks=col.breaks, border_color='black',
            cluster_rows=FALSE, cluster_cols=FALSE,
            annotation_col=col.meta, annotation_row=col.meta, annotation_colors=ann.colors,
            show_colnames=FALSE, show_rownames=FALSE,
            legend=FALSE, annotation_legend=FALSE, annotation_names_row=FALSE, annotation_names_col=FALSE,
            filename=tmp.file.name, heigh=5, width=5
        )
        # Raw data.
        tmp.file.name <- paste0(fig.tcr.reps.path, '/SimBetSites_', idx.type, '_Lin-', tmp.lin, '.csv')
        write.csv(file=tmp.file.name, x=to.plot, row.names=TRUE)
    }
}


############    -----------------------------------------    ############
### --------------------------- Figure on --------------------------- ###
### ---------------------- TRM TCR repertoires ---------------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
fig.trm.reps.path <- paste0(reports.path, '/figure_on_trm_reps')
if(!dir.exists(fig.trm.reps.path)) dir.create(fig.trm.reps.path)
# ---> Comments: Main dataset for this section, Lower-CD8 and Lower-CD4


### ------------------------- Text details -------------------------- ###


### -------------------------- Main Figure -------------------------- ###

# ----> Fraction of clonally expanded cells per donor.
for(data.set in names(aggr.set.main)){
    tmp.data <- meta.data.l[[data.set]][
        !is.na(clonotype.tag) & !is.na(donor.id.tag),
        .(
            clone.size=.N,
            in.trm=any(trm.tag=='TRM')
        ),
        by=.(clonotype.tag, donor.id.tag, gex.lin.tag)
    ]
    tmp.data <- merge(x=meta.data.l[[data.set]], y=tmp.data, by=c('clonotype.tag', 'donor.id.tag', 'gex.lin.tag'))
    tmp.data <- tmp.data[
        in.trm==TRUE,
        .(exp.cell.frac=.SD[clone.size>1, .N]/.SD[, .N]),
        by=.(donor.id.tag, gex.lin.tag)
    ]
    to.test <- spread(data=tmp.data, key=gex.lin.tag, value=exp.cell.frac)
    tmp.test <- wilcox.test(x=to.test[, 'CD4'], y=to.test[, 'CD8'], paired=TRUE)
    tmp.caption <- paste0('P-value is ', tmp.test$p.value)
    tmp.ggplot <- ggplot(data=tmp.data, aes(x=gex.lin.tag, y=exp.cell.frac)) +
        geom_boxplot(aes(color=gex.lin.tag), width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
        geom_jitter(shape=1, stroke=5, size=10, color='black', width=0, height=0) +
        geom_line(aes(group=donor.id.tag), linewidth=1) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=4)) +
        scale_color_manual(values=cell.type.cols) +
        labs(x='Cell type', y='Fraction of expanded cells', color='', caption=tmp.caption)
    tmp.lab <- paste0(
        '/', aggr.set.labs[data.set], '_ExpandedCellFract_CellType_Donor'
    )
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.trm.reps.path, file.name=tmp.lab, type='pdf',
        blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=7, height=14
    )
}

# ---> Cummulative relative frequency according to frequency rank.
tmp.groups <- c(names(gen.subsets), 'group.class.tag')
for(aggr.set in names(aggr.set.main)){
    for(tmp.subset in tmp.groups){
        if(tmp.subset=='group.class.tag'){
            subset.labs <- meta.data.l[[aggr.set]][
                !is.na(group.class.tag),
                sort(unique(as.character(group.class.tag)))
            ][1:4]
        }else{
            subset.labs <- gen.subsets[[tmp.subset]][1]
        }
        for(subset.lab in subset.labs){
            # ---> Retrieve clonotype's relative frequencies.
            tmp.data.1 <- meta.data.l[[aggr.set]][
                !is.na(clonotype.tag) & !is.na(donor.id.tag) & get(tmp.subset)==subset.lab,
                .(freq.abs=.N),
                by=.(
                    data.set, donor.id.tag, clonotype.tag
                )
            ]
            if(nrow(tmp.data.1)<100) next
            tmp.data.2 <- tmp.data.1[,
                .(freq.total=sum(freq.abs)),
                by=.(data.set, donor.id.tag)
            ]
            tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by=c('data.set', 'donor.id.tag'))
            tmp.data <- tmp.data[, freq.rel:=freq.abs/freq.total]
            setorderv(x=tmp.data, col='freq.rel', order=-1)
            # ---> For each sample, determine cummulative relative frequency & scaled rank.
            cell.lins <- tmp.data[, unique(data.set)]
            donor.ids <- tmp.data[, unique(donor.id.tag)]
            tmp.data <- lapply(X=cell.lins, FUN=function(cell.lin){
                # cat(cell.lin, '\n')
                tmp.data <- lapply(X=donor.ids, FUN=function(donor.id){
                    tmp.data <- tmp.data[data.set==cell.lin & donor.id.tag==donor.id]
                    tmp.data[, rank.raw:=1:.N]
                    tmp.data[, rank.scl:=rank.raw*100/.N]
                    tmp.data$cumm.freq.rel <- Reduce(x=tmp.data$freq.rel, f=sum, accumulate=TRUE)
                    if(nrow(tmp.data)==0) tmp.data <- NA
                    return(tmp.data)
                })
                tmp.data <- tmp.data[!is.na(tmp.data)]
                tmp.data <- rbindlist(l=tmp.data, use.names=TRUE)
                return(tmp.data)
            })
            tmp.data <- rbindlist(l=tmp.data, use.names=TRUE)
            # ---> Determine closes point to each percent value.
            # For plotting purposes. Otherwise, too many dots to be easily handled.
            tmp.vals <- seq(from=0.01, to=1, length.out=100)
            tmp.data <- lapply(X=tmp.vals, FUN=function(tmp.pct){
                tmp.data.1 <- tmp.data[,
                    .(rank.raw=.SD[
                        abs(cumm.freq.rel-tmp.pct)==min(abs(cumm.freq.rel-tmp.pct)),
                        rank.raw
                    ]),
                    by=.(data.set, donor.id.tag)
                ]
                tmp.data.1 <- merge(x=tmp.data.1, y=tmp.data, by=c('data.set', 'donor.id.tag', 'rank.raw'))
                return(tmp.data.1)
            })
            names(tmp.data) <- tmp.vals*100
            tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='rank.pct')
            # ---> Plot.
            for(tmp.lin in aggr.set.main[[aggr.set]]){
                to.plot <- tmp.data[data.set==tmp.lin]
                # @ Stats. Median cummulative frequency close to 10%.
                tmp.caption <- to.plot[,
                    .(cumm.freq.rel=.SD[
                        which.min(abs(rank.scl-10)),
                        cumm.freq.rel
                    ]),
                    by=donor.id.tag
                ][, median(cumm.freq.rel)*100]
                tmp.caption <- paste0('Median cummulative rel. freq. at 10% of scaled rank: ', round(x=tmp.caption, digits=4))
                # Plot.
                tmp.col <- if(tmp.subset=='group.class.tag') clusts.cols[[tmp.lin]][subset.lab] else gen.subset.cols[subset.lab]
                tmp.ggplot <- ggplot(data=to.plot, aes(x=rank.scl, y=cumm.freq.rel)) +
                    geom_line(aes(group=donor.id.tag), linewidth=0.8, color=tmp.col) +
                    # geom_point(size=1, shape=1, color=cell.type.cols[tmp.lin]) +
                    geom_smooth(
                        method='loess', formula=y~x, se=TRUE, level=0.95,
                        color='black', linewidth=4
                    ) +
                    geom_vline(xintercept=10, linewidth=2, color='#000080', linetype='dashed') +
                    scale_y_continuous(expand=expansion(add=c(0, 0)), breaks=scales::pretty_breaks(n=3)) +
                    scale_x_continuous(expand=expansion(add=c(0, 0)), breaks=scales::pretty_breaks(n=3)) +
                    labs(x='Scaled rank', y='Cummulative rel. freq.', color='', caption=tmp.caption)
                tmp.lab <- paste0('/', aggr.set.labs[aggr.set], '_CummRelFreq_ScaledRank_CellType-', gen.cell.types[tmp.lin], '_GenSs-', subset.lab)
                publish.plot(
                    tmp.ggplot=tmp.ggplot, output.path=fig.trm.reps.path, file.name=tmp.lab, type='pdf',
                    blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=14, height=14
                )
            }
        }
    }
}

# ---> Retrieve subset-bias rep. diversity metrics.
tmp.metrics <- c('gini')
bias.div.metrics <- lapply(X=names(aggr.set.main), FUN=function(data.set){
    div.metrics <- lapply(X=tmp.metrics, FUN=function(x){
        cat(x, '\n')
        tmp.data <- summ.div.metrics(
            tcr.meta=meta.data.l[[data.set]],
            donor.meta=donor.meta, div.metric=x,
            main.sets=aggr.set.main[[data.set]],
            tag.of.int='group.class.tag'
        )
        return(tmp.data)
    })
    names(div.metrics) <- tmp.metrics
    div.metrics <- rbindlist(l=div.metrics, use.names=TRUE, idcol='metric')
    return(div.metrics)
})
names(bias.div.metrics) <- names(aggr.set.main)

# ---> Compare diversity between subset-biased cells vs raw clusters.
groups.of.int <- c('CD4.0', 'CD8.0')
tmp.metrics <- c('gini')
for(data.set in names(aggr.set.main)){
    tmp.data <- list(
        `raw`=div.metrics[[data.set]],
        `bias`=bias.div.metrics[[data.set]]
    )
    tmp.data <- lapply(X=tmp.data, FUN=function(tmp.data){
        tmp.cols <- colnames(tmp.data)[str_detect(string=colnames(tmp.data), pattern='^CD[84]\\.\\d+$')]
        tmp.cols <- c('metric', 'donor.id.tag', tmp.cols)
        tmp.data <- tmp.data[, ..tmp.cols]
        tmp.data <- gather(data=tmp.data, key=subset, value=idx, -`metric`, -`donor.id.tag`)
        tmp.data <- as.data.table(tmp.data)
        return(tmp.data)
    })
    tmp.data <- merge(
        x=tmp.data[['raw']], y=tmp.data[['bias']],
        by=c('metric', 'donor.id.tag', 'subset'),
        suffixes=c('.raw', '.bias')
    )
    tmp.check <- all(groups.of.int %in% tmp.data[, subset])
    if(!tmp.check) stop('Unexpected error.\n')
    for(tmp.metric in tmp.metrics){
        tmp.limits <- if(tmp.metric=='gini') c(0, 1) else NULL
        for(tmp.group in groups.of.int){
            to.plot <- tmp.data[subset==tmp.group]
            to.plot[, `:=`(metric=NULL, subset=NULL)]
            # Statistical test. Paired Wilcoxon signed-rank test.
            tmp.test <- wilcox.test(
                x=to.plot[, idx.raw],
                y=to.plot[, idx.bias],
                alternative='two.sided', paired=TRUE, exact=TRUE
            )
            tmp.caption <- paste0('Wilcoxon signed-rank test nominal P-value is: ', tmp.test$p.value)
            # Plot.
            to.plot <- as.data.table(gather(data=to.plot, key='cell.subset', value='metric', -`donor.id.tag`))
            to.plot[, cell.subset:=str_extract(string=cell.subset, pattern='bias$|raw$')]
            tmp.cols <- c(
                `bias`=unname(gen.subset.cols['TRM']),
                `raw`='#CBCBCB'
            )
            tmp.ggplot <- ggplot(data=to.plot, aes(x=cell.subset, y=metric)) +
                geom_boxplot(aes(color=cell.subset), width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
                geom_jitter(shape=1, stroke=5, size=10, color='black', width=0, height=0) +
                geom_line(aes(group=donor.id.tag), linewidth=1) +
                scale_y_continuous(limits=tmp.limits, expand=expansion(add=c(0, 0)), breaks=scales::pretty_breaks(n=3)) +
                scale_color_manual(values=tmp.cols) +
                labs(x='Cell type', y='Metric', color='', caption=tmp.caption)
            tmp.lab <- paste0('/', aggr.set.labs[data.set], '_DivBiasRawComp_', '_Metric-', tmp.metric, '_', tmp.group)
            publish.plot(
                tmp.ggplot=tmp.ggplot, output.path=fig.trm.reps.path, file.name=tmp.lab, type='pdf',
                blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=7, height=14
            )
        }
    }
}

# ---> Repertoire overlap between anatomical sites, combined strategy. General function to apply process for specific sets of clonotypes as defined on phenotype bias.
plot.size.aware.clone.sharing.s1 <- function(
    subset.rules, y.log=FALSE,
    this.reports.path, rep.lab
){
    # Retrieve input data.
    tmp.data <- lapply(X=names(aggr.set.main), FUN=function(tmp.site){
        # Retrieve input TCR data.
        tmp.data <- copy(meta.data.l[[tmp.site]])
        if('tmp.tag' %in% colnames(tmp.data)) tmp.data[, tmp.tag:=NULL]
        tmp.data <- tmp.data[
            !(is.na(group.class.tag) & is.na(gen.clonotype.status))
        ]
        # Assign targets.
        tmp.data[
            !is.na(clonotype.tag) & !is.na(donor.id.tag),
            tmp.tag:='Rest'
        ]
        for(tmp.lin in aggr.set.main[[tmp.site]]){
            # Define target clonotype set.
            target.clusts <- names(clusts.defs[[tmp.lin]])[clusts.defs[[tmp.lin]] %in% subset.rules[[tmp.site]]]
            tmp.data[
                data.set==tmp.lin &
                !is.na(group.class.tag) &
                group.class.tag %in% target.clusts,
                tmp.tag:='Target'
            ]
        }
        return(tmp.data)
    })
    names(tmp.data) <- names(aggr.set.main)
    # Run process.
    tmp.reports.path <- paste0(this.reports.path, '/', rep.lab)
    if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
    plot.size.aware.clone.sharing(
        gen.tcr.data=tmp.data,
        clone.group.tags='tmp.tag',
        groups.to.filter=list('Target'),
        y.log=y.log,
        size.tholds=c(2, 5, 10, 15),
        this.reports.path=tmp.reports.path, file.prfx='FreqCompBetSites'
    )
}

# ---> Repertoire overlap between anatomical sites, combined strategy. Apply function for subsets defined below.

#       Forcing TRM on both sides.
tmp.rules <- list(
    `c40.lg`='TRM',
    # `c40.ln`=unique(unname(unlist(clusts.defs)))
    `c40.ln`='TRM'
)
plot.size.aware.clone.sharing.s1(
    subset.rules=tmp.rules,
    this.reports.path=fig.trm.reps.path, rep.lab='site_comp_lg-trm_vs_ln-trm'
)

#       Forcing TRM on the lung side only.
tmp.rules <- list(
    `c40.lg`='TRM',
    `c40.ln`=unique(unname(unlist(clusts.defs)))
)
plot.size.aware.clone.sharing.s1(
    subset.rules=tmp.rules,
    this.reports.path=fig.trm.reps.path, rep.lab='site_comp_lg-trm_vs_ln-all'
)
plot.size.aware.clone.sharing.s1(
    subset.rules=tmp.rules, y.log=TRUE,
    this.reports.path=fig.trm.reps.path, rep.lab='site_comp_lg-trm_vs_ln-all_log-y'
)

#       Forcing CTL on lung side only.
tmp.rules <- list(
    `c40.lg`=c('CD4-CTL', 'TEM'),
    `c40.ln`=unique(unname(unlist(clusts.defs)))
)
plot.size.aware.clone.sharing.s1(
    subset.rules=tmp.rules,
    this.reports.path=fig.trm.reps.path, rep.lab='site_comp_lg-ctl_vs_ln-all'
)

#       Everything but TRM on lung side.
tmp.rules <- list(
    `c40.lg`=setdiff(x=unique(unname(unlist(clusts.defs))), y='TRM'),
    `c40.ln`=unique(unname(unlist(clusts.defs)))
)
plot.size.aware.clone.sharing.s1(
    subset.rules=tmp.rules,
    this.reports.path=fig.trm.reps.path, rep.lab='site_comp_lg-no-trm_vs_ln-all'
)

#       Forcing GZMK+ on lung side only.
# Mock for CD4s. Note that results will not apply here for that lineage.
tmp.rules <- list(
    `c40.lg`=c('CD4-CTL', 'GZMKhi'),
    `c40.ln`=unique(unname(unlist(clusts.defs)))
)
plot.size.aware.clone.sharing.s1(
    subset.rules=tmp.rules,
    this.reports.path=fig.trm.reps.path, rep.lab='site_comp_lg-gzmk_vs_ln-all'
)


### ------------------------- Supp. figure -------------------------- ###



############    -----------------------------------------    ############
### --------------------------- Figure on --------------------------- ###
### ---------------------- proliferation assay ---------------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
fig.prol.path <- paste0(reports.path, '/figure_on_prol_assay')
if(!dir.exists(fig.prol.path)) dir.create(fig.prol.path)
# ---> Comments: Main dataset for this section, DICE-TCR cohort TCR data (the "resource").


### ------------------------- Text details -------------------------- ###

# ---> T-cell subset-specific median of the combined number of alpha and beta sequences per peptide pool and donor.
# @ Including cross-reactive TCRs.
tmp.data.1 <- prol.assay.meta[
  !is.na(specificity.class.tag),
  .(
    tcr.id, tcr.chain, specificity.class.tag,
    donor.id.tag=donor.id, cell.subset
  )
]
tmp.data.1 <- tmp.data.1[,
    .(clonotype.abs.freq=uniqueN(tcr.id)),
    by=.(donor.id.tag, cell.subset, specificity.class.tag)
]
tmp.data.2 <- tmp.data.1[,
  .(median=median(clonotype.abs.freq)),
  by=.(cell.subset)
]

# ---> Fraction of clonotypes and cells with an uncovered specificity solely based on the proliferation assay approach vs the integrative approach.
# @ Proliferation assay approach only.
tmp.list <- list(
  'CD8'=pred.ref.1,
  'CD4'=pred.ref.2
)
tmp.data <- lapply(X=names(tmp.list), FUN=function(t.subset){
  tmp.data <- tmp.list[[t.subset]]
  # Stacked clone info.
  tmp.data.1 <- tmp.data[,
    .(
        cell.fract=.SD[!is.na(exp.pred) & exp.pred!='Multiple', sum(full.size)] /
            sum(full.size),
        clone.fract=.SD[!is.na(exp.pred) & exp.pred!='Multiple', uniqueN(qry.clone.id)] /
            uniqueN(qry.clone.id)
    ),
      by=.(
          donor.id.tag
      )
  ]
  tmp.data.1 <- tmp.data.1[,
    .(
        cell.fract.min=min(cell.fract),
        cell.fract.max=max(cell.fract),
        cell.fract.median=median(cell.fract),
        clone.fract.min=min(clone.fract),
        clone.fract.max=max(clone.fract),
        clone.fract.median=median(clone.fract)
    )
  ]
  return(tmp.data.1)
})
names(tmp.data) <- names(tmp.list)
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='cell.subset')
# @ Integrative approach.
tmp.list <- list(
  'CD8'=pred.ref.1,
  'CD4'=pred.ref.2
)
tmp.data <- lapply(X=names(tmp.list), FUN=function(t.subset){
  tmp.data <- tmp.list[[t.subset]]
  # Stacked clone info.
  tmp.data.1 <- tmp.data[,
    .(
        cell.fract=.SD[!is.na(consensus.pred) & consensus.pred!='Multiple', sum(full.size)] /
            sum(full.size),
        clone.fract=.SD[!is.na(consensus.pred) & consensus.pred!='Multiple', uniqueN(qry.clone.id)] /
            uniqueN(qry.clone.id)
    ),
      by=.(
          donor.id.tag
      )
  ]
  tmp.data.1 <- tmp.data.1[,
    .(
        cell.fract.min=min(cell.fract),
        cell.fract.max=max(cell.fract),
        cell.fract.median=median(cell.fract),
        clone.fract.min=min(clone.fract),
        clone.fract.max=max(clone.fract),
        clone.fract.median=median(clone.fract)
    )
  ]
  return(tmp.data.1)
})
names(tmp.data) <- names(tmp.list)
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='cell.subset')

# ---> Number of donors w/ 2 or more specficities inferred solely based on the proliferation assay approach.
tmp.list <- list(
  'CD8'=pred.ref.1,
  'CD4'=pred.ref.2
)
tmp.data <- lapply(X=names(tmp.list), FUN=function(t.subset){
  tmp.data <- tmp.list[[t.subset]]
  tmp.data <- tmp.data[
    !is.na(exp.pred) & exp.pred!='Multiple' & exp.pred!='Cross-reactive',
    .(spc.count=uniqueN(exp.pred)),
    by=donor.id.tag
  ]
  tmp.data <- tmp.data[,
    .(
        donor.rel.freq=.SD[spc.count>1, .N] /
            tmp.list[[t.subset]][, uniqueN(donor.id.tag)],
        donor.abs.freq=.SD[spc.count>1, .N]
    )
  ]
  return(tmp.data)
})
names(tmp.data) <- names(tmp.list)
rbindlist(l=tmp.data, use.names=TRUE, idcol='cell.subset')

# ---> Number of both chain-matched clonotypes per specificity accordance type.
tmp.list <- list(
  'CD8'=pred.ref.1,
  'CD4'=pred.ref.2
)
tmp.data <- lapply(X=names(tmp.list), FUN=function(t.subset){
    tmp.data <- tmp.list[[t.subset]]
    # Clonotype count across donors.
    tmp.data <- tmp.data[
        !is.na(exp.pred) & exp.pred!='Multiple' & exp.chain.class=='BC',
        .(
            clonotype.rel.freq=.SD[exp.match.class=='Perfect', uniqueN(qry.clone.id)] /
                .N,
            clonotype.abs.freq=.SD[exp.match.class=='Perfect', uniqueN(qry.clone.id)]
        )
    ]
    return(tmp.data)
})
names(tmp.data) <- names(tmp.list)
rbindlist(l=tmp.data, use.names=TRUE, idcol='cell.subset')


### -------------------------- Main Figure -------------------------- ###

# ---> Combined number of alpha and beta sequences together per peptide pool and donor.
# @ Including cross-reactive TCRs.
tmp.data.1 <- prol.assay.meta[
  !is.na(specificity.class.tag),
  .(
    tcr.id, tcr.chain, specificity.class.tag,
    donor.id.tag=donor.id, cell.subset
  )
]
tmp.data.1 <- tmp.data.1[,
    .(clonotype.abs.freq=uniqueN(tcr.id)),
    by=.(donor.id.tag, cell.subset, specificity.class.tag)
]
tmp.data.2 <- tmp.data.1[,
  .(median=median(clonotype.abs.freq)),
  by=.(cell.subset)
]
tmp.vals <- tmp.data.1[, unique(cell.subset)]
for(t.subset in tmp.vals){
  tmp.data <- tmp.data.1[cell.subset==t.subset]
  tmp.ggplot <- ggplot(data=tmp.data, aes(x=specificity.class.tag, y=clonotype.abs.freq)) +
      geom_boxplot(aes(color=specificity.class.tag), width=0.7, linewidth=7, outlier.shape=NA) +
      geom_jitter(shape=1, stroke=2.5, size=5, color='black', width=0.22) +
      geom_hline(yintercept=tmp.data.2[cell.subset==t.subset, median], color=cell.type.cols[t.subset], linewidth=5) +
      scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
      scale_color_manual(values=pp.cols) +
      labs(x='Reactivity', y='Clonotype count', color='')
  tmp.lab <- paste0('ChainTotalCount_Reactivity_Donor_', t.subset)
  publish.plot(
      tmp.ggplot=tmp.ggplot, output.path=fig.prol.path, file.name=tmp.lab, type='pdf',
      blank.comp=blank.complement.3.1, do.legend=FALSE,
      width=16, height=8
  )
}
tmp.file.name <- paste0(fig.prol.path, '/ChainTotalCount_Subset.csv')
fwrite(file=tmp.file.name, x=tmp.data.2)

# ---> Combined number of alpha and beta sequences together per peptide pool across donors.
# @ Including cross-reactive TCRs.
tmp.data.1 <- prol.assay.meta[
  !is.na(specificity.class.tag),
  .(
    tcr.id, tcr.chain, specificity.class.tag,
    donor.id.tag=donor.id, cell.subset
  )
]
tmp.data.1 <- tmp.data.1[,
    .(clonotype.abs.freq=uniqueN(tcr.id)),
    by=.(cell.subset, specificity.class.tag)
]
tmp.vals <- tmp.data.1[, unique(cell.subset)]
for(t.subset in tmp.vals){
  tmp.data <- tmp.data.1[cell.subset==t.subset]
  tmp.ggplot <- ggplot(data=tmp.data, aes(x=specificity.class.tag, y=clonotype.abs.freq)) +
      geom_bar(aes(fill=specificity.class.tag), stat='identity', width=0.7, linewidth=3, color='black') +
      scale_y_continuous(expand=expansion(add=c(0, NA)), breaks=scales::pretty_breaks(n=3)) +
      scale_fill_manual(values=pp.cols) +
      labs(x='Reactivity', y='Clonotype count', color='')
  tmp.lab <- paste0('ChainTotalCount_Reactivity_AcrossDonors_', t.subset)
  publish.plot(
      tmp.ggplot=tmp.ggplot, output.path=fig.prol.path, file.name=tmp.lab, type='pdf',
      blank.comp=blank.complement.3.1, do.legend=FALSE,
      width=14, height=7
  )
}

# ---> Combined number of alpha and beta sequences separately per peptide pool and donor.
# @ Including cross-reactive TCRs.
tmp.data.1 <- prol.assay.meta[
  !is.na(specificity.class.tag),
  .(
    tcr.id, tcr.chain, specificity.class.tag,
    donor.id.tag=donor.id, cell.subset
  )
]
tmp.data.1 <- tmp.data.1[,
    .(clonotype.abs.freq=uniqueN(tcr.id)),
    by=.(donor.id.tag, cell.subset, tcr.chain, specificity.class.tag)
]
tmp.data.2 <- tmp.data.1[,
  .(median=median(clonotype.abs.freq)),
  by=.(cell.subset, tcr.chain)
]
# Set upper bounds.
tmp.tholds <- list(
    `CD4`=c(`TCRA`=NA, `TCRB`=600),
    `CD8`=c(`TCRA`=300, `TCRB`=250)
)
tmp.vals.1 <- tmp.data.1[, unique(cell.subset)]
tmp.vals.2 <- tmp.data.1[, unique(tcr.chain)]
for(t.subset in tmp.vals.1){
  for(t.chain in tmp.vals.2){
    # Retrieve data.
    tmp.thold <- tmp.tholds[[t.subset]][t.chain]
    tmp.data <- tmp.data.1[cell.subset==t.subset & tcr.chain==t.chain]
    # Spread to get real distributions.
    tmp.data <- tmp.data[, .(donor.id.tag, specificity.class.tag, clonotype.abs.freq)]
    tmp.data <- spread(data=tmp.data, key=specificity.class.tag, value=clonotype.abs.freq, fill=0)
    tmp.data <- as.data.table(gather(data=tmp.data, key=specificity.class.tag, value=clonotype.abs.freq, -`donor.id.tag`))
    # Set factors.
    tmp.vals <- factor(x=tmp.data[, specificity.class.tag], levels=names(pp.cols)[names(pp.cols)%in%tmp.data[, specificity.class.tag]])
    set(x=tmp.data, j='specificity.class.tag', value=tmp.vals)
    # Set upper bound.
    tmp.data[, bounded.freq:=clonotype.abs.freq]
    tmp.data[, point.status:='ori']
    if(!is.na(tmp.thold)) tmp.data[bounded.freq>tmp.thold, `:=`(bounded.ratio=tmp.thold, point.status='bounded')]
    # Plot.
    tmp.ggplot <- ggplot(data=tmp.data, aes(x=specificity.class.tag)) +
        geom_jitter(aes(y=bounded.freq, shape=point.status), stroke=5, size=10, color=gen.dot.col, width=0.2, height=0) +
        geom_boxplot(aes(y=clonotype.abs.freq, color=specificity.class.tag, fill=specificity.class.tag), alpha=0.4, width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
        coord_cartesian(ylim=c(NA, tmp.thold)) + # Custom close up.
        scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
        scale_shape_manual(values=tmp.shapes) +
        scale_color_manual(values=pp.cols) +
        scale_fill_manual(values=pp.cols) +
        labs(x='Reactivity', y='Clonotype count', color='')
    tmp.lab <- paste0('Chain', t.chain, 'Count_Reactivity_Donor_', t.subset)
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.prol.path, file.name=tmp.lab, type='pdf',
        blank.comp=blank.complement.3, do.legend=FALSE,
        width=28, height=14
    )
  }
}
tmp.file.name <- paste0(fig.prol.path, '/ChainSepCount_Subset.csv')
fwrite(file=tmp.file.name, x=tmp.data.2)


### ------------------------- Supp. figure -------------------------- ###

# ---> Number of clonotypes per general clonotype status.
tmp.data <- fread(file=prol.assay.meta.file)
tmp.data[, tmp.var:=as.character(gen.clonotype.status)]
tmp.data[total.abs.freq==1, tmp.var:='Singleton']
tmp.vals <- factor(x=tmp.data[, tmp.var], levels=c('Singleton', 'Unambiguous', 'Ambiguous'))
set(x=tmp.data, j='tmp.var', value=tmp.vals)
tmp.data.1 <- tmp.data[,
    .(clonotype.abs.freq=uniqueN(tcr.id)),
    by=.(donor.id.tag=donor.id, cell.subset, tcr.chain, gen.clonotype.status=tmp.var)
]
tmp.thold <- 3300
tmp.ggplot <- ggplot(data=tmp.data.1, aes(x=gen.clonotype.status, y=clonotype.abs.freq)) +
    # geom_violin(linewidth=4) +
    geom_boxplot(width=0.5, linewidth=4, fatten=4, outlier.shape=NA) +
    # geom_jitter(shape=1, stroke=5, size=5, color='black', width=0) +
    coord_cartesian(ylim=c(NA, tmp.thold)) + # Custom close up.
    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
    labs(x='Reactivity', y='Clonotype count', color='')
tmp.lab <- 'ChainTotalCount_CloneaStatus'
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.prol.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.3, do.legend=FALSE,
    width=14, height=14
)

# ---> Specificity ratio
# @ Linea scale.
tmp.data <- fread(file=prol.assay.meta.file)
tmp.data <- tmp.data[
  gen.clonotype.status=='Ambiguous',
  .(
    tcr.id, tcr.chain, specificity.class.tag,
    donor.id.tag=donor.id, cell.subset,
    gen.clonotype.status, spc.fc
  )
]
tmp.ggplot <- ggplot(data=tmp.data, aes(x=spc.fc)) +
    geom_density(linewidth=3.5, color='black', fill='gray', alpha=0) +
    geom_vline(xintercept=2, color='red', linetype='dashed', linewidth=2) +
    scale_y_continuous(expand=expansion(add=c(0, 0)), breaks=scales::pretty_breaks(n=3)) +
    scale_x_continuous(expand=expansion(add=c(0, NA)), breaks=scales::pretty_breaks(n=3)) +
    coord_cartesian(xlim=c(NA, 100)) + # Custom close up.
    labs(x='Specificity ratio', y='Emprirical density')
tmp.lab <- 'SpcFC_EmpDensity_Linear'
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.prol.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.3, do.legend=FALSE,
    width=10, height=10
)
# @ Log scale.
tmp.data <- fread(file=prol.assay.meta.file)
tmp.data <- tmp.data[
  gen.clonotype.status=='Ambiguous',
  .(
    tcr.id, tcr.chain, specificity.class.tag,
    donor.id.tag=donor.id, cell.subset,
    gen.clonotype.status, spc.fc=log2(spc.fc)
  )
]
# tmp.tholds <- tmp.data[, quantile(x=spc.fc, probs=c(0.6, 0.8))]
tmp.ggplot <- ggplot(data=tmp.data, aes(x=spc.fc)) +
    geom_density(linewidth=3.5, color='black', fill='gray', alpha=0) +
    geom_vline(xintercept=log2(2), color='red', linetype='dashed', linewidth=2) +
    scale_y_continuous(expand=expansion(add=c(0, 0)), breaks=scales::pretty_breaks(n=3)) +
    scale_x_continuous(expand=expansion(add=c(0, NA)), breaks=scales::pretty_breaks(n=3)) +
    labs(x='Specificity ratio (log10)', y='Emprirical density')
tmp.lab <- 'SpcFC_EmpDensity_Log'
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.prol.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.3, do.legend=FALSE,
    width=10, height=10
)

# ---> Number of both chain-matched clonotypes per specificity accordance type.
tmp.list <- list(
  'CD8'=pred.ref.1.l[['c40.lg']],
  'CD4'=pred.ref.2.l[['c40.lg']]
)
for(t.subset in names(tmp.list)){
  tmp.data <- tmp.list[[t.subset]]
  # Clonotype count across donors.
  tmp.data.1 <- tmp.data[
    !is.na(exp.pred) & exp.pred!='Multiple' & exp.chain.class=='BC',
    .(clonotype.abs.freq=uniqueN(qry.clone.id)),
    by=.(exp.match.class)
  ]
  tmp.ggplot <- ggplot(data=tmp.data.1, aes(x=exp.match.class, y=clonotype.abs.freq)) +
    geom_bar(stat='identity', color=cell.type.cols[t.subset], alpha=0, width=0.7, linewidth=5) +
    scale_y_continuous(expand=expansion(add=c(0, NA)), breaks=scales::pretty_breaks(n=3)) +
    labs(x='Match class', y='Clonotype count')
  tmp.lab <- paste0('ClonotypeCount_AccordanceClass_Across_CellType-', t.subset)
  publish.plot(
      tmp.ggplot=tmp.ggplot, output.path=fig.prol.path, file.name=tmp.lab, type='pdf',
      blank.comp=blank.complement.3, do.legend=FALSE,
      width=10, height=10
  )
}


############    -----------------------------------------    ############
### --------------------------- Figure on --------------------------- ###
### --------------- Class-II-restricted TCR resource ---------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
fig.resource.path <- paste0(reports.path, '/figure_on_resource')
if(!dir.exists(fig.resource.path)) dir.create(fig.resource.path)
# ---> Comments: Main dataset for this section, DICE-TCR cohort TCR data (the "resource").


### ------------------------- Text details -------------------------- ###

### ----------------> Class-I-restricted TCR resource
# ---> Number of clonotypes.
ref.1.meta[
    !is.na(beta.chain) &
    !is.na(specificity.class.tag),
    uniqueN(tcr.id)
]
# ---> Number of clonotypes w/ >1,000 entries.
tmp.thold <- 1000
ref.1.meta[
    !is.na(beta.chain) &
    !is.na(consensus.organism),
    .(freq.abs=uniqueN(tcr.id)),
    by=.(specificity.class.tag)
][freq.abs>tmp.thold]
# ---> Number of clonotypes w/ relevant specificities.
tmp.spcs <- c(
    'CMV', 'EBV',
    'SARS-CoV-2', 'IAV'
)
ref.1.meta[
    !is.na(beta.chain) &
    !is.na(consensus.organism) & consensus.organism%in%tmp.spcs,
    uniqueN(tcr.id)
]
# ---> Number of organisms covered by the compiled (class I-restricted TCR) reference overall.
ref.1.meta[, uniqueN(consensus.organism)]
# ---> 
tmp.data <- ref.1.meta[,
    .(freq.abs=uniqueN(data.base)),
    by=tcr.id
]
tmp.data[, .SD[freq.abs>1, .N]/.N]


### ----------------> Class-II-restricted TCR resource
# ---> Number of cells with both an available clonotype and a defined donor identity
ref.2.meta[!is.na(donor.id.tag), sum(total.abs.freq)]
# ---> Number of clonotypes.
ref.2.meta[!is.na(specificity.class.tag), uniqueN(clonotype.tag)]
# ---> Number and fraction of clonotypes for each general clonotype status.
ref.2.meta[,
    .(
        clone.abs.freq=uniqueN(clonotype.tag),
        clone.rel.freq=uniqueN(clonotype.tag) / ref.2.meta[, uniqueN(clonotype.tag)]
    ),
    by=gen.clonotype.status
]
# ---> Number and fraction of cells for each general clonotype status.
ref.2.meta[,
    .(
        clone.abs.freq=sum(total.abs.freq),
        clone.rel.freq=sum(total.abs.freq) / ref.2.meta[, sum(total.abs.freq)]
    ),
    by=gen.clonotype.status
]
# ---> Number of clonotypes per pathogen specificity for the median donor
tmp.data <- ref.2.meta[
    !is.na(specificity.class.tag),
    .(
        clone.abs.freq=uniqueN(clonotype.tag)
    ),
    by=.(donor.id.tag, specificity.class.tag)
]
tmp.data <- tmp.data[, .(clone.freq.stat=median(clone.abs.freq)), by=specificity.class.tag]
tmp.data[, median(clone.freq.stat)]
tmp.data[, range(clone.freq.stat)]
# ---> Resource overall size
ref.2.meta[
    !is.na(specificity.class.tag),
    uniqueN(clonotype.tag)
]
# ---> Stat for the number of different pathogen peptide pools covered per donor.
ref.2.meta[
    !is.na(specificity.class.tag),
    .(pp.abs.freq=uniqueN(specificity.class.tag)),
    by=donor.id.tag
][, .(
    min=min(pp.abs.freq),
    max=max(pp.abs.freq),
    median=median(pp.abs.freq)
)]


# ---> Overlap between individual resources.
merge.cols <- c('cdr3b.aa.seq', 'cdr3a.aa.seq', 'trb.v', 'trb.j', 'tra.v', 'tra.j')
tmp.data.1 <- ref.1.meta[
    specificity.class.tag %in% rel.spcs,
    .(
        cdr3b.aa.seq=beta.chain, cdr3a.aa.seq=alpha.chain,
        trb.v=trbv.base, trb.j=trbj.base,
        tra.v=trav.base, tra.j=traj.base,
        ag=as.character(specificity.class.tag)
    )
]
tmp.data.2 <- ref.2.meta[!is.na(specificity.class.tag), .(cdr3b.aa.seq, cdr3a.aa.seq, trb.v, trb.j, tra.v, tra.j, ag=as.character(specificity.class.tag))]
tmp.data <- merge(
    x=tmp.data.1, y=tmp.data.2,
    by=merge.cols,
    all=FALSE,
    suffixes=c('.1', '.2')
)
# Overlap. Fraction of new dataset.
tmp.data[, .N]*100/tmp.data.2[, .N]
# Accordance degree between references.
tmp.data[, .SD[ag.1==ag.2, .N]/.N]
tmp.data[, unique(ag.1)]
# Discordant overlapping TCRs.
tmp.data[ag.1!=ag.2]
# Overlap w/ validated TCRs.
tmp.data <- merge(
    x=tmp.data, y=val.res,
    by=merge.cols,
    all=FALSE,
    suffixes=c('.1', '.2')
)


### ---------------- Class-I-restricted TCR resource ---------------- ###

# ---> Define relevant specificities
# Now defined above. Might remove.
# rel.spcs <- c(
#     'CMV', 'EBV',
#     'IAV', 'SARS-CoV-2', 'MPV', 'PIV', 'RSV',
#     'B-pertussis-Vax', 'B-pertussis-Rest',
#     'Aspergillus'
# )

# ---> Number of clonotypes per reactivity.
# @ Including cross-reactive TCRs.
tmp.data <- ref.1.meta[
    !is.na(specificity.class.tag) &
    !is.na(beta.chain),
    .(tcr.id, specificity.class.tag, mhc.gene)
]
tmp.data <- unique(tmp.data)
tmp.ggplot <- ggplot(data=tmp.data, aes(x=specificity.class.tag)) +
    geom_bar(aes(fill=specificity.class.tag), color='black', width=0.7, linewidth=2.5) +
    scale_y_continuous(expand=expansion(add=c(0, 200)), breaks=scales::pretty_breaks(n=3)) +
    scale_fill_manual(values=pp.cols) +
    labs(x='Reactivity', y='Clonotype count', color='')
tmp.lab <- 'Ref-1_ClonotypeCount_Reactivity_Donor'
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.resource.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.3.3, do.legend=FALSE,
    width=15, height=10.5
)

# @ Plot only for relevant specificities.
tmp.data <- ref.1.meta[
    !is.na(specificity.class.tag) &
    !is.na(beta.chain),
    .(tcr.id, specificity.class.tag, mhc.gene)
]
tmp.data <- unique(tmp.data)
tmp.data <- tmp.data[specificity.class.tag %in% rel.spcs]
tmp.vals <- factor(x=tmp.data[, as.character(specificity.class.tag)], levels=rel.spcs)
set(x=tmp.data, j='specificity.class.tag', value=tmp.vals)
tmp.ggplot <- ggplot(data=tmp.data, aes(x=specificity.class.tag)) +
    geom_bar(aes(fill=specificity.class.tag), color='black', width=0.7, linewidth=3) +
    scale_x_discrete(drop=FALSE) +
    scale_y_continuous(expand=expansion(add=c(0, 200)), breaks=scales::pretty_breaks(n=3)) +
    scale_fill_manual(values=pp.cols) +
    labs(x='Reactivity', y='Clonotype count', color='')
tmp.lab <- 'Ref-1_ClonotypeCount_Reactivity_Relevant'
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.resource.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.3, do.legend=FALSE,
    width=21, height=10.15
)

# ---> Number of clonotypes per HLA gene.
# @ Including cross-reactive TCRs.
tmp.data <- ref.1.meta[
    !is.na(specificity.class.tag) & specificity.class.tag %in% rel.spcs &
    !is.na(beta.chain) &
    !is.na(mhc.gene),
    .(tcr.id, specificity.class.tag, mhc.gene)
]
tmp.data <- unique(tmp.data)
tmp.vals <- c(
  'HLA-A2', 'HLA-B7', 'HLA-E'
)
tmp.data[mhc.gene%in%tmp.vals, mhc.gene:='Class-I Rest']
# tmp.vals <- c(
#   'HLA-DQA1', 'HLA-DRB3', 'HLA-DPA1', 'HLA-DQ8', 'HLA-DRB4', 'HLA-DR11', 'HLA-DR3', 'HLA-DR5', 'HLA-DRB5'
# )
# tmp.data[mhc.gene%in%tmp.vals, mhc.gene:='Class-II Rest']
tmp.vals <- c(
  'HLA-A', 'HLA-B', 'HLA-C',
  'Class-I Rest'
)
tmp.data[!mhc.gene%in%tmp.vals, mhc.gene:='Class-II All']
tmp.lvls <- c(
  'HLA-A', 'HLA-B', 'HLA-C',
  'Class-I Rest',
  # 'HLA-DRA', 'HLA-DRB1', 'HLA-DPB1', 'HLA-DQB1',
  'Class-II All'
)
tmp.vals <- factor(x=tmp.data[, mhc.gene], levels=tmp.lvls)
set(x=tmp.data, j='mhc.gene', value=tmp.vals)
tmp.ggplot <- ggplot(data=tmp.data, aes(x=mhc.gene)) +
    geom_bar(color='black', fill='#d3d3d3', alpha=1, width=0.7, linewidth=2.5) +
    scale_y_continuous(expand=expansion(add=c(0, 200)), breaks=scales::pretty_breaks(n=4)) +
    labs(x='Reactivity', y='Clonotype count', color='')
tmp.lab <- 'Ref-1_ClonotypeCount_MHCGene'
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.resource.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.3, do.legend=FALSE,
    width=10, height=10
)

# ---> Comparison between class-I- and class-II-restricted DB entries
# @ Including cross-reactive TCRs.
tmp.data <- ref.1.meta[
    !is.na(specificity.class.tag) &
    !is.na(beta.chain) &
    specificity.class.tag %in% rel.spcs &
    !is.na(mhc.gene),
    .(tcr.id, mhc.gene)
]
tmp.data <- unique(tmp.data)
tmp.vals <- c(
  'HLA-A', 'HLA-B', 'HLA-C', 'HLA-A2', 'HLA-B7', 'HLA-E'
)
# Sanity check.
# tmp.data[
#     !mhc.gene %in% tmp.vals,
#     sort(unique(mhc.gene))
# ]]
tmp.data[mhc.gene%in%tmp.vals, mhc.rest:='Class-I']
tmp.data[!mhc.gene%in%tmp.vals, mhc.rest:='Class-II']
tmp.lvls <- c('Class-II', 'Class-I')
tmp.vals <- factor(x=tmp.data[, mhc.rest], levels=tmp.lvls)
set(x=tmp.data, j='mhc.rest', value=tmp.vals)
tmp.ggplot <- ggplot(data=tmp.data, aes(x=mhc.rest, fill=mhc.rest)) +
    geom_bar(color='black', alpha=1, width=0.7, linewidth=4) +
    scale_y_continuous(expand=expansion(add=c(0, 200)), breaks=scales::pretty_breaks(n=4)) +
    scale_fill_manual(values=mhc.rest.cols) +
    labs(x='Restriction', y='Clonotype count', color='')
tmp.lab <- 'Ref-1_ClonotypeCount_MHCRest'
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.resource.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.3.1, do.legend=FALSE,
    width=7, height=14
)


### --------------- Class-II-restricted TCR resource ---------------- ###

# ---> Number of clonotypes per reactivity across donors.
tmp.data <- unique(ref.2.meta[
    !is.na(specificity.class.tag),
    .(clonotype.tag, specificity.class.tag)
])
tmp.ggplot <- ggplot(data=tmp.data, aes(x=specificity.class.tag)) +
    geom_bar(aes(fill=specificity.class.tag), color='black', width=0.7, linewidth=3) +
    scale_y_continuous(expand=expansion(add=c(0, 200)), breaks=scales::pretty_breaks(n=3)) +
    scale_fill_manual(values=pp.cols) +
    labs(x='Reactivity', y='Clonotype count', color='')
tmp.lab <- 'Ref-2_ClonotypeCount_Reactivity_AcrossDonors'
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.resource.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.3, do.legend=FALSE,
    width=21, height=10.5
)

# ---> Number of clonotypes per donor per reactivity.
tmp.data <- ref.2.meta[!is.na(specificity.class.tag), .(donor.id.tag, clonotype.tag, specificity.class.tag)]
tmp.data <- tmp.data[,
    .(clonotype.abs.freq=.N),
    by=.(donor.id.tag, specificity.class.tag)
]
# Set upper bound.
tmp.thold <- 350
tmp.data[, bounded.abs.freq:=clonotype.abs.freq]
tmp.data[, point.status:='ori']
tmp.data[bounded.abs.freq>tmp.thold, `:=`(bounded.abs.freq=tmp.thold, point.status='bounded')]
# Reset levels.
tmp.vals <- factor(x=tmp.data[, as.character(specificity.class.tag)], levels=rel.spcs)
set(x=tmp.data, j='specificity.class.tag', value=tmp.vals)
# Plot.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=specificity.class.tag)) +
    geom_jitter(aes(y=bounded.abs.freq, shape=point.status), stroke=5, size=10, color=gen.dot.col, width=0.2, height=0) +
    geom_boxplot(aes(y=clonotype.abs.freq, color=specificity.class.tag, fill=specificity.class.tag), alpha=0.4, width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
    coord_cartesian(ylim=c(0, tmp.thold)) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
    scale_shape_manual(values=tmp.shapes) +
    scale_color_manual(values=pp.cols) +
    scale_fill_manual(values=pp.cols) +
    labs(x='Reactivity', y='Clonotype count', color='')
tmp.lab <- 'Ref-2_ClonotypeCount_Reactivity_Donor'
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.resource.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.3, do.legend=FALSE,
    width=28, height=14
)

# ---> TCR repertoire profile per donor and cell type. Plot that depicts clonal expansion info.
# Set order of clonotypes according to clone size.
tmp.data.1 <- ref.2.meta[
  !is.na(specificity.class.tag),
  .(
    donor.id.tag, clonotype.tag, clone.size=total.abs.freq, specificity.class.tag
  )
]
setorderv(x=tmp.data.1, cols=c('clone.size'), order=c(-1))
tmp.data.1[, order.tag:=paste(donor.id.tag, specificity.class.tag, clonotype.tag, sep=';')]
tmp.data.1$order.tag <- factor(x=tmp.data.1$order.tag, levels=tmp.data.1$order.tag)
# Plot.
tmp.ggplot <- ggplot(data=tmp.data.1, aes(x=donor.id.tag, y=clone.size)) +
    geom_bar(aes(group=order.tag, color=specificity.class.tag), stat='identity', position='stack', width=0.8, linewidth=1, alpha=0) +
    scale_y_continuous(expand=expansion(add=c(20, 0)), breaks=scales::pretty_breaks(n=4)) +
    scale_color_manual(values=pp.cols) +
    labs(x='Donor ID', y='Cell count', color='')
tmp.lab <- '/Ref-2_StackedSizedClones_Donor_Reactivity'
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.resource.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=30, height=15
)

# ---> Bioinformatics validation of reference.
# @ Perform statistics from SARS perspective.
tmp.data <- lapply(X=1:pub.motif[, .N], FUN=function(tmp.row){
    # Retrieve motif info.
    tmp.pttn <- pub.motif[tmp.row, motif]
    to.rep <- paste0('[RNDCQEGHILKMFPSTWYV]')
    tmp.pttn <- str_replace(string=tmp.pttn, pattern='\\*', replacement=to.rep)
    tmp.chain <- pub.motif[tmp.row, chain]
    tmp.col.1 <- ifelse(test=tmp.chain=='TRA', yes='cdr3a.aa.seq', no='cdr3b.aa.seq')
    tmp.col.2 <- ifelse(test=tmp.chain=='TRA', yes='tra.v', no='trb.v')
    tmp.col.3 <- ifelse(test=tmp.chain=='TRA', yes='tra.j', no='trb.j')
    # Retrieve vdj gene info.
    tmp.trv <- pub.motif[tmp.row, trv]
    if(is.na(tmp.trv)){
        tmp.trv <- ref.2.meta[, unique(get(tmp.col.2))]
    }else{
        tmp.trv <- str_split(string=tmp.trv, pattern=';')[[1]]
    }
    tmp.trj <- pub.motif[tmp.row, trj]
    if(is.na(tmp.trj)){
        tmp.trj <- ref.2.meta[, unique(get(tmp.col.3))]
    }else{
        tmp.trj <- str_split(string=tmp.trj, pattern=';')[[1]]
    }
    # Perform motif search per clonotype in reference.
    tmp.data.1 <- ref.2.meta[
        !is.na(specificity.class.tag) &
        get(tmp.col.1) %like% tmp.pttn &
        get(tmp.col.2) %in% tmp.trv &
        get(tmp.col.3) %in% tmp.trj,
        .N,
        # .(N=sum(total.abs.freq)),
        by=specificity.class.tag
    ]
    # Create contingency table for Fisher's test
    tmp.data.2 <- ref.2.meta[
        !is.na(specificity.class.tag) &
        !(
            get(tmp.col.1) %like% tmp.pttn &
            get(tmp.col.2) %in% tmp.trv &
            get(tmp.col.3) %in% tmp.trj)
        ,
        .N,
        by=specificity.class.tag
    ]
    tmp.data.1 <- tmp.data.1[,
        .(
            class='is.int.y',
            is.sars.y=.SD[specificity.class.tag=='SARS-CoV-2', sum(N)],
            is.sars.n=.SD[specificity.class.tag!='SARS-CoV-2', sum(N)]
        )
    ]
    tmp.data.2 <- tmp.data.2[,
        .(
            class='is.int.n',
            is.sars.y=.SD[specificity.class.tag=='SARS-CoV-2', sum(N)],
            is.sars.n=.SD[specificity.class.tag!='SARS-CoV-2', sum(N)]
        )
    ]
    tmp.data <- rbind(tmp.data.1, tmp.data.2)
    # Sanity check
    tmp.check <- tmp.data[, sum(is.sars.y)] + tmp.data[, sum(is.sars.n)] == ref.2.meta[!is.na(specificity.class.tag), .N]
    if(!tmp.check) stop('Unexpected error while creating contigency table.\n')
    # Perform test.
    for.test <- as.matrix(tmp.data[, .(is.sars.y, is.sars.n)])
    row.names(for.test) <- tmp.data$class
    tmp.p <- fisher.test(x=for.test)$p.value
    # Collect pieces of info.
    tmp.data <- data.table(
        int.sars=tmp.data[class=='is.int.y', is.sars.y],
        int.other=tmp.data[class=='is.int.y', is.sars.n],
        fisher.p=tmp.p
    )
    return(tmp.data)
})
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE)
tmp.data <- cbind(pub.motif, tmp.data)
tmp.data[, motif.id:=paste(motif.id, chain, sep=' ')]

# @ Set general motif order.
# Keep only motifs coming from specific study.
tmp.data <- tmp.data[pmid=='35841887']
# Set order
setorderv(x=tmp.data, cols=c('ref.specificity', 'motif.id'), na.last=TRUE)
# Set factor.
tmp.lvls <- tmp.data$motif.id
tmp.data[, motif.id:=factor(x=motif.id, levels=tmp.lvls)]

# @ Viz.
to.plot.1 <- tmp.data[,
    .(
        motif.id,
        `SARS-CoV-2`=int.sars, Other=int.other
    )
]
to.plot.1 <- as.data.table(gather(data=to.plot.1, key='int.type', value='int.size', -motif.id))
# to.plot.1 <- to.plot.1[int.size>0]
to.plot.2 <- tmp.data[, .(motif.id, fisher.p=round(fisher.p, digits=3))]
tmp.capt <- paste0('Respectively, P values are:', to.plot.2[, paste0(fisher.p, collapse=', ')])
# Plot
tmp.cols <- c(pp.cols['SARS-CoV-2'], `Other`='#CBCBCB')
tmp.ggplot <- ggplot(data=to.plot.1, aes(x=motif.id)) +
    geom_bar(aes(y=int.size, fill=int.type), stat='identity', color='black', width=0.7, linewidth=2.5) +
    # geom_text(data=to.plot.2, aes(y=to.plot.1[, max(int.size)+1], label=fisher.p)) +
    scale_y_continuous(expand=expansion(add=c(0, 3)), breaks=scales::pretty_breaks(n=4)) +
    scale_fill_manual(values=tmp.cols) +
    labs(x='Motif', y='Intersected clonotypes', fill='Type of\nintersect.', caption=tmp.capt)
tmp.lab <- 'Ref-2_Val-InSil-1'
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.resource.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.3.3, do.legend=FALSE,
    width=21, height=14
)

# ---> Comparison between class-I- and class-II-restricted clonotypes across resources
# @ Data from DBs
tmp.data.1 <- ref.1.meta[
    !is.na(specificity.class.tag) & specificity.class.tag %in% rel.spcs &
    !is.na(mhc.gene),
    .(tcr.id, mhc.gene)
]
tmp.data.1 <- unique(tmp.data.1)
tmp.vals <- c(
  'HLA-A', 'HLA-B', 'HLA-C', 'HLA-A2', 'HLA-B7', 'HLA-E'
)
# Sanity check.
# tmp.data[
#     !mhc.gene %in% tmp.vals,
#     sort(unique(mhc.gene))
# ]]
tmp.data.1[mhc.gene%in%tmp.vals, mhc.rest:='Class-I']
tmp.data.1[!mhc.gene%in%tmp.vals, mhc.rest:='Class-II']
# @ Data from in-house resource.
tmp.data.2 <- ref.2.meta[
    !is.na(specificity.class.tag),
    .(tcr.id=clonotype.tag, mhc.rest='Class-II')
]
tmp.data.2 <- unique(tmp.data.2)
# @ Summarize
tmp.data <- rbindlist(l=list(`Public`=tmp.data.1, `Generated`=tmp.data.2), use.names=TRUE, fill=TRUE, idcol='source')
tmp.lvls <- c('Class-II', 'Class-I')
tmp.vals <- factor(x=tmp.data[, mhc.rest], levels=tmp.lvls)
set(x=tmp.data, j='mhc.rest', value=tmp.vals)
# tmp.data[source!='Generated', source:=as.character(mhc.rest)]
# tmp.lvls <- c('Class-II', 'Class-I', 'Generated')
tmp.lvls <- c('Generated', 'Public')
tmp.vals <- factor(x=tmp.data[, source], levels=tmp.lvls)
set(x=tmp.data, j='source', value=tmp.vals)
# Plot
# tmp.cols <- c(mhc.rest.cols, `Generated`='#FFC500')
tmp.cols <- c(`Public`='#C3B1E1', `Generated`='#B5EAD7')
tmp.ggplot <- ggplot(data=tmp.data, aes(x=mhc.rest)) +
    geom_bar(aes(fill=source), color='black', alpha=1, width=0.7, linewidth=4) +
    scale_y_continuous(expand=expansion(add=c(0, 200)), breaks=scales::pretty_breaks(n=4)) +
    scale_fill_manual(values=tmp.cols) +
    labs(x='Restriction', y='Clonotype count', fill='Source')
tmp.lab <- 'Ref-1-2_ClonotypeCount_MHCRest'
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.resource.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.3, do.legend=FALSE,
    width=7, height=14
)


### -------------------------- Supp Figure -------------------------- ###

# ---> Number of clonotypes per general clonotype status.
# ref.2.meta[is.na(gen.clonotype.status), summary(total.abs.freq)]
tmp.data <- ref.2.meta[,
  .(
    tcr.id=clonotype.tag, specificity.class.tag,
    donor.id.tag,
    gen.clonotype.status, spc.fc
  )
]
tmp.data[is.na(gen.clonotype.status), gen.clonotype.status:='Singleton']
tmp.vals <- factor(x=tmp.data[, gen.clonotype.status], levels=c('Singleton', 'Unambiguous', 'Ambiguous'))
set(x=tmp.data, j='gen.clonotype.status', value=tmp.vals)
tmp.data.1 <- tmp.data[,
    .(clonotype.abs.freq=uniqueN(tcr.id)),
    by=.(donor.id.tag, gen.clonotype.status)
]
tmp.thold <- 13000
tmp.ggplot <- ggplot(data=tmp.data.1, aes(x=gen.clonotype.status, y=clonotype.abs.freq)) +
    geom_boxplot(width=0.5, linewidth=4, fatten=4, outlier.shape=NA) +
    # geom_jitter(size=5, color='black', width=0.25) +
    coord_cartesian(ylim=c(NA, tmp.thold)) + # Custom close up.
    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
    labs(x='Reactivity', y='Clonotype count', color='')
tmp.lab <- 'Ref-2_ChainTotalCount_CloneaStatus'
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.resource.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.3, do.legend=FALSE,
    width=14, height=14
)

# ---> Specificity ratio
# @ Linea scale.
tmp.data <- ref.2.meta[
  gen.clonotype.status=='Ambiguous',
  .(
    tcr.id=clonotype.tag, specificity.class.tag,
    donor.id.tag,
    gen.clonotype.status, spc.fc
  )
]
tmp.ggplot <- ggplot(data=tmp.data, aes(x=spc.fc)) +
    geom_density(linewidth=3.5, color='black', fill='gray', alpha=0) +
    geom_vline(xintercept=2, color='red', linetype='dashed', linewidth=2) +
    scale_y_continuous(expand=expansion(add=c(0, 0)), breaks=scales::pretty_breaks(n=3)) +
    scale_x_continuous(expand=expansion(add=c(0, NA)), breaks=scales::pretty_breaks(n=3)) +
    coord_cartesian(xlim=c(NA, tmp.tholds[3])) + # Custom close up.
    labs(x='Specificity ratio', y='Emprirical density')
tmp.lab <- 'Ref-2_SpcFC_EmpDensity_Linear'
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.resource.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.3, do.legend=FALSE,
    width=10, height=10
)
# @ Log scale.
tmp.data <- ref.2.meta[
  gen.clonotype.status=='Ambiguous',
  .(
    tcr.id=clonotype.tag, specificity.class.tag,
    donor.id.tag,
    gen.clonotype.status, spc.fc=log2(spc.fc)
  )
]
tmp.ggplot <- ggplot(data=tmp.data, aes(x=spc.fc)) +
    geom_density(linewidth=3.5, color='black', fill='gray', alpha=0) +
    geom_vline(xintercept=log2(2), color='red', linetype='dashed', linewidth=2) +
    scale_y_continuous(expand=expansion(add=c(0, 0)), breaks=scales::pretty_breaks(n=3)) +
    scale_x_continuous(expand=expansion(add=c(0, NA)), breaks=scales::pretty_breaks(n=3)) +
    labs(x='Specificity ratio (log10)', y='Emprirical density')
tmp.lab <- 'Ref-2_SpcFC_EmpDensity_Log'
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.resource.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.3, do.legend=FALSE,
    width=10, height=10
)
    geom_density(linewidth=3.5, color='black', fill='gray', alpha=0) +
    geom_vline(xintercept=log2(2), color='red', linetype='dashed', linewidth=2) +


############    -----------------------------------------    ############
### --------------------------- Figure on --------------------------- ###
### ----------------- TCR UCM systematic comparison ----------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
fig.ucm.path <- paste0(reports.path, '/figure_on_ucm_comp')
if(!dir.exists(fig.ucm.path)) dir.create(fig.ucm.path)
# ---> Comments: Main dataset for this section, TCR UCM comparison results, "ucm.comp.res"


### ------------------------- Text details -------------------------- ###
tmp.data <- fetch.data(full.report=ucm.comp.res, res.type='pred.results', item.name='auc')
tmp.data[,
    .(auc.median=median(auc)),
    by=.(tool.name, reactivity)
][, .(auc.median=median(auc.median)), by=tool.name]


### -------------------------- Main Figure -------------------------- ###

# ---> Brief summary on comparison results.
# Comparison between tools when considering the full dataset, even with singletons.
# Fetch AUC summaries for full dataset single trial per peptide pool.
tmp.data <- fetch.data(full.report=ucm.comp.res, res.type='pred.results', item.name='auc')
tmp.data <- tmp.data[size.thold==0]
tmp.data[, tool.name:=factor(x=tool.name, levels=tcr.ucms)]
# Plot.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=tool.name, y=auc)) +
    geom_boxplot(aes(fill=tool.name), color='black', alpha=0.5, width=0.7, linewidth=5, outlier.shape=NA) +
    geom_jitter(aes(color=reactivity), size=8, width=0.25) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
    scale_color_manual(values=pp.cols) +
    scale_fill_manual(values=ucm.cols) +
    coord_cartesian(ylim=c(0.3, 1)) + # Custom close up.
    labs(x='TCR TCM', y='Model AUC', color='Peptide pool', fill='Tool')
tmp.lab <- 'AUC_UCM_Spc_RefSet-Whole'
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.ucm.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.3.1, do.legend=FALSE,
    width=14, height=9.7
)
# Comparison between tools when disregarding singletons across different rel. freq. values.
# Fetch AUC summaries for data subset w/ curated/assigned antigen specificities.
tmp.data <- fetch.data(full.report=ucm.comp.res, res.type='pred.results', item.name='auc')
# NOTE: TEMPORARY FIX TO ALLOW TCRDIST3 TO BE SHOWN ON THE PLOT AND MAKE IT CLEAR THAT WE DON'T KNOW AS OF NOW WHY THERE AREN'T ANY RESULTS.
tmp.data[tool.name=='TCRdist3', `:=`(size.thold=100, auc=0)]
tmp.data[, tool.name:=factor(x=tool.name, levels=tcr.ucms)]
tmp.data <- tmp.data[size.thold!=0]
# Plot.
tmp.ggplot <- ggplot(data=tmp.data, aes(x=tool.name, y=auc)) +
    geom_boxplot(aes(fill=tool.name), color='black', alpha=0.5, width=0.7, linewidth=3, outlier.shape=NA) +
    geom_jitter(aes(color=reactivity), size=8, width=0.25) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
    scale_color_manual(values=pp.cols) +
    scale_fill_manual(values=ucm.cols) +
    coord_cartesian(ylim=c(0.3, 1)) + # Custom close up.
    labs(x='TCR TCM', y='Model AUC', color='Peptide pool', fill='Tool')
tmp.lab <- 'AUC_UCM_Spc_RefSet-Subset'
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.ucm.path, file.name=tmp.lab, type='pdf',
    blank.comp=blank.complement.3, do.legend=FALSE,
    width=14, height=10
)


### -------------------------- Supp Figure -------------------------- ###

# ---> Extended summary in terms of AUCs
# Fetch AUC summaries (general)
tmp.data <- fetch.data(full.report=ucm.comp.res, res.type='pred.results', item.name='auc')
tmp.data[, tool.name:=factor(x=tool.name, levels=tcr.ucms)]
# Plot per peptide pool.
tmp.vals <- tmp.data[, unique(reactivity)]
for(pp in tmp.vals){
    to.plot <- tmp.data[reactivity==pp]
    tmp.ggplot <- ggplot(data=to.plot, aes(x=size.thold, y=auc, color=tool.name)) +
        geom_point(size=14) +
        geom_line(aes(group=tool.name), linewidth=4.5, linetype='dashed') +
        scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
        scale_color_manual(values=ucm.cols) +
        coord_cartesian(ylim=c(0.3, 1)) + # Custom close up.
        labs(x='Threshold as the percentile of the rel. freq. distribution', y='Model AUC', color='TCR UCM')
    tmp.lab <- paste0('AUC_FreqThold_Spc-', pp)
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.ucm.path, file.name=tmp.lab, type='pdf',
        blank.comp=blank.complement.3, do.legend=FALSE,
        width=10, height=10
    )
}

# ---> Extended summary in terms of ROC curves
# Fetch AUC summaries only for size threshold of 1.
tmp.data <- fetch.data(full.report=ucm.comp.res, res.type='pred.results', item.name='metrics')
tmp.data[, tool.name:=factor(x=tool.name, levels=tcr.ucms)]
tmp.data <- tmp.data[size.thold==1]
# Plot per peptide pool.
tmp.vals <- tmp.data[, unique(reactivity)]
for(pp in tmp.vals){
    to.plot <- unique(tmp.data[reactivity==pp])
    tmp.ggplot <- ggplot(data=to.plot, aes(x=specificity, y=sensitivity, color=tool.name)) +
        geom_abline(slope=1, intercept=1, linewidth=4, color='black') +
        geom_line(linewidth=4) +
        scale_color_manual(values=ucm.cols) +
        scale_x_reverse(limits=c(1, 0), expand=c(0.01, 0.01), breaks=scales::pretty_breaks(n=3)) +
        # scale_x_reverse() +
        scale_y_continuous(limits=c(0, 1), expand=c(0.01, 0.01), breaks=scales::pretty_breaks(n=3)) +
        labs(x='Specificity', y='Sensitivity', color='TCR UCM') +
        theme_bw() + theme(legend.position='bottom')
    tmp.lab <- paste0('Sensitivity_Specificity_UCM_Spc-', pp)
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.ucm.path, file.name=tmp.lab, type='pdf',
        blank.comp=blank.complement.3, do.legend=FALSE,
        width=10, height=10
    )
}   


# @ Fetch AUC summaries.
tmp.data <- fetch.data(full.report=ucm.comp.res, res.type='pred.results', item.name='auc')
# @ Fetch ROC curve summaries.
tmp.data <- fetch.data(full.report=ucm.comp.res, res.type='pred.results', item.name='metrics')


############    -----------------------------------------    ############
### --------------------------- Figure on --------------------------- ###
### ---------------------- Integrative strategy --------------------- ###
############    -----------------------------------------    ############

# ----> Output directory.
fig.int.strat.path <- paste0(reports.path, '/figure_on_int_strat')
if(!dir.exists(fig.int.strat.path)) dir.create(fig.int.strat.path)


### ------------------------- Text details -------------------------- ###

for(data.set in names(aggr.set.labs)){
    cat(aggr.set.labs[data.set], '\n')
    tmp.list <- list(
        'CD8'=pred.ref.1.l[[data.set]],
        'CD4'=pred.ref.2.l[[data.set]]
    )
    tmp.data <- lapply(X=names(tmp.list), FUN=function(t.subset){
        tmp.data <- tmp.list[[t.subset]]
        tmp.data.1 <- tmp.data[
            !is.na(match.pred) & match.pred!='Multiple',
            (
                match.fract=.SD[
                    match.chain.class=='IA' | match.chain.class=='IB' | match.chain.class=='BC',
                    uniqueN(qry.clone.id)
                ] /
                uniqueN(qry.clone.id)
            )
        ]
        return(tmp.data.1)
    })
    tmp.data <- unlist(tmp.data)
    names(tmp.data) <- names(tmp.list)
    cat(paste(names(tmp.data), tmp.data, sep=': ', collapse='\n'), '\n')
}

# ----> Fraction of cells with an assigned specificity.
for(data.set in names(aggr.set.labs)){
    cat(aggr.set.labs[data.set], '\n')
    tmp.data <- meta.data.l[[data.set]][
        !is.na(clonotype.tag),
        .(
            spc.cell.no=.SD[!is.na(ag.spc.tag), .N],
            spc.cell.fract=.SD[!is.na(ag.spc.tag), .N]/.N
        ),
        by=.(gex.lin.tag, donor.id.tag)
    ]
    to.output <- tmp.data[, .(stat=sum(spc.cell.no)), by=gex.lin.tag]
    cat('Total no. of cells w/ assigned specificity\n')
    cat(paste(to.output[, gex.lin.tag], to.output[, stat], sep=': ', collapse='\n'), '\n')
    to.output <- tmp.data[, .(stat=median(spc.cell.fract)), by=gex.lin.tag]
    cat('Median\n')
    cat(paste(to.output[, gex.lin.tag], to.output[, stat], sep=': ', collapse='\n'), '\n')
    cat('General range\n')
    cat(tmp.data[, range(spc.cell.fract)], '\n')
    cat('\n')
}


### -------------------------- Main Figure -------------------------- ###

# ----> Direct comparison between strategies.
# For experimental prediction only and for the integrative approach
tmp.apps <- c(`Integrative`='consensus.pred', `Experimental`='exp.pred')
# tmp.tholds <- list(
#     `c40.lg`=c(`Integrative`=0.35, `Experimental`=0.27),
#     `c40.ln`=c(`Integrative`=0.15, `Experimental`=0.05)
# )
tmp.tholds <- list(
    `c40.lg`=c(`CD4`=0.22, `CD8`=0.35),
    `c40.ln`=c(`CD4`=0.11, `CD8`=0.13)
)
for(data.set in names(aggr.set.labs)){
    for(cell.lin in names(cell.type.cols)){
        # app.col <- tmp.apps[tmp.app]
        tmp.thold <- tmp.tholds[[data.set]][cell.lin]
        tmp.data <- if(cell.lin=='CD8'){
            pred.ref.1.l[[data.set]][, .(donor.id.tag, full.size, consensus.pred, exp.pred)]
        }else{
            pred.ref.2.l[[data.set]][, .(donor.id.tag, full.size, consensus.pred, exp.pred)]
        }
        tmp.data[consensus.pred=='Multiple', consensus.pred:=NA]
        tmp.data[exp.pred=='Multiple', exp.pred:=NA]
        tmp.data <- tmp.data[,
            .(
                freq.total=sum(full.size),
                freq.abs.exp=.SD[!is.na(exp.pred), sum(full.size)],
                freq.abs.int=.SD[!is.na(consensus.pred), sum(full.size)],
                freq.rel.exp=.SD[!is.na(exp.pred), sum(full.size)]/sum(full.size),
                freq.rel.int=.SD[!is.na(consensus.pred), sum(full.size)]/sum(full.size)
            ),
            by=.(donor.id.tag)
        ]
        tmp.data <- tmp.data[, .(donor.id.tag, `Experimental`=freq.rel.exp, `Integrative`=freq.rel.int)]
        tmp.data <- as.data.table(gather(data=tmp.data, key='approach', value='rel.freq', -`donor.id.tag`))
        # Set upper bound.
        tmp.data[, bound.freq:=rel.freq]
        tmp.data[, point.status:='ori']
        tmp.data[bound.freq>tmp.thold, `:=`(bound.freq=tmp.thold, point.status='bounded')]
        # @ Plot.
        tmp.capt <- tmp.data[,
            .(
                median=round(x=median(rel.freq), digits=5),
                min=round(x=min(rel.freq), digits=5),
                max=round(x=max(rel.freq), digits=5)
            ),
            by=.(approach)
        ][, paste(
            paste0('Median ', paste(approach, median, sep=': '), collapse='. '),
            paste0('Min ', paste(approach, min, sep=': '), collapse='. '),
            paste0('Max ', paste(approach, max, sep=': '), collapse='. '),
            sep='\n'
        )]
        tmp.ggplot <- ggplot(data=tmp.data, aes(x=approach)) +
            geom_line(aes(group=donor.id.tag, y=bound.freq), color=gen.link.col, linewidth=gen.link.width) +
            geom_jitter(aes(y=bound.freq, shape=point.status), stroke=5, size=10, color=gen.dot.col, width=0, height=0) +
            geom_boxplot(aes(y=rel.freq, color=approach, fill=approach), alpha=0.4, width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
            coord_cartesian(ylim=c(0, tmp.thold)) +
            scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
            scale_shape_manual(values=tmp.shapes) +
            scale_color_manual(values=app.cols) +
            scale_fill_manual(values=app.cols) +
            labs(x='Strategy', y='Fraction of cells w/ assigned specificity', color='', caption=tmp.capt)
        tmp.lab <- paste0(
            '/', aggr.set.labs[data.set], '_StrategyComp_Lin-', cell.lin
        )
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.int.strat.path, file.name=tmp.lab, type='pdf',
            blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=7, height=14
        )
    }
}


### -------------------------- Supp Figure -------------------------- ###

# ---> Number of clonotypes per chain class of match between TCR reference and query TCRs.
tmp.list <- list(
  'CD8'=pred.ref.1.l[['c40.lg']],
  'CD4'=pred.ref.2.l[['c40.lg']]
)
tmp.breaks <- list(
    `CD4`=c(600, 7000),
    `CD8`=c(500, 1400)
)
for(t.subset in names(tmp.list)){
  tmp.data <- tmp.list[[t.subset]]
  # Clonotype count per donor.
  tmp.data.1 <- tmp.data[
    !is.na(match.pred) & match.pred!='Multiple',
    .(clonotype.abs.freq=uniqueN(qry.clone.id)),
    by=.(donor.id.tag, match.chain.class)
  ]
  tmp.ggplot <- ggplot(data=tmp.data.1, aes(x=match.chain.class, y=clonotype.abs.freq)) +
    geom_boxplot(color=cell.type.cols[t.subset], width=0.7, linewidth=5, outlier.shape=NA) +
    geom_jitter(size=5, color='black', width=0.25) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
    labs(x='Match class', y='Clonotype count')
  tmp.lab <- paste0('ClonotypeCount_MatchClass_Donor_CellType-', t.subset)
  publish.plot(
      tmp.ggplot=tmp.ggplot, output.path=fig.int.strat.path, file.name=tmp.lab, type='pdf',
      blank.comp=blank.complement.3, do.legend=FALSE,
      width=21, height=15
  )
  # Clonotype count across donors.
  tmp.data.1 <- tmp.data[
    !is.na(match.pred) & match.pred!='Multiple',
    .(clonotype.abs.freq=uniqueN(qry.clone.id)),
    by=.(match.chain.class)
  ]
  tmp.ggplot <- ggplot(data=tmp.data.1, aes(x=match.chain.class, y=clonotype.abs.freq)) +
    geom_bar(stat='identity', color=cell.type.cols[t.subset], alpha=0, width=0.7, linewidth=5) +
    # geom_col(color=cell.type.cols[t.subset], alpha=0, width=0.7, linewidth=5) +
    scale_y_continuous(expand=expansion(add=c(0, NA)), breaks=scales::pretty_breaks(n=3)) +
    # scale_y_break(expand=expansion(add=c(0, NA)), breaks=tmp.breaks[[t.subset]]) +
    labs(x='Match class', y='Clonotype count')
  tmp.lab <- paste0('ClonotypeCount_MatchClass_Across_CellType-', t.subset)
  publish.plot(
      tmp.ggplot=tmp.ggplot, output.path=fig.int.strat.path, file.name=tmp.lab, type='pdf',
      blank.comp=blank.complement.3, do.legend=FALSE,
      width=21, height=15
  )
}

# ---> Number of both chain-matched clonotypes per match specificity accordance class.
tmp.list <- list(
  'CD8'=pred.ref.1.l[['c40.lg']],
  'CD4'=pred.ref.2.l[['c40.lg']]
)
for(t.subset in names(tmp.list)){
  tmp.data <- tmp.list[[t.subset]]
  # Clonotype count across donors.
  tmp.data.1 <- tmp.data[
    !is.na(match.pred) & match.pred!='Multiple' & (match.chain.class=='IA' | match.chain.class=='IB' | match.chain.class=='BC'),
    .(clonotype.abs.freq=uniqueN(qry.clone.id)),
    by=.(match.match.class)
  ]
  tmp.ggplot <- ggplot(data=tmp.data.1, aes(x=match.match.class, y=clonotype.abs.freq)) +
    geom_bar(stat='identity', color=cell.type.cols[t.subset], alpha=0, width=0.7, linewidth=5) +
    scale_y_continuous(expand=expansion(add=c(0, NA)), breaks=scales::pretty_breaks(n=3)) +
    labs(x='Match class', y='Clonotype count')
  tmp.lab <- paste0('ClonotypeCount_MatchMatchClass_Across_CellType-', t.subset)
  publish.plot(
      tmp.ggplot=tmp.ggplot, output.path=fig.int.strat.path, file.name=tmp.lab, type='pdf',
      blank.comp=blank.complement.3, do.legend=FALSE,
      width=14, height=10
  )
}

# ---> Number of clonotypes per general approach.
tmp.list <- list(
  'CD8'=pred.ref.1.l[['c40.lg']],
  'CD4'=pred.ref.2.l[['c40.lg']]
)
for(t.subset in names(tmp.list)){
  tmp.data <- tmp.list[[t.subset]]
  tmp.data <- tmp.data[
    !is.na(consensus.approach.class), 
    .(
        qry.clone.id, donor.id.tag, full.size, consensus.gen.class, consensus.approach.class, clust.pred
    )
  ]
  # Clonotype count across donors.
  tmp.data.1 <- tmp.data[,
    .(
        clonotype.abs.freq=.N,
        cell.abs.freq=sum(full.size)
    ),
    by=.(approach=consensus.approach.class)
  ]
  tmp.cols <- c(
    `Cell`='cell.abs.freq',
    `Clonotype`='clonotype.abs.freq'
  )
  for(tmp.lab in names(tmp.cols)){
    tmp.var <- tmp.cols[tmp.lab]
    tmp.ggplot <- ggplot(data=tmp.data.1, aes(x=approach, y=get(tmp.var))) +
        geom_bar(stat='identity', color=cell.type.cols[t.subset], alpha=0, width=0.7, linewidth=5) +
        scale_y_continuous(expand=expansion(add=c(0, NA)), breaks=scales::pretty_breaks(n=3)) +
        labs(x='Approach', y='Clonotype count')
    tmp.lab <- paste0(tmp.lab, 'Count_GenApproach_Across_CellType-', t.subset)
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.int.strat.path, file.name=tmp.lab, type='pdf',
        blank.comp=blank.complement.3, do.legend=FALSE,
        width=15, height=15
    )
  }
}


############    -----------------------------------------    ############
### --------------------------- Figure on --------------------------- ###
### ------------------- Antigen-specific T cells -------------------- ###
### ------------ Repertoire comparisons between lineages ------------ ###
############    -----------------------------------------    ############

# ----> Output directory.
fig.as.rep.comps.path <- paste0(reports.path, '/figure_on_ag-spc_rep_comps')
if(!dir.exists(fig.as.rep.comps.path)) dir.create(fig.as.rep.comps.path)
# ---> Comments: Main dataset for this section, CD4 T cell set.
main.obj.name <- NULL

### ------------------------- Text details -------------------------- ###

# ---> Fetch data.
tmp.data.1 <- meta.data.l[['c40.lg']][
    !is.na(clonotype.tag) & !is.na(ag.spc.tag),
    .(
        cell.count=.N,
        clone.count=uniqueN(clonotype.tag)
    ),
    by=.(gex.lin.tag, donor.id.tag, ag.spc.tag)
]
tmp.data.2 <- meta.data.l[['c40.lg']][
    !is.na(clonotype.tag),
    .(total.count=.N),
    by=.(gex.lin.tag, donor.id.tag)
]
tmp.data.1 <- merge(x=tmp.data.1, y=tmp.data.2, by=c('gex.lin.tag', 'donor.id.tag'), all.x=TRUE, all.y=FALSE)
tmp.data.1[, cell.fract:=cell.count/total.count]
# @ Spread to capture a donor-wide distribution.
tmp.data.2 <- tmp.data.1[, .(gex.lin.tag, donor.id.tag, ag.spc.tag, cell.fract)]
tmp.data.2 <- spread(data=tmp.data.2, key=ag.spc.tag, value=cell.fract, fill=0)
tmp.data.2 <- as.data.table(gather(data=tmp.data.2, key=ag.spc.tag, value=cell.fract, -`gex.lin.tag`, -`donor.id.tag`))
# @ Range of pathogen-specific cell fractions.
tmp.data.2[, round(x=range(cell.fract)*100, digits=3)]
# @ Number of pathogens w/ uncovered memory responses to different numbers of donors. A response is counted as 10 cells.
to.check <- tmp.data.1[
    cell.count>10,
    .(donor.count=uniqueN(donor.id.tag)),
    by=.(gex.lin.tag, ag.spc.tag)
]
setorderv(x=to.check, col='donor.count', order=-1)
# @ Median frequency across donors.
to.check <- tmp.data.2[,
    .(median.cell.fract=mean(cell.fract)),
    by=.(gex.lin.tag, ag.spc.tag)
]
to.check[gex.lin.tag=='CD8' & ag.spc.tag=='CMV']
tmp.vals <- c('EBV', 'IAV', 'SARS-CoV-2')
tmp.data.2[
    gex.lin.tag=='CD8' & ag.spc.tag %in% tmp.vals,
    sum(cell.fract),
    by=donor.id.tag
][, mean(V1)]
# @ Clonotype frequency per pool.
tmp.vals <- c('CMV', 'EBV', 'IAV', 'SARS-CoV-2')
to.check.1 <- tmp.data.1[
    ag.spc.tag %in% tmp.vals,
    .(count=sum(clone.count)),
    by=.(gex.lin.tag, ag.spc.tag)
]
to.check.2 <- tmp.data.1[
    ag.spc.tag %in% tmp.vals,
    .(total=sum(clone.count)),
    by=.(gex.lin.tag)
]
to.check <- merge(x=to.check.1, y=to.check.2, by='gex.lin.tag')
to.check[gex.lin.tag=='CD4', count/total]
tmp.vals <- c('CMV', 'EBV', 'IAV', 'SARS-CoV-2')
tmp.data.1[
    ag.spc.tag %in% tmp.vals,
    .SD[ag.spc.tag=='CMV', sum(clone.count)]/sum(clone.count),
    by=gex.lin.tag
]
# @ 
tmp.data.1[
    ag.spc.tag=='CMV' & 
    cell.fract<0.01,
    uniqueN(donor.id.tag)
]

# ---> Clonality details.
tmp.data <- meta.data.l[['c40.lg']][
    !is.na(ag.spc.tag),
    .(abs.freq=uniqueN(barcode)),
    by=.(gex.lin.tag, donor.id.tag, clonotype.tag)
]
tmp.data[abs.freq>2, summary(abs.freq)]


### -------------------------- Main Figure -------------------------- ###

# ---> Preflights.
subset.rel.spcs <- list(
    `c40.lg.cd4.d40`=rel.spcs,
    `c40.lg.cd8.d40`=c('CMV', 'EBV', 'SARS-CoV-2', 'IAV'),
    `c40.lg.cd8.trm`=c('CMV', 'EBV', 'SARS-CoV-2', 'IAV'),
    `c40.lg.cd8.tem`=c('CMV', 'EBV', 'SARS-CoV-2', 'IAV'),
    `c40.ln.cd4.d08`=rel.spcs,
    `c40.ln.cd8.d08`=c('CMV', 'EBV', 'SARS-CoV-2', 'IAV')
)

# ---> Prevalence of pathogen-specific T cells in the lung. Donor-specific info.
# Define lineage-specific color breaks.
tmp.col.breaks <- list(
    `c40.lg.cd4.d40`=c(0, 0.02, 0.04),
    `c40.lg.cd8.d40`=c(0.05, 0.1, 0.15),
    `c40.ln.cd4.d08`=c(0.01, 0.02),
    `c40.ln.cd8.d08`=c(0, 0.04, 0.08)
)
tmp.donor.ords <- list(
    `c40.lg.cd4.d40`=custom.donor.ord.2,
    `c40.lg.cd8.d40`=custom.donor.ord.3,
    `c40.ln.cd4.d08`=custom.donor.ord.1,
    `c40.ln.cd8.d08`=custom.donor.ord.1
)
for(data.set in names(aggr.set.labs)){
    for(t.subset in meta.data.l[[data.set]][, unique(data.set)]){
        # @ Data fetch
        # Fetch data to plot.
        tmp.data <- meta.data.l[[data.set]][data.set==t.subset]
        tmp.data.1 <- tmp.data[
            !is.na(donor.id.tag) & !is.na(clonotype.tag) & !is.na(ag.spc.tag),
            .(cell.count=.N),
            by=.(donor.id.tag, ag.spc.tag)
        ]
        tmp.data.2 <- tmp.data[
            !is.na(donor.id.tag) & !is.na(clonotype.tag),
            .(total.count=.N),
            by=.(donor.id.tag)
        ]
        tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='donor.id.tag', all.x=TRUE, all.y=FALSE)
        tmp.data[, cell.fract:=cell.count/total.count]
        tmp.data <- tmp.data[, .(donor.id.tag, ag.spc.tag, cell.fract)]
        tmp.thold <- tmp.data[, quantile(x=cell.fract, probs=0.95)]
        if(tmp.thold<0.05) tmp.thold <- 0.05
        tmp.data[cell.fract>tmp.thold, cell.fract:=tmp.thold]
        tmp.data <- spread(data=tmp.data, key=ag.spc.tag, value=cell.fract, drop=FALSE)
        row.names(tmp.data) <- tmp.data$donor.id.tag; tmp.data$donor.id.tag <- NULL
        tmp.data <- as.matrix(tmp.data)
        # @ Relevant specificities.
        tmp.data <- tmp.data[,
            subset.rel.spcs[[t.subset]]
        ]
        # tmp.data <- tmp.data[rev(row.names(tmp.data)), ]
        # @ Set donor order
        custom.donor.ord <- tmp.donor.ords[[t.subset]]
        tmp.vals <- custom.donor.ord[custom.donor.ord %in% row.names(tmp.data)]
        tmp.data <- tmp.data[tmp.vals, ]
        # Set metadata.
        col.metadata <- data.frame(
            row.names=colnames(tmp.data),
            `Specificity`=colnames(tmp.data)
        )
        # Define colors for metadata tracks.
        tmp.cols <- pp.cols[names(pp.cols) %in% colnames(tmp.data)]
        ann.colors <- list(
            `Specificity`=tmp.cols
        )
        # Set color scale and breaks for heatmap.
        break.no <- 200
        col.breaks <- seq(from=min(range(tmp.data, na.rm=TRUE)), to=max(range(tmp.data, na.rm=TRUE)), length.out=break.no)
        mid.point <- .005
        mid.point <- which.min(abs(col.breaks-mid.point))
        # hmap.col.scale.1 <- colorRampPalette(c('#ffffff', '#ffff9a'))(mid.point)
        # hmap.col.scale.2 <- colorRampPalette(c('#ffff9a', '#ffd34d', '#ffc04d', '#ffa500', '#ff6500', '#ff4d4d', '#ff0000', '#b30000', '#670000'))(break.no-(mid.point+1))
        hmap.col.scale.1 <- colorRampPalette(c('#414487', '#2A788E'))(mid.point)
        hmap.col.scale.2 <- colorRampPalette(c('#2A788E', '#22A884', '#FDE725', '#FFEA00'))(break.no-(mid.point+1))
        hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
        # File specs.
        tmp.height <- (nrow(tmp.data) * 6.68/40) + 0.32
        tmp.width.c <- (ncol(tmp.data) * 7)/10
        tmp.width.b <- (ncol(tmp.data) * 1.8)/10
        # @ Version with hierarchical clustering.
        # To complete only if possible.
        tmp.check <- any(c(
            any(rowSums(is.na(tmp.data))==ncol(tmp.data)),
            any(colSums(is.na(tmp.data))==nrow(tmp.data))
        ))
        if(!tmp.check){
            # Complete version.
            tmp.file.name <- paste0(fig.as.rep.comps.path, '/', aggr.set.labs[data.set], '_', obj.extended.names[t.subset], '_Donor_Reactivity_CellFract_Clust.C.pdf')
            pheatmap(
                mat=tmp.data, scale='none',
                color=hmap.col.scale, breaks=col.breaks, na_col='#FFFFFF', border_color='black',
                cluster_rows=TRUE, cluster_cols=TRUE,
                annotation_col=col.metadata, annotation_colors=ann.colors,
                show_colnames=FALSE, show_rownames=TRUE,
                legend=TRUE, annotation_legend=TRUE, annotation_names_row=TRUE, annotation_names_col=TRUE,
                filename=tmp.file.name, heigh=tmp.height, width=tmp.width.c
            )
            # Blank version.
            tmp.file.name <- paste0(fig.as.rep.comps.path, '/', aggr.set.labs[data.set], '_', obj.extended.names[t.subset],'_Donor_Reactivity_CellFract_Clust.B.pdf')
            pheatmap(
                mat=tmp.data, scale='none',
                color=hmap.col.scale, breaks=col.breaks, na_col='#FFFFFF', border_color='black',
                cluster_rows=TRUE, cluster_cols=TRUE,
                annotation_col=col.metadata, annotation_colors=ann.colors,
                show_colnames=FALSE, show_rownames=FALSE,
                legend=FALSE, annotation_legend=FALSE, annotation_names_row=FALSE, annotation_names_col=FALSE,
                filename=tmp.file.name, heigh=tmp.height, width=tmp.width.b+1
            )
        }

        # @ Version w/out hierarchical clustering.
        # Complete version.
        tmp.file.name <- paste0(fig.as.rep.comps.path, '/', aggr.set.labs[data.set], '_', obj.extended.names[t.subset], '_Donor_Reactivity_CellFract_NoClust.C.pdf')
        pheatmap(
            mat=tmp.data, scale='none',
            color=hmap.col.scale, breaks=col.breaks, na_col='#FFFFFF', border_color='black',
            cluster_rows=FALSE, cluster_cols=FALSE,
            annotation_col=col.metadata, annotation_colors=ann.colors,
            show_colnames=FALSE, show_rownames=TRUE,
            legend=TRUE, annotation_legend=TRUE, annotation_names_row=TRUE, annotation_names_col=TRUE,
            filename=tmp.file.name, heigh=tmp.height, width=tmp.width.c
        )
        # Blank version.
        tmp.file.name <- paste0(fig.as.rep.comps.path, '/', aggr.set.labs[data.set], '_', obj.extended.names[t.subset], '_Donor_Reactivity_CellFract_NoClust.B.pdf')
        pheatmap(
            mat=tmp.data, scale='none',
            color=hmap.col.scale, breaks=col.breaks, na_col='#FFFFFF', border_color='black',
            cluster_rows=FALSE, cluster_cols=FALSE,
            annotation_col=col.metadata, annotation_colors=ann.colors,
            show_colnames=FALSE, show_rownames=FALSE,
            legend=FALSE, annotation_legend=FALSE, annotation_names_row=FALSE, annotation_names_col=FALSE,
            filename=tmp.file.name, heigh=tmp.height, width=tmp.width.b
        )

        # @ Color legend only.
        tmp.data <- as.data.frame(tmp.data)
        tmp.data$donor.id.tag <- row.names(tmp.data)
        tmp.data <- gather(data=tmp.data, key='spc', value='freq.rel', -`donor.id.tag`)
        tmp.data <- tmp.data[!is.na(tmp.data$freq.rel), ]
        # W/ labels
        tmp.ggplot <- ggplot(data=tmp.data, aes(x=freq.rel, y=freq.rel, col=freq.rel)) + 
        # tmp.ggplot <- ggplot(data=as.data.frame(tmp.data), aes(x=`SARS-CoV-2`, y=`SARS-CoV-2`, col=`SARS-CoV-2`)) + 
            geom_point() +
            # scale_color_gradientn(colors=hmap.col.scale, breaks=col.breaks[c(65, 130, 195)]) +
            scale_color_gradientn(colors=hmap.col.scale, breaks=tmp.col.breaks[[t.subset]]) +
            theme(
                legend.ticks=element_line(color='black', linewidth=0.6),
                legend.ticks.length=unit(0.22, "cm"),
                legend.frame=element_rect(color='black', linewidth=0.6)
            )
        tmp.ggplot <- get_legend(p=tmp.ggplot)
        tmp.file.name <- paste0(fig.as.rep.comps.path, '/', aggr.set.labs[data.set], '_', obj.extended.names[t.subset],'_Donor_Reactivity_CellFract.L1.pdf')
        pdf(file=tmp.file.name, height=2, width=2)
        print(as_ggplot(tmp.ggplot))
        dev.off()
        # W/out labels
        tmp.ggplot <- ggplot(data=as.data.frame(tmp.data), aes(x=freq.rel, y=freq.rel, col=freq.rel)) + 
        # tmp.ggplot <- ggplot(data=as.data.frame(tmp.data), aes(x=`SARS-CoV-2`, y=`SARS-CoV-2`, col=`SARS-CoV-2`)) + 
            geom_point() +
            # scale_color_gradientn(colors=hmap.col.scale, breaks=col.breaks[c(65, 130, 195)]) +
            scale_color_gradientn(colors=hmap.col.scale, breaks=tmp.col.breaks[[t.subset]], name=NULL, labels=NULL) +
            theme(
                legend.ticks=element_line(color='black', linewidth=0.6),
                legend.ticks.length=unit(0.22, "cm"),
                legend.frame=element_rect(color='black', linewidth=0.6)
            )
        tmp.ggplot <- get_legend(p=tmp.ggplot)
        tmp.file.name <- paste0(fig.as.rep.comps.path, '/', aggr.set.labs[data.set], '_', obj.extended.names[t.subset],'_Donor_Reactivity_CellFract.L2.pdf')
        pdf(file=tmp.file.name, height=1.25, width=0.3)
        print(as_ggplot(tmp.ggplot))
        dev.off()
    }
}

# ---> Prevalence of pathogen-specific T cells in the lung. Summary across donors
# tmp.bounds <- c(
#     `c40.lg.cd4.d40`=0.042,
#     `c40.lg.cd8.d40`=0.175,
#     `c40.ln.cd4.d08`=0.02,
#     `c40.ln.cd8.d08`=0.065
# )
tmp.bounds <- c(
    `c40.lg.cd4.d40`=0.06,
    `c40.lg.cd8.d40`=0.21,
    `c40.ln.cd4.d08`=0.02,
    `c40.ln.cd8.d08`=0.065
)
for(data.set in names(aggr.set.labs)){
    for(t.subset in meta.data.l[[data.set]][, unique(data.set)]){
        # @ Data fetch
        # Fetch data to plot.
        tmp.data <- meta.data.l[[data.set]][data.set==t.subset]
        tmp.data.1 <- tmp.data[
            !is.na(ag.spc.tag),
            .(cell.count=.N),
            by=.(donor.id.tag, ag.spc.tag=as.character(ag.spc.tag))
        ]
        tmp.data.2 <- tmp.data[
            !is.na(clonotype.tag),
            .(total.count=.N),
            by=.(donor.id.tag)
        ]
        tmp.data.1 <- merge(x=tmp.data.1, y=tmp.data.2, by='donor.id.tag', all.x=TRUE, all.y=FALSE)
        tmp.data.1[, cell.fract:=cell.count/total.count]
        # Spread to capture a donor-wide distribution.
        tmp.data.1 <- tmp.data.1[, .(donor.id.tag, ag.spc.tag, cell.fract)]
        tmp.data.1 <- spread(data=tmp.data.1, key=ag.spc.tag, value=cell.fract, fill=0)
        tmp.data.1 <- as.data.table(gather(data=tmp.data.1, key=ag.spc.tag, value=cell.fract, -`donor.id.tag`))
        # Relevant specificities.
        tmp.data.1 <- tmp.data.1[ag.spc.tag %in% subset.rel.spcs[[t.subset]]]
        # Determine basic stats.
        tmp.capt.1 <- tmp.data.1[, paste0(
            'Median: ', round(median(cell.fract), 5), '. ',
            'Min: ', round(min(cell.fract), 5), '. ',
            'Max: ', round(max(cell.fract), 5)
        )]
        tmp.pps <- c('IAV', 'SARS-CoV-2', 'RSV', 'PIV', 'MPV')
        tmp.capt.2 <- tmp.data.1[
            ag.spc.tag %in% tmp.pps,
            paste0(
                'Median: ', round(median(cell.fract), 5), '. ',
                'Min: ', round(min(cell.fract), 5), '. ',
                'Max: ', round(max(cell.fract), 5)
            )
        ]
        tmp.pps <- c('CMV', 'EBV')
        tmp.capt.3 <- tmp.data.1[
            ag.spc.tag %in% tmp.pps,
            paste0(
                'Median: ', round(median(cell.fract), 5), '. ',
                'Min: ', round(min(cell.fract), 5), '. ',
                'Max: ', round(max(cell.fract), 5)
            )
        ]
        tmp.pps <- c('B-pertussis-Vax', 'B-pertussis-Rest')
        tmp.capt.4 <- tmp.data.1[
            ag.spc.tag %in% tmp.pps,
            paste0(
                'Median: ', round(median(cell.fract), 5), '. ',
                'Min: ', round(min(cell.fract), 5), '. ',
                'Max: ', round(max(cell.fract), 5)
            )
        ]
        tmp.pps <- c('Aspergillus')
        tmp.capt.5 <- tmp.data.1[
            ag.spc.tag %in% tmp.pps,
            paste0(
                'Median: ', round(median(cell.fract), 5), '. ',
                'Min: ', round(min(cell.fract), 5), '. ',
                'Max: ', round(max(cell.fract), 5)
            )
        ]
        tmp.caption <- paste0(
            'All. ', tmp.capt.1, '\n',
            'Resp. viruses. ', tmp.capt.2, '\n',
            'Herpesviruses. ', tmp.capt.3, '\n',
            'Pertussis. ', tmp.capt.4, '\n',
            'Aspergillus. ', tmp.capt.5
        )
        # Add CMV status.
        tmp.data.2 <- tmp.data.1[ag.spc.tag=='CMV']
        tmp.data.2 <- merge(
            x=tmp.data.2, y=ser.cmv.data[, .(donor.id.tag, cmv.status)],
            by='donor.id.tag', all.x=TRUE, all.y=FALSE
        )
        tmp.data.2[, `:=`(
            ag.spc.tag=paste(ag.spc.tag, cmv.status, sep='.'),
            cmv.status=NULL, sero=TRUE
        )]
        tmp.data <- rbind(tmp.data.1, tmp.data.2, fill=TRUE)
        tmp.data[is.na(sero), sero:=FALSE]
        # Set upper bound.
        tmp.data[, bounded.cell.fract:=cell.fract]
        tmp.data[, point.status:='ori']
        tmp.bound <- tmp.bounds[t.subset]
        tmp.data[cell.fract>tmp.bound, `:=`(bounded.cell.fract=tmp.bound, point.status='bounded')]
        # Set factors.
        tmp.vals <- factor(x=tmp.data[, ag.spc.tag], levels=c(rel.spcs, 'CMV.Positive', 'CMV.Negative'))
        set(x=tmp.data, j='ag.spc.tag', value=tmp.vals)
        # ---> Across pathogens, only pathogens.
        tmp.vals <- c('CMV.Positive', 'CMV.Negative')
        tmp.width <- ((0.8*14/0.53) * length(subset.rel.spcs[[t.subset]]))/10
        # @ Plot.
        tmp.ggplot <- ggplot(data=tmp.data[!ag.spc.tag %in% tmp.vals], aes(x=ag.spc.tag)) +
            geom_jitter(aes(y=bounded.cell.fract, shape=point.status), stroke=4, size=8, color=gen.dot.col, width=0.2, height=0) +
            geom_boxplot(aes(color=ag.spc.tag, fill=ag.spc.tag, y=cell.fract), alpha=0.4, width=0.7, linewidth=3, fatten=4, outlier.shape=NA) +
            # geom_hline(yintercept=0.01, linewidth=2, linetype='dashed', color='black') +
            coord_cartesian(ylim=c(0, tmp.bound)) +
            # scale_x_discrete(drop=FALSE) +
            scale_x_discrete(drop=TRUE) +
            scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
            scale_shape_manual(values=tmp.shapes) +
            scale_color_manual(values=pp.cols) +
            scale_fill_manual(values=pp.cols) +
            labs(x='Reactivity', y='% of T cells', color='', caption=tmp.caption)
        tmp.lab <- paste0('/', aggr.set.labs[data.set], '_', obj.extended.names[t.subset], '_CellFract_Reactivity_Donor_Summ')
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.as.rep.comps.path, file.name=tmp.lab, type='pdf',
            blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=tmp.width, height=14
        )
        # ---> Across pathogens, including CMV serostatus.
        tmp.cols <- c(pp.cols, `CMV.Positive`=unname(pp.cols['CMV']), `CMV.Negative`=unname(pp.cols['CMV']))
        tmp.width <- (30 * length(subset.rel.spcs[[t.subset]]))/10
        if(t.subset=='c40.lg.cd8.d40') tmp.width <- tmp.width + 2
        # @ Plot.
        tmp.ggplot <- ggplot(data=tmp.data, aes(x=ag.spc.tag)) +
            geom_jitter(aes(y=bounded.cell.fract, shape=point.status), stroke=4, size=8, color=gen.dot.col, width=0.2, height=0) +
            geom_boxplot(aes(color=ag.spc.tag, fill=ag.spc.tag, y=cell.fract), alpha=0.4, width=0.7, linewidth=3, fatten=4, outlier.shape=NA) +
            geom_vline(xintercept=-Inf, color='black', linewidth=5) +
            # geom_hline(yintercept=0.01, linewidth=2, linetype='dashed', color='black') +
            facet_grid(facets=~sero, scales='free_x', space='free_x') +
            coord_cartesian(ylim=c(0, tmp.bound)) +
            scale_x_discrete(drop=TRUE) +
            scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
            scale_shape_manual(values=tmp.shapes) +
            scale_color_manual(values=tmp.cols) +
            scale_fill_manual(values=tmp.cols) +
            labs(x='Reactivity', y='% of T cells', color='', caption=tmp.caption) +
            theme(panel.spacing=unit(5, "lines"))   # Increase space between grids.
        tmp.lab <- paste0('/', aggr.set.labs[data.set], '_', obj.extended.names[t.subset], '_CellFract_Reactivity_Donor_Summ_CSS')
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.as.rep.comps.path, file.name=tmp.lab, type='pdf',
            blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=tmp.width, height=14
        )
        # ---> Between pertussis vaccine and non-vaccine antigens.
        tmp.vals <- c('B-pertussis-Vax', 'B-pertussis-Rest')
        to.plot <- tmp.data[ag.spc.tag %in% tmp.vals]
        if(nrow(to.plot)==0) next
        to.test <- spread(data=to.plot[, .(donor.id.tag, ag.spc.tag, cell.fract)], key=ag.spc.tag, value=cell.fract)
        tmp.test <- wilcox.test(x=to.test[, tmp.vals[1]], y=to.test[, tmp.vals[2]], paired=TRUE)
        tmp.caption <- paste0('P-value is ', tmp.test$p.value)
        tmp.ggplot <- ggplot(data=to.plot, aes(x=ag.spc.tag)) +
            geom_line(aes(y=cell.fract, group=donor.id.tag), linewidth=gen.link.width, color=gen.link.col) +
            geom_jitter(aes(y=bounded.cell.fract, shape=point.status), stroke=4, size=8, color=gen.dot.col, width=0, height=0) +
            geom_boxplot(aes(color=ag.spc.tag, fill=ag.spc.tag, y=cell.fract), alpha=0.4, width=0.7, linewidth=3, fatten=4, outlier.shape=NA) +
            # coord_cartesian(ylim=c(0, tmp.bound)) +
            scale_x_discrete(drop=TRUE) +
            scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
            scale_shape_manual(values=tmp.shapes) +
            scale_color_manual(values=pp.cols) +
            scale_fill_manual(values=pp.cols) +
            labs(x='Reactivity', y='% of T cells', color='', caption=tmp.caption)
        tmp.lab <- paste0('/', aggr.set.labs[data.set], '_', obj.extended.names[t.subset], '_CellFract_Reactivity-Sel-1_Donor_Summ')
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.as.rep.comps.path, file.name=tmp.lab, type='pdf',
            blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=7, height=14
        )
    }
}

# ---> Prevalence of pathogen-specific T cells in the lung. Comparison for enrichment between CMV-reactive cells and cells reactive to other major viruses.
# Define specificities that were well-represented in both resources.
tmp.spcs <- c(
    'CMV', 'EBV',
    'SARS-CoV-2', 'IAV'
)
for(data.set in names(aggr.set.labs)){
    for(t.subset in meta.data.l[[data.set]][, unique(data.set)]){
        tmp.data <- unique(meta.data.l[[data.set]][
            data.set==t.subset &
            ag.spc.tag %in% tmp.spcs,
            .(clonotype.tag, ag.spc.tag)
        ])
        tmp.ggplot <- ggplot(data=tmp.data, aes(x='', fill=ag.spc.tag)) +
            geom_bar(width=1, color='black', linewidth=1) +
            scale_fill_manual(values=pp.cols) +
            coord_polar('y', start=0) +
            labs(x='', y='', fill='Reactivity') +
            theme_void()
        tmp.lab <- paste0('/', aggr.set.labs[data.set], '_', obj.extended.names[t.subset], '_CellFractOfAllPreds_Reactivity_Summ')
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.as.rep.comps.path, file.name=tmp.lab, type='pdf',
            blank.comp=blank.complement.4, do.legend=FALSE, do.rotate=TRUE#,
            # width=15, height=5
        )
    }
}

# ---> Number of cells per donor per clonotype, bins filles according to reactivity.
for(data.set in names(aggr.set.labs)){
    for(t.subset in meta.data.l[[data.set]][, unique(data.set)]){
        tmp.data <- meta.data.l[[data.set]][
            data.set==t.subset
        ]
        # Stacked clone info.
        tmp.data.1 <- tmp.data[
            !is.na(donor.id.tag) & !is.na(ag.spc.tag) &
            ag.spc.tag %in% subset.rel.spcs[[t.subset]],
            .(clone.size=.N),
            by=.(
                donor.id.tag, clonotype.tag, consensus.pred=ag.spc.tag
            )
        ]
        to.order <- 1:length(pp.cols); names(to.order) <- names(pp.cols)
        tmp.data.1[, order.tag:=to.order[consensus.pred]]
        setorderv(x=tmp.data.1, cols=c('order.tag', 'clone.size'), order=c(1, -1))
        tmp.data.1[, order.tag:=paste(donor.id.tag, clonotype.tag, consensus.pred, sep=';')]
        tmp.data.1$order.tag <- factor(x=tmp.data.1$order.tag, levels=tmp.data.1$order.tag)
        # Number of unique clonotypes per donor and cell type.
        tmp.data.2 <- tmp.data[
            !is.na(donor.id.tag) & !is.na(ag.spc.tag) &
            ag.spc.tag %in% subset.rel.spcs[[t.subset]],
            .(clone.count=uniqueN(clonotype.tag)),
            by=.(donor.id.tag)
        ]
        tmp.data.2[clone.count<100, clone.count:=NA]
        # To force all donors to appear on the plot.
        if(data.set=='c40.lg'){
            tmp.data.3 <- data.table(
                donor.id.tag=levels(tmp.data.2$donor.id.tag)[!levels(tmp.data.2$donor.id.tag) %in% tmp.data.2$donor.id.tag]
            )
            tmp.data.1 <- rbindlist(
                l=list(tmp.data.1, tmp.data.3),
                use.names=TRUE, fill=TRUE
            )
        }
        # @ Set donor order
        custom.donor.ord <- tmp.donor.ords[[t.subset]]
        tmp.lvls <- rev(custom.donor.ord[custom.donor.ord %in% tmp.data.1$donor.id.tag])
        tmp.vals <- factor(x=as.character(tmp.data.1$donor.id.tag), levels=tmp.lvls)
        set(x=tmp.data.1, j='donor.id.tag', value=tmp.vals)
        # Define file width according to donor amount.
        # tmp.width <- tmp.data.1[, uniqueN(donor.id.tag)]*14.62/40 + 0.38
        # tmp.width <- tmp.data[, uniqueN(donor.id.tag)] * 15/40
        tmp.width <- (tmp.data.1[, uniqueN(donor.id.tag)] * 6.68/40) + 0.32
        tmp.width <- (tmp.width * 5)/1.8
        # Plot.
        tmp.ggplot <- ggplot(data=tmp.data.1, aes(x=donor.id.tag)) +
            geom_bar(aes(y=clone.size, group=order.tag, color=consensus.pred), stat='identity', position='stack', width=0.8, linewidth=1.4, alpha=0) +
            # geom_point(data=tmp.data.2, aes(y=clone.count), size=7, shape=21, fill='white', color='black') +
            scale_color_manual(values=pp.cols) +
            scale_y_continuous(expand=c(0, 0), breaks=scales::pretty_breaks(n=3)) +
            labs(x='Donor ID', y='Cell count', fill='Reactivity')
        tmp.lab <- paste0('/', aggr.set.labs[data.set], '_', obj.extended.names[t.subset], '_CellCount_Donor_ConsReactivity-StackedClones')
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.as.rep.comps.path, file.name=tmp.lab, type='pdf',
            blank.comp=blank.complement.7.2, do.legend=FALSE, do.rotate=TRUE, width=tmp.width, height=5
        )
    }
}

# ---> Percentage of clones accounted for by each pathogen among all cells w/ specificity assignments. Summary across donors.
subset.rel.spcs <- list(
    `c40.lg.cd4.d40`=rel.spcs,
    `c40.lg.cd8.d40`=c('CMV', 'EBV', 'IAV', 'SARS-CoV-2'),
    `c40.lg.cd8.trm`=c('CMV', 'EBV', 'IAV', 'SARS-CoV-2'),
    `c40.lg.cd8.tem`=c('CMV', 'EBV', 'IAV', 'SARS-CoV-2'),
    `c40.ln.cd4.d08`=rel.spcs,
    `c40.ln.cd8.d08`=c('CMV', 'EBV', 'IAV', 'SARS-CoV-2')
)
for(data.set in names(aggr.set.labs)){
    tmp.vals <- meta.data.l[[data.set]][, unique(data.set)]
    for(t.subset in tmp.vals){
        tmp.data <- meta.data.l[[data.set]][
            data.set==t.subset &
            !is.na(ag.spc.tag) &
            ag.spc.tag %in% subset.rel.spcs[[t.subset]],
        ]
        tmp.data[, clone.id:=paste(clonotype.tag, donor.id.tag, ';')]
        tmp.vals <- c('Cell', 'Clone')
        for(tmp.item in tmp.vals){
            if(tmp.item=='Cell'){
                to.plot <- tmp.data[,
                    .(freq.abs=.N),
                    by=.(cluster=clusters.tag)
                ]
            }else{
                to.plot <- tmp.data[,
                    .(freq.abs=uniqueN(clone.id)),
                    by=.(cluster=clusters.tag)
                ]
            }
            tmp.ggplot <- ggplot(data=to.plot, aes(x='', y=freq.abs, fill=cluster)) +
                geom_bar(stat='identity', linewidth=0.5, color='black') +
                coord_polar(theta='y') +
                scale_fill_manual(values=clusts.cols[[t.subset]]) +
                labs(x='', y='', fill='Specificity') +
                theme_void()
            tmp.lab <- paste0('/', obj.extended.names[t.subset], '_', tmp.item, 'Fract_Cluster-Summ')
            publish.plot(
                tmp.ggplot=tmp.ggplot, output.path=fig.as.rep.comps.path, file.name=tmp.lab, type='pdf',
                blank.comp=blank.complement.4, do.legend=FALSE, do.rotate=FALSE, width=14, height=14
            )
        }
    }
}

# ---> Percentage of clones accounted for by each pathogen among selected subset's cells w/ specificity assignments. Summary across donors.
# Preflights.
tmp.reports.path <- paste0(fig.as.rep.comps.path, '/pop-wise_pat_fract_summ')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
subset.rel.spcs <- list(
    `c40.lg.cd4.d40`=rel.spcs,
    `c40.lg.cd8.d40`=c('CMV', 'EBV', 'IAV', 'SARS-CoV-2'),
    `c40.lg.cd8.trm`=c('CMV', 'EBV', 'IAV', 'SARS-CoV-2'),
    `c40.lg.cd8.tem`=c('CMV', 'EBV', 'IAV', 'SARS-CoV-2'),
    `c40.ln.cd4.d08`=rel.spcs,
    `c40.ln.cd8.d08`=c('CMV', 'EBV', 'IAV', 'SARS-CoV-2')
)
# Identify ourliers in CMV-specific CD4 T cell phenotypic landscape.
tmp.data <- meta.data.l[['c40.lg']][
    gex.lin.tag=='CD4' & clusters.tag=='2' & !is.na(ag.spc.tag)
]
tmp.data.1 <- tmp.data[,
    .(freq.abs=.N),
    by=.(donor.id.tag, ag.spc.tag)
]
tmp.data.2 <- tmp.data[,
    .(freq.total=.N),
    by=.(donor.id.tag)
]
tmp.data <- merge(
    x=tmp.data.1, y=tmp.data.2,
    by='donor.id.tag'
)
tmp.data[, freq.rel:=freq.abs/freq.total]
tmp.donors <- tmp.data[ag.spc.tag=='MPV' & freq.rel>0.5, donor.id.tag]
all.donors <- meta.data.l[['c40.lg']][, unique(donor.id.tag)]
tmp.donors <- setdiff(x=all.donors, y=tmp.donors)
donor.subsets <- lapply(X=subset.rel.spcs, FUN=function(x) return(list(`All`=all.donors)))
donor.subsets[['c40.lg.cd4.d40']] <- c(donor.subsets[['c40.lg.cd4.d40']], list(`S1`=tmp.donors))
# Process per T-cell dataset
tmp.subsets <- c(gen.subsets, list(`gzmk.tag`=c('GZMKhi')))
for(data.set in names(aggr.set.labs)){
    tmp.vals <- meta.data.l[[data.set]][, unique(data.set)]
    for(t.subset in tmp.vals){
        for(tmp.subset in names(tmp.subsets)){
            for(donor.subset in names(donor.subsets[[t.subset]])){
                tmp.val <- tmp.subsets[[tmp.subset]][1]
                tmp.data <- meta.data.l[[data.set]][
                    data.set==t.subset &
                    !is.na(ag.spc.tag) &
                    ag.spc.tag %in% subset.rel.spcs[[t.subset]] &
                    donor.id.tag %in% donor.subsets[[t.subset]][[donor.subset]] &
                    get(tmp.subset)==tmp.val,
                ]
                tmp.data[, clone.id:=paste(clonotype.tag, donor.id.tag, ';')]
                if(tmp.data[, .N]==0) next
                tmp.vals <- c('Cell', 'Clone')
                for(tmp.item in tmp.vals){
                    if(tmp.item=='Cell'){
                        to.plot <- tmp.data[,
                            .(freq.abs=.N),
                            by=.(ag=ag.spc.tag)
                        ]
                    }else{
                        to.plot <- tmp.data[,
                            .(freq.abs=uniqueN(clone.id)),
                            by=.(ag=ag.spc.tag)
                        ]
                    }
                    tmp.ggplot <- ggplot(data=to.plot, aes(x='', y=freq.abs, fill=ag)) +
                        geom_bar(stat='identity', linewidth=0.5, color='black') +
                        coord_polar(theta='y') +
                        scale_fill_manual(values=pp.cols) +
                        labs(x='', y='', fill='Specificity') +
                        theme_void()
                    tmp.lab <- paste0('/', obj.extended.names[t.subset], '_', tmp.item, 'Fract_Pat-Summ_', tmp.val, '_DonorSet-', donor.subset)
                    publish.plot(
                        tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
                        blank.comp=blank.complement.4, do.legend=FALSE, do.rotate=FALSE, width=14, height=14
                    )
                }
            }
        }
    }
}

# ---> Prevalence associations between anatomical sites.
# Retrieve cell fractions per lineage and anatomical site.
tmp.data <- lapply(X=names(aggr.set.labs), FUN=function(data.set){
    tmp.vals <- meta.data.l[[data.set]][, unique(data.set)]
    tmp.data <- lapply(X=tmp.vals, FUN=function(t.subset){
        # @ Data fetch
        # Fetch data to plot.
        tmp.data <- meta.data.l[[data.set]][data.set==t.subset]
        tmp.data.1 <- tmp.data[
            !is.na(donor.id.tag) & !is.na(clonotype.tag) & !is.na(ag.spc.tag),
            .(cell.count=.N),
            by=.(donor.id.tag, ag.spc.tag)
        ]
        tmp.data.2 <- tmp.data[
            !is.na(donor.id.tag) & !is.na(clonotype.tag),
            .(total.count=.N),
            by=.(donor.id.tag)
        ]
        tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='donor.id.tag', all.x=TRUE, all.y=FALSE)
        tmp.data[, cell.fract:=cell.count/total.count]
        tmp.data <- tmp.data[, .(donor.id.tag, ag.spc.tag, cell.fract)]
        tmp.data <- tmp.data[donor.id.tag %in% site.donor.ovlp]
        return(tmp.data)
    })
    names(tmp.data) <- gen.cell.types[tmp.vals]
    tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='cell.lin')
    return(tmp.data)

})
names(tmp.data) <- aggr.simp.labs[names(aggr.set.labs)]
tmp.data <- merge(
    x=tmp.data[['LG']], y=tmp.data[['LN']],
    by=c('cell.lin', 'donor.id.tag', 'ag.spc.tag'),
    all=TRUE, suffixes=c('.lg', '.ln')
)
tmp.data[is.na(cell.fract.lg), cell.fract.lg:=0]
tmp.data[is.na(cell.fract.ln), cell.fract.ln:=0]
tmp.data[, tmp.var:='All']
for(tmp.lin in names(cell.type.cols)){
    to.plot <- tmp.data[cell.lin==tmp.lin]
    tmp.ggplot <- ggplot(data=to.plot, aes(x=cell.fract.lg, y=cell.fract.ln)) +
        geom_point(aes(color=ag.spc.tag), shape=1, stroke=5, size=10) +
        geom_smooth(method='lm', formula=y~x, se=FALSE, color='black', linewidth=4) +
        scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
        scale_color_manual(values=pp.cols) +
        labs(
            x='Lung cell fraction',
            y='LN cell fraction'
        )
    tmp.lab <- paste0('/SimBetSites_CellFract_Lin-', tmp.lin)
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.as.rep.comps.path, file.name=tmp.lab, type='pdf',
        stat.cor=TRUE, cor.group='tmp.var',
        blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=14, height=14
    )
}

# ---> Retrieve ag-specific rep. diversity metrics.
# tmp.metrics <- c('gini', 'inv.simp', 'gini.simp', 'chao1')
tmp.metrics <- c('gini', 'inv.simp')
spc.div.metrics <- lapply(X=names(aggr.set.main), FUN=function(data.set){
    div.metrics <- lapply(X=tmp.metrics, FUN=function(x){
        cat(x, '\n')
        tmp.data <- summ.div.metrics(
            tcr.meta=meta.data.l[[data.set]],
            donor.meta=donor.meta, div.metric=x,
            main.sets=aggr.set.main[[data.set]],
            tag.of.int='ag.spc.tag'
        )
        return(tmp.data)
    })
    names(div.metrics) <- tmp.metrics
    div.metrics <- rbindlist(l=div.metrics, use.names=TRUE, idcol='metric')
    return(div.metrics)
})
names(spc.div.metrics) <- names(aggr.set.main)

# ---> Donor-specific repertoire diversity per T-cell lineage (CD4 and CD8).
# Process per metric.
for(data.set in names(aggr.set.labs)){
    for(tmp.metric in tmp.metrics){
        # Retrieve general data
        tmp.data <- spc.div.metrics[[data.set]][metric==tmp.metric]
        tmp.cols <- colnames(tmp.data)[colnames(tmp.data) %like% 'CD[48]']
        tmp.cols <- tmp.cols[!tmp.cols %like% '\\.cell$|\\.clone$|\\.All$']
        tmp.cols <- c('donor.id.tag', tmp.cols)
        tmp.data <- tmp.data[, ..tmp.cols]
        tmp.data <- as.data.table(gather(data=tmp.data, key='condition', value='metric', -`donor.id.tag`))
        tmp.data <- tmp.data[!is.na(metric)]
        tmp.data <- as.data.table(separate(data=tmp.data, col=condition, into=c('cell.type', 'ag.spc'), sep='\\.'))
        # tmp.lvls <- names(pp.cols)[names(pp.cols) %in% tmp.data[, ag.spc]]
        tmp.lvls <- rel.spcs
        tmp.vals <- factor(x=tmp.data[, ag.spc], levels=tmp.lvls)
        set(x=tmp.data, j='ag.spc', value=tmp.vals)
        # Plot per cell type.
        for(t.subset in names(obj.extended.names)){
            tmp.val <- gen.cell.types[t.subset]
            to.plot <- tmp.data[cell.type==tmp.val]
            tmp.ggplot <- ggplot(data=to.plot, aes(x=ag.spc)) +
                geom_boxplot(aes(color=ag.spc, y=metric), width=0.7, linewidth=3, fatten=4, outlier.shape=NA) +
                geom_jitter(aes(y=metric), shape=1, stroke=4, size=8, color='black', width=0) +
                scale_x_discrete(drop=FALSE) +
                scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
                scale_color_manual(values=pp.cols) +
                labs(x='Reactivity', y='Metric', color='')
            tmp.lab <- paste0('/', aggr.set.labs[[data.set]], '_', obj.extended.names[t.subset], '_Metric-', tmp.metric, '_Spc_Donor')
            publish.plot(
                tmp.ggplot=tmp.ggplot, output.path=fig.as.rep.comps.path, file.name=tmp.lab, type='pdf',
                blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=24, height=14
            )
        }
    }
}
# Output table w/ metric values.
for(data.set in names(aggr.set.labs)){
    tmp.file.name <- paste0(fig.as.rep.comps.path, '/', aggr.set.labs[[data.set]], '_DivMetricPerPopAndDonor.csv')
    fwrite(file=tmp.file.name, x=spc.div.metrics[[data.set]], na=NA)
}

# ---> Sharing status of pathogen-specific clonotypes between sites
# Process per set of vars.
for(data.set in names(aggr.set.labs)){
    tmp.vals <- meta.data.l[[data.set]][, unique(data.set)]
    for(t.subset in tmp.vals){
        # @ Data fetch
        # Fetch data to plot.
        tmp.data.1 <- meta.data.l[[data.set]][
            data.set==t.subset & !is.na(clonotype.tag) & !is.na(ag.spc.tag) &
            !is.na(donor.id.tag) & donor.id.tag %in% site.donor.ovlp,
            .(barcode, clonotype.tag, donor.id.tag, ag.spc.tag)
        ]
        tmp.data.2 <- gen.tcr.data.1[[data.set]][
            gex.lin.tag==gen.cell.types[t.subset] &
            !is.na(donor.id.tag) & donor.id.tag %in% site.donor.ovlp,
            .(clonotype.tag, donor.id.tag, comb.ids, site.ovlp)
        ]
        tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by=c('clonotype.tag', 'donor.id.tag'), all.x=TRUE, all.y=FALSE)
        tmp.data.1 <- tmp.data[,
            .(freq.abs=.N),
            by=.(donor.id.tag, ag.spc.tag, site.ovlp)
        ]
        tmp.data.2 <- tmp.data[
            !is.na(donor.id.tag) & !is.na(clonotype.tag),
            .(freq.total=.N),
            by=.(donor.id.tag, ag.spc.tag)
        ]
        tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by=c('donor.id.tag', 'ag.spc.tag'))
        tmp.data[, freq.rel:=freq.abs/freq.total]
        # Apply limit of detection.
        tmp.data <- tmp.data[freq.abs>=limit.of.detect]
        # Force unique values.
        tmp.data[, freq.abs:=NULL]
        tmp.data <- as.data.table(spread(data=tmp.data, key=site.ovlp, value=freq.rel, fill=0))
        # Define file width according to donor amount.
        tmp.width <- (tmp.data[, uniqueN(donor.id.tag)] * 6.68/40) + 0.32
        tmp.width <- (tmp.width * 5)/1.8
        # Plot
        tmp.ggplot <- ggplot(data=tmp.data, aes(x=donor.id.tag)) +
            geom_bar(data=unique(tmp.data[, .(donor.id.tag)]), aes(y=1), stat='identity', color='black', width=0.6, alpha=0, linewidth=0.3) +
            geom_jitter(aes(y=Overlapped, color=ag.spc.tag), shape=1, stroke=2.2, size=3, height=0, width=0.15) +
            scale_color_manual(values=pp.cols) +
            scale_y_continuous(expand=c(0.05, 0), limits=c(0, 1), breaks=scales::pretty_breaks(n=3)) +
            labs(x='Donor ID', y='% of ag-spc cells w/ overlapped clonotypes', color='Reactivity')
        tmp.lab <- paste0('/SimBetSites_Donor_OvlpCellFract_', aggr.set.labs[data.set], '_', obj.extended.names[t.subset])
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.as.rep.comps.path, file.name=tmp.lab, type='pdf',
            blank.comp=blank.complement.7, do.legend=FALSE, do.rotate=TRUE, width=tmp.width, height=5
        )
    }
}

# ---> Clonotype sharing between sites, summary across donors.
for(tmp.lin in names(cell.type.cols)){
    to.plot <- lapply(X=names(gen.tcr.data.1), FUN=function(tmp.site){
        tmp.data.1 <- gen.tcr.data.1[[tmp.site]][gex.lin.tag==tmp.lin]
        tmp.data.2 <- meta.data.l[[tmp.site]][
            gex.lin.tag==tmp.lin &
            !is.na(clonotype.tag) & !is.na(donor.id.tag),
            .(clonotype.tag, donor.id.tag, ag.spc.tag, check=TRUE)
        ]
        tmp.check <- tmp.data.2[, .(to.check=uniqueN(ag.spc.tag)), by=.(clonotype.tag, donor.id.tag)][, all(to.check==1)]
        if(!tmp.check) stop('Unexpected error!\n')
        tmp.data.2 <- unique(tmp.data.2)
        tmp.data <- merge(
            x=tmp.data.1, y=tmp.data.2,
            by=c('clonotype.tag', 'donor.id.tag'),
            all.x=TRUE, all.y=FALSE
        )
        tmp.check <- tmp.data[, !any(is.na(tmp.check)) & all(tmp.check)]
        if(!tmp.check) stop('Unexpected error!\n')
        x <- tmp.data[!is.na(ag.spc.tag), comb.ids]
        x <- unlist(str_split(string=x, pattern='\\/'))
        return(x)
    })
    names(to.plot) <- names(gen.tcr.data.1)
    tmp.vals <- c(`c40.lg`='LG', `c40.ln`='LN')
    names(to.plot) <- tmp.vals[names(to.plot)]
    # Complete v
    tmp.file.name <- paste0(fig.as.rep.comps.path, '/SimBetSites_Ovlp-Clone', '_Lin-', tmp.lin, '.C.tiff')
    VennDiagram::venn.diagram(
        x=to.plot,
        hyper.test=TRUE, total.population=uniqueN(unlist(to.plot)),
        filename=tmp.file.name, lwd=1.5, fill=an.site.cols, col='black',
        disable.logging=TRUE
    )
    # Blank v
    tmp.file.name <- paste0(fig.as.rep.comps.path, '/SimBetSites_Ovlp-Clone', '_Lin-', tmp.lin, '.B.tiff')
    VennDiagram::venn.diagram(
        x=to.plot,
        hyper.test=FALSE, total.population=uniqueN(unlist(to.plot)),
        filename=tmp.file.name, lwd=1.5, fill=an.site.cols, col='black',
        cat.cex=FALSE, cex=FALSE,
        disable.logging=TRUE
    )
}

# ---> Percentage of T cell clonally overlapped between anatomic sites per subset.
tmp.tholds <- c(
    `CD4`=0.002,
    `CD8`=0.08
)
for(tmp.lin in unique(gen.cell.types)){
    # Identify unique sets from each aggr.
    lin.sets <- unlist(lapply(X=aggr.set.main, FUN=function(x){
        names(gen.cell.types[x])[gen.cell.types[x]==tmp.lin]
    }))
    # Collect details for each site as reference.
    tmp.data <- list(
        `LG`=c('c40.lg', 'c40.ln'),
        `LN`=c('c40.ln', 'c40.lg')
    )
    tmp.data <- lapply(X=tmp.data, FUN=function(tmp.ids){
        tmp.data.1 <- get.multi.site.info(
            data.set=tmp.ids[1], t.subset=lin.sets[tmp.ids[1]],
            anti.set=tmp.ids[2], anti.subset=lin.sets[tmp.ids[2]],
            rel.cols=NULL
        )
        tmp.data.2 <- unique(meta.data.l[[tmp.ids[1]]][, 
            .(clonotype.tag, donor.id.tag, ag.spc.tag, check=TRUE)
        ])
        tmp.data <- merge(
            x=tmp.data.1, y=tmp.data.2,
            by=c('clonotype.tag', 'donor.id.tag'),
            all.x=TRUE, all.y=FALSE
        )
        tmp.check <- tmp.data[, !any(is.na(check)) & all(check)]
        if(!tmp.check) stop('Unexpected error!\n')
        tmp.data.1 <- tmp.data[,
            .(freq.abs=.N),
            by=.(comb.ids, donor.id.tag, ag.spc.tag)
        ]
        tmp.data.2 <- tmp.data[,
            .(freq.total=.N),
            by=.(donor.id.tag)
        ]
        tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='donor.id.tag')
        tmp.data[, freq.rel:=freq.abs/freq.total]
        tmp.data[, `:=`(freq.abs=NULL, freq.total=NULL)]
        tmp.data <- tmp.data[!is.na(ag.spc.tag)]
        return(tmp.data)
    })
    tmp.data <- merge(
        x=tmp.data[[1]], y=tmp.data[[2]],
        by=c('donor.id.tag', 'comb.ids'), # Note that considering the specificity ensures we're comparing the exact same clonotypes.
        suffixes=c('.LG', '.LN'),
        all=TRUE
    )
    # ---> Assessment of clonotypes unique to lung.
    to.plot <- tmp.data[
        !is.na(ag.spc.tag.LG) & is.na(ag.spc.tag.LN)
    ]
    to.plot[, `:=`(site='LG', freq.rel=freq.rel.LG)]
    tmp.caption <- paste0(
        'Total shared clonotypes: ', to.plot[, .N], '.'
    )
    # Set upper bound.
    tmp.thold <- tmp.tholds[tmp.lin]
    to.plot[, bounded.freq.rel:=freq.rel]
    to.plot[, point.status:='ori']
    to.plot[bounded.freq.rel>tmp.thold, `:=`(bounded.freq.rel=tmp.thold, point.status='bounded')]
    # Plot.
    tmp.ggplot <- ggplot(data=to.plot, aes(x=site)) +
        geom_boxplot(aes(y=freq.rel, color=site), width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
        geom_jitter(aes(y=bounded.freq.rel, shape=point.status), stroke=3, size=8, color='black', width=0) +
        coord_cartesian(ylim=c(0, tmp.thold)) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=4)) +
        scale_shape_manual(values=tmp.shapes) +
        scale_color_manual(values=an.site.cols) +
        labs(x='Anatomic site', y='Clonotype relative freq. within site', color='', caption=tmp.caption)
    tmp.lab <- paste0(
        '/', tmp.lin, '_Site-LG-Only_AgSpc-CloneRelFreq'
    )
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.as.rep.comps.path, file.name=tmp.lab, type='pdf',
        blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=3.5, height=14
    )
    # ---> Assessment of overlapped clonotypes.
    to.plot <- tmp.data[
        !is.na(ag.spc.tag.LG) & !is.na(ag.spc.tag.LN) &
        ag.spc.tag.LG==ag.spc.tag.LN,
    ]
    to.plot[, `:=`(
        ag.spc.tag=ag.spc.tag.LG,
        ag.spc.tag.LG=NULL, ag.spc.tag.LN=NULL
    )]
    tmp.test <- wilcox.test(x=to.plot[, freq.rel.LG], y=to.plot[, freq.rel.LN], paired=TRUE)
    tmp.caption <- paste0(
        'P-value is ', tmp.test$p.value, '\n',
        'Total shared clonotypes: ', to.plot[, .N], '.'
    )
    to.plot <- as.data.table(gather(data=to.plot, key=site, value='freq.rel', -`donor.id.tag`, -`comb.ids`, -`ag.spc.tag`))
    to.plot[, site:=str_replace(string=site, pattern='freq.rel.', replacement='')]
    to.plot[, group:=paste(donor.id.tag, comb.ids, sep='.')]
    # Set upper bound.
    tmp.thold <- tmp.tholds[tmp.lin]
    to.plot[, bounded.freq.rel:=freq.rel]
    to.plot[, point.status:='ori']
    to.plot[bounded.freq.rel>tmp.thold, `:=`(bounded.freq.rel=tmp.thold, point.status='bounded')]
    # Plot
    tmp.ggplot <- ggplot(data=to.plot, aes(x=site)) +
        geom_boxplot(aes(y=freq.rel, color=site), width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
        geom_jitter(aes(y=bounded.freq.rel, shape=point.status), stroke=3, size=8, color='black', width=0) +
        geom_line(aes(y=bounded.freq.rel, group=group), linewidth=0.2) +
        coord_cartesian(ylim=c(0, tmp.thold)) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=4)) +
        scale_shape_manual(values=tmp.shapes) +
        scale_color_manual(values=an.site.cols) +
        labs(x='Anatomic site', y='Clonotype relative freq. within site', color='', caption=tmp.caption)
    tmp.lab <- paste0(
        '/', tmp.lin, '_Site-Ovlped_AgSpc-CloneRelFreq'
    )
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.as.rep.comps.path, file.name=tmp.lab, type='pdf',
        blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=7, height=14
    )
}


### -------------------------- Supp Figure -------------------------- ###

# ----> Number of query clonotypes w/ an assigned specificity according to each individual strategy.
#      Across donor.
# tmp.data.1 <- pred.ref.2[
#     !is.na(consensus.pred),
#     .(clonotype.abs.freq=uniqueN(qry.clone.id)),
#     by=.(consensus.approach.class)
# ]
# tmp.ggplot <- ggplot(data=tmp.data.1, aes(x=consensus.approach.class, y=clonotype.abs.freq)) +
#     geom_bar(stat='identity', color=cell.type.cols['CD4'], alpha=0, width=0.7, linewidth=5) +
#     scale_y_continuous(expand=expansion(add=c(0, NA)), breaks=scales::pretty_breaks(n=3)) +
#     labs(x='Consensus approach', y='Clonotype count')
# tmp.lab <- paste0('ClonotypeCount_ConsensusApproach_Across')
# publish.plot(
#     tmp.ggplot=tmp.ggplot, output.path=fig.cd4.rep.path, file.name=tmp.lab, type='pdf',
#     blank.comp=blank.complement.3, do.legend=FALSE,
#     width=10, height=10
# )
# #      Per donor.
# tmp.data.1 <- pred.ref.2[
#     !is.na(consensus.pred),
#     .(clonotype.abs.freq=uniqueN(qry.clone.id)),
#     by=.(donor.id.tag, consensus.approach.class)
# ]
# tmp.ggplot <- ggplot(data=tmp.data.1, aes(x=consensus.approach.class, y=clonotype.abs.freq)) +
#     geom_boxplot(color=cell.type.cols['CD4'], width=0.7, linewidth=5, outlier.shape=NA) +
#     geom_jitter(size=5, color='black', width=0.25) +
#     scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
#     labs(x='Consensus approach', y='Clonotype count')
# tmp.lab <- paste0('ClonotypeCount_ConsensusApproach_Donor')
# publish.plot(
#     tmp.ggplot=tmp.ggplot, output.path=fig.cd4.rep.path, file.name=tmp.lab, type='pdf',
#     blank.comp=blank.complement.3, do.legend=FALSE,
#     width=10, height=10
# )



############    -----------------------------------------    ############
### --------------------------- Figure on --------------------------- ###
### ------------- Prediction validations and exploration ------------ ###
############    -----------------------------------------    ############

# ----> Output directory.
fig.pred.val.path <- paste0(reports.path, '/figure_on_pred_vals')
if(!dir.exists(fig.pred.val.path)) dir.create(fig.pred.val.path)
# ---> Comments: Main dataset for this section, CD4 T cell set.
main.obj.name <- NULL


### ------------------------- Text details -------------------------- ###

# ---> CMV signal when startifying donors according to CMV serology status.
tmp.data <- merge(
    x=meta.data.l[['c40.lg']][, .(gex.lin.tag, barcode, donor.id.tag, ag.spc.tag)],
    y=ser.cmv.data,
    by='donor.id.tag',
    all.x=TRUE, all.y=FALSE
)
tmp.data <- tmp.data[,
    .(cmv.cell.fract=.SD[ag.spc.tag=='CMV', .N]/.N),
    by=.(gex.lin.tag, donor.id.tag, cmv.status)
]
tmp.data[,
    .(mean(cmv.cell.fract)*100),
    by=.(gex.lin.tag, cmv.status)
]


### -------------------------- Main Figure -------------------------- ###

# ---> General data retieval.
# Retrieve donor-specific fractions.
tmp.data.1 <- meta.data.l[['c40.lg']][ # It's correct to use this object as the source of cell counts for ag-specific cells. ;)
    !is.na(gex.lin.tag) & !is.na(donor.id.tag) & !is.na(clonotype.tag) & !is.na(ag.spc.tag),
    .(
        cell.count=.N
    ),
    by=.(
        gex.lin.tag,
        donor.id.tag,
        ag.spc=ag.spc.tag
    )
]
tmp.data.2 <- meta.data.l[['c40.lg']][
    !is.na(gex.lin.tag) & !is.na(donor.id.tag) & !is.na(clonotype.tag),
    .(
        total.count=.N
    ),
    by=.(
        gex.lin.tag,
        donor.id.tag
    )
]
merge.cols <- c('gex.lin.tag', 'donor.id.tag')
fig.data <- merge(x=tmp.data.1, y=tmp.data.2, by=merge.cols)
fig.data[, cell.fract:=cell.count/total.count]
fig.rel.spcs <- c('CMV', 'EBV', 'SARS-CoV-2', 'IAV', 'PIV', 'RSV', 'MPV')
fig.data <- fig.data[ag.spc %in% fig.rel.spcs]
fig.data <- fig.data[, .(donor.id.tag, ag.spc, gex.lin.tag, cell.fract)]
# tmp.data[, cell.fract:=car::logit(p=cell.fract)]
fig.data <- as.data.table(spread(data=fig.data, key=gex.lin.tag, value=cell.fract, fill=0))
# Merge w/ donor clinics and demographics
#       CMV serology.
fig.data <- merge(x=fig.data, y=ser.cmv.data, by='donor.id.tag')
# setorderv(x=tmp.data, cols='cmv.igg')
fig.data[is.na(cmv.igg), cmv.igg:=0]
#       Other clinical an demographical variables.
fig.data <- merge(x=fig.data, y=donor.meta, by='donor.id.tag', all.x=TRUE)


# ---> Correlation of donor-specific cell fractions between lineages on a pathogen specific basis.
# Plot.
tmp.ggplot <- ggplot(data=fig.data, aes(x=CD4, y=CD8, color=ag.spc)) +
    geom_point(shape=1, stroke=5, size=10) +
    # geom_smooth(method='lm', formula=y~x, se=FALSE, linewidth=4) +
    scale_color_manual(values=pp.cols) +
    scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
    labs(x='CD4 donor-specific cell fractions', y='CD8 donor-specific cell fractions', color='Pathogen\nspecificity')
tmp.lab <- paste0(
    '/TLineageCellFract_Donor_Specificity'
)
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.pred.val.path, file.name=tmp.lab, type='pdf',
    stat.cor=TRUE, cor.group='ag.spc',
    blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=14, height=14
)

# ---> Association between predictions and continuous variables.
# Define continuous variables.
cont.vars <- c(
    `CMV-IgG`='cmv.igg',
    `Age`='age',
    `BMI`='bmi',
    `AlcoholUnits`='alcohol.units',
    `VaccSpan-SARS`='covid.shot.span',
    `VaccSpan-IAV`='flu.shot.span'
)
# Process per cell type.
for(t.subset in names(cell.type.cols)){
    for(tmp.spc in fig.rel.spcs){
        for(cont.var in names(cont.vars)){
            tmp.data <- fig.data[
                ag.spc==tmp.spc,
                .(
                    cont.var=get(cont.vars[cont.var]),
                    freq.var=get(t.subset),
                    tmp.var='All'
                )
            ]
            tmp.ggplot <- ggplot(data=tmp.data, aes(x=cont.var, y=freq.var)) +
                geom_point(shape=1, stroke=5, size=12) +
                geom_smooth(method='lm', formula=y~x, se=FALSE, fullrange=TRUE, color='black', linewidth=5) +
                scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
                scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
                labs(x='Continuous variable', y='% of ag-specific T cells')
            tmp.lab <- paste0(
                '/Association_', tmp.spc, '_',  cont.var, '_CellType-', t.subset
            )
            publish.plot(
                tmp.ggplot=tmp.ggplot, output.path=fig.pred.val.path, file.name=tmp.lab, type='pdf',
                stat.cor=TRUE, cor.group='tmp.var',
                blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=14, height=14
            )
        }
    }
}

# ---> Association between predictions and discrete variables.
# Define continuous variables.
disc.vars <- c(
    `CMV-Serology`='cmv.status',
    `Gender`='gender',
    `SmokingStatus`='smoking.status'
)
# Plotting settings.
tmp.bounds <- c(
    `CD4.CMV.CMV-Serology`=0.06,
    `CD8.CMV.CMV-Serology`=0.20
)
disc.var.cols <- list(
    `CMV-Serology`=c(`Positive`=unname(pp.cols['CMV']), `Negative`='#CBCBCB')
)
# Process per cell type.
for(t.subset in names(cell.type.cols)){
    for(tmp.spc in fig.rel.spcs){
        for(disc.var in names(disc.vars)){
            tmp.data <- fig.data[
                ag.spc==tmp.spc,
                .(
                    disc.var=get(disc.vars[disc.var]),
                    freq.var=get(t.subset),
                    tmp.var='All'
                )
            ]
            # Stats.
            if(tmp.data[, uniqueN(disc.var)==2]){
                tmp.groups <- tmp.data[, unique(disc.var)]
                tmp.test <- wilcox.test(
                    x=tmp.data[disc.var==tmp.groups[1], freq.var],
                    y=tmp.data[disc.var==tmp.groups[2], freq.var],
                    paired=FALSE
                )
                tmp.caption <- paste0('P-value is ', tmp.test$p.value)
            }else{
                tmp.caption <- ''
            }
            # Set upper bound.
            tmp.lab <- paste(t.subset, tmp.spc, disc.var, sep='.')
            tmp.bound <- if(tmp.lab %in% names(tmp.bounds)) tmp.bounds[tmp.lab] else tmp.data[, max(freq.var)]
            tmp.data[, bounded.var:=freq.var]
            tmp.data[, point.status:='ori']
            tmp.data[freq.var>tmp.bound, `:=`(bounded.var=tmp.bound, point.status='bounded')]
            # Plot.
            tmp.ggplot <- ggplot(data=tmp.data, aes(x=disc.var)) +
                geom_boxplot(aes(y=freq.var, color=disc.var), width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
                geom_jitter(aes(y=bounded.var, shape=point.status), stroke=5, size=10, color='black', width=0) +
                scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
                scale_shape_manual(values=tmp.shapes) +
                # scale_color_manual(values=tmp.cols) +
                coord_cartesian(ylim=c(0, tmp.bound)) +
                labs(x='', y='% of ag-specific T cells', caption=tmp.caption, color='') +
                theme(legend.position='none')
            if(disc.var %in% names(disc.var.cols)){
                tmp.ggplot <- tmp.ggplot +
                    scale_color_manual(values=disc.var.cols[[disc.var]])
            }
            tmp.lab <- paste0(
                '/Comparison_', tmp.spc, '_',  disc.var, '_CellType-', t.subset
            )
            publish.plot(
                tmp.ggplot=tmp.ggplot, output.path=fig.pred.val.path, file.name=tmp.lab, type='pdf',
                blank.comp=blank.complement.3.1, do.legend=FALSE,
                width=7, height=14
            )
            # tmp.ggplot <- ggplot(data=tmp.data, aes(x=gex.lin.tag, y=cell.count)) +
            #     geom_boxplot(aes(color=gex.lin.tag), width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
            #     geom_jitter(shape=1, stroke=5, size=10, color='black', width=0) +
            #     geom_line(aes(group=donor.id.tag), linewidth=1) +
            #     scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
            #     scale_color_manual(values=cell.type.cols) +
            #     labs(x='Cell type', y='Number of cells', color='', caption=tmp.caption)
            # tmp.lab <- paste0(
            #     '/', aggr.set.labs[data.set], '_CellCount_CellType_Donor'
            # )
            # publish.plot(
            #     tmp.ggplot=tmp.ggplot, output.path=fig.tcr.reps.path, file.name=tmp.lab, type='pdf',
            #     blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=7, height=14
            # )
        }
    }
}


############    -----------------------------------------    ############
### --------------------------- Figure on --------------------------- ###
### ------------- Phenotype comparisons between lineages ------------ ###
############    -----------------------------------------    ############

# ----> Output directory.
fig.phen.comps.path <- paste0(reports.path, '/figure_on_phen_comps')
if(!dir.exists(fig.phen.comps.path)) dir.create(fig.phen.comps.path)

### ------------------------- Text details -------------------------- ###

# Number of cells per lineage.
lapply(X=meta.data.l, FUN=function(meta.data)
    meta.data[, .N, by=gex.lin.tag]
)
lapply(X=srt.objs.list, FUN=function(x) length(Cells(x)))

### -------------------------- Main Figure -------------------------- ###

# ---> UMAP plots depicting clusters.
for(t.subset in names(obj.extended.names)){
    cat(t.subset, '\n')
    tmp.vals <- c(
        'c40.lg.cd8.d40',
        'c40.lg.cd4.d40'
    )
    dot.size <- if(t.subset %in% tmp.vals) 0.05 else 0.1
    # Whole dataset.
    tmp.ggplot <- get.umap.gg(obj.name=t.subset, attempted.format='tiff', cell.subset=NULL, dot.size=dot.size)
    tmp.lab <- paste0(obj.extended.names[t.subset], '_PopsOnUMAP_Whole')
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.phen.comps.path, file.name=tmp.lab, type='tiff',
        blank.comp=blank.complement.1, do.legend=FALSE,
        height=10, width=10
    )
    # Subset
    tmp.ggplot <- get.umap.gg(obj.name=t.subset, attempted.format='tiff', cell.subset=subset.list[[t.subset]], dot.size=dot.size)
    tmp.lab <- paste0(obj.extended.names[t.subset], '_PopsOnUMAP_Subset')
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.phen.comps.path, file.name=tmp.lab, type='tiff',
        blank.comp=blank.complement.1, do.legend=FALSE,
        height=10, width=10
    )
}

# ---> Proportion of cells per cluster.
for(t.subset in names(obj.extended.names)){
    tmp.ggplot <- get.cell.fracs.gg(
        seurat.obj=srt.objs.list[[t.subset]],
        x.var=NULL, fill.var=clust.labs[t.subset], fill.cols=clusts.cols[[t.subset]],
        break.no=3
    )
    tmp.ggplot
    tmp.lab <- paste0(obj.extended.names[t.subset], '_CellFracs_AllCells_Clusters')
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.phen.comps.path, file.name=tmp.lab, type='pdf',
        blank.comp=blank.complement.3.1, do.legend=FALSE,
        width=2, height=21
    )
}

# ---> Feature plots.
tmp.reports.path <- paste0(fig.phen.comps.path, '/feature_plots')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
# Preflights
this.color.scale <- c("#fff8a6", "#f5e74e", "#f0dd0c", "#FFA500", "#FF3500", "#670000")
markers.of.int <- c(
  'HLA-DRB1', 'HLA-DRB5', 'HLA-DRA', 'HLA-DPA1', 'HLA-DQA1',
  'ZNF683', 'CXCR6', 'FABP5', 'BATF', 'RBPJ',
  'IFNG', 'GZMB', 'IL7R', 'CD27', 'TCF7'
)
# these.subsets <- list(
#   # `Full`=lapply(X=srt.objs.list, FUN=function(x) Cells(x)), # If ever necessary.
#   `50K`=gpr25.subset.list,
#   `30K`=lapply(X=gpr25.subset.list, FUN=function(x) x[1:30000])
# )
# Process per setting option.
for(tmp.obj in names(srt.objs.list)[3]){
  # Retrieve expression and dimentionality values
  tmp.data.1 <- as.data.table(cbind(
    data.frame(barcode=Cells(srt.objs.list[[tmp.obj]])),
    # srt.objs.list[[tmp.obj]]@meta.data,
    srt.objs.list[[tmp.obj]]@reductions$umap@cell.embeddings
  ))
  tmp.data.2 <- as.data.table(t(GetAssayData(object=srt.objs.list[[tmp.obj]], slot='data')[
    markers.of.int, ,
    drop=F
  ]))
#   colnames(tmp.data.2) <- markers.of.int
  tmp.data <- cbind(tmp.data.1, tmp.data.2)
#   for(subset.name in names(these.subsets)){
#     tmp.subset <- these.subsets[[subset.name]][[tmp.obj]]
    # to.plot <- tmp.data[barcode %in% tmp.subset]
    to.plot <- tmp.data
    # for(att.format in c('pdf', 'tiff')){
    for(att.format in c('pdf')){
      for(tmp.marker in markers.of.int){
        if('tmp.marker' %in% colnames(to.plot)) to.plot[, tmp.marker:=NULL]
        to.plot[, tmp.marker:=get(tmp.marker)]
        # to.plot[, size.col:=ifelse(test=tmp.marker>0, yes='+', no='-')]
        # setorderv(x=to.plot, cols=tmp.marker, order=1)
        # size.vals <- if(subset.name=='Full') c('-'=0.3, '+'=0.5) else c('-'=0.5, '+'=1.2)
        # size.vals <- c('-'=0.3, '+'=0.5)
        tmp.ggplot <- ggplot(data=to.plot, aes(x=UMAP_1, y=UMAP_2, col=tmp.marker)) +
            # geom_point(aes(size=size.col), alpha=1) +
            geom_point(alpha=1) +
            scale_color_gradientn(colors=this.color.scale, breaks=scales::pretty_breaks(n=3)) +
            # scale_size_manual(values=size.vals) +
            labs(x='UMAP 1', y='UMAP 2', col='Expression\n(Seurat-\nnormalized)')
        # Output.
        # tmp.lab <- paste0(obj.extended.names[tmp.obj], '_FeaturePlot_Feat-', tmp.marker, '_CellSet-', subset.name)
        tmp.lab <- paste0(obj.extended.names[tmp.obj], '_FeaturePlot_Feat-', tmp.marker)
        publish.plot(
          tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type=att.format,
          blank.comp=blank.complement.1, do.legend=TRUE,
          height=10, width=10
        )
      }
    # }
  }
}

# ---> Correlation of donor-specific cell fractions of phenotype counterparts between lineages.
pop.corrs <- list(
    `c40.lg`=c(
        `CD4-0`='TRM',
        `CD4-1`='TCM',
        `CD4-2`='TEM',
        `CD8-0`='TRM',
        `CD8-1`='TEM',
        `CD8-3`='TCM'
    ),
    `c40.ln`=c(
        `CD4-6`='TRM',
        `CD4-0`='TCM',
        # `CD4-4`='GPR183hi',
        `CD4-2`='IFNR',
        `CD8-3`='TRM',
        `CD8-2`='TCM',
        # `CD8-4`='GPR183hi',
        `CD8-1`='IFNR'
    )
)
for(data.set in names(aggr.set.main)){
    # Retrieve metadata
    tmp.data <- meta.data.l[[data.set]]
    # Retrieve donor-specific fractions.
    tmp.data <- tmp.data[, .(gex.lin.tag, barcode, clonotype.tag, donor.id.tag, clusters.tag)]
    # Potentially merge some clusters.
    # meta.data[gex.lin.tag=='CD8' & tmp.cluster.tag=='4', tmp.cluster.tag:='0']
    tmp.data.1 <- tmp.data[
        !is.na(gex.lin.tag) & !is.na(donor.id.tag) & !is.na(clonotype.tag),
        .(
            cell.count=.N
        ),
        by=.(
            gex.lin.tag,
            donor.id.tag,
            clusters.tag
        )
    ]
    tmp.data.2 <- tmp.data[
        !is.na(gex.lin.tag) & !is.na(donor.id.tag) & !is.na(clonotype.tag),
        .(
            total.count=.N
        ),
        by=.(
            gex.lin.tag,
            donor.id.tag
        )
    ]
    merge.cols <- c('gex.lin.tag', 'donor.id.tag')
    tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by=merge.cols)
    tmp.data[, cell.fract:=cell.count/total.count]
    # Assign general phenotype disregarding lineage label.
    tmp.data[, pop.tag:=paste(gex.lin.tag, clusters.tag, sep='-')]
    tmp.vals <- pop.corrs[[data.set]]
    tmp.data <- tmp.data[pop.tag %in% names(tmp.vals)]
    tmp.data[, pop.tag:=tmp.vals[pop.tag]]
    tmp.data <- tmp.data[, .(donor.id.tag, pop.tag, gex.lin.tag, cell.fract)]
    tmp.data <- spread(data=tmp.data, key=gex.lin.tag, value=cell.fract)
    # Plot.
    tmp.ggplot <- ggplot(data=tmp.data, aes(x=CD4, y=CD8, color=pop.tag)) +
        geom_point(shape=1, stroke=5, size=10, alpha=0.6) +
        geom_smooth(method='lm', formula=y~x, se=FALSE, fullrange=TRUE, linewidth=2.5) +
        scale_color_manual(values=gen.subset.cols) +
        scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
        labs(x='CD4 donor-specific cell fractions', y='CD8 donor-specific cell fractions', color='Phenotype')
    tmp.lab <- paste0(
        '/', aggr.set.labs[[data.set]], '_TLineageCellFract_Donor_Phenotype'
    )
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.phen.comps.path, file.name=tmp.lab, type='pdf',
        stat.cor=TRUE, cor.group='pop.tag',
        blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=14, height=14
    )
}

# ---> Expanded cell fraction, general function for different group-related columns in metadata.
set.rel.pops <- list(
    `c40.lg`=list(
        `c40.lg.cd4.d40`=c('0', '1'),
        `c40.lg.cd8.d40`=c('0', '3')
    ),
    `c40.ln`=list(
        `c40.ln.cd4.d08`=c('6', '4'),
        `c40.ln.cd8.d08`=c('3', '4')
    )
) 
comp.group.expansion <- function(
    group.tag='clusters.tag', group.lab='Cluster', group.cols,
    set.rel.pops, freq.thold=1,
    tmp.reports.path
){
    for(data.set in names(aggr.set.main)){
        for(t.subset in aggr.set.main[[data.set]]){
            tmp.data.1 <- meta.data.l[[data.set]][
                data.set==t.subset &
                !is.na(clonotype.tag) & !is.na(donor.id.tag) &
                !is.na(get(group.tag)),
                .(freq.total=.N),
                by=.(donor.id.tag, group=get(group.tag))
            ]
            tmp.data.2 <- meta.data.l[[data.set]][
                data.set==t.subset &
                !is.na(clonotype.tag) & !is.na(donor.id.tag) &
                !is.na(get(group.tag)),
                .(freq.abs=.N),
                by=.(clonotype.tag, donor.id.tag, group=get(group.tag))
            ]
            tmp.data.2 <- tmp.data.2[
                freq.abs>freq.thold,
                .(freq.abs=sum(freq.abs)),
                by=.(donor.id.tag, group)
            ]
            tmp.data <- merge(
                x=tmp.data.1, y=tmp.data.2,
                by=c('donor.id.tag', 'group'),
                all=TRUE
            )
            tmp.data[is.na(freq.abs), freq.abs:=0]
            tmp.data[, freq.rel:=freq.abs/freq.total]
            # ---> Plot w/ all groups.
            tmp.data$group <- factor(x=tmp.data$group, levels=names(group.cols[[t.subset]]))
            tmp.ggplot <- ggplot(data=tmp.data, aes(x=group, y=freq.rel)) +
                geom_jitter(shape=1, stroke=4, size=10, color=gen.dot.col, width=0, height=0) +
                geom_boxplot(aes(color=group, fill=group), alpha=0.4, width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
                scale_y_continuous(limits=c(0, 1), expand=c(0, 0), breaks=c(0, 0.5, 1)) +
                scale_color_manual(values=group.cols[[t.subset]]) +
                scale_fill_manual(values=group.cols[[t.subset]]) +
                labs(x='Cell type', y='Expanded clones (% of T cells)', color='')
            tmp.lab <- paste0(obj.extended.names[t.subset], '_ExpandedCellFract', '_', group.lab, '_DonorID')
            publish.plot(
                tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
                height=14, width=(tmp.data[, uniqueN(group)]*2.8),
                blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE
            )
            # ---> Plot w/ selected populations.
            if(!is.null(set.rel.pops)){
                tmp.data[, group:=as.character(group)]
                if(length(set.rel.pops[[data.set]][[t.subset]])!=2) stop('Unexpected error. Expected comparison between TRM and TCM, specified in that order.\n')
                tmp.vals <- c('TRM', 'TCM')
                names(tmp.vals) <- set.rel.pops[[data.set]][[t.subset]]
                to.plot <- tmp.data[group %in% names(tmp.vals)]
                to.plot[, group:=tmp.vals[group]]
                to.plot$group <- factor(x=to.plot$group, levels=tmp.vals)
                # Stats
                tmp.stats <- to.plot[, .(donor.id.tag, group, freq.rel)]
                tmp.stats <- as.data.table(spread(data=tmp.stats, key=group, value=freq.rel))
                tmp.test <- wilcox.test(
                    x=tmp.stats[, TRM],
                    y=tmp.stats[, TCM],
                    alternative='two.sided', paired=TRUE, exact=TRUE
                )
                tmp.caption <- paste0('Wilcoxon signed-rank test nominal P-value is: ', tmp.test$p.value)
                # Plot
                tmp.ggplot <- ggplot(data=to.plot, aes(x=group, y=freq.rel)) +
                    geom_line(aes(group=donor.id.tag), color=gen.link.col, linewidth=gen.link.width) +
                    geom_jitter(shape=1, stroke=5, size=10, color=gen.dot.col, width=0, height=0) +
                    geom_boxplot(aes(color=group, fill=group), alpha=0.4, width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
                    scale_y_continuous(limits=c(0, 1), expand=c(0, 0), breaks=c(0, 0.5, 1)) +
                    scale_color_manual(values=gen.subset.cols) +
                    scale_fill_manual(values=gen.subset.cols) +
                    labs(x='Cell type', y='Expanded clones (% of T cells)', color='', caption=tmp.caption)
                tmp.lab <- paste0(obj.extended.names[t.subset], '_ExpandedCellFract', '_', group.lab, '-Summ_DonorID')
                publish.plot(
                    tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
                    blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=7, height=14
                )
            }
        }
    }
}
# ---> Expanded cell fraction, for phenotype-biased cells.
comp.group.expansion(
    group.tag='group.class.tag', group.lab='Bias', group.cols=clusts.cols,
    set.rel.pops=set.rel.pops, freq.thold=10,
    tmp.reports.path=fig.phen.comps.path
)
# ---> Expanded cell fraction, for raw clusters.
comp.group.expansion(
    group.tag='clusters.tag', group.lab='Cluster', group.cols=clusts.cols,
    set.rel.pops=set.rel.pops, freq.thold=1,
    tmp.reports.path=fig.phen.comps.path
)

# ---> Diversity metrics.
comp.group.div <- function(
    these.metrics, group.tag='clusters.tag',
    group.lab='Cluster', group.cols,
    set.rel.pops,
    tmp.reports.path
){
    tmp.metrics <- these.metrics[[1]][, unique(metric)][1]
    names(tmp.metrics) <- str_to_title(string=str_replace_all(string=tmp.metrics, pattern='\\.', replacement=' '))
    for(data.set in names(aggr.set.main)){
        for(tmp.metric in names(tmp.metrics)){
            for(t.subset in aggr.set.main[[data.set]]){
                # Retrieve metrics data.
                tmp.cols <- c(
                    'donor.id.tag',
                    paste(gen.cell.types[t.subset], names(group.cols[[t.subset]]), sep='.')
                )
                tmp.check <- all(tmp.cols %in% colnames(these.metrics[[data.set]]))
                if(!tmp.check){
                    tmp.warn <- setdiff(x=tmp.cols, y=colnames(these.metrics[[data.set]]))
                    tmp.warn <- paste0('Following columns could not be found in diversity metrics table:\n', paste0(tmp.warn, collapse='\n'))
                    warning(tmp.warn)
                    tmp.cols <- tmp.cols[tmp.cols %in% colnames(these.metrics[[data.set]])]
                }
                tmp.data.1 <- these.metrics[[data.set]][
                    metric %in% tmp.metrics[tmp.metric],
                    ..tmp.cols
                ]
                colnames(tmp.data.1) <- str_replace(string=colnames(tmp.data.1), pattern=paste0('^', gen.cell.types[t.subset], '\\.'), replacement='')
                tmp.data.1 <- as.data.table(gather(data=tmp.data.1, key='group', value='metric', -`donor.id.tag`))
                tmp.data.2 <- meta.data.l[[data.set]][
                    data.set==t.subset & !is.na(clonotype.tag) &
                    !is.na(get(group.tag)),
                    .(cell.count=.N),
                    by=.(
                        donor.id.tag,
                        group=get(group.tag)
                    )
                ]
                tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by=c('donor.id.tag', 'group'), all=TRUE)
                # Remove entries w/ small number of cells.
                cell.thold <- 50 # Actually, this is mandated from the beginning when metrics are calculated
                tmp.data <- tmp.data[cell.count>cell.thold]
                tmp.data <- tmp.data[!is.na(metric)]
                # Factors
                tmp.data$group <- factor(x=tmp.data$group, levels=names(group.cols[[t.subset]]))
                # ---> Plot w/ all populations.
                tmp.lims <- if(tmp.metric=='Gini') c(0, 1) else c(NA, NA)
                # tmp.breaks <- if(tmp.metric=='Gini') c(0, 0.3, 0.6, 0.9) else scales::pretty_breaks(n=3)
                tmp.breaks <- if(tmp.metric=='Gini') c(0, 0.5, 1) else scales::pretty_breaks(n=3)
                tmp.ggplot <- ggplot(data=tmp.data, aes(x=group, y=metric)) +
                    geom_jitter(shape=1, stroke=4, size=10, color=gen.dot.col, width=0, height=0) +
                    geom_boxplot(aes(color=group, fill=group), alpha=0.4, width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
                    scale_y_continuous(limits=tmp.lims, expand=c(0, 0), breaks=tmp.breaks) +
                    scale_color_manual(values=group.cols[[t.subset]]) +
                    scale_fill_manual(values=group.cols[[t.subset]]) +
                    labs(x='Cell type', y='Metric', color='')
                tmp.lab <- paste0(obj.extended.names[t.subset], '_Idx-', tmp.metric, '_', group.lab, '_DonorID')
                publish.plot(
                    tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
                    height=14, width=(tmp.data[, uniqueN(group)]*2.8),
                    blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE
                )
                # ---> Plot w/ selected populations.
                tmp.data[, group:=as.character(group)]
                if(length(set.rel.pops[[data.set]][[t.subset]])!=2) stop('Unexpected error. Expected comparison between TRM and TCM, specified in that order.\n')
                tmp.vals <- c('TRM', 'TCM')
                names(tmp.vals) <- set.rel.pops[[data.set]][[t.subset]]
                to.plot <- tmp.data[group %in% names(tmp.vals)]
                to.plot[, group:=tmp.vals[group]]
                to.plot$group <- factor(x=to.plot$group, levels=tmp.vals)
                # Stats
                tmp.stats <- to.plot[, .(donor.id.tag, group, metric)]
                tmp.stats <- as.data.table(spread(data=tmp.stats, key=group, value=metric))
                tmp.test <- wilcox.test(
                    x=tmp.stats[, TRM],
                    y=tmp.stats[, TCM],
                    alternative='two.sided', paired=TRUE, exact=TRUE
                )
                tmp.caption <- paste0('Wilcoxon signed-rank test nominal P-value is: ', tmp.test$p.value)
                # Plot
                tmp.ggplot <- ggplot(data=to.plot, aes(x=group, y=metric)) +
                    geom_line(aes(group=donor.id.tag), color=gen.link.col, linewidth=gen.link.width) +
                    geom_jitter(shape=1, stroke=5, size=10, color=gen.dot.col, width=0, height=0) +
                    geom_boxplot(aes(color=group, fill=group), alpha=0.4, width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
                    scale_y_continuous(limits=tmp.lims, expand=c(0, 0), breaks=tmp.breaks) +
                    scale_color_manual(values=gen.subset.cols) +
                    scale_fill_manual(values=gen.subset.cols) +
                    labs(x='Cell type', y='Metric', color='', caption=tmp.caption)
                tmp.lab <- paste0(obj.extended.names[t.subset], '_Idx-', tmp.metric, '_', group.lab, '-Summ_DonorID')
                publish.plot(
                    tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
                    blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=7, height=14
                )
            }
        }
    }
}
# ---> Diversity metrics, for phenotype-biased cells.
comp.group.div(
    these.metrics=bias.div.metrics, group.tag='group.class.tag',
    group.lab='Bias', group.cols=clusts.cols,
    set.rel.pops=set.rel.pops,
    tmp.reports.path=fig.phen.comps.path
)
# ---> Diversity metrics, for raw clusters.
comp.group.div(
    these.metrics=div.metrics, group.tag='clusters.tag',
    group.lab='Cluster', group.cols=clusts.cols,
    set.rel.pops=set.rel.pops,
    tmp.reports.path=fig.phen.comps.path
)

# ---> Diversity metrics, association between lineages.
pop.corrs <- list(
    `c40.lg`=c(
        `CD4.0`='TRM',
        `CD4.1`='TCM',
        `CD4.2`='TEM',
        `CD8.0`='TRM',
        `CD8.1`='TEM',
        `CD8.3`='TCM'
    ),
    `c40.ln`=c(
        `CD4.6`='TRM',
        `CD4.0`='TCM',
        # `CD4-4`='GPR183hi',
        `CD4.2`='IFNR',
        `CD8.3`='TRM',
        `CD8.2`='TCM',
        # `CD8-4`='GPR183hi',
        `CD8.1`='IFNR'
    )
)
tmp.metrics <- div.metrics[, unique(metric)]
names(tmp.metrics) <- str_to_title(string=str_replace_all(string=tmp.metrics, pattern='\\.', replacement=' '))
for(data.set in names(aggr.set.main)){
    for(tmp.metric in names(tmp.metrics)){
        # Retrieve metrics data.
        tmp.vals <- pop.corrs[[data.set]]
        tmp.cols <- c(
            'donor.id.tag',
            names(tmp.vals)
        )
        # tmp.check <- all(tmp.cols %in% colnames(div.metrics[[data.set]]))
        # if(!tmp.check){
        #     tmp.warn <- setdiff(x=tmp.cols, y=colnames(div.metrics[[data.set]]))
        #     tmp.warn <- paste0('Following columns could not be found in diversity metrics table:\n', paste0(tmp.warn, collapse='\n'))
        #     warning(tmp.warn)
        #     tmp.cols <- tmp.cols[tmp.cols %in% colnames(div.metrics[[data.set]])]
        # }
        tmp.data.1 <- div.metrics[[data.set]][
            metric %in% tmp.metrics[tmp.metric],
            ..tmp.cols
        ]
        tmp.data.1 <- as.data.table(gather(data=tmp.data.1, key='pop', value='metric', -`donor.id.tag`))
        tmp.data.2 <- meta.data.l[[data.set]][
            !is.na(clonotype.tag),
            .(cell.count=.N),
            by=.(
                donor.id.tag,
                gex.lin.tag, clusters.tag
                # pop=paste(gex.lin.tag, clusters.tag, sep='.')
            )
        ]
        tmp.data.2[,
            `:=`(
                pop=paste(gex.lin.tag, clusters.tag, sep='.'),
                gex.lin.tag=NULL, clusters.tag=NULL
            )
        ]
        tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by=c('donor.id.tag', 'pop'), all.x=TRUE, all.y=FALSE)
        # Remove entries w/ small number of cells.
        cell.thold <- 50 # Actually, this is mandated from the beginning when metrics are calculated
        tmp.data <- tmp.data[cell.count>cell.thold]
        tmp.data <- tmp.data[!is.na(metric)]
        tmp.data[, cell.type:=str_extract(string=pop, pattern='CD4|CD8')]
        tmp.data[, pop:=tmp.vals[pop]]
        # Final format.
        tmp.data[, cell.count:=NULL]
        tmp.data <- spread(data=tmp.data, key=cell.type, value=metric)
        # Factors
        # tmp.data$cluster <- factor(x=tmp.data$cluster, levels=names(clusts.cols[[t.subset]]))
        # Plot.
        tmp.ggplot <- ggplot(data=tmp.data, aes(x=CD4, y=CD8, color=pop)) +
            geom_point(shape=1, stroke=5, size=10, alpha=0.6) +
            geom_smooth(method='lm', formula=y~x, se=FALSE, fullrange=TRUE, linewidth=2.5) +
            scale_color_manual(values=gen.subset.cols) +
            scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
            scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
            labs(x='Index in CD4', y='Index in CD8', color='Phenotype')
        tmp.lab <- paste0(
            '/', aggr.set.labs[[data.set]], '_Idx-', tmp.metric, '_TLineage_Donor_Phenotype'
        )
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.phen.comps.path, file.name=tmp.lab, type='pdf',
            stat.cor=TRUE, cor.group='pop',
            blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=14, height=14
        )
    }
}

# ---> Raw clonotype sharing between TRM and other subsets.
# Only considering TCRs w/ high frequency.
freq.thold <- 10
clusts.to.rm <- list(
    `c40.lg.cd4.d40`=c('5'),
    `c40.lg.cd8.d40`=c('5', '6')
)
for(data.set in names(aggr.set.labs)){
    tmp.vals <- meta.data.l[[data.set]][, unique(data.set)]
    for(t.subset in tmp.vals){
        tmp.data.1 <- meta.data.l[[data.set]][
            data.set==t.subset &
            !is.na(donor.id.tag) & !is.na(clonotype.tag),
            .(
                barcode,
                trm.tag, clusters.tag=get(clust.labs[t.subset]),
                clone.id=paste(clonotype.tag, donor.id.tag, sep=';')
            )
        ]
        if(t.subset %in% names(clusts.to.rm)) tmp.data.1 <- tmp.data.1[!clusters.tag%in%clusts.to.rm[[t.subset]]]
        tmp.data.2 <- tmp.data.1[, .(freq.abs=.N), by=clone.id]
        tmp.data <- merge(
            x=tmp.data.1, y=tmp.data.2,
            by='clone.id'
        )
        tmp.data <- tmp.data[freq.abs>=freq.thold]
        to.plot <- list(
            `TRM`=tmp.data[trm.tag=='TRM', unique(clone.id)],
            `Rest`=tmp.data[trm.tag=='Rest', unique(clone.id)]
        )
        tmp.cluster <- names(which(clusts.defs[[t.subset]]=='TRM'))
        tmp.cols <- c(clusts.cols[[t.subset]][tmp.cluster], `Rest`='#CBCBCB'); names(tmp.cols)[1] <- 'TRM'
        # Complete v
        tmp.file.name <- paste0(fig.phen.comps.path, '/', obj.extended.names[t.subset], '_TRM-vs-Rest_Ovlp-Clone', '.C.tiff')
        VennDiagram::venn.diagram(
            x=to.plot,
            hyper.test=TRUE, total.population=uniqueN(unlist(to.plot)),
            filename=tmp.file.name, lwd=3, fill=tmp.cols, col='black',
            disable.logging=TRUE
        )
        # Blank v
        tmp.file.name <- paste0(fig.phen.comps.path, '/', obj.extended.names[t.subset], '_TRM-vs-Rest_Ovlp-Clone', '.B.tiff')
        VennDiagram::venn.diagram(
            x=to.plot,
            hyper.test=FALSE, total.population=uniqueN(unlist(to.plot)),
            filename=tmp.file.name, lwd=3, fill=tmp.cols, col='black',
            cat.cex=FALSE, cex=FALSE,
            disable.logging=TRUE
        )
    }
}

# ---> Plastic vs subset-biased clonotypes, attempts 1 and 2, general function.
gen.pvb.perform <- function(
    this.att, tmp.tholds, data.set,
    set.rel.pops, set.up.bounds=NULL,
    gen.row.clust=TRUE, gen.col.clust=TRUE,
    tmp.reports.path, height=7, width=7
){
    for(t.subset in aggr.set.main[[data.set]]){
        cat(t.subset, '\n')
        for(tmp.thold in tmp.tholds){
            cat(tmp.thold, '\n')
            # @ Data fetch
            tmp.data <- meta.data.l[[data.set]][gex.lin.tag==gen.cell.types[t.subset]]
            # Fetch data to plot.
            tmp.data.1 <- tmp.data[
                !is.na(clonotype.tag) & !is.na(donor.id.tag),
                .(cell.count=.N),
                by=.(donor.id.tag, clonotype.tag, cluster=clusters.tag)
            ]
            tmp.data.2 <- tmp.data[
                !is.na(clonotype.tag) & !is.na(donor.id.tag),
                .(total.count=.N),
                by=.(donor.id.tag)
            ]
            tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='donor.id.tag', all.x=TRUE, all.y=FALSE)
            tmp.data[, cell.fract:=cell.count/total.count]
            # Subsets of interest
            tmp.pops <- set.rel.pops[[t.subset]]
            tmp.data <- tmp.data[
                cluster %in% tmp.pops
            ]
            # Clonotype subset.
            if(this.att=='att.1'){
                to.filter <- tmp.data[
                    cell.fract>=tmp.thold,
                    .(donor.id.tag, clonotype.tag)
                ]
                to.filter <- unique(to.filter)
            }else{
                if(this.att=='att.2'){
                    to.filter <- lapply(X=tmp.pops, FUN=function(tmp.pop){
                        to.filter <- tmp.data[cluster==tmp.pop]
                        setorderv(x=to.filter, cols='cell.fract', order=-1)
                        if(tmp.thold>to.filter[, .N]) tmp.thold <- to.filter[, .N]
                        to.filter <- to.filter[1:tmp.thold, .(donor.id.tag, clonotype.tag)]
                        return(to.filter)
                    })
                    to.filter <- rbindlist(l=to.filter, use.names=TRUE)
                    to.filter <- unique(to.filter)
                }else{
                    if(this.att=='att.3'){
                        to.filter <- tmp.data[,
                            .(freq.total=sum(cell.count)),
                            by=.(donor.id.tag, clonotype.tag)
                        ]
                        to.filter <- to.filter[freq.total>=tmp.thold, .(donor.id.tag, clonotype.tag)]
                    }else{
                        stop('Unvalid attempt label. Valid labels are att.1 and att2.\n')
                    }
                }
            }
            tmp.check <- to.filter[, .N<10]
            if(tmp.check){
                warning('Less than 10 clonotypes were recovered for set filters.\n')
                next
            }else{
                cat(paste0('Unique number of TCRs: ', to.filter[, .N], '\n'))
            }
            tmp.data <- merge(
                x=to.filter, y=tmp.data,
                by=c('donor.id.tag', 'clonotype.tag')
            )
            # Final input format.
            tmp.data <- tmp.data[, .(donor.id.tag, clonotype.tag, cluster, cell.fract)]
            tmp.data <- spread(data=tmp.data, key=cluster, value=cell.fract, fill=0)
            row.names(tmp.data) <- paste(tmp.data$donor.id.tag, tmp.data$clonotype.tag, sep='.')
            tmp.data$donor.id.tag <- NULL; tmp.data$clonotype.tag <- NULL
            tmp.data <- as.matrix(tmp.data)
            # tmp.data <- tmp.data[, rev(x=colnames(tmp.data))]
            # Hierarchical clustering.
            #   Rows
            if(gen.row.clust){
                dist.mat <- vegan::vegdist(x=tmp.data, method='bray')
                row.clust <- hclust(d=dist.mat, method='average')
            }else{
                row.clust <- FALSE
            }
            #   Columns
            if(gen.col.clust){
                dist.mat <- vegan::vegdist(x=t(tmp.data), method='bray')
                col.clust <- hclust(d=dist.mat, method='average')
            }else{
                col.clust <- FALSE
            }
            # Upper fraction bound
            if(is.null(set.up.bounds) | !(t.subset %in% names(set.up.bounds))){
                tmp.bound <- as.vector(tmp.data)
                tmp.bound <- tmp.bound[tmp.bound>0]
                tmp.bound <- quantile(x=tmp.bound, probs=0.99)
            }else{
                if(is.numeric(set.up.bounds)){
                    tmp.bound <- set.up.bounds
                }else{
                    tmp.bound <- set.up.bounds[t.subset]
                }
            }
            tmp.data[tmp.data>tmp.bound] <- tmp.bound
            # Set metadata.
            col.metadata <- data.frame(
                row.names=colnames(tmp.data),
                `Cluster`=colnames(tmp.data)
            )
            # Define colors for metadata tracks.
            tmp.cols <- clusts.cols[[t.subset]]
            ann.colors <- list(
                `Cluster`=tmp.cols
            )
            # Set color scale and breaks for heatmap.
            break.no <- 200
            col.breaks <- seq(from=min(range(tmp.data, na.rm=TRUE)), to=max(range(tmp.data, na.rm=TRUE)), length.out=break.no)
            # mid.point <- .005
            # mid.point <- which.min(abs(col.breaks-mid.point))
            # hmap.col.scale.1 <- colorRampPalette(c('#FFFFFF', '#414487', '#2A788E'))(mid.point)
            # hmap.col.scale.2 <- colorRampPalette(c('#2A788E', '#22A884', '#FDE725', '#FFEA00'))(break.no-(mid.point+1))
            # hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
            hmap.col.scale <- colorRampPalette(c('#FFFFFF', '#414487', '#2A788E', '#22A884', '#FDE725', '#FFEA00'))(break.no+1)
            # @ Version with hierarchical clustering.
            # Complete version.
            tmp.file.name <- paste0(
                tmp.reports.path,
                '/', obj.extended.names[t.subset], '_Clonotype_Cluster_CellFract_Thold-', tmp.thold, '.C.pdf'
            )
            tmp.plot <- pheatmap(
                mat=tmp.data, scale='none',
                color=hmap.col.scale, breaks=col.breaks, na_col='#FFFFFF', border_color='black',
                cluster_rows=row.clust, cluster_cols=col.clust,
                annotation_col=col.metadata, annotation_colors=ann.colors,
                show_colnames=FALSE, show_rownames=FALSE,
                legend=TRUE, annotation_legend=TRUE, annotation_names_row=TRUE,
                filename=NA
            )
            pdf(file=tmp.file.name, height=height, width=width)
            grid.draw(rectGrob(gp=gpar(fill='white', lwd=0)))
            grid.draw(tmp.plot)
            dev.off()
            # grid.gedit("layout", gp = gpar(col = "white", text = ""))
            # Blank version.
            tmp.file.name <- paste0(
                tmp.reports.path,
                '/', obj.extended.names[t.subset], '_Clonotype_Cluster_CellFract_Thold-', tmp.thold, '.B.pdf'
            )
            tmp.plot <- pheatmap(
                mat=tmp.data, scale='none',
                color=hmap.col.scale, breaks=col.breaks, na_col='#FFFFFF', border_color='black',
                cluster_rows=row.clust, cluster_cols=col.clust,
                annotation_col=col.metadata, annotation_colors=ann.colors,
                show_colnames=FALSE, show_rownames=FALSE,
                legend=FALSE, annotation_legend=FALSE, annotation_names_row=FALSE, annotation_names_col=FALSE,
                filename=NA
            )
            pdf(file=tmp.file.name, height=height, width=width)
            grid.draw(rectGrob(gp=gpar(fill='white', lwd=0)))
            grid.draw(tmp.plot)
            dev.off()
            # @ Color legend only.
            # W/ labels
            tmp.ggplot <- ggplot(data=tmp.data, aes(x=`0`, y=`0`, col=`0`)) + 
                geom_point() +
                # scale_color_gradientn(colors=hmap.col.scale, breaks=col.breaks[c(65, 130, 195)]) +
                scale_color_gradientn(colors=hmap.col.scale) +
                theme(
                    legend.ticks=element_line(color='black', linewidth=0.6),
                    legend.ticks.length=unit(0.22, "cm"),
                    legend.frame=element_rect(color='black', linewidth=0.6)
                )
            tmp.ggplot <- get_legend(p=tmp.ggplot)
            tmp.file.name <- paste0(
                tmp.reports.path,
                '/', obj.extended.names[t.subset], '_Clonotype_Cluster_CellFract_Thold-', tmp.thold, '.L1.pdf'
            )
            pdf(file=tmp.file.name, height=2, width=2)
            print(as_ggplot(tmp.ggplot))
            dev.off()
            # W/out labels
            tmp.ggplot <- ggplot(data=tmp.data, aes(x=`0`, y=`0`, col=`0`)) + 
                geom_point() +
                # scale_color_gradientn(colors=hmap.col.scale, breaks=col.breaks[c(65, 130, 195)]) +
                scale_color_gradientn(colors=hmap.col.scale, name=NULL, labels=NULL) +
                theme(
                    legend.ticks=element_line(color='black', linewidth=0.6),
                    legend.ticks.length=unit(0.22, "cm"),
                    legend.frame=element_rect(color='black', linewidth=0.6)
                )
            tmp.ggplot <- get_legend(p=tmp.ggplot)
            tmp.file.name <- paste0(
                tmp.reports.path,
                '/', obj.extended.names[t.subset], '_Clonotype_Cluster_CellFract_Thold-', tmp.thold, '.L2.pdf'
            )
            pdf(file=tmp.file.name, height=1.25, width=0.3)
            print(as_ggplot(tmp.ggplot))
            dev.off()
        }
    }
    return(NULL)
}

# ---> General definition of relevant populations per dataset.
set.rel.pops <- list(
    `c40.lg`=list(
        `c40.lg.cd4.d40`=c('0', '1', '2', '3', '4'),
        `c40.lg.cd8.d40`=c('0', '1', '2', '3', '4')
    ),
    `c40.ln`=list(
        `c40.ln.cd4.d08`=c('0', '1', '2', '3', '4', '5', '6'),
        `c40.ln.cd8.d08`=c('0', '1', '2', '3', '4', '5')
    )
)

# ---> Plastic vs subset-biased clonotypes, attempt 1.
# Attempt description. Consider all clonotypes with freq. above certain threshold.
# Temporary directory
tmp.reports.path <- paste0(fig.phen.comps.path, '/clone_plast_concept_1')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
# Predefinitions.
# set.up.bounds <- list(
#     `c40.lg`=c(
#         `c40.lg.cd4.d40`=0.01,
#         `c40.lg.cd8.d40`=0.01
#     ),
#     `c40.ln`=c(
#         `c40.ln.cd4.d08`=0.01,
#         `c40.ln.cd8.d08`=0.01
#     )
# )
freq.tholds <- c(0.0005, 0.001, 0.005)
# Perform.
for(data.set in names(aggr.set.main)){
    gen.pvb.perform(
        this.att='att.1', tmp.tholds=freq.tholds, data.set=data.set,
        set.rel.pops=set.rel.pops[[data.set]], set.up.bounds=0.90,
        tmp.reports.path=tmp.reports.path
    )
}

# ---> Plastic vs subset-biased clonotypes, attempt 2.
# Attempt description. Consider top n clonotypes for each subset.
# Temporary directory
tmp.reports.path <- paste0(fig.phen.comps.path, '/clone_plast_concept_2')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
# Predefinitions.
freq.tholds <- c(50, 100, 150, 200, 250)
# Perform.
for(data.set in names(aggr.set.main)){
    gen.pvb.perform(
        this.att='att.2', tmp.tholds=freq.tholds, data.set=data.set,
        set.rel.pops=set.rel.pops[[data.set]], set.up.bounds=0.90,
        tmp.reports.path=tmp.reports.path
    )
}

# ---> Plastic vs subset-biased clonotypes, attempt 3.
# Attempt description. Consider all clonotypes with freq. above certain threshold.
# Temporary directory
tmp.reports.path <- paste0(fig.phen.comps.path, '/clone_plast_concept_3')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
# Predefinitions.
freq.tholds <- c(2, 10)
# Perform.
for(data.set in names(aggr.set.main)){
    gen.pvb.perform(
        this.att='att.3', tmp.tholds=freq.tholds, data.set=data.set,
        set.rel.pops=set.rel.pops[[data.set]], set.up.bounds=0.90,
        gen.row.clust=TRUE, gen.col.clust=FALSE,
        tmp.reports.path=tmp.reports.path, height=6, width=7
    )
}


# ---> Plastic vs subset-biased clonotypes, attempt 4

# @ General functions
# For circos initialization
pvb.init.circos <- function(segment.data, freq.rel.thold){ # pvb stands for plastic vs biased
    # Initialize the Circos plot
    tmp.data <- segment.data[, .(width=sum(freq.rel)), by=.(t.pop)]
    tmp.vals <- Reduce(x=tmp.data$width, f=sum, accumulate=TRUE)
    tmp.vals <- c(0, tmp.vals)
    tmp.data[, `:=`(
        pos.start=tmp.vals[1:length(tmp.vals)-1],
        pos.end=tmp.vals[2:length(tmp.vals)]
    )]
    circos.par('track.height'=0.3)
    circos.initialize(
        sectors=tmp.data[, t.pop],
        sector.width=tmp.data[, width],
        xlim=as.matrix(tmp.data[, .(pos.start, pos.end)])
    )
    # Add sector tracks
    circos.track(ylim=c(0, 1), panel.fun=function(x, y){
        sector <- CELL_META$sector.index
        xlim <- CELL_META$xlim
        ylim <- CELL_META$ylim
        circos.rect(
            xlim[1], 0, xlim[2], 1,
            col=clusts.cols[[tmp.lin]][sector]
        )
        # Lines showing clonotype relative frequency.
        tmp.data <- segment.data[t.pop==sector & freq.rel>freq.rel.thold, .N]>0
        if(tmp.data){
            circos.lines(
                x=segment.data[t.pop==sector & freq.rel>freq.rel.thold, freq.cumm],
                y=rep(x=0.4, times=segment.data[t.pop==sector & freq.rel>freq.rel.thold, .N]),
                type='h', col='black'
            )
            # Rectangle distinguishing between singletons and clonotypes (expanded TCRs).
            tmp.lim <- segment.data[t.pop==sector & freq.rel>freq.rel.thold, max(freq.cumm)]
            circos.rect(
                xleft=xlim[1], ybottom=0, xright=tmp.lim, ytop=0.4,
                col=adjustcolor("#FFFFFF", alpha.f=0)
            )
            circos.rect(
                xleft=tmp.lim, ybottom=0, xright=xlim[2], ytop=0.4,
                col=adjustcolor("#FFFFFF", alpha.f=1)
            )
        }else{
            # Rectangle indicating singletons across segment.
            circos.rect(
                xleft=xlim[1], ybottom=0, xright=xlim[2], ytop=0.4,
                col=adjustcolor("#FFFFFF", alpha.f=1)
            )
        }
    }, track.height=0.22, bg.border=NA)
}
# To add links.
pvb.add.links <- function(link.data){
    for(idx in 1:nrow(link.data)){
        tmp.link <- unlist(link.data[idx, ])
        coords.1 <- unlist(segment.data[
            donor.id.tag==tmp.link['donor.id.tag'] &
            clonotype.tag==tmp.link['clonotype.tag'] &
            t.pop==tmp.link['t.pop.1'],
            .(start=freq.cumm-freq.rel, end=freq.cumm)
        ])
        coords.2 <- unlist(segment.data[
            donor.id.tag==tmp.link['donor.id.tag'] &
            clonotype.tag==tmp.link['clonotype.tag'] &
            t.pop==tmp.link['t.pop.2'],
            .(start=freq.cumm-freq.rel, end=freq.cumm)
        ])
        circos.link(
            sector.index1=tmp.link['t.pop.1'], point1=coords.1,
            sector.index2=tmp.link['t.pop.2'], point2=coords.2,
            col=adjustcolor("#FDFD96", alpha.f=0.3),
            border = "black"
        )
    }
}
# @ Process per lineage and frequency threshold.
tmp.reports.path <- paste0(fig.phen.comps.path, '/clone_plast_concept_4')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
data.set <- 'c40.lg'
# tmp.tholds <- c(0.0001, 0.0005, 0.001, 0.005, 0.01)
tmp.tholds <- c(0.0005, 0.001, 0.005, 0.01)
for(tmp.lin in aggr.set.main[[data.set]]){
    # ---> Segment-specific info
    tmp.data <- meta.data.l[[data.set]][data.set==tmp.lin]
    tmp.data.1 <- tmp.data[
        !is.na(clonotype.tag) & !is.na(donor.id.tag),
        .(freq.abs=uniqueN(barcode)),
        by=.(
            clonotype.tag, donor.id.tag,
            t.pop=as.character(get(clust.labs[tmp.lin]))
        )
    ]
    tmp.data.2 <- tmp.data[
        !is.na(clonotype.tag) & !is.na(donor.id.tag),
        .(freq.total=uniqueN(barcode)),
        by=.(donor.id.tag)
    ]
    segment.data <- merge(x=tmp.data.1, y=tmp.data.2, by='donor.id.tag')
    segment.data[, freq.rel:=freq.abs/freq.total]
    setorderv(x=segment.data, cols=c('t.pop', 'freq.rel'), order=c(1, -1))
    # Determine cummulative relative frequency.
    segment.data$freq.cumm <- Reduce(x=segment.data$freq.rel, f=sum, accumulate=TRUE)
    # ---> Link-specific info.
    tmp.vals <- segment.data[, unique(t.pop)]
    link.data <- lapply(X=1:(length(tmp.vals)-1), FUN=function(idx.1){
        tmp.data <- lapply(X=(idx.1+1):length(tmp.vals), FUN=function(idx.2){
            tmp.cols <- c('donor.id.tag', 'clonotype.tag', 't.pop', 'freq.rel')
            tmp.data.1 <- segment.data[t.pop==tmp.vals[idx.1], ..tmp.cols]
            tmp.data.2 <- segment.data[t.pop==tmp.vals[idx.2], ..tmp.cols]
            tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by=c('donor.id.tag', 'clonotype.tag'), all=FALSE, suffixes=c('.1', '.2'))
            return(tmp.data)
        })
        tmp.data <- rbindlist(l=tmp.data, use.names=TRUE)
        return(tmp.data)
    })
    link.data <- rbindlist(l=link.data, use.names=TRUE)
    # ---> Process per relative frequency threshold.
    for(freq.rel.thold in tmp.tholds){
        # ---> Circos plot, option 1
        tmp.file.name <- paste0(
            tmp.reports.path, '/',
            aggr.set.labs[data.set], '_SimBetSites_Circos-1', '_Lin-', gen.cell.types[tmp.lin], '_PctThold-', freq.rel.thold*100, '.B.pdf'
        )
        pdf(file=tmp.file.name, width=14, height=14)
        # Initialize the Circos plot
        pvb.att.3.init.circos(segment.data=segment.data, freq.rel.thold=freq.rel.thold)
        # Create links between sectors (representing shared TCRs)
        tmp.data <- link.data[freq.rel.1>freq.rel.thold & freq.rel.2>freq.rel.thold]
        pvb.att.3.add.links(link.data=tmp.data)
        # Clear the plot
        dev.off()
        circos.clear()
        # ---> Circos plot, option 2
        # tmp.file.name <- paste0(
        #     tmp.reports.path, '/',
        #     aggr.set.labs[data.set], '_SimBetSites_Circos-2', '_Lin-', gen.cell.types[tmp.lin], '_PctThold-', freq.rel.thold*100, '.B.pdf'
        # )
        # pdf(file=tmp.file.name, width=14, height=14)
        # # Initialize the Circos plot
        # pvb.att.3.init.circos(segment.data=segment.data, freq.rel.thold=freq.rel.thold)
        # # Create links between sectors (representing shared TCRs)
        # tmp.data <- link.data[freq.rel.1>freq.rel.thold | freq.rel.2>freq.rel.thold]
        # pvb.att.3.add.links(link.data=tmp.data)
        # # Clear the plot
        # dev.off()
        # circos.clear()
    }
}

# ---> Formal determination of subset-biased and plastic clonotypes. Cell counts per subset.
# Potentially repetitive. Remove if necessary.
set.rel.pops <- list(
    `c40.lg`=list(
        `c40.lg.cd4.d40`=c('0', '1', '2', '3', '4'),
        `c40.lg.cd8.d40`=c('0', '1', '2', '3', '4')
    ),
    `c40.ln`=list(
        `c40.ln.cd4.d08`=c('0', '1', '2', '3', '4', '5', '6'),
        `c40.ln.cd8.d08`=c('0', '1', '2', '3', '4', '5')
    )
)
for(data.set in names(aggr.set.labs)){
    tmp.vals <- meta.data.l[[data.set]][, unique(data.set)]
    for(t.subset in tmp.vals){
        # @ Retrieve data.
        tmp.data <- meta.data.l[[data.set]][gex.lin.tag==gen.cell.types[t.subset]]
        tmp.data <- tmp.data[
            !is.na(clonotype.tag) & !is.na(donor.id.tag),
            .(
                barcode, donor.id.tag, clonotype.tag,
                cluster=as.character(get(clust.labs[t.subset])),
                gen.clonotype.status=gen.clonotype.status,
                group.bias=as.character(group.class.tag)
            )
        ]
        tmp.data <- tmp.data[cluster %in% set.rel.pops[[data.set]][[t.subset]]]
        # @ Group assignment.
        tmp.data[
            cluster==group.bias,
            group.bias:='Biased to same cluster'
        ]
        tmp.data[
            !is.na(group.bias) & cluster!=group.bias & group.bias!='Biased to same cluster',
            group.bias:='Biased to diff. cluster'
        ]
        tmp.data[
            gen.clonotype.status=='Singleton',
            group.bias:='Singleton'
        ]
        tmp.data[
            is.na(group.bias),
            group.bias:='Plastic'
        ]
        # @ Subset - Highly expanded clones.
        tmp.data.1 <- tmp.data[,
            .(freq.abs=uniqueN(barcode)),
            by=.(donor.id.tag, clonotype.tag)
        ]
        tmp.data <- merge(
            x=tmp.data, y=tmp.data.1,
            by=c('donor.id.tag', 'clonotype.tag')
        )
        tmp.data <- tmp.data[freq.abs>=hec.freq.thold]
        # @ Set cluster order
        # tmp.lvls <- rev(names(clusts.cols[[t.subset]]))
        tmp.lvls <- names(clusts.cols[[t.subset]])
        tmp.vals <- factor(x=as.character(tmp.data$cluster), levels=tmp.lvls)
        set(x=tmp.data, j='cluster', value=tmp.vals)
        # @ Set cell types (according to bias) order
        tmp.lvls <- names(bias.1.cols)
        tmp.vals <- factor(x=as.character(tmp.data$group.bias), levels=tmp.lvls)
        set(x=tmp.data, j='group.bias', value=tmp.vals)
        # Define file width according to cluster count.
        # tmp.width <- (tmp.data[, uniqueN(cluster)] * 6.68/40) + 0.32
        # tmp.width <- (tmp.width * 5)/1.8
        tmp.width <- 3.672222 * (5 / tmp.data[, uniqueN(cluster)])
        # Plot, v1.
        tmp.ggplot <- ggplot(data=tmp.data, aes(x=cluster, fill=group.bias)) +
            # geom_bar(position='stack', width=0.8, linewidth=0, alpha=1) +
            geom_bar(position='fill', width=0.8, linewidth=0, alpha=1) +
            scale_fill_manual(values=bias.1.cols) +
            scale_y_continuous(expand=c(0, 0), breaks=scales::pretty_breaks(n=3)) +
            labs(x='Cluster', y='Cell fraction', fill='Clonotype status')
        tmp.lab <- paste0('/', obj.extended.names[t.subset], '_CellFract_Cluster_BiasStatus-Dir')
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.phen.comps.path, file.name=tmp.lab, type='pdf',
            blank.comp=blank.complement.7.1, do.legend=FALSE, do.rotate=TRUE, width=tmp.width*3, height=3.5*3
        )
        # Plot, v2.
        tmp.data[, clonotype.id:=paste(clonotype.tag, donor.id.tag, sep=';')]
        tmp.check <- tmp.data[, uniqueN(clonotype.id)]
        cat(paste0(obj.extended.names[t.subset], ': ', tmp.check, '\n'))
        tmp.vals <- c('Cell', 'Clone')
        for(tmp.item in tmp.vals){
            if(tmp.item=='Cell'){
                to.plot <- tmp.data[,
                    .(
                        bias.fract=.SD[group.bias=='Biased to same cluster', .N]/.N
                    ),
                    by=cluster
                ]
            }else{
                to.plot <- tmp.data[,
                    .(
                        bias.fract=.SD[group.bias=='Biased to same cluster', uniqueN(clonotype.id)]/uniqueN(clonotype.id)
                    ),
                    by=cluster
                ]
            }
            tmp.lvls <- rev(levels(tmp.data$cluster))
            tmp.vals <- factor(x=as.character(tmp.data$cluster), levels=tmp.lvls)
            set(x=tmp.data, j='cluster', value=tmp.vals)
            # Define file width according to cluster count.
            tmp.width <- 3.672222 * (5 / to.plot[, uniqueN(cluster)])
            tmp.ggplot <- ggplot(data=to.plot, aes(x=cluster, y=bias.fract)) +
                # geom_bar(position='stack', width=0.8, linewidth=0, alpha=1) +
                geom_bar(aes(fill=cluster), stat='identity', position='dodge', width=0.8, linewidth=0, alpha=1) +
                scale_fill_manual(values=clusts.cols[[t.subset]]) +
                scale_y_continuous(limits=c(0, 1), expand=c(0, 0), breaks=scales::pretty_breaks(n=3)) +
                labs(x='Cluster', y=paste0(tmp.item, ' fraction'), fill='Clonotype status')
            tmp.lab <- paste0('/', obj.extended.names[t.subset], '_', tmp.item, 'Fract_Cluster_BiasStatus-Dir-2')
            publish.plot(
                tmp.ggplot=tmp.ggplot, output.path=fig.phen.comps.path, file.name=tmp.lab, type='pdf',
                blank.comp=blank.complement.7, do.legend=FALSE, do.rotate=TRUE, width=tmp.width*3, height=2.8*3
            )
        }
    }
}

# ---> Percentage of T cell clonally overlapped between anatomic sites per subset.
for(data.set in names(aggr.set.labs)){
    tmp.vals <- meta.data.l[[data.set]][, unique(data.set)]
    for(t.subset in tmp.vals){
        anti.set <- if(data.set=='c40.lg') 'c40.ln' else 'c40.lg'
        anti.subset <- aggr.set.main[[anti.set]]
        anti.subset <- names(gen.cell.types[anti.subset])[gen.cell.types[anti.subset]==gen.cell.types[t.subset]]
        # Collect details on anti dataset, including bias info.
        tmp.data <- get.multi.site.info(
            data.set=data.set, t.subset=t.subset,
            anti.set=anti.set, anti.subset=anti.subset,
            rel.cols=NULL
        )
        # @ Set cluster order
        tmp.lvls <- rev(names(clusts.cols[[t.subset]]))
        tmp.vals <- factor(x=as.character(tmp.data$cluster), levels=tmp.lvls)
        set(x=tmp.data, j='cluster', value=tmp.vals)
        # Define file width according to donor amount.
        tmp.width <- (tmp.data[, uniqueN(cluster)] * 6.68/40) + 0.32
        tmp.width <- (tmp.width * 5)/1.8
        # ---> Donor-wise information.
        tmp.vals <- c('Cell', 'Clone')
        for(tmp.item in tmp.vals){
            # Overlap fraction per subset.
            if(tmp.item=='Cell'){
                to.plot <- tmp.data[,
                    .(ovlp.freq.rel=.SD[site.ovlp=='Overlapped', .N]/.N),
                    by=.(donor.id.tag, cluster)
                ]
            }else{
                to.plot <- tmp.data[,
                    .(ovlp.freq.rel=.SD[site.ovlp=='Overlapped', uniqueN(clonotype.tag)]/uniqueN(clonotype.tag)),
                    by=.(donor.id.tag, cluster)
                ]
            }
            # Plot.
            tmp.ggplot <- ggplot(data=to.plot, aes(x=cluster, y=ovlp.freq.rel)) +
                geom_jitter(shape=1, stroke=5, size=10, color=gen.dot.col, width=0, height=0) +
                geom_boxplot(aes(color=cluster, fill=cluster), alpha=0.4, width=0.7, linewidth=4, outlier.shape=NA) +
                scale_color_manual(values=clusts.cols[[t.subset]]) +
                scale_fill_manual(values=clusts.cols[[t.subset]]) +
                scale_y_continuous(expand=c(0, 0), limits=c(0, 1), breaks=scales::pretty_breaks(n=3)) +
                labs(x='Cluster', y=paste0(tmp.item, ' fraction'), fill='Cluster')
            tmp.lab <- paste0('/', obj.extended.names[t.subset], '_Ovlped', tmp.item, 'Fract_Cluster-Donor')
            publish.plot(
                tmp.ggplot=tmp.ggplot, output.path=fig.phen.comps.path, file.name=tmp.lab, type='pdf',
                blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=tmp.width*3, height=14
            )
        }
        # ---> Summary across donors.
        tmp.vals <- c('Cell', 'Clone')
        for(tmp.item in tmp.vals){
            # Overlap fraction per subset.
            if(tmp.item=='Cell'){
                to.plot <- tmp.data[,
                    .(ovlp.freq.rel=.SD[site.ovlp=='Overlapped', .N]/.N),
                    by=cluster
                ]
            }else{
                to.plot <- tmp.data[,
                    .(ovlp.freq.rel=.SD[site.ovlp=='Overlapped', uniqueN(clonotype.tag)]/uniqueN(clonotype.tag)),
                    by=cluster
                ]
            }
            # Plot.
            tmp.ggplot <- ggplot(data=to.plot, aes(x=cluster, y=ovlp.freq.rel, fill=cluster)) +
                # geom_bar(stat='identity', width=0.8, linewidth=1.2, alpha=1, color='black') +
                geom_bar(stat='identity', width=0.8, linewidth=0, alpha=1, color='black') +
                scale_fill_manual(values=clusts.cols[[t.subset]]) +
                scale_y_continuous(expand=c(0, 0), limits=c(0, 1), breaks=scales::pretty_breaks(n=3)) +
                labs(x='Cluster', y=paste0(tmp.item, ' fraction'), fill='Cluster')
            tmp.lab <- paste0('/', obj.extended.names[t.subset], '_Ovlped', tmp.item, 'Fract_Cluster')
            publish.plot(
                tmp.ggplot=tmp.ggplot, output.path=fig.phen.comps.path, file.name=tmp.lab, type='pdf',
                blank.comp=blank.complement.7, do.legend=FALSE, do.rotate=TRUE, width=tmp.width, height=5
            )
        }
    }
}

# ---> Percentage of T cell clonally overlapped between anatomic sites per subset.
tmp.tholds <- c(
    `CD4`=0.008,
    `CD8`=0.02
)
for(tmp.lin in unique(gen.cell.types)){
    # Identify unique sets from each aggr.
    lin.sets <- unlist(lapply(X=aggr.set.main, FUN=function(x){
        names(gen.cell.types[x])[gen.cell.types[x]==tmp.lin]
    }))
    # Collect details for each site as reference.
    tmp.data <- list(
        `LG`=c('c40.lg', 'c40.ln'),
        `LN`=c('c40.ln', 'c40.lg')
    )
    tmp.data <- lapply(X=tmp.data, FUN=function(tmp.ids){
        tmp.data <- get.multi.site.info(
            data.set=tmp.ids[1], t.subset=lin.sets[tmp.ids[1]],
            anti.set=tmp.ids[2], anti.subset=lin.sets[tmp.ids[2]],
            rel.cols=NULL
        )
        tmp.data.1 <- tmp.data[,
            .(freq.abs=.N),
            by=.(comb.ids, donor.id.tag)
        ]
        tmp.data.2 <- tmp.data[,
            .(freq.total=.N),
            by=.(donor.id.tag)
        ]
        tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='donor.id.tag')
        tmp.data[, freq.rel:=freq.abs/freq.total]
        tmp.data[, `:=`(freq.abs=NULL, freq.total=NULL)]
        return(tmp.data)
    })
    tmp.data <- merge(
        x=tmp.data[[1]], y=tmp.data[[2]],
        by=c('donor.id.tag', 'comb.ids'),
        suffixes=c('.LG', '.LN')
    )
    tmp.test <- wilcox.test(x=tmp.data[, freq.rel.LG], y=tmp.data[, freq.rel.LN], paired=TRUE)
    tmp.caption <- paste0('P-value is ', tmp.test$p.value)
    to.plot <- as.data.table(gather(data=tmp.data, key=site, value='freq.rel', -`donor.id.tag`, -`comb.ids`))
    to.plot[, site:=str_replace(string=site, pattern='freq.rel.', replacement='')]
    to.plot[, group:=paste(donor.id.tag, comb.ids, sep='.')]
    # Set upper bound.
    tmp.thold <- tmp.tholds[tmp.lin]
    to.plot[, bounded.freq.rel:=freq.rel]
    to.plot[, point.status:='ori']
    to.plot[bounded.freq.rel>tmp.thold, `:=`(bounded.freq.rel=tmp.thold, point.status='bounded')]
    # Plot.
    tmp.ggplot <- ggplot(data=to.plot, aes(x=site)) +
        geom_boxplot(aes(y=freq.rel, color=site), width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
        geom_jitter(aes(y=bounded.freq.rel, shape=point.status), stroke=3, size=8, color='black', width=0) +
        geom_line(aes(y=bounded.freq.rel, group=group), linewidth=0.2) +
        coord_cartesian(ylim=c(0, tmp.thold)) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=4)) +
        scale_shape_manual(values=tmp.shapes) +
        scale_color_manual(values=an.site.cols) +
        labs(x='Anatomic site', y='Clonotype relative freq. within site', color='', caption=tmp.caption)
    tmp.lab <- paste0(
        '/', tmp.lin, '_Site_CloneRelFreq'
    )
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.phen.comps.path, file.name=tmp.lab, type='pdf',
        blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=7, height=14
    )
}

# ---> External phenotype of clonotypes overlapped between anatomical sites.
for(data.set in names(aggr.set.labs)){
    tmp.vals <- meta.data.l[[data.set]][, unique(data.set)]
    for(t.subset in tmp.vals){
        anti.set <- if(data.set=='c40.lg') 'c40.ln' else 'c40.lg'
        anti.subset <- aggr.set.main[[anti.set]]
        anti.subset <- names(gen.cell.types[anti.subset])[gen.cell.types[anti.subset]==gen.cell.types[t.subset]]
        # Collect details on anti dataset, including bias info.
        tmp.data <- get.multi.site.info(
            data.set=data.set, t.subset=t.subset,
            anti.set=anti.set, anti.subset=anti.subset,
            rel.cols=NULL
        )
        # To plot.
        to.plot <- tmp.data[site.ovlp=='Overlapped' & !is.na(group.bias)]
        anti.trm.clust <- names(clusts.defs[[anti.subset]])[clusts.defs[[anti.subset]]=='TRM']
        to.plot[, summ.group.bias:='Other']
        to.plot[group.bias==anti.trm.clust, summ.group.bias:='TRM']
        # @ Set cluster order
        tmp.lvls <- rev(names(clusts.cols[[t.subset]]))
        tmp.vals <- factor(x=as.character(to.plot$cluster), levels=tmp.lvls)
        set(x=to.plot, j='cluster', value=tmp.vals)
        # @ Set cell types (according to bias) order
        tmp.cols <- c(
            clusts.cols[[anti.subset]],
            bias.1.cols[names(bias.1.cols)=='Plastic']
        )
        tmp.lvls <- names(tmp.cols)
        tmp.vals <- factor(x=as.character(to.plot$group.bias), levels=tmp.lvls)
        set(x=to.plot, j='group.bias', value=tmp.vals)
        # Define file width according to donor amount.
        tmp.width <- (to.plot[, uniqueN(cluster)] * 6.68/40) + 0.32
        tmp.width <- (tmp.width * 5)/1.8
        # Plot.
        tmp.ggplot <- ggplot(data=to.plot, aes(x=cluster, fill=group.bias)) +
            geom_bar(position='stack', width=0.8, linewidth=0, alpha=1) +
            scale_fill_manual(values=tmp.cols) +
            scale_y_continuous(expand=c(0, 0), breaks=scales::pretty_breaks(n=3)) +
            labs(x='Cluster', y='Cell count', fill='Clonotype status')
        tmp.lab <- paste0('/', obj.extended.names[t.subset], '_CellFract_Cluster_BiasStatus-Anti')
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.phen.comps.path, file.name=tmp.lab, type='pdf',
            blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=tmp.width*3, height=14
        )
        tmp.cols <- c(gen.subset.cols, `Other`='#CBCBCB')
        tmp.ggplot <- ggplot(data=to.plot, aes(x=cluster, fill=summ.group.bias)) +
            geom_bar(position='stack', width=0.8, linewidth=0, alpha=1) +
            scale_fill_manual(values=tmp.cols) +
            scale_y_continuous(expand=c(0, 0), breaks=scales::pretty_breaks(n=3)) +
            labs(x='Cluster', y='Cell count', fill='Clonotype status')
        tmp.lab <- paste0('/', obj.extended.names[t.subset], '_CellFract_Cluster_BiasStatus-Anti-Summ')
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.phen.comps.path, file.name=tmp.lab, type='pdf',
            blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=tmp.width*3, height=14
        )
    }
}


### -------------------------- Supp Figure -------------------------- ###

# ---> Dot plots.
# Define genes of interest.
these.markers <- list(
    `c40.lg.cd4.d40`=c(
        # ---> THIFNR
        # 'IFI6', 'MX1',
        # ---> Proliferating
        'MKI67', 'TOP2A',
        # ---> TRM
        'CXCR6', 'ITGAE', 'ZNF683', 'IFNG', 'HOPX',
        # ---> Classical CD4-CTL
        'PRF1', 'GZMB',
        # ---> TCM
        'KLF2', 'TCF7', 'SELL',
        # ---> TREG
        'FOXP3', 'IL2RA', 'TIGIT', 'CTLA4',
        # ---> TFH
        'PDCD1', 'CXCL13'
    ),
    `c40.lg.cd8.d40`=c(
        # ---> Proliferating
        'MKI67', 'TOP2A',
        # ---> TRM
        'CXCR6', 'ITGAE', 'IFNG', 'ZNF683',
        # ---> NKG2C+ TRM
        'KLRC2', 'AREG', 'KIR2DL4',
        # ---> GZMKhi
        'CRTAM', 'TNFSF9', 'GZMK',
        # ---> Effectors
        'FGFBP2', 'PRF1', 'GZMB',
        # ---> TCM
        'KLF2', 'SELL', 'TCF7',
        # ---> MAIT
        'KLRB1', 'TRAV1-2', 'TRBV20-1'
    ),
    `c40.lg.cd8.trm`=c(
        # ---> BATFhi
        'HLA-DRB1', 'HLA-DRB5', 'HLA-DRA', 'HLA-DPA1', 'HLA-DQA1',
        # 'HLA-DQB1', 'HLA-DQA2', 'HLA-DMA',
        'ZNF683', 'CXCR6', 'FABP5',
        # 'ITGB1', 'ITGB2', 'APOBEC3G',
        'BATF', 'RBPJ',
        # 'LGALS1','LGALS3',
        'IFNG', 'GZMB',
        # 'GZMA', 'GZMH', 'MIR155HG'
        # ---> IL7Rhi
        'IL7R', 'CD27', 'TCF7'
        # 'CD55', 'KLF2', 'KLRC1', 'KLRC2', 'GPR183'
    ),
    `c40.lg.cd8.tem`=c(
        # ---> GZMKpos
        'GZMK', # 'CD74', 'DUSP2',
        'IL7R', 'TCF7', # 'CD27',
        # 'CMC1', 'LTB', 'RBFOX2', 'CAMK4',
        # ---> ZNF683pos
        'GZMB', 'GNLY', # 'FGFBP2', 'GZMH', 'CX3CR1', 'BHLHE40', 'CD6',
        'IFITM1', 'ISG15', 'IFI6', # 'OASL',
        'ZNF683', # 'ITGB1', 'LGALS1',
        # ---> NKG2Cpos
        'KIR2DL3', 'KLRC3', 'KLRC2',
        'IKZF2', 'TIGIT', # 'TYROBP',
        # 'CCL3', 'XCL2',
        # 'FCGR3A', 'KIR3DL2', 'KLRF1',
        # 'TRGV7', 'TRGV3', 'TRGC2',
        # 'CCL4L2', 'CEP78',
        # ---> TEM 3
        'TRBV7-6', 'TRBV7-4', #'TRAV14DV4',
        # ---> TEM 4
        'TRBV20-1', 'TRAV8-4' #'TRGV9'
    ),
    `c40.ln.cd4.d08`=c(
        # ---> TRM
        # 'LGALS1', 'LGALS3',
        'ITGAE', 'ITGB1', # 'ITGA1',
        # 'ZNF683',
        'CXCR6', 'CXCR3',
        # 'IFNG', 'TNF',
        # ---> TCM / TN
        'S1PR1', 'KLF2', 'SELL',
        'TCF7', # 'LEF1', 'IL7R', 'CCR7',
        # 'IL2', 'KLRB1',
        'GPR183', 'CD40LG',
        # ---> IFNR
        'ISG15', 'MX1', # 'IFIT1', 'IFIT2',
        # ---> TFH
        # 'MS4A6A', 'CD79A', 'FABP5',
        'CD200', 'BTLA',# 'BCL6',
        # 'BCL6',
        'PDCD1', 'CTLA4', # 'ICOS',
        # ---> TREG
        'FOXP3', 'IL2RA' # 'IKZF2',
    ),
    `c40.ln.cd8.d08`=c(
        # ---> TRM
        # 'LGALS1', 'LGALS3', 
        'ZNF683', 'CXCR6', 'CXCR3',
        'ITGAE', # 'ITGA1', 'ITGB1',
        # 'IFNG', 'TNF',
        # ---> TCM
        'S1PR1', 'KLF2',
        'TCF7', # 'LEF1', 'IL7R', # 'CCR7',
        'SELL',
        'GPR183',
        # 'IL2', 'CD40LG',
        # ---> GZMKhi
        'GZMK', 'EOMES', # 'CRTAM',
        # ---> HLApos
        # 'HLA-DMA', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1',
        # ---> IFNR
        'ISG15', 'MX1', # 'IFIT1', 'IFIT2',
        # ---> Elongation factors.
        # 'EEF1G', 'EEF2', 'EEF1B2',
        # ---> MAIT
        'TRAV1-2', 'KLRB1' #, 'TRBV6-1', 'TRBV20-1'
        # ---> TEM
        # 'PRF1', 'GZMB', #'GNLY',
        # 'KLRG1' # 'KLRD1' 'FCGR3A', 'FGFBP2',
    )
)
# Define clusters' order.
these.pops <- list(
    `c40.lg.cd4.d40`=c('5', '0', '2', '1', '3', '4'),
    `c40.lg.cd8.d40`=c('6', '0', '4', '2', '1', '3', '5'),
    `c40.lg.cd8.trm`=c('0', '1'),
    `c40.lg.cd8.tem`=c('0', '1', '2', '3', '4'),
    `c40.ln.cd4.d08`=c('6', '1', '0', '4', '2', '5', '3'),
    `c40.ln.cd8.d08`=c('3', '2', '5', '4', '0', '1', '7')
)
# Process per lineage.
for(t.subset in names(obj.extended.names)){
    # Specific options.
    do.norm <- if(t.subset=='c40.lg.cd8.trm') FALSE else TRUE
    tmp.color.scale <- if(t.subset=='c40.lg.cd8.trm') signatures.col.scale else 'hot.and.cold'
    # Plot
    tmp.ggplot <- dot.plot(
        seurat.obj=srt.objs.list[[t.subset]], features=these.markers[[t.subset]], slot='data', ensembl=FALSE,
        # do.norm=FALSE,
        do.norm=do.norm,
        groups.tag=clust.labs[t.subset], groups.order=these.pops[[t.subset]], groups.of.int=NULL, filter.tags=NULL, groups.to.filter=NULL, keep=FALSE, na.rm=TRUE, feature.thold=NULL,
        this.color.scale=tmp.color.scale,
        col.min=NULL, col.max=NULL,
        scale.by='radius', dot.scale=12,
        size.min=0, size.max=NA,
        file.name=NULL
    )
    tmp.ggplot <- tmp.ggplot + theme(legend.position='bottom')
    # Output.
    tmp.width <- (length(these.pops[[t.subset]]) * 0.5) + 0.3
    # tmp.width <- 5
    tmp.height <- (length(these.markers[[t.subset]]) * 0.4) + 0.3
    tmp.lab <- paste0(obj.extended.names[t.subset], '_Markers_DotPlot_Opt-A')
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.phen.comps.path, file.name=tmp.lab, type='pdf', blank.comp=blank.complement.2, do.legend=TRUE,
        width=tmp.width, height=tmp.height, c.width=tmp.width+2, c.height=tmp.height, legend.height=5, legend.width=0.5
    )
}

# ---> Module scoring calculation (if necessary).
mods.of.int <- names(modules.feats)
for(data.set in names(srt.objs.list)){
    if(!all(mods.of.int %in% colnames(srt.objs.list[[data.set]]@meta.data))){
        mod.in.path <- paste0(module.feats.dir, '/base_signatures')
        # mod.in.path <- paste0(module.feats.dir, '/base_signatures/tmp') # Replace w/ this if a specific subset of signatures needs to be added and they were included as symbolic links in this folder.
        mod.out.path <- paste0(module.feats.dir, '/tmp_dir')
        if(!dir.exists(mod.out.path)) dir.create(mod.out.path)
        srt.objs.list[[data.set]] <- get.module.scores(
            seurat.obj=srt.objs.list[[data.set]], module.feats.dir=mod.in.path, reports.path=mod.out.path,
            signal.thold=0, is.ensembl=FALSE
        )
        # Remove tmp dir.
        tmp.cmmd <- paste0('rm -r ', mod.out.path)
        system(command=tmp.cmmd, intern=FALSE)
    }
}

# ---> Module scoring as UMAP heatmaps.
tmp.reports.path <- paste0(fig.phen.comps.path, '/module_scoring')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
# Define relevant signature per T-cell lineage.
mods.of.int <- list(
    `CD4`=c(
        `Clarke-TRM`='trm.signature.clarke', `Pauken-TCM`='TCM.Pauken.Tumor', `Scott-TRM-Up`='trm.up.scott', `Scott-TRM-Down`='trm.down.scott',
        `Patil-CTL`='cell.cytotoxicity.patil', `Arlehamn-TH1`='th1.signature1.arlehamn', `Arlehamn-TH17`='th17.signature2.arlehamndan', `Locci-TFH`='tfh.signature.locci', `Best-Cylcing`='cell.cycling.best',
        `MSIGDB-IFNR-1-2`='type1n2.interferon.broadreact',
        `DICE-TREGN`='R24.TREGN.RESTING', `DICE-TREGM`='R24.TREGMEM.RESTING', `DICE-TFH`='R24.TFH.RESTING', `DICE-TH2`='R24.TH2.RESTING', `DICE-TH1`='R24.TH1.RESTING', `DICE-TH17`='R24.TH17.RESTING', `DICE-THSTAR`='R24.THSTAR.RESTING', `DICE-CD4N`='R24.CD4N.RESTING'
        # `Clarke-TRM`='trm.signature.clarke', `Pauken-TCM`='TCM.Pauken.Tumor', `Patil-CTL`='cell.cytotoxicity.patil', `DICE-TREGN`='R24.TREGN.RESTING', `DICE-TREGM`='R24.TREGMEM.RESTING'
    ),
    `CD8`=c(
        `Guo-TEM`='tcell.cytotoxicity.guo', `Clarke-TRM`='trm.signature.clarke', `Pauken-TCM`='TCM.Pauken.Tumor',
        `Jonsson-GZMK`='gzmkpos-cd8.jonsson', `Lan-GZMK`='gzmkpos-cd8.lan',
        `Scott-TRM-Up`='trm.up.scott', `Scott-TRM-Down`='trm.down.scott',
        `Cullen-Unhelped`='unhelped.cullen', `Consens-Exhaust`='dixhaust.consensus2.vj', `Best-Cylcing`='cell.cycling.best',
        `MSIGDB-IFNR-1-2`='type1n2.interferon.broadreact', `MSIGDB-Glycolysis`='hallmark.glycolysis'
        `Guo-TEM`='tcell.cytotoxicity.guo', `Clarke-TRM`='trm.signature.clarke', `Pauken-TCM`='TCM.Pauken.Tumor'
    )
)
# Custom score bounds.
gen.scale.tholds <- list(
    `c40.lg.cd4.d40`=list(
        `Clarke-TRM`=0.2, `Pauken-TCM`=0.4, `Scott-TRM-Up`=NA, `Scott-TRM-Down`=NA,
        `Patil-CTL`=0.2, `Arlehamn-TH1`=NA, `Arlehamn-TH17`=NA, `Locci-TFH`=NA, `Best-Cylcing`=NA,
        `MSIGDB-IFNR-1-2`=NA,
        `DICE-TREGN`=0.15, `DICE-TREGM`=0.15, `DICE-TFH`=NA, `DICE-TH2`=NA, `DICE-TH1`=NA, `DICE-TH17`=NA, `DICE-THSTAR`=NA, `DICE-CD4N`=NA
        # TRM=0.3, TREG=0.25
    ),
    `c40.lg.cd8.d40`=list(
        `Guo-TEM`=NA, `Clarke-TRM`=0.25, `Pauken-TCM`=0.4,
        `Jonsson-GZMK`=NA, `Lan-GZMK`=0.6,
        `Scott-TRM-Up`=NA, `Scott-TRM-Down`=NA,
        `Cullen-Unhelped`=NA, `Consens-Exhaust`=NA, `Best-Cylcing`=NA,
        `MSIGDB-IFNR-1-2`=NA, `MSIGDB-Glycolysis`=NA
        # TRM=0.3
    ),
    `c40.ln.cd4.d08`=list(
        `Clarke-TRM`=0.2, `Pauken-TCM`=NA, `Scott-TRM-Up`=NA, `Scott-TRM-Down`=NA,
        `Patil-CTL`=NA, `Arlehamn-TH1`=NA, `Arlehamn-TH17`=NA, `Locci-TFH`=NA, `Best-Cylcing`=NA,
        `MSIGDB-IFNR-1-2`=0.4,
        `DICE-TREGN`=0.2, `DICE-TREGM`=0.2, `DICE-TFH`=NA, `DICE-TH2`=NA, `DICE-TH1`=NA, `DICE-TH17`=NA, `DICE-THSTAR`=NA, `DICE-CD4N`=NA
    ),
    `c40.ln.cd8.d08`=list(
        `Guo-TEM`=NA, `Clarke-TRM`=0.2, `Pauken-TCM`=0.4,
        `Jonsson-GZMK`=NA, `Lan-GZMK`=0.8,
        `Scott-TRM-Up`=NA, `Scott-TRM-Down`=NA,
        `Cullen-Unhelped`=NA, `Consens-Exhaust`=NA, `Best-Cylcing`=NA,
        `MSIGDB-IFNR-1-2`=0.4, `MSIGDB-Glycolysis`=NA
        # TRM=0.3
    )
)
# Process per aggregated dataset.
for(data.set in names(obj.extended.names)){
    cat(paste0(obj.extended.names[data.set], '\n'))
    this.reports.path <- paste0(tmp.reports.path, '/', obj.extended.names[data.set])
    if(!dir.exists(this.reports.path)) dir.create(this.reports.path)
    tmp.tholds <- if(data.set %in% names(gen.scale.tholds)) gen.scale.tholds[[data.set]] else NULL
    tmp.vals <- c(
        'c40.lg.cd8.d40',
        'c40.lg.cd4.d40'
    )
    dot.size <- if(t.subset %in% tmp.vals) 0.05 else 0.8
    get.mod.score.plots(
        seurat.obj=srt.objs.list[[data.set]], 
        these.mods=mods.of.int[[gen.cell.types[data.set]]],
        scale.tholds=tmp.tholds, dot.size=dot.size,
        output.path=this.reports.path, file.pffx=obj.extended.names[data.set]
    )
}

# ---> FGSEA
tmp.reports.path <- paste0(fig.phen.comps.path, '/fgsea')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
# @ Seeting definition
# Reference group
# Cluster definition pattern
# Signature
gsea.settings <- list(
    `c40.lg.cd4.d40`=list(
        `TRM`=c('TRM', '^TRM$', 'trm.signature.clarke'),
        `TCM`=c('TCM', '^TCM$', 'TCM.Pauken.Tumor'),
        `CTL`=c('CTL', '^CD4-CTL$', 'cell.cytotoxicity.patil'),
        `TREG`=c('TREG', '^TREG$', 'R24.TREGMEM.RESTING')
    ),
    `c40.lg.cd8.d40`=list(
        # `TRM`=c('TRM', '^TRM$', 'trm.signature.clarke'),
        # `TCM`=c('TCM', '^TCM$', 'TCM.Pauken.Tumor'),
        # `TEM`=c('TEM', '^TEM$', 'tcell.cytotoxicity.guo'),
        `GZMK`=c('GZMK', '^GZMK', 'gzmkpos-cd8.jonsson')
    ),
    `c40.lg.cd8.trm`=list(
        `TRM`=c('TRM', '^BATFhi$', 'trm.signature.clarke')
    ),
    `c40.ln.cd4.d08`=list(
        `TRM`=c('TRM', '^TRM$', 'trm.signature.clarke'),
        `TCM`=c('TCM', '^TCM$', 'TCM.Pauken.Tumor'),
        `TREG`=c('TREG', '^TREG$', 'R24.TREGMEM.RESTING'),
        `IFNR`=c('IFNR', '^IFNR$', 'type1n2.interferon.broadreact')
    ),
    `c40.ln.cd8.d08`=list(
        # `TRM`=c('TRM', '^TRM$', 'trm.signature.clarke'),
        # `TCM`=c('TCM', '^TCM$', 'TCM.Pauken.Tumor'),
        `GZMK`=c('GZMK', '^GZMK', 'gzmkpos-cd8.jonsson')#,
        # `IFNR`=c('IFNR', '^IFNR$', 'type1n2.interferon.broadreact')
    )
)
for(data.set in names(gsea.settings)[3]){
    cat(paste0(obj.extended.names[data.set], '\n'))
    this.reports.path <- paste0(tmp.reports.path, '/', obj.extended.names[data.set])
    if(!dir.exists(this.reports.path)) dir.create(this.reports.path)
    for(tmp.eval in names(gsea.settings[[data.set]])){
        tmp.ref <- gsea.settings[[data.set]][[tmp.eval]][1]
        tmp.pttn <- gsea.settings[[data.set]][[tmp.eval]][2]
        tmp.mod <- gsea.settings[[data.set]][[tmp.eval]][3]
        srt.objs.list[[data.set]]@meta.data[, 'tmp.tag'] <- ifelse(
            test=str_detect(
                string=clusts.defs[[data.set]][
                    as.character(srt.objs.list[[data.set]]@meta.data[, clust.labs[data.set]])
                ],
                pattern=tmp.pttn
            ),
            yes=tmp.ref,
            no='Rest'
        )
        do.fgsea(
            metrics.src=srt.objs.list[[data.set]],
            modules.feats=modules.feats[tmp.mod], 
            metric='signal.to.noise', cell.subset=subset.list[[data.set]],
            tag.of.int='tmp.tag',
            output.path=this.reports.path,
            vals.to.depict=tmp.ref,
            files.preffix=paste0(obj.extended.names[data.set])
        )
        srt.objs.list[[data.set]]@meta.data[, 'tmp.tag'] <- NULL
    }
}

# ---> Fraction of cells per cluster per donor.
for(t.subset in names(obj.extended.names)){
    tmp.blank.comp <- if(t.subset=='c40.lg.cd8.d40' | t.subset=='c40.ln.cd8.d08') blank.complement.3.2 else blank.complement.3.1
    tmp.ggplot <- get.cell.fracs.gg(
        seurat.obj=srt.objs.list[[t.subset]],
        x.var='donor.id.tag', fill.var=clust.labs[t.subset], fill.cols=clusts.cols[[t.subset]],
        break.no=3, bar.width=0.6, line.width=1.3
    )
    tmp.width <- length(unique(srt.objs.list[[t.subset]]@meta.data[, 'donor.id.tag']))
    tmp.width <- ((tmp.width * 28)/40) + 2
    tmp.lab <- paste0(obj.extended.names[t.subset], '_CellFracs_DonorID_Clusters')
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.phen.comps.path, file.name=tmp.lab, type='pdf',
        blank.comp=tmp.blank.comp, do.legend=FALSE,
        width=tmp.width, height=10
    )
}

# ---> Phenotype of all T cells and T cells clonally overlapped between anatomic sites.
for(data.set in names(aggr.set.labs)){
    tmp.vals <- meta.data.l[[data.set]][, unique(data.set)]
    for(t.subset in tmp.vals){
        anti.set <- if(data.set=='c40.lg') 'c40.ln' else 'c40.lg'
        anti.subset <- aggr.set.main[[anti.set]]
        anti.subset <- names(gen.cell.types[anti.subset])[gen.cell.types[anti.subset]==gen.cell.types[t.subset]]
        # Collect details on anti dataset, including bias info.
        tmp.data <- get.multi.site.info(
            data.set=data.set, t.subset=t.subset,
            anti.set=anti.set, anti.subset=anti.subset,
            rel.cols=NULL
        )
        to.plot <- list(
            `Full dataset`=tmp.data,
            `Overlapped only`=tmp.data[site.ovlp=='Overlapped']
        )
        to.plot <- rbindlist(l=to.plot, use.names=TRUE, idcol='data.subset')
        # Plot.
        tmp.ggplot <- ggplot(data=to.plot, aes(x=data.subset, fill=cluster)) +
            # geom_bar(stat='identity', width=0.8, linewidth=1.2, alpha=1, color='black') +
            geom_bar(position='fill', width=0.8, linewidth=0, alpha=1, color='black') +
            scale_fill_manual(values=clusts.cols[[t.subset]]) +
            scale_y_continuous(expand=c(0, 0), limits=c(0, 1), breaks=scales::pretty_breaks(n=3)) +
            labs(x='Data subset', y='Cell fraction', fill='Cluster')
        tmp.lab <- paste0('/', obj.extended.names[t.subset], '_ClusterCellFract_DataSubset')
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.phen.comps.path, file.name=tmp.lab, type='pdf',
            blank.comp=blank.complement.7, do.legend=FALSE, do.rotate=TRUE, width=3, height=5
        )
    }
}


############    -----------------------------------------    ############
### --------------------------- Figure on --------------------------- ###
### ------------------- Antigen-specific T cells -------------------- ###
### ------------- Phenotype comparisons between lineages ------------ ###
############    -----------------------------------------    ############

# ----> Output directory.
fig.as.phen.comps.path <- paste0(reports.path, '/figure_on_ag-spc_phen_comps')
if(!dir.exists(fig.as.phen.comps.path)) dir.create(fig.as.phen.comps.path)
# ---> Comments: Main dataset for this section, CD4 T cell set.
main.obj.name <- NULL

### -------------------------- Main Figure -------------------------- ###

# ---> Predicted TCRs on UMAP
for(t.subset in names(main.extended.names)){
    # @ Fetch data.
    tmp.data.1 <- meta.data.l[['c40.lg']][, .(barcode, clusters.tag, ag.spc.tag)]
    tmp.data.2 <- as.data.frame(srt.objs.list[[t.subset]]@reductions$umap@cell.embeddings)
    tmp.data.2$barcode <- row.names(tmp.data.2)
    tmp.data <- merge(
        x=tmp.data.1, y=tmp.data.2,
        by='barcode'
    )
    # @ Plot.
    # Define scale values.
    alpha.scale <- c(
        `no-spc`=0,
        `spc`=1,
        `den`=0.3
    )
    # Process per specificity
    tmp.spcs <- names(pp.cols)[names(pp.cols) %in% tmp.data[, ag.spc.tag]]
    tmp.spcs <- c('All', tmp.spcs)
    for(tmp.spc in tmp.spcs){
        # Set color variable
        tmp.data[, col.val:=clusters.tag]
        # Set scaling variable
        tmp.data[, scale.val:=ifelse(
            test=!is.na(ag.spc.tag) & ag.spc.tag==tmp.spc,
            yes='spc',
            no='no-spc'
        )]
        # Set row order, allowing for dots w/ predictions to come on top.
        setorderv(x=tmp.data, col='scale.val', order=-1)
        # Density data.
        plot.data <- if(tmp.spc!='All') tmp.data[!is.na(ag.spc.tag) & ag.spc.tag==tmp.spc] else tmp.data
        # Plot
        tmp.ggplot <- ggplot(data=tmp.data, aes(x=UMAP_1, y=UMAP_2, color=col.val)) +
            geom_point(aes(alpha=scale.val), size=0.3) +
            stat_density_2d(
                data=plot.data,
                aes(fill=after_stat(level), alpha='den'),
                geom='polygon', linewidth=1,
                contour=TRUE,
                color='black', bins=5
            ) +
            scale_color_manual(values=clusts.cols[[t.subset]]) +
            scale_fill_gradient(low="#FFFFD4", high="#F6407F") +
            scale_alpha_manual(values=alpha.scale) +
            labs(x='UMAP 1', y='UMAP 2', col='Cluster +\nReactivity', size='Pred?', alpha='Pred?')
        tmp.lab <- paste0(obj.extended.names[t.subset], '_ImpClonesOnUMAP_Spc-', tmp.spc)
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=fig.as.phen.comps.path, file.name=tmp.lab, type='tiff',
            blank.comp=blank.complement.1, do.legend=FALSE,
            height=10, width=10
        )
    }
}

# ---> Validated TCRs on UMAP
for(t.subset in names(main.extended.names)){
    # @ Fetch data.
    tmp.data.1 <- meta.data.l[['c40.lg']][, .(barcode, clusters.tag)]
    tmp.data.2 <- as.data.frame(srt.objs.list[[t.subset]]@reductions$umap@cell.embeddings)
    tmp.data.2$barcode <- row.names(tmp.data.2)
    tmp.data.1 <- merge(
        x=tmp.data.1, y=tmp.data.2,
        by='barcode'
    )
    tmp.data.2 <- meta.data.l[['c40.lg']][
        gex.lin.tag==gen.cell.types[t.subset],
        .(barcode, val.data.set, val.clonotype.tag, val.specificity.tag)
    ]
    tmp.data <- merge(
        x=tmp.data.1, y=tmp.data.2,
        by='barcode',
        all.x=TRUE, all.y=FALSE
    )
    # Set color variable
    tmp.data[, col.val:=val.specificity.tag]
    tmp.data[is.na(col.val), col.val:=clusters.tag]
    # Set scaling variable
    tmp.data[, scale.val:=ifelse(
        test=is.na(val.specificity.tag),
        yes='no',
        no='yes'
    )]
    # Set row order, allowing for dots w/ predictions to come on top.
    setorderv(x=tmp.data, col='scale.val', order=-1)
    # @ Plot.
    # Define scale values.
    alpha.scale <- c(
        `no`=0.05,
        `yes`=1
    )
    size.scale <- c(
        `no`=0.1,
        `yes`=5
    )
    # Colors
    tmp.cols <- c(
        pp.cols[names(pp.cols) %in% tmp.data[, val.specificity.tag]],
        clusts.cols[[t.subset]]
    )
    # Plot
    tmp.ggplot <- ggplot(data=tmp.data, aes(x=UMAP_1, y=UMAP_2, color=col.val, size=scale.val, alpha=scale.val)) +
        geom_point() +
        scale_color_manual(values=tmp.cols) +
        scale_size_manual(values=size.scale) +
        scale_alpha_manual(values=alpha.scale) +
        labs(x='UMAP 1', y='UMAP 2', col='Cluster +\nReactivity', size='Pred?', alpha='Pred?')
    tmp.lab <- paste0(obj.extended.names[t.subset], '_ValClonesOnUMAP')
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.as.phen.comps.path, file.name=tmp.lab, type='tiff',
        blank.comp=blank.complement.1, do.legend=FALSE,
        height=10, width=10
    )
}

# ---> Validated TCR summary
size.tholds <- c(
    `c40.lg.cd4.d40`=1,
    `c40.lg.cd8.d40`=50
)
for(t.subset in names(main.extended.names)){
    tmp.data <- meta.data.l[['c40.lg']][
        gex.lin.tag==gen.cell.types[t.subset] &
        !is.na(val.clonotype.tag),
        .(clonotype.tag=paste(val.clonotype.tag, donor.id.tag, sep='.'), clusters.tag)
    ]
    tmp.lvls <- tmp.data[, .N, by=.(clonotype.tag)]
    to.keep <- tmp.lvls[N>size.tholds[t.subset], clonotype.tag]
    tmp.lvls <- tmp.lvls[clonotype.tag %in% to.keep]
    setorderv(x=tmp.lvls, cols='N', order=1)
    tmp.lvls <- tmp.lvls$clonotype.tag
    tmp.data <- tmp.data[clonotype.tag %in% to.keep]
    tmp.vals <- factor(x=tmp.data$clonotype.tag, levels=tmp.lvls)
    set(x=tmp.data, j='clonotype.tag', value=tmp.vals)
    tmp.ggplot <- ggplot(data=tmp.data, aes(x=clonotype.tag, fill=clusters.tag)) +
        geom_bar(position='stack', width=0.8, linewidth=1, color='black') +
        scale_fill_manual(values=clusts.cols[[t.subset]]) +
        scale_y_continuous(expand=expansion(add=c(0, 2)), breaks=scales::pretty_breaks(n=3)) +
        labs(x='Clonotype ID', y='Number of cells', fill='Cluster')
    tmp.lab <- paste0(main.extended.names[t.subset], '_ValClone_CellCount_Cluster_Stacked')
    tmp.width <- (tmp.data[, uniqueN(clonotype.tag)]*1)+1.2
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.as.phen.comps.path, file.name=tmp.lab, type='pdf',
        blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=tmp.width, height=7
    )
}


### ---------- Pathogen-specific, subset-biased responses ----------- ###
### --------------------------- attempt 1 --------------------------- ###

tmp.reports.path <- paste0(fig.as.phen.comps.path, '/subset_bias_attempt-1')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

# ---> Pathogen frequency according to limit of detection.
size.tholds <- 2:3
clone.no.tholds <- 1:10
for(data.set in names(aggr.set.labs)){
    tmp.vals <- meta.data.l[[data.set]][, unique(data.set)]
    for(t.subset in tmp.vals){
        tmp.data <- meta.data.l[[data.set]][data.set==t.subset]
        tmp.data.1 <- tmp.data[
            !is.na(clonotype.tag) & !is.na(donor.id.tag) & !is.na(group.class.tag),
            .(freq.abs=.N),
            by=.(donor.id.tag, clonotype.tag, group.class.tag, ag.spc.tag)
        ]
        tmp.data.2 <- tmp.data[
            !is.na(clonotype.tag) & !is.na(donor.id.tag) & !is.na(group.class.tag),
            .(freq.total=.N),
            by=.(donor.id.tag)
        ]
        tmp.data <- merge(
            x=tmp.data.1, y=tmp.data.2,
            by=c('donor.id.tag')
        )
        tmp.data[, freq.rel:=freq.abs/freq.total]
        tmp.check <- tmp.data[, .N]==tmp.data[, .N, by=.(donor.id.tag, clonotype.tag)][, .N]
        if(!tmp.check) stop('Unexpected error!\n')
        bias.groups <- tmp.data[, sort(unique(as.character(group.class.tag)))]
        for(tmp.group in bias.groups){
            for(tmp.size.thold in size.tholds){
                to.plot <- tmp.data[
                    group.class.tag==tmp.group &
                    freq.abs>=tmp.size.thold
                ]
                to.plot <- to.plot[
                    !is.na(ag.spc.tag),
                    .(spc.freq.abs=.N),
                    by=.(donor.id.tag, ag.spc.tag)
                ]
                to.plot <- lapply(X=clone.no.tholds, FUN=function(x){
                    tmp.data <- to.plot[spc.freq.abs>=x, uniqueN(ag.spc.tag), by=donor.id.tag]
                    colnames(tmp.data)[2] <- as.character(x)
                    return(tmp.data)
                })
                to.plot <- Reduce(x=to.plot, f=function(x, y) merge(x=x, y=y, by='donor.id.tag', all=TRUE))
                to.plot <- as.data.table(gather(data=to.plot, key=clone.no.thold, value=spc.freq.abs, -`donor.id.tag`))
                # to.plot[, clone.no.thold:=as.integer(clone.no.thold)]
                tmp.lvls <- as.character(clone.no.tholds)
                to.plot[, clone.no.thold:=factor(x=clone.no.thold, levels=tmp.lvls)]
                to.plot[is.na(spc.freq.abs), spc.freq.abs:=0]
                # @ Stats. Median across donors for each specificity number threshold.
                tmp.stats <- to.plot[, .(stat=median(spc.freq.abs)), by=clone.no.thold]
                tmp.caption <- tmp.stats[, paste(clone.no.thold, stat, sep=':')]
                tmp.caption <-  paste0(tmp.caption, collapse=' ')
                # Plot.
                tmp.ggplot <- ggplot(data=to.plot, aes(x=clone.no.thold, y=spc.freq.abs)) +
                    geom_boxplot(aes(color=clone.no.thold), width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
                    geom_jitter(shape=1, stroke=5, size=10, color='black', width=0.25, height=0) +
                    # geom_line(data=tmp.stats, aes(y=stat, group='stat'), linewidth=4) +
                    scale_y_continuous(expand=expansion(add=c(0.5, 0.5)), breaks=scales::pretty_breaks(n=3)) +
                    scale_color_viridis_d(direction=-1) +
                    labs(x='Ag-specific clonotype threshold', y='No. of pathogens w/ positive response', color='', caption=tmp.caption)
                tmp.lab <- paste0('/', obj.extended.names[t.subset], '_BiasedClonotypeResp_SizeThold-', tmp.size.thold, '_Bias-', tmp.group)
                publish.plot(
                    tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
                    blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=3*to.plot[, uniqueN(clone.no.thold)], height=14
                )
            }
        }
    }
}


### -------------------- Quantitative assessment -------------------- ###
### ------------------------- First attempt ------------------------- ###

### -------------------- Quantitative assessment -------------------- ###

# ---> Quantitative enrichment for ag groups within T-cell subsets.
# @ Define relevant populations per T-cell lineage
pops.of.int <- list(
    `c40.lg.cd4.d40`=paste('c40.lg.cd4.d40', c('0', '2', '1', '3', '4'), sep='.'),
    `c40.lg.cd8.d40`=paste('c40.lg.cd8.d40', c('0', '1', '3', '2', '4'), sep='.'),
    `c40.lg.cd8.trm`=paste('c40.lg.cd8.trm', c('0', '1', '2', '3'), sep='.'),
    `c40.lg.cd8.tem`=paste('c40.lg.cd8.tem', c('0', '1', '2', '3', '4'), sep='.'),
    `c40.ln.cd4.d08`=paste('c40.ln.cd4.d08', names(clusts.cols[['c40.ln.cd4.d08']]), sep='.'),
    `c40.ln.cd8.d08`=paste('c40.ln.cd8.d08', names(clusts.cols[['c40.ln.cd8.d08']]), sep='.')
)
# @ Subset-specific relative specificities
subset.rel.spcs <- list(
    `c40.lg.cd4.d40`=rel.spcs,
    `c40.lg.cd8.d40`=c('CMV', 'EBV', 'SARS-CoV-2', 'IAV'),
    `c40.lg.cd8.trm`=c('CMV', 'EBV', 'SARS-CoV-2', 'IAV'),
    `c40.lg.cd8.tem`=c('CMV', 'EBV', 'SARS-CoV-2', 'IAV'),
    `c40.ln.cd4.d08`=rel.spcs,
    `c40.ln.cd8.d08`=c('CMV', 'EBV', 'SARS-CoV-2', 'IAV')
)
# @ Define dataset-specific limits of detection.
subset.lods <- c(
    `c40.lg.cd4.d40`=limit.of.detect,
    `c40.lg.cd8.d40`=limit.of.detect,
    `c40.lg.cd8.trm`=limit.of.detect,
    `c40.lg.cd8.tem`=limit.of.detect,
    `c40.ln.cd4.d08`=limit.of.detect,
    `c40.ln.cd8.d08`=limit.of.detect
)
subset.lods <- data.table(cell.type=names(subset.lods), lod=subset.lods)

# ---> Preflights
# General definitions.
# lg.aggr <- names(aggr.set.corr)[aggr.set.corr=='c40.lg']
# tmp.data <- lapply(X=lg.aggr, FUN=function(t.subset){
tmp.data <- lapply(X=names(obj.extended.names), FUN=function(t.subset){
    tmp.data <- as.data.table(srt.objs.list[[t.subset]]@meta.data)
    tmp.data[, barcode:=Cells(srt.objs.list[[t.subset]])]
    tmp.data <- tmp.data[, .(
        barcode, donor.id.tag,
        clusters.tag=get(clust.labs[t.subset]),
        clonotype.tag, consensus.pred=ag.spc.tag
    )]
})
# names(tmp.data) <- lg.aggr
names(tmp.data) <- names(obj.extended.names)
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='cell.type')
tmp.data[, pop.tag:=paste(cell.type, clusters.tag, sep='.')]
tmp.vals <- factor(x=as.character(tmp.data$donor.id.tag), levels=sort(unique(tmp.data$donor.id.tag)))
set(x=tmp.data, j='donor.id.tag', value=tmp.vals)
tmp.groups <- rel.spcs
names(tmp.groups) <- tmp.groups
tmp.groups <- c(
    as.list(tmp.groups),
    ag.groups
)
tmp.pops <- unlist(pops.of.int)

#   @ Previous attempt discarded.

#   @ Fractions of ag set-specific cells/clonotypes in the group (cluster) out of all ag set-specific cells/clonotypes.
# Preflights
cat.var <- 'clusters.tag'
fig.data.2 <- lapply(X=names(tmp.groups), FUN=function(ag.group){
    tmp.vals <- tmp.groups[[ag.group]]
    tmp.data.1 <- tmp.data[
        !is.na(donor.id.tag) &
        pop.tag %in% tmp.pops,
        .(
            cell.count=.SD[consensus.pred %in% tmp.vals, .N],
            clone.count=.SD[consensus.pred %in% tmp.vals, uniqueN(clonotype.tag)]
        ),
        by=.(cell.type, donor.id.tag, cat.tag=get(cat.var))
    ]
    tmp.data.2 <- tmp.data[
        !is.na(donor.id.tag) &
        pop.tag %in% tmp.pops &
        consensus.pred %in% tmp.vals,
        .(
            cell.total=.N,
            clone.total=uniqueN(clonotype.tag)
        ),
        by=.(cell.type, donor.id.tag)
    ]
    # Apply limit of detection filtering step.
    tmp.data.2 <- merge(x=tmp.data.2, y=subset.lods, by='cell.type', all.x=TRUE, all.y=FALSE)
    tmp.check <- tmp.data.2[, !any(is.na(lod))]
    if(!tmp.check) stop('Unexpected error!\n')
    tmp.data.2 <- tmp.data.2[cell.total>=lod]
    tmp.data.2[, lod:=NULL]
    tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by=c('cell.type', 'donor.id.tag'), all=FALSE)
    tmp.data[, `:=`(
        cell.fract=cell.count/cell.total,
        clone.fract=clone.count/clone.total
    )]
    # tmp.cols <- setdiff(x=colnames(tmp.data), y=c('cell.count', 'clone.count', 'cell.total', 'clone.total'))
    tmp.cols <- setdiff(x=colnames(tmp.data), y=c('cell.total', 'clone.total'))
    tmp.data <- tmp.data[, ..tmp.cols]
    tmp.check <- tmp.data[, .(check=round(x=sum(cell.fract))), by=.(cell.type, donor.id.tag)][, all(check==1)]
    if(!tmp.check) stop('Unexpected error!\n')
    return(tmp.data)
})
names(fig.data.2) <- names(tmp.groups)
fig.data.2 <- rbindlist(l=fig.data.2, use.names=TRUE, idcol='ag.group')

#   @ Put different types of fractions together.
fig.data <- list(
    # `1`=fig.data.1,
    `2`=fig.data.2#,
    # `3`=fig.data.3
)

# @ Tidy donor meta.
tmp.vals <- c('flu.shot.span', 'covid.shot.span', 'age', 'gender')
donor.meta.choice <- donor.meta.hmap[donor.var %in% tmp.vals]
donor.meta.choice[, fill.val:=str_replace(string=fill.val, pattern=paste0(donor.var, '\\.'), replacement='')]
donor.meta.choice[, donor.var:=str_to_title(string=str_replace_all(
        string=donor.var,
        pattern='\\.', replacement=' '
    ))
]
donor.meta.choice <- spread(data=donor.meta.choice[, .(donor.id.tag, donor.var, fill.val)], key=donor.var, value=fill.val)
donor.meta.choice <- merge(
    x=donor.meta.choice,
    y=ser.cmv.data[, .(donor.id.tag, `CMV Ser.`=cmv.status)],
    by='donor.id.tag'
)
donor.meta.choice <- as.data.frame(donor.meta.choice, stringsAsFactors=FALSE)
row.names(donor.meta.choice) <- donor.meta.choice$donor.id.tag; donor.meta.choice$donor.id.tag <- NULL
tmp.cols <- c(
    'CMV Ser.', 'Gender', 'Age', 'Covid Shot Span', 'Flu Shot Span'
)
donor.meta.choice <- donor.meta.choice[, tmp.cols]
donor.cols.choice <- list(
    `Gender`=c(
        'Female'='#B63EB8',
        'Male'='#78B852'
    ),
    `Age`=c(
        '<66'='#26EDFF',
        '66-70'='#23A4E8',
        '71-75'='#3383FF',
        '>75'='#233DE8'
    ),
    `Covid Shot Span`=c(
        '<1m'='#52A858',
        '>1m & <3m'='#66D16E',
        '>3m'='#79F581'
    ),
    `Flu Shot Span`=c(
        '<1m'='#52A858',
        '>1m & < 3m'='#66D16E',
        '>3m'='#79F581'
    ),
    `CMV Ser.`=c(
        `Positive`=unname(pp.cols['CMV']),
        `Negative`='#CBCBCB'
    )
)

# ---> Apply full assessment for each type of fractions.
# @ General function
do.quant.assess <- function(
    plot.data, att.lab,
    subset.lods,
    ag.group.tholds.1,
    tmp.rel.spcs=NULL,
    hmap.rel.pops=NULL,
    order.clust.opts=NULL,
    order.functs=c(`Max`=max),
    ref.groups=NULL
){
    # ---> Retrieve fractions.
    # plot.data <- fig.data[[att]]

    # ---> Donor-specific information.
    tmp.reports.path <- paste0(fig.as.phen.comps.path, '/donor-spc_info_attempt-', att.lab)
    if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
    # Set general thresholds.
    # for(t.subset in names(main.extended.names)){
    for(t.subset in plot.data[, unique(cell.type)]){
        cat(t.subset, '\n')
        # Fetch data to plot (general).
        gen.data <- plot.data[
            cell.type==t.subset
        ]
        if(!is.null(tmp.rel.spcs)) gen.data <- gen.data[ag.group %in% tmp.rel.spcs[[t.subset]]]
        if(!is.null(hmap.rel.pops) & !is.null(hmap.rel.pops[[t.subset]])){
            gen.data <- gen.data[cat.tag %in% hmap.rel.pops[[t.subset]]]
        }
        uniq.pops <- gen.data[, sort(unique(as.character(cat.tag)))]
        uniq.pps <- names(pp.cols)[names(pp.cols) %in% gen.data$ag.group]
        # Fetch data to plot (bar plots, absolute cell frequency).
        brpt.data <- gen.data[ag.group %in% uniq.pps]
        brpt.data <- brpt.data[, .(cell.count=sum(cell.count)), by=.(ag.group, donor.id.tag)]
        tmp.vals <- factor(x=brpt.data$ag.group, levels=uniq.pps)
        set(x=brpt.data, j='ag.group', value=tmp.vals)
        # Fetch data to plot (heatmap).
        hmap.data <- gen.data[, .(donor.id.tag=as.character(donor.id.tag), group=paste(cat.tag, ag.group, sep='.'), cell.fract)]
        hmap.data <- spread(data=hmap.data, key=group, value=cell.fract, drop=FALSE)
        row.names(hmap.data) <- hmap.data$donor.id.tag; hmap.data$donor.id.tag <- NULL
        hmap.data <- as.matrix(hmap.data)
        tmp.cols <- paste(rep(x=uniq.pops, each=length(uniq.pps)), uniq.pps, sep='.')
        hmap.data <- hmap.data[, tmp.cols]
        # Hierarchical clustering (if possible)
        clust.obj <- dist(x=hmap.data, method='euclidean')
        clust.obj <- tryCatch(
            expr={hclust(d=clust.obj, method='ward.D')},
            error=function(e) NULL
        )
        # Set metadata.
        col.meta <- data.frame(
            row.names=colnames(hmap.data),
            `Cluster`=str_extract(
                string=colnames(hmap.data),
                pattern='^\\d+'
            ),
            `Specificity`=str_replace(
                string=colnames(hmap.data),
                pattern='^\\d+\\.', replacement=''
            )
        )
        # Define colors for metadata tracks.
        ann.colors <- list(
            `Cluster`=clusts.cols[[t.subset]][names(clusts.cols[[t.subset]]) %in% col.meta$Cluster],
            `Specificity`=pp.cols[names(pp.cols) %in% col.meta$Specificity]
        )
        ann.colors <- c(
            ann.colors,
            donor.cols.choice
        )
        # Set color scale and breaks for heatmap.
        break.no <- 200
        col.breaks <- seq(from=min(range(hmap.data, na.rm=TRUE)), to=max(range(hmap.data, na.rm=TRUE)), length.out=break.no)
        mid.point <- (max(col.breaks) - min(col.breaks)) / 2
        mid.point <- which.min(abs(col.breaks-mid.point))
        # hmap.col.scale.1 <- colorRampPalette(c('#414487', '#2A788E'))(mid.point)
        # hmap.col.scale.2 <- colorRampPalette(c('#2A788E', '#22A884', '#FDE725', '#FFEA00'))(break.no-(mid.point+1))
        hmap.col.scale.1 <- colorRampPalette(c('#414487', '#355F8A', '#2A788E', '#22A884'))(mid.point)
        hmap.col.scale.2 <- colorRampPalette(c('#22A884', '#7AD151', '#FDE725', '#FFEA00'))(break.no-(mid.point+1))
        hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
        na.col <- '#000000'
        # File specs.
        tmp.height <- (nrow(hmap.data) * 0.18) + 0.3
        tmp.width <- ncol(hmap.data) * 0.18
        # @ Version with hierarchical clustering.
        if(!is.null(clust.obj)){
            # Complete version.
            tmp.file.name <- paste0(
                tmp.reports.path, '/',
                obj.extended.names[t.subset], '_Donor_Reactivity_CellFract_WClust.C.pdf'
            )
            pheatmap(
                mat=hmap.data, scale='none',
                color=hmap.col.scale, breaks=col.breaks, na_col=na.col, border_color='black',
                cluster_rows=clust.obj, cluster_cols=FALSE,
                # clustering_method='ward.D',
                annotation_row=donor.meta.choice, annotation_col=col.meta, annotation_colors=ann.colors,
                show_colnames=FALSE, show_rownames=TRUE, legend=TRUE, annotation_legend=TRUE,
                filename=tmp.file.name, height=tmp.height+3, width=tmp.width+3
            )
            # Blank version.
            tmp.file.name <- paste0(
                tmp.reports.path, '/',
                obj.extended.names[t.subset], '_Donor_Reactivity_CellFract_WClust.B.pdf'
            )
            pheatmap(
                mat=hmap.data, scale='none',
                color=hmap.col.scale, breaks=col.breaks, na_col=na.col, border_color='black',
                cluster_rows=clust.obj, cluster_cols=FALSE,
                # clustering_method='ward.D',
                # annotation_row=donor.meta.choice,
                annotation_col=col.meta, annotation_colors=ann.colors,
                show_colnames=FALSE, show_rownames=FALSE,
                legend=FALSE, annotation_legend=FALSE, annotation_names_row=FALSE, annotation_names_col=FALSE,
                filename=tmp.file.name, height=tmp.height+1, width=tmp.width+1
            )
            # @ Number of cells per specificity per donor.
            tmp.lvls <- rev(x=clust.obj$labels[clust.obj$order])
            tmp.vals <- factor(x=as.character(brpt.data[, donor.id.tag]), levels=tmp.lvls)
            set(x=brpt.data, j='donor.id.tag', value=tmp.vals)
            tmp.ggplot <- ggplot(data=brpt.data, aes(x=donor.id.tag, y=cell.count, fill=ag.group)) +
                geom_bar(stat='identity', position='stack', width=0.6, color='black', linewidth=1.3) +
                scale_y_continuous(expand=expansion(add=c(0, 0)), breaks=scales::pretty_breaks(n=3)) +
                scale_fill_manual(values=pp.cols) +
                labs(x='Donor ID', y='Cell count', fill='Specificity')
            tmp.lab <- paste0(obj.extended.names[t.subset], '_SpcCounts_DonorID_Clusters_HClustOrder')
            publish.plot(
                tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
                blank.comp=blank.complement.3.1, do.legend=FALSE,
                width=tmp.height*5, height=7
            )
            # @ Fraction of cells per cluster per donor.
            tmp.lvls <- rev(x=clust.obj$labels[clust.obj$order])
            tmp.vals <- factor(x=as.character(srt.objs.list[[t.subset]]@meta.data[, 'donor.id.tag']), levels=tmp.lvls)
            srt.objs.list[[t.subset]]@meta.data[, 'tmp.tag'] <- tmp.vals
            tmp.ggplot <- get.cell.fracs.gg(
                seurat.obj=srt.objs.list[[t.subset]],
                x.var='tmp.tag', fill.var=clust.labs[t.subset], fill.cols=clusts.cols[[t.subset]],
                break.no=3, bar.width=0.6, line.width=1.3
            )
            tmp.lab <- paste0(obj.extended.names[t.subset], '_PopFracs_DonorID_Clusters_HClustOrder')
            publish.plot(
                tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
                blank.comp=blank.complement.3.1, do.legend=FALSE,
                width=tmp.height*5, height=7
            )
        }
        # @ Version w/out hierarchical clustering.
        # Complete version.
        tmp.file.name <- paste0(
            tmp.reports.path, '/',
            obj.extended.names[t.subset], '_Donor_Reactivity_CellFract_NoClust.C.pdf'
        )
        pheatmap(
            mat=hmap.data, scale='none',
            color=hmap.col.scale, breaks=col.breaks, na_col=na.col, border_color='black',
            cluster_rows=FALSE, cluster_cols=FALSE,
            annotation_row=donor.meta.choice, annotation_col=col.meta, annotation_colors=ann.colors,
            show_colnames=FALSE, show_rownames=TRUE, legend=TRUE, annotation_legend=TRUE,
            filename=tmp.file.name, height=tmp.height+2, width=tmp.width+2
        )
        # Blank version.
        tmp.file.name <- paste0(
            tmp.reports.path, '/',
            obj.extended.names[t.subset], '_Donor_Reactivity_CellFract_NoClust.B.pdf'
        )
        pheatmap(
            mat=hmap.data, scale='none',
            color=hmap.col.scale, breaks=col.breaks, na_col=na.col, border_color='black',
            cluster_rows=FALSE, cluster_cols=FALSE,
            # annotation_row=donor.meta.choice,
            annotation_col=col.meta, annotation_colors=ann.colors,
            show_colnames=FALSE, show_rownames=FALSE,
            legend=FALSE, annotation_legend=FALSE, annotation_names_row=FALSE, annotation_names_col=FALSE,
            filename=tmp.file.name, height=tmp.height, width=tmp.width
        )
        # @ Version w/out hierarchical clustering w/ specific order.
        if(!is.null(order.clust.opts)){
            if(!is.null(order.clust.opts[[t.subset]])){
                for(order.clust.opt in order.clust.opts[[t.subset]]){
                    for(order.funct.name in names(order.functs)){
                        # @ Heatmap.
                        order.funct <- order.functs[[order.funct.name]]
                        tmp.cols <- colnames(hmap.data)[str_detect(string=colnames(hmap.data), pattern=order.clust.opt)]
                        to.order <- hmap.data[, tmp.cols]
                        to.order <- apply(X=to.order, MARGIN=1, FUN=order.funct, na.rm=TRUE)
                        to.order <- sort(x=to.order, decreasing=TRUE)
                        to.order <- names(to.order)
                        hmap.data <- hmap.data[to.order, ]
                        # Complete version.
                        tmp.file.name <- paste0(
                            tmp.reports.path, '/',
                            obj.extended.names[t.subset], '_Donor_Reactivity_CellFract_CustOdr-', order.funct.name, '-', order.clust.opt, '.C.pdf'
                        )
                        pheatmap(
                            mat=hmap.data, scale='none',
                            color=hmap.col.scale, breaks=col.breaks, na_col=na.col, border_color='black',
                            cluster_rows=FALSE, cluster_cols=FALSE,
                            annotation_row=donor.meta.choice, annotation_col=col.meta, annotation_colors=ann.colors,
                            show_colnames=FALSE, show_rownames=TRUE, legend=TRUE, annotation_legend=TRUE,
                            filename=tmp.file.name, height=tmp.height+2, width=tmp.width+2
                        )
                        # Blank version.
                        tmp.file.name <- paste0(
                            tmp.reports.path, '/',
                            obj.extended.names[t.subset], '_Donor_Reactivity_CellFract_CustOdr-', order.funct.name, '-', order.clust.opt, '.B.pdf'
                        )
                        pheatmap(
                            mat=hmap.data, scale='none',
                            color=hmap.col.scale, breaks=col.breaks, na_col=na.col, border_color='black',
                            cluster_rows=FALSE, cluster_cols=FALSE,
                            # annotation_row=donor.meta.choice,
                            annotation_col=col.meta, annotation_colors=ann.colors,
                            show_colnames=FALSE, show_rownames=FALSE,
                            legend=FALSE, annotation_legend=FALSE, annotation_names_row=FALSE, annotation_names_col=FALSE,
                            filename=tmp.file.name, height=tmp.height, width=tmp.width
                        )
                        # @ Number of cells per specificity per donor.
                        tmp.lvls <- rev(x=to.order)
                        tmp.vals <- factor(x=as.character(brpt.data[, donor.id.tag]), levels=tmp.lvls)
                        set(x=brpt.data, j='donor.id.tag', value=tmp.vals)
                        tmp.ggplot <- ggplot(data=brpt.data, aes(x=donor.id.tag, y=cell.count, fill=ag.group)) +
                            geom_bar(stat='identity', position='stack', width=0.6, color='black', linewidth=1.3) +
                            scale_y_continuous(expand=expansion(add=c(0, 0)), breaks=scales::pretty_breaks(n=3)) +
                            scale_fill_manual(values=pp.cols) +
                            labs(x='Donor ID', y='Cell count', fill='Specificity')
                        tmp.lab <- paste0(obj.extended.names[t.subset], '_SpcCounts_DonorID_Clusters_CustOdr-', order.funct.name, '-', order.clust.opt)
                        publish.plot(
                            tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
                            blank.comp=blank.complement.3.1, do.legend=FALSE,
                            width=tmp.height*5, height=7
                        )
                        # @ Fraction of cells per cluster per donor.
                        tmp.lvls <- rev(x=to.order)
                        tmp.vals <- factor(x=as.character(srt.objs.list[[t.subset]]@meta.data[, 'donor.id.tag']), levels=tmp.lvls)
                        if('tmp.tag' %in% colnames(srt.objs.list[[t.subset]]@meta.data)) srt.objs.list[[t.subset]]@meta.data[, 'tmp.tag'] <- NULL
                        srt.objs.list[[t.subset]]@meta.data[, 'tmp.tag'] <- tmp.vals
                        tmp.check <- (!is.null(hmap.rel.pops) & !is.null(hmap.rel.pops[[t.subset]]) & length(hmap.rel.pops[[t.subset]])>1) | is.null(hmap.rel.pops[[t.subset]])
                        if(tmp.check){
                            tmp.ggplot <- get.cell.fracs.gg(
                                seurat.obj=srt.objs.list[[t.subset]],
                                x.var='tmp.tag', fill.var=clust.labs[t.subset], fill.cols=clusts.cols[[t.subset]],
                                break.no=3, bar.width=0.6, line.width=1.3
                            )
                        }else{
                            tmp.data <- as.data.table(srt.objs.list[[t.subset]]@meta.data)
                            tmp.data <- tmp.data[
                                !is.na(tmp.tag),
                                .(freq.rel=.SD[
                                    get(clust.labs[t.subset])==hmap.rel.pops[[t.subset]], .N
                                    ] / .N
                                ),
                                by=.(donor.id.tag=tmp.tag)
                            ]
                            tmp.cols <- clusts.cols[[t.subset]]
                            tmp.cols <- tmp.cols[names(tmp.cols)==hmap.rel.pops[[t.subset]]]
                            tmp.ggplot <- tmp.ggplot <- ggplot(data=tmp.data, aes(x=donor.id.tag, y=freq.rel)) +
                                geom_bar(stat='identity', color='black', fill=tmp.cols, width=0.6, linewidth=1.3) +
                                scale_y_continuous(expand=expansion(add=c(0, 0)), limits=c(0, 1), breaks=scales::pretty_breaks(n=3)) +
                                labs(x='', y='Cell fraction')
                        }
                        tmp.lab <- paste0(obj.extended.names[t.subset], '_PopFracs_DonorID_Clusters_CustOdr-', order.funct.name, '-', order.clust.opt)
                        publish.plot(
                            tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
                            blank.comp=blank.complement.3.1, do.legend=FALSE,
                            width=tmp.height*5, height=7
                        )
                    }
                }
            }
        }
        # @ Color legend only.
        # W/ labels
        # tmp.breaks <- if(tmp.group %in% names(ag.group.breaks)) NULL else ag.group.breaks[[tmp.group]]
        tmp.breaks <- NULL
        tmp.data <- as.data.table(gather(data=as.data.frame(hmap.data), key=key, value=value))
        tmp.data <- tmp.data[!is.na(value)]
        if(is.null(tmp.breaks)){
            tmp.ggplot <- ggplot(data=tmp.data, aes(x=value, y=value, col=value)) + 
                geom_point() +
                scale_color_gradientn(colors=hmap.col.scale) +
                theme(
                    legend.ticks=element_line(color='black', linewidth=0.6),
                    legend.ticks.length=unit(0.22, "cm"),
                    legend.frame=element_rect(color='black', linewidth=0.6)
                )
        }else{
            tmp.ggplot <- ggplot(data=tmp.data, aes(x=value, y=value, col=value)) + 
                geom_point() +
                scale_color_gradientn(colors=hmap.col.scale) +
                scale_color_gradientn(colors=hmap.col.scale, breaks=tmp.breaks) +
                theme(
                    legend.ticks=element_line(color='black', linewidth=0.6),
                    legend.ticks.length=unit(0.22, "cm"),
                    legend.frame=element_rect(color='black', linewidth=0.6)
                )
        }
        tmp.ggplot <- get_legend(p=tmp.ggplot)
        tmp.file.name <- paste0(
            tmp.reports.path, '/',
            obj.extended.names[t.subset], '_Donor_Reactivity_CellFract_Legend.L1.pdf'
        )
        pdf(file=tmp.file.name, height=2, width=2)
        print(as_ggplot(tmp.ggplot))
        dev.off()
        # W/out labels
        if(is.null(tmp.breaks)){
            tmp.ggplot <- ggplot(data=tmp.data, aes(x=value, y=value, col=value)) + 
                geom_point() +
                scale_color_gradientn(colors=hmap.col.scale, name=NULL, labels=NULL) +
                theme(
                    legend.ticks=element_line(color='black', linewidth=0.6),
                    legend.ticks.length=unit(0.22, "cm"),
                    legend.frame=element_rect(color='black', linewidth=0.6)
                )
        }else{
            tmp.ggplot <- ggplot(data=tmp.data, aes(x=value, y=value, col=value)) + 
                geom_point() +
                scale_color_gradientn(colors=hmap.col.scale) +
                scale_color_gradientn(colors=hmap.col.scale, breaks=tmp.breaks, name=NULL, labels=NULL) +
                theme(
                    legend.ticks=element_line(color='black', linewidth=0.6),
                    legend.ticks.length=unit(0.22, "cm"),
                    legend.frame=element_rect(color='black', linewidth=0.6)
                )
        }
        tmp.ggplot <- get_legend(p=tmp.ggplot)
        tmp.file.name <- paste0(
            tmp.reports.path, '/',
            obj.extended.names[t.subset], '_Donor_Reactivity_CellFract_Legend.L2.pdf'
        )
        pdf(file=tmp.file.name, height=1.25, width=0.3)
        print(as_ggplot(tmp.ggplot))
        dev.off()
    }

    # ---> Summarization output dir.
    tmp.reports.path <- paste0(fig.as.phen.comps.path, '/summ_info_attempt-', att.lab)
    if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

    # ---> Summarization across donors per population.
    # Process per aggregated dataset.
    for(t.subset in names(obj.extended.names)){
        # General cluster-wise fraction of all T cells.
        tmp.data <- as.data.table(srt.objs.list[[t.subset]]@meta.data)
        tmp.data.1 <- tmp.data[
            !is.na(donor.id.tag),
            .(freq.abs=.N),
            by=.(donor.id.tag, cat.tag=get(clust.labs[t.subset]))
        ]
        tmp.data.2 <- tmp.data[
            !is.na(donor.id.tag),
            .(freq.total=.N),
            by=.(donor.id.tag)
        ]
        gen.plot.data <- merge(x=tmp.data.1, y=tmp.data.2, by='donor.id.tag')
        gen.plot.data[, cell.fract:=freq.abs/freq.total]
        # Apply consistent absolute frequency threshold.
        gen.plot.data <- gen.plot.data[freq.total>subset.lods[cell.type==t.subset, lod]]
        gen.plot.data[, `:=`(freq.abs=NULL, freq.total=NULL)]
        gen.plot.data[, ag.group:='ALL T']
        # Process per cluster.
        tmp.pops <- as.character(plot.data[cell.type==t.subset, unique(cat.tag)])
        for(tmp.pop in tmp.pops){
            # Fetch data to plot.
            tmp.data <- plot.data[cell.type==t.subset & cat.tag==tmp.pop, ]
            tmp.data <- tmp.data[ag.group %in% names(pp.cols)]
            if(!is.null(tmp.rel.spcs)) tmp.data <- tmp.data[ag.group %in% tmp.rel.spcs[[t.subset]]]
            # ---> Stats.
            # Kruskal Wallis test.
            tmp.test.1 <- kruskal.test(
                x=tmp.data,
                formula=cell.fract~ag.group
            )
            tmp.test.1 <- tmp.test.1$p.value
            # Multiple pairwise Wilcoxon signed-rank test.
            tmp.test.2 <- pairwise.wilcox.test(
                x=tmp.data$cell.fract, g=tmp.data$ag.group,
                p.adjust.method='BH'
            )
            tmp.test.2 <- tmp.test.2$p.value
            tmp.lab <- paste0(obj.extended.names[t.subset], '_Pop-', tmp.pop, '_CellFract_Pat_Donor_MultPairwiseComps.csv')
            tmp.file.name <- paste0(tmp.reports.path, '/', tmp.lab)
            write.csv(file=tmp.file.name, x=tmp.test.2, quote=TRUE, row.names=TRUE)
            # Continue w/ preparation.
            tmp.thold <- if(tmp.pop %in% names(ag.group.tholds.1[[t.subset]])) ag.group.tholds.1[[t.subset]][as.character(tmp.pop)] else 1
            # Add cell fractions for all T cells.
            tmp.data <- list(
                tmp.data,
                gen.plot.data
            )
            tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, fill=TRUE)
            # Set factors
            tmp.lvls <- names(pp.cols)[names(pp.cols) %in% tmp.data[, ag.group]]
            tmp.lvls <- c('ALL T', tmp.lvls)
            tmp.vals <- factor(x=tmp.data[, ag.group], levels=tmp.lvls)
            set(x=tmp.data, j='ag.group', value=tmp.vals)
            # Set upper bound.
            tmp.data[, bounded.cell.fract:=cell.fract]
            tmp.data[, point.status:='ori']
            tmp.data[bounded.cell.fract>tmp.thold, `:=`(bounded.cell.fract=tmp.thold, point.status='bounded')]
            # Set bound.
            if(tmp.thold==1) tmp.thold <- NA
            # Dimensions and theme
            tmp.height <- 14
            tmp.width <- ((tmp.data[, uniqueN(ag.group)]-1)*2)
            if(!is.null(ref.groups) & tmp.pop %in% ref.groups[[t.subset]]){
                tmp.blank.comp <- blank.complement.3.1
                tmp.width <- tmp.width + 0.4
            }else{
                tmp.blank.comp <- blank.complement.3.2
            }
            # ---> Version w/ ag-specific cell fractions only.
            tmp.ggplot <- ggplot(data=tmp.data[ag.group!='ALL T'], aes(x=ag.group)) +
                geom_jitter(aes(y=bounded.cell.fract, shape=point.status), stroke=5, size=8, color=gen.dot.col, width=0, height=0) +
                geom_boxplot(aes(y=cell.fract, color=ag.group, fill=ag.group), alpha=0.4, width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
                # coord_cartesian(ylim=c(0, tmp.thold)) +
                scale_y_continuous(
                    limits=c(0, 1),
                    breaks=scales::pretty_breaks(n=3)
                ) +
                scale_shape_manual(values=tmp.shapes) +
                scale_color_manual(values=pp.cols) +
                scale_fill_manual(values=pp.cols) +
                labs(
                    x='Group',
                    y=paste0('Percentage ag-specific cells in indicated groups'),
                    color='',
                    caption=paste0(
                        'Each dot represents an individual donor.\n',
                        'Kruskal-Wallis Rank Sum Test P-value: ',
                        tmp.test.1
                    )
                )
            tmp.lab <- paste0(obj.extended.names[t.subset], '_Pop-', tmp.pop, '_CellFract_Pat_Donor')
            publish.plot(
                tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
                height=tmp.height, width=tmp.width,
                blank.comp=tmp.blank.comp, do.legend=FALSE, do.rotate=TRUE
            )
            # ---> Version including ALL T cell group.
            tmp.cols <- c(
                `ALL T`='#000000',
                pp.cols
            )
            tmp.width <- tmp.width + 2
            tmp.ggplot <- ggplot(data=tmp.data, aes(x=ag.group)) +
                geom_jitter(aes(y=bounded.cell.fract, shape=point.status), stroke=5, size=8, color=gen.dot.col, width=0, height=0) +
                geom_boxplot(aes(y=cell.fract, color=ag.group, fill=ag.group), alpha=0.4, width=0.7, linewidth=4, fatten=4, outlier.shape=NA) +
                # coord_cartesian(ylim=c(0, tmp.thold)) +
                scale_y_continuous(
                    limits=c(0, 1),
                    breaks=scales::pretty_breaks(n=3)
                ) +
                scale_shape_manual(values=tmp.shapes) +
                scale_color_manual(values=tmp.cols) +
                scale_fill_manual(values=tmp.cols) +
                labs(
                    x='Group',
                    y=paste0('Percentage ag-specific cells in indicated groups'),
                    color='',
                    caption=paste0(
                        'Each dot represents an individual donor.\n',
                        'Kruskal-Wallis Rank Sum Test P-value: ',
                        tmp.test.1
                    )
                )
            tmp.lab <- paste0(obj.extended.names[t.subset], '_Pop-', tmp.pop, '_CellFract_Pat-and-ALL_Donor')
            publish.plot(
                tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
                height=tmp.height, width=tmp.width,
                blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE
            )
        }
    }
    return(NULL)
}

# @ Previous attempts discarded.

# @ Process for fractions of ag set-specific cells/clonotypes in the group (cluster) out of all cells in the group (cluster).
ag.group.tholds.1 <- list(
    `c40.lg.cd4.d40`=c(),
    `c40.lg.cd8.d40`=c(),
    `c40.lg.cd8.trm`=c(),
    `c40.lg.cd8.tem`=c()
)
hmap.rel.pops <- list(
    `c40.lg.cd4.d40`=c('0'),
    `c40.lg.cd8.d40`=c('0', '1', '2'),
    `c40.lg.cd8.trm`=NULL,
    `c40.lg.cd8.tem`=NULL,
    `c40.ln.cd4.d08`=c('6'),
    `c40.ln.cd8.d08`=c('0', '1', '3')
)
order.clust.opts <- list(
    `c40.lg.cd4.d40`=c('0'),
    `c40.lg.cd8.d40`='0',
    `c40.ln.cd4.d08`='6',
    `c40.ln.cd8.d08`='3'
)
ref.groups <- list(
    `c40.lg.cd4.d40`=c('0', '2'),
    `c40.lg.cd8.d40`='0',
    `c40.lg.cd8.trm`='0',
    `c40.ln.cd4.d08`='6',
    `c40.ln.cd8.d08`='3'
)
do.quant.assess(
    plot.data=fig.data[['2']], att.lab='2',
    subset.lods=subset.lods,
    ag.group.tholds.1=ag.group.tholds.1,
    tmp.rel.spcs=subset.rel.spcs,
    hmap.rel.pops=hmap.rel.pops,
    order.clust.opts=order.clust.opts,
    order.functs=c(`Max`=max, `Mean`=mean, `Median`=median),
    ref.groups=ref.groups
)

# ---> Pathogen-specific fractions across donors.
bg.cell.set <- rbindlist(l=list(pred.bc.ref.1, pred.bc.ref.2), use.names=TRUE)
bg.cell.set <- bg.cell.set[ext.consensus.pred=='Background', barcode]
for(t.subset in names(obj.extended.names)){
    # Retrieve specificities.
    data.set <- aggr.set.corr[t.subset]
    tmp.data <- as.data.frame(meta.data.l[[data.set]][, .(barcode, tmp.tag=as.character(ag.spc.tag))])
    # Keep only relevant specificities.
    tmp.data[!tmp.data$tmp.tag %in% subset.rel.spcs[[t.subset]], 'tmp.tag'] <- NA
    # Add background set.
    # tmp.data[tmp.data$barcode %in% bg.cell.set, 'tmp.tag'] <- 'Background'
    # Add info to seurat object and plot.
    row.names(tmp.data) <- tmp.data$barcode; tmp.data$barcode <- NULL
    tmp.check <- all(Cells(srt.objs.list[[t.subset]]) %in% row.names(tmp.data))
    if(!tmp.check) stop('Unexpected error.\n')
    tmp.data$tmp.tag <- factor(x=tmp.data$tmp.tag, levels=c(names(pp.cols), 'Background'))
    srt.objs.list[[t.subset]]@meta.data[, 'tmp.tag'] <- tmp.data[Cells(srt.objs.list[[t.subset]]), 'tmp.tag']
    tmp.ggplot <- get.cell.fracs.gg(
        seurat.obj=srt.objs.list[[t.subset]],
        x.var='tmp.tag', fill.var=clust.labs[t.subset], fill.cols=clusts.cols[[t.subset]],
        # x.var='ag.spc.tag', fill.var='clusters.tag', fill.cols=clusts.cols[[t.subset]]
        break.no=3
    )
    tmp.width <- uniqueN(tmp.data$tmp.tag) * 1.5
    tmp.lab <- paste0(obj.extended.names[t.subset], '_CellFracs_Spcs_Clusters')
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=fig.as.phen.comps.path, file.name=tmp.lab, type='pdf',
        blank.comp=blank.complement.3, do.legend=FALSE,
        width=tmp.width, height=7
    )
    srt.objs.list[[t.subset]]@meta.data[, 'tmp.tag'] <- NULL
}

# ----> Associations, process oer variable.
# @ General function
do.assoc.assess <- function(
    plot.data, att.lab,
    tmp.rel.spcs=NULL
){
    # ---> Population-specific frequencies
    tmp.data <- lapply(X=aggr.set.main[['c40.lg']], FUN=function(tmp.lin){
        tmp.data.1 <- meta.data.l[['c40.lg']][
            data.set==tmp.lin,
            .(freq.abs=.N),
            by=.(
                donor.id.tag,
                pop.tag=paste0(
                    gen.cell.types[tmp.lin], '.',
                    clusts.defs[[tmp.lin]][get(clust.labs[tmp.lin])]
                )
            )
        ]
        tmp.data.2 <- meta.data.l[['c40.lg']][
            data.set==tmp.lin,
            .(freq.total=.N),
            by=.(donor.id.tag)
        ]
        tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='donor.id.tag')
        tmp.data[, freq.rel:=freq.abs/freq.total]
        tmp.data <- tmp.data[, .(donor.id.tag, pop.tag, freq.rel)]
        tmp.data <- as.data.table(spread(data=tmp.data, key=pop.tag, value=freq.rel, fill=0))
        return(tmp.data)
    })
    tmp.data <- merge(x=tmp.data[[1]], y=tmp.data[[2]], by='donor.id.tag', all=TRUE)

    # ---> Full data.
    plot.data <- plot.data[cell.type %in% aggr.set.main[['c40.lg']]]
    # Merge w/ antigen-specific fractions.
    plot.data <- merge(
        x=plot.data, y=tmp.data,
        by='donor.id.tag', all.x=TRUE, all.y=FALSE
    )
    # Merge w/ donor clinics and demographics
    #       CMV serology.
    plot.data <- merge(x=plot.data, y=ser.cmv.data, by='donor.id.tag')
    # setorderv(x=tmp.data, cols='cmv.igg')
    plot.data[is.na(cmv.igg), cmv.igg:=0]
    #       Other clinical an demographical variables.
    plot.data <- merge(x=plot.data, y=donor.meta, by='donor.id.tag', all.x=TRUE)

    # ---> Association between predictions and continuous variables.
    tmp.reports.path <- paste0(fig.as.phen.comps.path, '/assoc_vars-cont_attemp-', att.lab)
    if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
    # Define continuous variables.
    gen.vars <- c(
        `CMV-IgG`='cmv.igg',
        `Age`='age',
        `BMI`='bmi',
        `AlcoholUnits`='alcohol.units',
        `VaccSpan-SARS`='covid.shot.span',
        `VaccSpan-IAV`='flu.shot.span'
    )
    cont.vars <- lapply(X=aggr.set.main[['c40.lg']], FUN=function(t.subset){
        tmp.vals <- colnames(plot.data)[colnames(plot.data) %like% paste0('^', gen.cell.types[t.subset])]
        names(tmp.vals) <- str_replace(string=tmp.vals, pattern='\\.', replacement='-')
        tmp.vals <- c(gen.vars, tmp.vals)
        return(tmp.vals)
    })
    names(cont.vars) <- aggr.set.main[['c40.lg']]
    # Process per cell type.
    for(t.subset in aggr.set.main[['c40.lg']]){
        for(cont.var in names(cont.vars[[t.subset]])){
            for(tmp.spc in fig.rel.spcs){
                if(!tmp.spc %in% tmp.rel.spcs[[t.subset]]) next
                tmp.data <- plot.data[
                    ag.group==tmp.spc,
                    .(
                        cont.var=get(cont.vars[[t.subset]][cont.var]),
                        freq.var=cell.fract,
                        tmp.var=cat.tag
                    )
                ]
                tmp.ggplot <- ggplot(data=tmp.data, aes(x=cont.var, y=freq.var, color=tmp.var)) +
                    geom_point(shape=1, stroke=5, size=12, color='black') +
                    geom_smooth(method='lm', formula=y~x, se=FALSE, fullrange=TRUE, linewidth=5) +
                    scale_x_continuous(breaks=scales::pretty_breaks(n=3)) +
                    scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
                    scale_color_manual(values=clusts.cols[[t.subset]]) +
                    labs(x='Continuous variable', y='% of ag-specific T cells')
                tmp.lab <- paste0(
                    '/Association_', tmp.spc, '_',  cont.var, '_CellType-', gen.cell.types[t.subset]
                )
                publish.plot(
                    tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
                    stat.cor=TRUE, cor.group='tmp.var',
                    blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=14, height=14
                )
            }
        }
    }

    # ---> Association between predictions and discrete variables.
    tmp.reports.path <- paste0(fig.as.phen.comps.path, '/assoc_vars-disc_attemp-', att.lab)
    if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
    # Define continuous variables.
    disc.vars <- c(
        `CMV-Serology`='cmv.status',
        `Gender`='gender',
        `SmokingStatus`='smoking.status'
    )
    # Plotting settings.
    tmp.bounds <- c(
        # `CD4.CMV.CMV-Serology`=0.06,
        # `CD8.CMV.CMV-Serology`=0.20
    )
    disc.var.cols <- list(
        `CMV-Serology`=c(`Positive`=unname(pp.cols['CMV']), `Negative`='#CBCBCB')
    )
    # Process per cell type.
    for(t.subset in aggr.set.main[['c40.lg']]){
        for(tmp.spc in fig.rel.spcs){
            if(!tmp.spc %in% tmp.rel.spcs[[t.subset]]) next
            for(disc.var in names(disc.vars)){
                tmp.vals <- sort(as.character(plot.data[, unique(cat.tag)]))
                tmp.file.name <- paste0(
                    tmp.reports.path, '/Comparison_', tmp.spc, '_',  disc.var, '_CellType-', gen.cell.types[t.subset], '.pdf'
                )
                pdf(file=tmp.file.name)
                for(tmp.val in tmp.vals){
                    tmp.data <- plot.data[
                        ag.group==tmp.spc &
                        cat.tag==tmp.val,
                        .(
                            disc.var=get(disc.vars[disc.var]),
                            freq.var=cell.fract
                        )
                    ]
                    # Stats.
                    if(tmp.data[, uniqueN(disc.var)==2]){
                        tmp.groups <- tmp.data[, unique(disc.var)]
                        tmp.test <- wilcox.test(
                            x=tmp.data[disc.var==tmp.groups[1], freq.var],
                            y=tmp.data[disc.var==tmp.groups[2], freq.var],
                            paired=FALSE
                        )
                        tmp.caption <- paste0('P-value is ', round(x=tmp.test$p.value, digits=6))
                    }else{
                        tmp.caption <- ''
                    }
                    # Set upper bound.
                    tmp.lab <- paste(gen.cell.types[t.subset], tmp.spc, disc.var, sep='.')
                    tmp.bound <- if(tmp.lab %in% names(tmp.bounds)) tmp.bounds[tmp.lab] else tmp.data[, max(freq.var)]
                    tmp.data[, bounded.var:=freq.var]
                    tmp.data[, point.status:='ori']
                    tmp.data[freq.var>tmp.bound, `:=`(bounded.var=tmp.bound, point.status='bounded')]
                    # Plot.
                    tmp.ggplot <- ggplot(data=tmp.data, aes(x=disc.var)) +
                        geom_boxplot(aes(y=freq.var, color=disc.var), width=0.7, linewidth=4, outlier.shape=NA) +
                        geom_jitter(aes(y=bounded.var, shape=point.status), stroke=5, size=10, color='black', width=0.22) +
                        scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
                        scale_shape_manual(values=tmp.shapes) +
                        # scale_color_manual(values=tmp.cols) +
                        coord_cartesian(ylim=c(0, tmp.bound)) +
                        labs(
                            title=paste0(tmp.val, ': ', clusts.defs[[t.subset]][tmp.val]),
                            x='', y='% of ag-specific T cells',
                            color='', caption=tmp.caption
                        ) +
                        theme(legend.position='none') +
                        theme_bw()
                    if(disc.var %in% names(disc.var.cols)){
                        tmp.ggplot <- tmp.ggplot +
                            scale_color_manual(values=disc.var.cols[[disc.var]])
                    }
                    print(tmp.ggplot)
                    # publish.plot(
                    #     tmp.ggplot=tmp.ggplot, output.path=fig.pred.val.path, file.name=tmp.lab, type='pdf',
                    #     blank.comp=blank.complement.3.1, do.legend=FALSE,
                    #     width=7, height=14
                    # )
                }
                dev.off()
            }
        }
    }

    return(NULL)
}

# @ Previous attempts discarded.

# @ Process for fractions of ag set-specific cells/clonotypes in the group (cluster) out of all cells in the group (cluster).
do.assoc.assess(
    plot.data=fig.data[['2']], att.lab='2',
    tmp.rel.spcs=subset.rel.spcs
)

# ----> Associations, summary.
att.lab <- '2'
tmp.reports.path <- paste0(fig.as.phen.comps.path, '/assoc_vars-cont-summ_attemp-', att.lab)
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
# Data retrieval
#       Response variables.
tmp.data.1 <- copy(fig.data[[att.lab]])
tmp.data.1 <- tmp.data.1[
    # cell.type %in% aggr.set.main[['c40.lg']] &
    (cell.type %in% aggr.set.main[['c40.lg']] |  cell.type=='c40.lg.cd8.trm') &
    ag.group %in% names(pp.cols),
    .(data.set=cell.type, ag.group, donor.id.tag, cluster=cat.tag, cell.fract)
]
tmp.data.2 <- lapply(X=clusts.defs, FUN=function(x){
    tmp.data <- data.table(
        cluster=names(x),
        cell.subset=x
    )
    return(tmp.data)
})
tmp.data.2 <- rbindlist(l=tmp.data.2, use.names=TRUE, idcol='data.set')
tmp.data.1 <- merge(
    x=tmp.data.1, y=tmp.data.2,
    by=c('data.set', 'cluster'), all.x=TRUE
)
tmp.data.1[data.set=='c40.lg.cd8.trm', data.set:='c40.lg.cd8.d40']
#       Explanatory variables.
tmp.data.2 <- lapply(X=names(aggr.set.main), FUN=function(data.set){
    tmp.lins <- aggr.set.main[[data.set]]
    if(data.set=='c40.lg') tmp.lins <- c(tmp.lins, 'c40.lg.cd8.trm')
    tmp.data.2 <- lapply(X=tmp.lins, FUN=function(tmp.lin){
        tmp.data <- as.data.table(srt.objs.list[[tmp.lin]]@meta.data)
        tmp.data[, barcode:=Cells(srt.objs.list[[tmp.lin]])]
        tmp.data.1 <- tmp.data[
            !is.na(donor.id.tag),
            .(freq.abs=.N),
            by=.(
                donor.id.tag,
                pop.tag=paste0(
                    gen.cell.types[tmp.lin], '.',
                    clusts.defs[[tmp.lin]][as.character(get(clust.labs[tmp.lin]))]
                )
            )
        ]
        tmp.data.2 <- tmp.data[
            !is.na(donor.id.tag),
            .(freq.total=.N),
            by=.(donor.id.tag)
        ]
        tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='donor.id.tag')
        tmp.data[, freq.rel:=freq.abs/freq.total]
        tmp.data <- tmp.data[, .(
            donor.id.tag,
            pop.tag=paste(str_replace(string=aggr.set.labs[data.set], pattern='C40-', replacement=''), pop.tag, sep='.'),
            freq.rel
        )]
        tmp.data <- as.data.table(spread(data=tmp.data, key=pop.tag, value=freq.rel, fill=0))
        return(tmp.data)
    })
    names(tmp.data.2) <- tmp.lins
    return(tmp.data.2)
})
tmp.data.2 <- c(tmp.data.2[[1]], tmp.data.2[[2]])
# lapply(X=tmp.data.2, FUN=colnames)
# Preflight.
rel.res.vars <- list( # Relevant response variables.
    `c40.lg.cd4.d40`=list(
        `Rel`=c('TRM', 'TCM', 'CD4-CTL'),
        `TRM`='TRM',
        `CTL`='CD4-CTL'
    ),
    `c40.lg.cd8.d40`=list(
        `Rel`=c('TRM', 'TEM', 'GZMKhi', 'TCM', 'BATFhi')
    )
)
rem.exp.vars <- list( # Non-relevant explanatory variables.
    `c40.lg.cd4.d40`=list(
        `Rel`=c('LG.CD4.Prolif'),
        `TRM`=c('LG.CD4.Prolif'),
        `CTL`=c('LG.CD4.Prolif')
    ),
    `c40.lg.cd8.d40`=list(
        `Rel`=c('LG.CD8.Prolif', 'LG.CD8.IL7Rhi')
    )
)
gen.vars <- c(
    `Age`='age',
    `BMI`='bmi',
    `AlcoholUnits`='alcohol.units',
    # `CMV-IgG`='cmv.igg',
    `VaccSpan-SARS`='covid.shot.span',
    `VaccSpan-IAV`='flu.shot.span'
)
# General process per T-cell lineage
for(tmp.lin in names(cell.type.cols)){
    # @ Aggr names for each site.
    lg.set <- aggr.set.main[['c40.lg']]
    lg.set <- lg.set[gen.cell.types[lg.set]==tmp.lin]
    # if(tmp.lin=='CD8') lg.set <- c(lg.set, 'c40.lg.cd8.trm')
    ln.set <- aggr.set.main[['c40.ln']]
    ln.set <- ln.set[gen.cell.types[ln.set]==tmp.lin]
    for(rel.pop.set in names(rel.res.vars[[lg.set]])){
        # @ "Response" variables
        to.plot.1 <- tmp.data.1[
            data.set==lg.set &
            ag.group %in% subset.rel.spcs[[lg.set]]
        ]
        # to.plot.1[, cell.subset:=clusts.defs[[lg.set]][as.character(cell.subset)]]
        to.plot.1 <- to.plot.1[cell.subset %in% rel.res.vars[[lg.set]][[rel.pop.set]]]
        to.plot.1[, key:=paste(ag.group, cell.subset, sep='.')]
        to.plot.1 <- to.plot.1[, .(donor.id.tag, key, cell.fract)]
        to.plot.1 <- as.data.table(spread(data=to.plot.1, key=key, value=cell.fract, fill=NA))
        var.set.1 <- setdiff(x=colnames(to.plot.1), y='donor.id.tag')
        tmp.cols <- paste(
            subset.rel.spcs[[lg.set]],
            rep(x=rel.res.vars[[lg.set]][[rel.pop.set]], each=length(subset.rel.spcs[[lg.set]])),
            sep='.'
        )
        var.set.1 <- tmp.cols[tmp.cols %in% var.set.1]
        # @ "Explanatory" variables
        to.plot.2 <- merge(
            x=tmp.data.2[[lg.set]], y=tmp.data.2[[ln.set]],
            by='donor.id.tag',
            all=TRUE
        )
        if(tmp.lin=='CD8') to.plot.2 <- merge(
            x=to.plot.2, y=tmp.data.2[['c40.lg.cd8.trm']],
            by='donor.id.tag',
            all=TRUE
        )
        #       Merge w/ donor clinics and demographics
        #       CMV serology.
        to.plot.2 <- merge(
            x=to.plot.2, y=ser.cmv.data[, .(donor.id.tag, cmv.igg)],
            by='donor.id.tag', all=TRUE
        )
        to.plot.2[is.na(cmv.igg), cmv.igg:=0]
        #       Other clinical an demographical variables.
        tmp.cols <- c('donor.id.tag', gen.vars[gen.vars %in% colnames(donor.meta)])
        to.plot.2 <- merge(
            x=to.plot.2, y=donor.meta[, ..tmp.cols],
            by='donor.id.tag', all=TRUE
        )
        var.set.2 <- setdiff(x=colnames(to.plot.2), y='donor.id.tag')
        tmp.cols <- paste('LG', gen.cell.types[lg.set], clusts.defs[[lg.set]], sep='.')
        if(tmp.lin=='CD8') tmp.cols <- c(
            tmp.cols,
            paste('LG', gen.cell.types['c40.lg.cd8.trm'], clusts.defs[['c40.lg.cd8.trm']], sep='.')
        )
        tmp.cols <- c(
            tmp.cols,
            # paste('LN', gen.cell.types[ln.set], clusts.defs[[ln.set]], sep='.'),
            gen.vars
        )
        var.set.2 <- tmp.cols[tmp.cols %in% var.set.2]
        tmp.cols <- rem.exp.vars[[lg.set]][[rel.pop.set]]
        if(all(!is.na(tmp.cols))) var.set.2 <- setdiff(x=var.set.2, y=tmp.cols)
        # @ Full set of vars.
        to.plot <- merge(x=to.plot.2, y=to.plot.1, by='donor.id.tag')
        to.plot <- as.data.frame(to.plot)
        row.names(to.plot) <- to.plot$donor.id.tag; to.plot$donor.id.tag <- NULL
        to.plot <- as.matrix(to.plot)
        # @ Pairwise distance and significance calculation
        corr.data <- rcorr(to.plot, type='spearman')
        r.mat <- corr.data$r
        r.mat <- r.mat[
            var.set.2,
            var.set.1
        ]
        p.mat <- corr.data$P
        p.mat <- p.mat[
            var.set.2,
            var.set.1
        ]
        # @ Correlation plot-like plot
        tmp.width <- ncol(r.mat)/2
        tmp.height <- nrow(r.mat)/2
        tmp.file.name <- paste0(
            tmp.reports.path, '/',
            obj.extended.names[lg.set], '_CorrPlot_Dist-Spearman_Set-', rel.pop.set,
            '.C.pdf'
        )
        pdf(file=tmp.file.name, width=tmp.width+2, height=tmp.height+2)
        corrplot(
            corr=r.mat,
            method = 'circle', type='full',
            col=rev(x=COL2('RdBu', 200)),
            diag=TRUE,
            p.mat=p.mat, sig.level=0.01, insig='blank',
            outline=TRUE,
            addgrid.col='black',
            tl.pos='lt',
            tl.cex=1, tl.col='black',
            cl.pos='r'
        )
        dev.off()
        tmp.file.name <- paste0(
            tmp.reports.path, '/',
            obj.extended.names[lg.set], '_CorrPlot_Dist-Spearman_Set-', rel.pop.set,
            '.B.pdf'
        )
        pdf(file=tmp.file.name, width=tmp.width, height=tmp.height)
        corrplot(
            corr=r.mat,
            method = 'circle', type='full',
            col=rev(x=COL2('RdBu', 200)),
            diag=TRUE,
            p.mat=p.mat, sig.level=0.01, insig='blank',
            outline=TRUE,
            addgrid.col='black',
            tl.pos='n',
            tl.cex=0.6, tl.col='black',
            cl.pos='n'
        )
        dev.off()
        # @ Color legend only.
        break.no <- 200
        col.breaks <- seq(from=-1, to=1, length.out=break.no)
        hmap.col.scale <- rev(COL2('RdBu', break.no+1))
        to.plot <- data.table(freq.total=col.breaks)
        # W/ labels
        tmp.ggplot <- ggplot(data=to.plot, aes(x=`freq.total`, y=`freq.total`, col=`freq.total`)) + 
            geom_point() +
            scale_color_gradientn(colors=hmap.col.scale, limits=c(-1, 1), breaks=scales::pretty_breaks(n=3)) +
            theme(
                legend.ticks=element_line(color='black', linewidth=0.6),
                legend.ticks.length=unit(0.22, "cm"),
                legend.frame=element_rect(color='black', linewidth=0.6)
            )
        tmp.ggplot <- get_legend(p=tmp.ggplot)
        tmp.file.name <- paste0(
            tmp.reports.path, '/',
            obj.extended.names[lg.set], '_CorrPlot_Dist-Spearman_Set-', rel.pop.set,
            '.L1.pdf'
        )
        pdf(file=tmp.file.name, height=2, width=2)
        print(as_ggplot(tmp.ggplot))
        dev.off()
        # W/out labels
        tmp.ggplot <- ggplot(data=to.plot, aes(x=`freq.total`, y=`freq.total`, col=`freq.total`)) + 
            geom_point() +
            scale_color_gradientn(colors=hmap.col.scale, limits=c(-1, 1), breaks=scales::pretty_breaks(n=3), name=NULL, labels=NULL) +
            theme(
                legend.ticks=element_line(color='black', linewidth=0.6),
                legend.ticks.length=unit(0.22, "cm"),
                legend.frame=element_rect(color='black', linewidth=0.6)
            )
        tmp.ggplot <- get_legend(p=tmp.ggplot)
        tmp.file.name <- paste0(
            tmp.reports.path, '/',
            obj.extended.names[lg.set], '_CorrPlot_Dist-Spearman_Set-', rel.pop.set,
            '.L2.pdf'
        )
        pdf(file=tmp.file.name, height=1.25, width=0.3)
        print(as_ggplot(tmp.ggplot))
        dev.off()
    }
}


### --------------------- Qualitative assessment -------------------- ###
### --------------------------- attempt 1 --------------------------- ###

tmp.reports.path <- paste0(fig.as.phen.comps.path, '/qual_gen_virus_attempt-1')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

# ---> Comparison among viruses for DEGs
subset.rel.spcs <- list(
    `c40.lg.cd4.d40`=list(
        `0`=rel.spcs,
        `2`=c('Aspergillus', 'CMV', 'EBV', 'IAV', 'SARS-CoV-2', 'MPV')
    ),
    `c40.lg.cd8.d40`=list(
        `0`=c('CMV', 'EBV', 'IAV', 'SARS-CoV-2')
    )
)
gene.list.1 <- list(
    `c40.lg.cd4.d40`=list(
        `0`=c(
            # ---> CMV
            # Up
            'GZMK', 'FKBP5', 'ZNF683',
            # Down
            # 'KLRB1', 'XCL1',
            # ---> IAV 
            # Up
            'LGALS3', 'GZMH', 'BHLHE40', 'LGALS1', 'CTSH',
            # Down
            # 'IL7R', 'GPR183',
            # ---> MPV
            # Up
            'BATF', 'IL7R', 'KLRB1', 'GPR25', 'CCL4', 'RUNX2', 'LGALS3BP', 'LCP1', 'IKZF3', 'CCR6', 'TXNIP', 'CRYBG1',
            # 'FKBP5',
            # 'GZMH',
            # Down
            # 'JUNB', 'FOS', 'NR4A2', 'NR4A1', 'XCL1', 'ZNF331', 'IFITM1', 'ISG15',
            # ---> PIV
            # Up
            'HLA-DQA2', 'HLA-DRB5', 'XCL2', 'XCL1', 'CD101', 'IL32', 'KLRG1', 'TIMP1', 'CD52', 'IL13', 'IL9R',
            # 'LGALS1', 'ZNF683',
            # 'BATF', 
            # Down
            # 'CCR6', 'IFI44', 'IFITM1',
            # ---> B. pertussis
            # Up
            'GZMB', 'CD7',
            # ---> Aspergillus
            # Up
            'IFITM1', 'MX1', 'TUBA1B', 'IRF1', 'IFI6', 'IFI44L',
            'CXCR4', 'TNFAIP3', 'CD9'
            # 'CCL4', 'GZMH',
        ),
        `2`=c(
            # ---> CMV
            # Upregulated
            # 'RGS1', 'PRDM1', 'BCL2', 'STAT1', 'NPDC1',
            # Downregulated
            'GNLY', 'PRF1', 'GZMB', 'FCGR3A',
            # 'GZMK',
            'SELL', #'IL7R',
            'HLA-DQA1', 'HLA-DRB1', 'BHLHE40', # 'HLA-DQB1',
            # 'PDCD1', 'TIGIT', 'CTLA4',
            # 'ISG20', 'IRF1',
            'KLRB1', 'KLRD1'
            # 'ZNF683', 'LGALS1'
        )
    )
)
gene.list.2 <- list(
    `c40.lg.cd4.d40`=list(
        `0`=c('TNF', 'IFNG')
    ),
    `c40.lg.cd8.d40`=list(
        `0`=c('TNF', 'IFNG')
    )
)
list.of.gene.lists <- list(
    `dea_res`=gene.list.1,
    `cytokines`=gene.list.2
)
do.norm.opts <- c(
    `dea_res`=TRUE,
    `cytokines`=FALSE
)
# Process per subset.
for(list.lab in names(list.of.gene.lists)){
    this.reports.path <- paste0(tmp.reports.path, '/', list.lab)
    if(!dir.exists(this.reports.path)) dir.create(this.reports.path)
    gene.lists <- list.of.gene.lists[[list.lab]]
    for(t.subset in names(gene.lists)){
        for(tmp.cluster in names(gene.lists[[t.subset]])){
            these.markers <- gene.lists[[t.subset]][[tmp.cluster]]
            do.norm.opt <- do.norm.opts[list.lab]
            if(do.norm.opt){
                scale.opt <- 'hot.and.cold'
            }else{
                scale.opt <- signatures.col.scale
            }
            # Plot
            tmp.ggplot <- dot.plot(
                seurat.obj=srt.objs.list[[t.subset]], features=these.markers, slot='data', ensembl=FALSE,
                do.norm=do.norm.opt,
                groups.tag='ag.spc.tag', groups.order=subset.rel.spcs[[t.subset]][[tmp.cluster]], groups.of.int=NULL,
                filter.tags=c(
                    'ag.spc.tag',
                    clust.labs[[t.subset]]
                ),
                groups.to.filter=list(
                    subset.rel.spcs[[t.subset]][[tmp.cluster]],
                    tmp.cluster
                ),
                keep=TRUE, na.rm=TRUE, feature.thold=NULL,
                this.color.scale=scale.opt,
                col.min=NULL, col.max=NULL,
                scale.by='radius', dot.scale=12,
                size.min=0, size.max=NA,
                file.name=NULL
            )
            tmp.ggplot <- tmp.ggplot + theme(legend.position='bottom')
            # Output.
            tmp.width <- (length(subset.rel.spcs[[t.subset]][[tmp.cluster]]) * 0.3) + 2
            tmp.height <- (length(these.markers) * 0.4) + 0.3
            tmp.lab <- paste0(obj.extended.names[t.subset], '_DotPlot_Comp-Pathogen_Pop-', tmp.cluster)
            publish.plot(
                tmp.ggplot=tmp.ggplot, output.path=this.reports.path, file.name=tmp.lab, type='pdf', blank.comp=blank.complement.2, do.rotate=TRUE, do.legend=TRUE,
                width=tmp.width, height=tmp.height, c.height=NULL, c.width=tmp.width+1.5, legend.height=3.2, legend.width=0.5
            )
        }
    }
}


### --------------------- Qualitative assessment -------------------- ###
### --------------------------- attempt 2 --------------------------- ###

tmp.reports.path <- paste0(fig.as.phen.comps.path, '/qual_gen_virus_attempt-2')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)


# # ---> Violin plots.
subset.rel.spcs <- list(
    `c40.lg.cd4.d40`=rel.spcs,
    # `c40.lg.cd4.d40`=c('IAV'),
    `c40.lg.cd8.d40`=c('CMV', 'EBV', 'SARS-CoV-2', 'IAV')
)
t.subset <- 'c40.lg.cd8.d40'
tmp.clusts <- '0'
tmp.genes <- c(
    'PRF1', 'GZMB',
    'LAG3', 'CXCR3', 'CXCR6', 'ITGAE',
    'TNF', 'IFNG',
    'HLA-DRA', 'HLA-DRB5',
    'BATF', 'BACH2',
    'CD40LG',
    'FABP5', 'RBPJ'
)
genes.w.ticks <- c('BATF', 'BACH2')
for(tmp.gene in tmp.genes){
    tmp.vars <- unname(c('ag.spc.tag', clust.labs[t.subset]))
    tmp.ggplot <- vln.plot(
        seurat.obj=srt.objs.list[[t.subset]], feature=tmp.gene, slot='data',
        groups.tag='ag.spc.tag', na.rm=TRUE, groups.of.int=NULL, groups.label=NULL,
        filter.tags=tmp.vars,
        groups.to.filter=list(
            subset.rel.spcs[[t.subset]],
            tmp.clusts
        ),
        keep=TRUE,
        feature.thold=NULL, vln.type='violin', color='burst.frequency', size.thold=0,
        adjust.val=1, trim.val=TRUE, box.color='white', plot.limits=c(0, 60), add.points=FALSE, this.color.scale=signatures.col.scale, do.fill=TRUE, file.name=NULL
    )
    tmp.ggplot <- tmp.ggplot +
        scale_y_continuous(expand=expansion(add=c(0, 0)), breaks=scales::pretty_breaks(n=3))
    tmp.width <- length(subset.rel.spcs[[t.subset]])*3
    tmp.lab <- paste0(
        obj.extended.names[t.subset], '_Vln-', tmp.gene, '_Clusts-', paste0(tmp.clusts, collapse='-')
    )
    tmp.blank.comp <- if(tmp.gene %in% genes.w.ticks) blank.complement.3 else blank.complement.3.1
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
        blank.comp=tmp.blank.comp, do.legend=TRUE,
        width=tmp.width, height=14
    )
}


### --------------------- Qualitative assessment -------------------- ###
### --------------------------- attempt 3 --------------------------- ###

tmp.reports.path <- paste0(fig.as.phen.comps.path, '/qual_gen_virus_attempt-3')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

# ---> DEG count, summary agnostic to directionality.
subset.rel.pops <- list(
    `c40.lg.cd8.d40`=c('0', '1', '2', '3'),
    `c40.lg.cd4.d40`=c('0', '1', '2')
)
comp.no.bounds <- list(
    `c40.lg.cd8.d40`=c(3, 3),
    `c40.lg.cd4.d40`=c(2, 4)
)
for(t.subset in names(main.extended.names)){
    tmp.data <- apc.res[[t.subset]]
    tmp.pops <- subset.rel.pops[[t.subset]]
    tmp.data <- tmp.data[
        sig.comp.no>=comp.no.bounds[[t.subset]][1] &
        cluster.tag %in% tmp.pops,
        .(deg.count=.N),
        by=.(cluster.tag, consensus.pred, sig.comp.no)
    ]
    tmp.data[
        sig.comp.no>=comp.no.bounds[[t.subset]][2],
        sig.comp.no:=comp.no.bounds[[t.subset]][2]
    ]
    tmp.data <- tmp.data[,
        .(deg.count=sum(deg.count)),
        by=.(cluster.tag, consensus.pred, sig.comp.no)
    ]
    # tmp.data[,
    #     group:=paste(cluster.tag, consensus.pred, sep='.')
    # ]
    # Set factors.
    tmp.lvls <- as.character(sort(tmp.data[, unique(sig.comp.no)]))
    tmp.vals <- factor(x=as.character(tmp.data[, sig.comp.no]), levels=tmp.lvls)
    set(x=tmp.data, j='sig.comp.no', value=tmp.vals)
    tmp.lvls <- tmp.data[, gtools::mixedsort(unique(cluster.tag))]
    if(t.subset=='c40.lg.cd8.d40') tmp.lvls <- rev(tmp.lvls)
    tmp.vals <- factor(x=as.character(tmp.data[, cluster.tag]), levels=tmp.lvls)
    set(x=tmp.data, j='cluster.tag', value=tmp.vals)
    tmp.lvls <- names(pp.cols)[names(pp.cols) %in% tmp.data[, consensus.pred]]
    if(t.subset=='c40.lg.cd8.d40') tmp.lvls <- rev(tmp.lvls)
    tmp.vals <- factor(x=as.character(tmp.data[, consensus.pred]), levels=tmp.lvls)
    set(x=tmp.data, j='consensus.pred', value=tmp.vals)
    # Plot
    tmp.ggplot <- ggplot(data=tmp.data, aes(x=consensus.pred, y=deg.count, fill=as.character(sig.comp.no))) +
        geom_bar(stat='identity', position='stack', width=0.7, linewidth=3, color='black') +
        facet_wrap(facets=~cluster.tag, nrow=1) +
        scale_fill_brewer() +
        # scale_y_continuous(expand=expansion(add=c(10, 10)), breaks=scales::pretty_breaks(n=3)) +
        scale_y_continuous(expand=expansion(add=c(0, 5)), breaks=scales::pretty_breaks(n=3)) +
        labs(x='Group', y='Number of DEGs', fill='Comp. no')
    tmp.lab <- paste0(main.extended.names[t.subset], '_DEGSummCount_DirAgnostic')
    tmp.width <- if(t.subset=='c40.lg.cd4.d40') 3.4 else 1.8
    tmp.width <- tmp.data[, uniqueN(cluster.tag)]*tmp.data[, uniqueN(cluster.tag)]*tmp.width
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
        blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=tmp.width, height=12
    )
}

# ---> DEA heatmap, full set of DEGs.
subset.rel.pops <- list(
    `c40.lg.cd8.d40`=c('0', '1', '2'),
    `c40.lg.cd4.d40`=c('0', '2')
)
subset.rel.spcs <- list(
    `c40.lg.cd4.d40`=rel.spcs,
    `c40.lg.cd8.d40`=c('CMV', 'EBV', 'IAV', 'SARS-CoV-2')
)
comp.no.tholds <- list(
    `c40.lg.cd4.d40`=c(1, 2, 3, 4),
    `c40.lg.cd8.d40`=c(1, 2, 3)
)
for(t.subset in names(main.extended.names)){
    tmp.data <- apc.res[[t.subset]]
    tmp.pops <- subset.rel.pops[[t.subset]]
    tmp.tholds <- comp.no.tholds[[t.subset]]
    for(rel.pop in tmp.pops){
        for(tmp.thold in tmp.tholds){
            # Retrieve gene order according to LFC and general significance
            data.to.order <- tmp.data[cluster.tag==rel.pop]
            data.to.order <- data.to.order[, 
                .(
                    comp.subset=consensus.pred,
                    gene.id=gene.id,
                    p.adj=ifelse(
                        test=!is.na(p.adj..max.sig) & sig.comp.no>=tmp.thold,
                        yes=p.adj..max.sig,
                        no=p.adj..max.gen
                    ),
                    lfc=ifelse(
                        test=!is.na(lfc..min.sig) & sig.comp.no>=tmp.thold,
                        yes=lfc..min.sig,
                        no=lfc..min.gen
                    )
                )
            ]
            tmp.groups <- names(pp.cols)[names(pp.cols) %in% data.to.order[, comp.subset]]
            data.to.order <- set.scalonated.order(
                dea.results=data.to.order,
                feat.id.col='gene.id', element.col='comp.subset', group.order=tmp.groups,
                p.col='p.adj', p.thold=0.05,
                lfc.col='lfc', lfc.thold=0.25,
                lfc.type='upregulated',
                sort.metric='lfc'
            )
            data.to.order <- data.to.order[['metric.mat']]
            # Retrieve expression stats.
            exp.stats <- get.exp.stats(
                seurat.obj=srt.objs.list[[t.subset]], norm.method='cpms',
                group.tag='ag.spc.tag', na.rm=TRUE,
                filter.tags=c(
                    clust.labs[t.subset],
                    'ag.spc.tag'
                ),
                groups.to.filter=list(
                    c(rel.pop),
                    subset.rel.spcs[[t.subset]]
                ),
                to.keep=TRUE,
                min.cell.count=10, do.pos.only=TRUE
            )
            exp.stats <- exp.stats[['exp.stats']]
            # Tidy data to get plot input
            data.to.plot <- exp.stats[, .(feature, group, stat=mean)]
            data.to.plot <- spread(data=data.to.plot, key=group, value=stat)
            row.names(data.to.plot) <- data.to.plot$feature; data.to.plot$feature <- NULL
            data.to.plot <- as.matrix(data.to.plot)
            data.to.plot <- data.to.plot[row.names(data.to.order), colnames(data.to.order)]
            data.to.plot <- scale(x=t(data.to.plot))
            data.to.plot <- t(data.to.plot)
            # Set color scale and breaks for heatmap.
            col.breaks <- seq(
                from=min(range(data.to.plot, na.rm=TRUE)),
                to=max(range(data.to.plot, na.rm=TRUE)),
                length.out=100
            )
            mid.point <- which.min(abs(col.breaks - 0))
            hmap.col.scale.1 <- colorRampPalette(c('blue', 'mediumblue', 'black'))(mid.point)
            hmap.col.scale.2 <- colorRampPalette(c('black', 'gold', 'yellow'))(100-(mid.point+1))
            hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
            # Set columns metadata.
            col.meta <- data.frame(
                row.names=colnames(data.to.plot),
                `Specificity`=colnames(data.to.plot),
                stringsAsFactors=FALSE
            )
            meta.cols <- list(
                `Specificity`=pp.cols[names(pp.cols) %in% row.names(col.meta)]
            )
            # Plot.
            tmp.width <- ncol(data.to.plot)*0.6
            tmp.file.name <- paste0(
                tmp.reports.path, '/',
                obj.extended.names[t.subset], '_DEGExp_Pop-', rel.pop, '_DEGCountThold-', tmp.thold,
                '.C.pdf'
            )
            pheatmap(
                mat=data.to.plot,
                color=hmap.col.scale, breaks=col.breaks, scale='none',
                na_col='#b3b3b3', border_color=NA,
                cluster_rows=FALSE, cluster_cols=FALSE,
                annotation_col=col.meta, annotation_names_col=FALSE, annotation_colors=meta.cols,
                show_rownames=FALSE, show_colnames=FALSE,
                filename=tmp.file.name, height=7, width=tmp.width
            )
            tmp.file.name <- paste0(
                tmp.reports.path, '/',
                obj.extended.names[t.subset], '_DEGExp_Pop-', rel.pop, '_DEGCountThold-', tmp.thold,
                '.B.pdf'
            )
            pheatmap(
                mat=data.to.plot,
                color=hmap.col.scale, breaks=col.breaks, scale='none',
                na_col='#b3b3b3', border_color=NA,
                cluster_rows=FALSE, cluster_cols=FALSE,
                annotation_col=col.meta, annotation_names_col=FALSE, annotation_colors=meta.cols,
                show_rownames=FALSE, show_colnames=FALSE, legend=FALSE, annotation_legend=FALSE,
                filename=tmp.file.name, height=7, width=tmp.width
            )
            # @ Color legend only.
            # W/ labels
            # tmp.breaks <- if(tmp.group %in% names(ag.group.breaks)) NULL else ag.group.breaks[[tmp.group]]
            tmp.breaks <- NULL
            data.to.plot <- as.data.table(gather(data=as.data.frame(data.to.plot), key=key, value=value))
            data.to.plot <- data.to.plot[!is.na(value)]
            if(is.null(tmp.breaks)){
                tmp.ggplot <- ggplot(data=data.to.plot, aes(x=value, y=value, col=value)) + 
                    geom_point() +
                    scale_color_gradientn(colors=hmap.col.scale) +
                    theme(
                        legend.ticks=element_line(color='black', linewidth=0.6),
                        legend.ticks.length=unit(0.22, "cm"),
                        legend.frame=element_rect(color='black', linewidth=0.6)
                    )
            }else{
                tmp.ggplot <- ggplot(data=data.to.plot, aes(x=value, y=value, col=value)) + 
                    geom_point() +
                    scale_color_gradientn(colors=hmap.col.scale) +
                    scale_color_gradientn(colors=hmap.col.scale, breaks=tmp.breaks) +
                    theme(
                        legend.ticks=element_line(color='black', linewidth=0.6),
                        legend.ticks.length=unit(0.22, "cm"),
                        legend.frame=element_rect(color='black', linewidth=0.6)
                    )
            }
            tmp.ggplot <- get_legend(p=tmp.ggplot)
            tmp.file.name <- paste0(
                tmp.reports.path, '/',
                obj.extended.names[t.subset], '_DEGExp_Pop-', rel.pop, '_DEGCountThold-', tmp.thold,
                '.L1.pdf'
            )
            pdf(file=tmp.file.name, height=2, width=2)
            print(as_ggplot(tmp.ggplot))
            dev.off()
            # W/out labels
            tmp.ggplot <- ggplot(data=data.to.plot, aes(x=value, y=value, col=value)) + 
                geom_point() +
                scale_color_gradientn(colors=hmap.col.scale, name=NULL, labels=NULL) +
                theme(
                    legend.ticks=element_line(color='black', linewidth=0.6),
                    legend.ticks.length=unit(0.22, "cm"),
                    legend.frame=element_rect(color='black', linewidth=0.6)
                )
            tmp.ggplot <- get_legend(p=tmp.ggplot)
            tmp.file.name <- paste0(
                tmp.reports.path, '/',
                obj.extended.names[t.subset], '_DEGExp_Pop-', rel.pop, '_DEGCountThold-', tmp.thold,
                '.L2.pdf'
            )
            pdf(file=tmp.file.name, height=1.25, width=0.3)
            print(as_ggplot(tmp.ggplot))
            dev.off()
        }
    }
}


### --------------------- Qualitative assessment -------------------- ###
### --------------------------- attempt 4 --------------------------- ###

tmp.reports.path <- paste0(fig.as.phen.comps.path, '/qual_gen_virus_attempt-4')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

# ---> Preflights
# @ Function to create heatmap depending of stat choice
plot.hmap.1 <- function(
    exp.stats, cell.counts,
    gene.list, process.lab,
    donor.meta.choice, donor.cols.choice=NULL,
    stat.choice=c('mean', 'mean.pos', 'prop'),
    cell.thold=10,
    cluster.cols=FALSE, cluster.rows=FALSE,
    hmap.height=7, hmap.width=14
){
    # Preflights
    stat.choice <- match.arg(arg=stat.choice)
    tmp.vals <- c(
        `mean`='AllMean',
        `mean.pos`='ExpMean',
        `prop`='ExpProp'
    )
    stat.lab <- tmp.vals[stat.choice]
    # Retrieve data to plot.
    tmp.groups <- cell.counts[cell.count>=cell.thold, group.tag]
    to.plot <- exp.stats[group %in% tmp.groups]
    to.filter <- to.plot[, .(to.check=any(prop>0) & any(mean>0)), by=feature]
    to.filter <- to.filter[to.check==TRUE, feature]
    gene.list <- tmp.list[gene.list %in% to.filter]
    to.plot <- to.plot[
        feature %in% gene.list,
        .(group, feature, stat=get(stat.choice))
    ]
    # Tidy data
    to.plot <- spread(data=to.plot, key=group, value=stat, fill=NA)
    row.names(to.plot) <- to.plot$feature; to.plot$feature <- NULL
    to.plot <- as.matrix(to.plot)
    to.plot <- t(scale(x=t(to.plot), center=TRUE, scale=TRUE))
    to.plot <- to.plot[gene.list, ]
    # Set metadata.
    col.meta <- data.frame(
        Row=colnames(to.plot),
        `Specificity`=str_extract(string=colnames(to.plot), pattern='^[^\\.]+'),
        `Donor`=str_extract(string=colnames(to.plot), pattern='L\\d+$')
    )
    if(any(is.na(col.meta$Donor))){
        if(is.null(donor.meta.choice)){
            col.meta$Donor <- NULL
        }else{
            stop('Faulty donor info.\n')
        }
    }
    if(!is.null(donor.meta.choice)){
        col.meta <- merge(
            x=col.meta, y=donor.meta.choice,
            by.x='Donor', by.y='donor.id.tag',
            all.x=TRUE, all.y=FALSE
        )
    }
    row.names(col.meta) <- col.meta$Row; col.meta$Row <- NULL
    # Define colors for metadata tracks.
    tmp.cols <- pp.cols[names(pp.cols) %in% col.meta$Specificity]
    ann.colors <- list(
        `Specificity`=tmp.cols
    )
    if(!is.null(donor.cols.choice) & is.list(donor.cols.choice)){
        ann.colors <- c(ann.colors, donor.cols.choice)
    }
    # Set color scale and breaks for heatmap.
    break.no <- 200
    col.breaks <- seq(from=min(range(to.plot, na.rm=TRUE)), to=max(range(to.plot, na.rm=TRUE)), length.out=break.no)
    # mid.point <- (max(col.breaks) - min(col.breaks)) / 2
    mid.point <- 0
    mid.point <- which.min(abs(col.breaks-mid.point))
    # hmap.col.scale.1 <- colorRampPalette(c('#ffffff', '#ffff9a', '#ffd34d', '#ffc04d', '#ffa500', '#ff6500'))(mid.point)
    # hmap.col.scale.2 <- colorRampPalette(c('#ff6500', '#ff4d4d', '#ff0000', '#b30000', '#670000'))(break.no-(mid.point+1))
    hmap.col.scale.1 <- colorRampPalette(c('blue', 'mediumblue', 'black'))(mid.point)
    hmap.col.scale.2 <- colorRampPalette(c('black', 'gold', 'yellow'))(break.no-(mid.point+1))
    hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
    # Complete version.
    tmp.file.name <- paste0(this.reports.path, '/', process.lab, '_Thold-', cell.thold, '_', stat.lab, '.C.pdf')
    pheatmap(
        mat=to.plot, scale='row',
        color=hmap.col.scale, breaks=col.breaks, border_color='black',
        cluster_rows=cluster.rows, cluster_cols=cluster.cols,
        annotation_col=col.meta, annotation_colors=ann.colors,
        show_colnames=FALSE, show_rownames=TRUE, legend=TRUE, annotation_legend=TRUE,
        filename=tmp.file.name, height=hmap.height, width=hmap.width
    )
    tmp.file.name <- paste0(this.reports.path, '/', process.lab, '_Thold-', cell.thold, '_', stat.lab, '.B.pdf')
    pheatmap(
        mat=to.plot, scale='row',
        color=hmap.col.scale, breaks=col.breaks, border_color='black',
        cluster_rows=cluster.rows, cluster_cols=cluster.cols,
        annotation_col=col.meta, annotation_colors=ann.colors,
        show_colnames=FALSE, show_rownames=FALSE, legend=FALSE, annotation_legend=FALSE,
        filename=tmp.file.name, height=hmap.height, width=hmap.width*0.4
    )
    # Cell counts per group
    tmp.groups <- colnames(to.plot)
    to.plot <- cell.counts[group.tag %in% tmp.groups]
    tmp.vals <- factor(x=to.plot$group.tag, levels=tmp.groups)
    set(x=to.plot, j='group.tag', value=tmp.vals)
    tmp.ggplot <- ggplot(data=to.plot, aes(x=group.tag, y=cell.count)) +
        geom_bar(stat='identity', position='dodge', width=0.7, linewidth=1, color='black', fill='#d3d3d3') +
        geom_hline(yintercept=30, color='red', linetype='dashed', linewidth=1.5) +
        # scale_y_continuous(expand=expansion(add=c(0, 15)), breaks=scales::pretty_breaks(n=3)) +
        labs(x='Group', y='Number of cells')
    tmp.lab <- paste0(
        '/', process.lab, '_Thold-', cell.thold, '_CellCounts'
    )
    publish.plot(
        tmp.ggplot=tmp.ggplot, output.path=this.reports.path, file.name=tmp.lab, type='pdf',
        blank.comp=blank.complement.3.1, do.legend=FALSE, do.rotate=TRUE, width=10, height=5
    )
    return(NULL)
}
# @ Function to run the process for all available stats given a cell threshold.
full.hmap.1 <- function(
    exp.stats, cell.counts,
    gene.list, process.lab, cell.thold,
    donor.meta.choice=donor.meta.choice, donor.cols.choice=donor.cols.choice,
    cluster.rows=FALSE, cluster.cols=FALSE,
    hmap.height, hmap.width
){
    plot.hmap.1(
        exp.stats=exp.stats, cell.counts=cell.counts,
        gene.list=gene.list, process.lab=process.lab,
        donor.meta.choice=donor.meta.choice, donor.cols.choice=donor.cols.choice,
        stat.choice='mean',
        cell.thold=cell.thold,
        cluster.rows=cluster.rows, cluster.cols=cluster.cols,
        hmap.height=hmap.height, hmap.width=hmap.width
    )
    plot.hmap.1(
        exp.stats=exp.stats, cell.counts=cell.counts,
        gene.list=gene.list, process.lab=process.lab,
        donor.meta.choice=donor.meta.choice, donor.cols.choice=donor.cols.choice,
        stat.choice='mean.pos',
        cell.thold=cell.thold,
        cluster.rows=FALSE, cluster.cols=cluster.cols,
        hmap.height=hmap.height, hmap.width=hmap.width
    )
    plot.hmap.1(
        exp.stats=exp.stats, cell.counts=cell.counts,
        gene.list=gene.list, process.lab=process.lab,
        donor.meta.choice=donor.meta.choice, donor.cols.choice=donor.cols.choice,
        stat.choice='prop',
        cell.thold=cell.thold,
        cluster.rows=cluster.rows, cluster.cols=cluster.cols,
        hmap.height=hmap.height, hmap.width=hmap.width
    )
    return(NULL)
}
# @ DEG lists
gene.list.1 <- list(
    `c40.lg.cd8.d40`=list(
        `0`=c(
            # ---> CMV
            'XCL1', 'CD7', 'NFKBIZ', 'NKG7', 'GZMK',
            # ---> EBV
            'PRF1', 'GNLY', 'CST7', 'PDCD1', 'RHOC', 'CISH',
            # ---> IAV 
            'CD40LG', 'RBPJ',
            'BATF', 'TNF', 'IFNG', 'CTSH',
            'ZNF683', 'FABP5', 'CXCR6',
            # 'ITGAE', 'LGALS3', 'HOPX', 'LGALS1',
            # 'GZMH', 'GZMB',
            'HLA-DRB1', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DRB5', 'HLA-DRA', 'HLA-DQA1',
            # 'HLA-DQA2', 'HLA-DMA',
            'TNFSF14', 'TOX2', 'RORC',
            # 'MIR155HG', 'IL5RA', 'CAPG', 'CPNE7', 'MT2A', 'CD74', 'FKBP5', 'STMN1', 'IL6ST', 'MAP4K1',
            # ---> SARS-CoV-2
            'IL7R', 'CD55', 'AREG', 'IL2',
            'KLRD1', 'KLRC1', 'KLRC2', 'KLRB1', 'ICOS', 'CXCR4'
            # 'GPR183','ERN1', 'SPRY1', 'PPP1R2C'
        ),
        `1`=c(
            # ---> CMV
            'GNLY', 'TNF', 'IFNG', 'CX3CR1', 'ICOS', 'BHLHE40', 'LAG3',
            # 'TNFSF9', 'FCRL3', 'FCRL6', 'CD69',
            'IFI6', 'MX1', 'MX2', 'ISG15', 'IRF7', 'IRF9', 'IFI35', 'IFI44', 'IFI44L', 'IFIT3', 'OASL', 'OAS1', 'OAS2', 'STAT1', 'STAT4',
            # 'IRF1', 'IFI16',
            'EIF2AK2', 'EIF4A1', 'EIF4A2', 'EIF4A3', 'EIF5', 'TUBA1A', 'TUBB2A', 'TUBB4B', 'MCL1',
            'XAF1', 'TNFAIP3',
            # 'CCL3', 'LGR6', 'SLAMF7', 'DDX60',
            # ---> EBV
            # 'GZMK', 'CRTAM', 'KIR3DL1', 'TPD52', 'ITGA4','HLA-DRB1',
            'BCL2', 'PDCD1', 'CXCR3', 'KLF2', 'IKZF2',
            'ZNF683',
            'CD55', 'CD27', 'GPR183', 'LEF1', 'COTL1',
            # ---> EBV & IAV
            'XCL1',
            # ---> IAV 
            'IKZF2', 'SELL', 'IL7R', 'S1PR1',
            'KLRC2', 'KLRC3', 'KLRC4', 'KLRF1', 'KLRG1', 'KIR2DL1', 'KIR2DL3',
            'TIGIT', 'CD300A', 'IL18RAP', 'COL6A2',
            # ---> SARS-CoV-2
            'PRF1', 'GZMB', 'GZMH', 'TBX21', 'FCGR3A', 'RUNX3', 'PRDM1', 'TNFRSF1B', 'FKBP5', 'FGFBP2', 'IL10RA',
            'ITGB1', 'LGALS1', 'LGALS3', 'CXCR4'
            # 'PRMT2', 'IGFBP7', 'ID2', 'GSTM2', 'IKZF1', 'AKR1C3', 'CEBPD'
        )
    ),
    `c40.lg.cd4.d40`=list(
        `0`=c(
            # ---> CMV
            # Up
            # 'GZMK', 'FKBP5', 'ZNF683',
            # Down
            # 'KLRB1', 'XCL1',
            # ---> IAV 
            # Up
            # 'LGALS3', 'GZMH', 'BHLHE40', 'LGALS1', 'CTSH',
            # Down
            # 'IL7R', 'GPR183',
            # ---> MPV
            # Up
            'BATF', 'IL7R', 'KLRB1', 'GPR25', 'CCL4', 'RUNX2', 'LGALS3BP', 'LCP1', 'FKBP5', 'GZMH', 'IKZF3', 'CCR6', 'TXNIP', 'CRYBG1',
            # Down
            # 'JUNB', 'FOS', 'NR4A2', 'NR4A1', 'XCL1', 'ZNF331', 'IFITM1', 'ISG15',
            # ---> PIV
            # Up
            'HLA-DQA2', 'HLA-DRB5', 'BATF', 'ZNF683', 'LGALS1', 'XCL2', 'XCL1', 'CD101', 'IL32', 'KLRG1', 'TIMP1', 'CD52', 'IL13', 'IL9R',
            # Down
            # 'CCR6', 'IFI44', 'IFITM1',
            # ---> B. pertussis
            # Up
            'GZMB', 'CD7',
            # ---> Aspergillus
            # Up
            'IFITM1', 'MX1', 'TUBA1B', 'IRF1', 'IFI6', 'IFI44L',
            'GZMH', 'CCL4', 'CXCR4', 'TNFAIP3', 'CD9'
        ),
        `2`=c(
            # ---> CMV
            # Upregulated
            'RGS1', 'PRDM1', 'BCL2', 'STAT1', 'NPDC1',
            # Downregulated
            'SELL', 'IL7R',
            'GNLY', 'GZMK', 'PRF1', 'GZMB', 'FCGR3A',
            'HLA-DQA1', 'HLA-DQB1', 'HLA-DRB1',
            'BHLHE40', 'PDCD1', 'TIGIT', 'CTLA4',
            'ISG20', 'IRF1',
            'KLRB1', 'KLRD1'
            # 'ZNF683', 'LGALS1'
        )
    )
)
gene.list.2 <- list(
    `c40.lg.cd4.d40`=list(
        `0`=c(
            row.names(srt.objs.list[['c40.lg.cd4.d40']])[str_detect(string=row.names(srt.objs.list[['c40.lg.cd4.d40']]), pattern='^IL\\d+')],
            'TNF', 'IFNG'
        )
    ),
    `c40.lg.cd8.d40`=list(
        `0`=c(
            row.names(srt.objs.list[['c40.lg.cd8.d40']])[str_detect(string=row.names(srt.objs.list[['c40.lg.cd8.d40']]), pattern='^IL\\d+')],
            'TNF', 'IFNG'
        )
    )
)
list.of.gene.lists <- list(
    `dea_res`=gene.list.1,
    `cytokines`=gene.list.2
)
clust.row.opts <- c(
    `dea_res`=FALSE,
    `cytokines`=TRUE
)
# @ Relevant specificities per cell subset
subset.rel.spcs <- list(
    `c40.lg.cd4.d40`=rel.spcs,
    `c40.lg.cd8.d40`=c('CMV', 'EBV', 'IAV', 'SARS-CoV-2')
)
# @ Tidy donor meta.
tmp.vals <- c('flu.shot.span', 'covid.shot.span', 'age', 'gender')
donor.meta.choice <- donor.meta.hmap[donor.var %in% tmp.vals]
donor.meta.choice <- spread(data=donor.meta.choice[, .(donor.id.tag, donor.var, fill.val)], key=donor.var, value=fill.val)
donor.cols.choice <- list(
    `gender`=c(
        'gender.Female'='#B63EB8',
        'gender.Male'='#78B852'
    ),
    `age`=c(
        'age.<66'='#26EDFF',
        'age.66-70'='#23A4E8',
        'age.71-75'='#3383FF',
        'age.>75'='#233DE8'
    ),
    `covid.shot.span`=c(
        'covid.shot.span.<1m'='#52A858',
        'covid.shot.span.>1m & <3m'='#66D16E',
        'covid.shot.span.>3m'='#79F581'
    ),
    `flu.shot.span`=c(
        'flu.shot.span.<1m'='#52A858',
        'flu.shot.span.>1m & < 3m'='#66D16E',
        'flu.shot.span.>3m'='#79F581'
    )

)

# @ Process per subset
for(list.lab in names(list.of.gene.lists)){
    this.reports.path <- paste0(tmp.reports.path, '/', list.lab)
    if(!dir.exists(this.reports.path)) dir.create(this.reports.path)
    gene.lists <- list.of.gene.lists[[list.lab]]
    for(t.subset in names(gene.lists)){
        for(tmp.cluster in names(gene.lists[[t.subset]])){
            # ---> Pathogen-specific DEGs in a cluster, donor-specific stats.
            # @ Preflights.
            # Define specificity/donor tag in seurat object.
            tmp.vals <- as.character(srt.objs.list[[t.subset]]@meta.data[, 'ag.spc.tag'])
            tmp.vals[!tmp.vals %in% subset.rel.spcs[[t.subset]]] <- NA
            tmp.vals[!is.na(tmp.vals)] <- paste(
                tmp.vals[!is.na(tmp.vals)],
                srt.objs.list[[t.subset]]@meta.data[
                    !is.na(tmp.vals), 'donor.id.tag'
                ],
                sep='.'
            )
            srt.objs.list[[t.subset]]@meta.data[, 'tmp.tag'] <- tmp.vals
            # Retrieve expression stats
            tmp.data <- get.exp.stats(
                seurat.obj=srt.objs.list[[t.subset]], norm.method='cpms',
                group.tag='tmp.tag', na.rm=TRUE,
                filter.tags=c(clust.labs[t.subset]), groups.to.filter=list(c(tmp.cluster)), to.keep=TRUE,
                min.cell.count=10, do.pos.only=TRUE
            )
            exp.stats <- tmp.data[['exp.stats']]
            cell.counts <- tmp.data[['group.counts']]
            srt.objs.list[[t.subset]]@meta.data[, 'tmp.tag'] <- NULL
            # Filter unexpressed cells from gene list
            tmp.list <- gene.lists[[t.subset]][[tmp.cluster]]
            # to.filter <- exp.stats[, .(to.check=any(prop>0)), by=feature]
            # to.filter <- to.filter[to.check==TRUE, feature]
            # tmp.list <- tmp.list[tmp.list %in% to.filter]

            # @ Processes
            # Process for 10 cell threshold
            tmp.lab <- paste0(obj.extended.names[t.subset], '_', tmp.cluster, '-DEGs_DonorSpc')
            full.hmap.1(
                exp.stats=exp.stats, cell.counts=cell.counts,
                gene.list=tmp.list, process.lab=tmp.lab, cell.thold=10,
                donor.meta.choice=donor.meta.choice, donor.cols.choice=donor.cols.choice,
                cluster.rows=clust.row.opts[list.lab], cluster.cols=FALSE,
                # cluster.rows=FALSE, cluster.cols=FALSE,
                hmap.height=length(tmp.list)*0.3,
                hmap.width=length(subset.rel.spcs[[t.subset]])*2
            )
            # Process for 30 cell threshold
            full.hmap.1(
                exp.stats=exp.stats, cell.counts=cell.counts,
                gene.list=tmp.list, process.lab=tmp.lab, cell.thold=30,
                donor.meta.choice=donor.meta.choice, donor.cols.choice=donor.cols.choice,
                cluster.rows=clust.row.opts[list.lab], cluster.cols=FALSE,
                # cluster.rows=FALSE, cluster.cols=FALSE,
                hmap.height=length(tmp.list)*0.3,
                hmap.width=length(subset.rel.spcs[[t.subset]])*2
            )
            # Process for 50 cell threshold
            full.hmap.1(
                exp.stats=exp.stats, cell.counts=cell.counts,
                gene.list=tmp.list, process.lab=tmp.lab, cell.thold=50,
                donor.meta.choice=donor.meta.choice, donor.cols.choice=donor.cols.choice,
                cluster.rows=clust.row.opts[list.lab], cluster.cols=FALSE,
                # cluster.rows=FALSE, cluster.cols=FALSE,
                hmap.height=length(tmp.list)*0.3,
                hmap.width=length(subset.rel.spcs[[t.subset]])*2
            )


            # ---> Pathogen-specific DEGs in a cluster, stats across donors.

            # @ Preflights.
            # Define specificity tag in seurat object.
            tmp.vals <- as.character(srt.objs.list[[t.subset]]@meta.data[, 'ag.spc.tag'])
            tmp.vals[!tmp.vals %in% subset.rel.spcs[[t.subset]]] <- NA
            srt.objs.list[[t.subset]]@meta.data[, 'tmp.tag'] <- tmp.vals
            # Retrieve expression stats
            tmp.data <- get.exp.stats(
                seurat.obj=srt.objs.list[[t.subset]], norm.method='cpms',
                group.tag='tmp.tag', na.rm=TRUE,
                filter.tags=c(clust.labs[t.subset]), groups.to.filter=list(c(tmp.cluster)), to.keep=TRUE,
                min.cell.count=10, do.pos.only=TRUE
            )
            exp.stats <- tmp.data[['exp.stats']]
            cell.counts <- tmp.data[['group.counts']]
            srt.objs.list[[t.subset]]@meta.data[, 'tmp.tag'] <- NULL
            # @ Processes
            # Process for 10 cell threshold
            tmp.lab <- paste0(obj.extended.names[t.subset], '_', tmp.cluster, '-DEGs_Summ')
            full.hmap.1(
                exp.stats=exp.stats, cell.counts=cell.counts,
                gene.list=tmp.list, process.lab=tmp.lab, cell.thold=10,
                donor.meta.choice=NULL, donor.cols.choice=NULL,
                cluster.rows=clust.row.opts[list.lab], cluster.cols=TRUE,
                hmap.height=length(tmp.list)*0.3,
                hmap.width=length(subset.rel.spcs[[t.subset]])
            )
        }
    }
}


### --------------------- Qualitative assessment -------------------- ###
### --------------------------- attempt 5 --------------------------- ###

tmp.reports.path <- paste0(fig.as.phen.comps.path, '/qual_gen_virus_attempt-5')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

# @ Module list per cell subset.
mod.list <- list(
    `c40.lg.cd8.d40`=list(
        `0`=c(
            'trm.signature.clarke',
            # 'trm.up.scott', 'trm.down.scott',
            'cell.cycling.best'#,
            # 'tcell.cytotoxicity.guo'
            # 'dixhaust.consensus2.vj'
        )
    )#,
    # `c40.lg.cd4.d40`=list()
)
# @ Relevant specificities per cell subset
subset.rel.spcs <- list(
    `c40.lg.cd4.d40`=rel.spcs,
    `c40.lg.cd8.d40`=c('CMV', 'EBV', 'SARS-CoV-2', 'IAV')
)
# @ Process per cell subset.
for(t.subset in names(mod.list)){
    for(tmp.cluster in names(mod.list[[t.subset]])){
        # Define cell subset to act on.
        cells.to.keep <- as.data.table(srt.objs.list[[t.subset]]@meta.data)
        cells.to.keep$barcode <- Cells(srt.objs.list[[t.subset]])
        cells.to.keep <- cells.to.keep[
            ag.spc.tag %in% subset.rel.spcs[[t.subset]] &
            get(clust.labs[t.subset])==tmp.cluster,
            barcode
        ]
        # Define annotations
        tmp.meta <- data.frame(
            row.names=subset.rel.spcs[[t.subset]],
            spc=subset.rel.spcs[[t.subset]]
        )
        tmp.cols <- list(
            spc=pp.cols
        )
        # Create heatmap
        tmp.pffx <- paste0(obj.extended.names[t.subset], '_', 'C-', tmp.cluster)
        do.fgsea(
            metrics.src=srt.objs.list[[t.subset]], 
            modules.feats=modules.feats[mod.list[[t.subset]][[tmp.cluster]]],
            metric='signal.to.noise', cell.subset=cells.to.keep,
            tag.of.int='ag.spc.tag', vals.to.depict=NULL, 
            output.path=tmp.reports.path, files.preffix=tmp.pffx,
            full.reports=FALSE,
            do.heatmap=TRUE, signature.order=NULL, ref.order=subset.rel.spcs[[t.subset]],
            col.metadata=tmp.meta, ann.colors=tmp.cols, just.get.res=FALSE
        )
    }
}


### ----------------- Comparison between lung & LN ------------------ ###

tmp.reports.path <- paste0(fig.as.phen.comps.path, '/comp_against_LN')
if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)

# ---> External phenotype of clonotypes overlapped between anatomical sites.
subset.rel.spcs <- list(
    `c40.lg.cd4.d40`=rel.spcs,
    `c40.lg.cd8.d40`=c('CMV', 'EBV', 'IAV', 'SARS-CoV-2')
)
for(data.set in names(aggr.set.labs)){
    tmp.vals <- meta.data.l[[data.set]][, unique(data.set)]
    for(t.subset in tmp.vals){
        anti.set <- if(data.set=='c40.lg') 'c40.ln' else 'c40.lg'
        anti.subset <- aggr.set.main[[anti.set]]
        anti.subset <- names(gen.cell.types[anti.subset])[gen.cell.types[anti.subset]==gen.cell.types[t.subset]]
        # Collect details on anti dataset, including bias info.
        tmp.data <- get.multi.site.info(
            data.set=data.set, t.subset=t.subset,
            anti.set=anti.set, anti.subset=anti.subset,
            rel.cols='ag.spc.tag'
        )
        # Define relevant specificities.
        if(t.subset %in% names(subset.rel.spcs)){
            tmp.lvls <- subset.rel.spcs[[t.subset]]
            tmp.vals <- factor(x=as.character(tmp.data$ag.spc.tag), levels=tmp.lvls)
            set(x=tmp.data, j='ag.spc.tag', value=tmp.vals)
        }
        # Plot for each cluster.
        tmp.clusts <- names(clusts.cols[[t.subset]])
        for(this.clust in tmp.clusts){
            # To plot.
            to.plot <- tmp.data[
                site.ovlp=='Overlapped' & !is.na(group.bias) &
                !is.na(ag.spc.tag) &
                cluster==this.clust
            ]
            # @ Set cell types (according to bias) order
            tmp.cols <- c(
                clusts.cols[[anti.subset]],
                bias.1.cols[names(bias.1.cols)=='Plastic']
            )
            tmp.lvls <- names(tmp.cols)
            tmp.vals <- factor(x=as.character(to.plot$group.bias), levels=tmp.lvls)
            set(x=to.plot, j='group.bias', value=tmp.vals)
            # Define file width according to donor amount.
            tmp.width <- (length(to.plot[, levels(ag.spc.tag)]) * 6.68/40) + 0.32
            tmp.width <- (tmp.width * 5)/1.8
            # Plot.
            tmp.ggplot <- ggplot(data=to.plot, aes(x=ag.spc.tag, fill=group.bias)) +
                geom_bar(position='stack', width=0.8, linewidth=0, alpha=1) +
                scale_fill_manual(values=tmp.cols) +
                scale_x_discrete(drop=FALSE) +
                scale_y_continuous(expand=c(0, 0), breaks=scales::pretty_breaks(n=3)) +
                labs(x='Cluster', y='Cell count', fill='Clonotype status')
            tmp.lab <- paste0('/', obj.extended.names[t.subset], '_CellCount_Cluster_BiasStatus-Anti_RefCluster-', this.clust)
            publish.plot(
                tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
                blank.comp=blank.complement.2, do.legend=FALSE, do.rotate=TRUE, width=tmp.width, height=5
            )
        }
    }
}

# @ General functions
# For circos initialization
lvl.att.1.init.circos <- function(segment.data, freq.rel.thold){ # lvl stands for lung vs LN
    # Initialize the Circos plot
    tmp.data <- segment.data[, .(width=sum(freq.rel)), by=.(t.pop)]
    tmp.vals <- Reduce(x=tmp.data$width, f=sum, accumulate=TRUE)
    tmp.vals <- c(0, tmp.vals)
    tmp.data[, `:=`(
        pos.start=tmp.vals[1:length(tmp.vals)-1],
        pos.end=tmp.vals[2:length(tmp.vals)]
    )]
    circos.par('track.height'=0.3)
    circos.initialize(
        sectors=tmp.data[, t.pop],
        sector.width=tmp.data[, width],
        xlim=as.matrix(tmp.data[, .(pos.start, pos.end)])
    )
    # Add sector tracks
    circos.track(ylim=c(0, 1), panel.fun=function(x, y){
        sector <- CELL_META$sector.index
        xlim <- CELL_META$xlim
        ylim <- CELL_META$ylim
        if(sector %like% 'LG'){
            tmp.set <- if(tmp.lin=='CD4') 'c40.lg.cd4.d40' else 'c40.lg.cd8.d40'
        }else{
            tmp.set <- if(tmp.lin=='CD4') 'c40.ln.cd4.d08' else 'c40.ln.cd8.d08'
        }
        tmp.site <- str_extract(string=sector, pattern='^LG|^LN')
        tmp.clust <- str_replace(string=sector, pattern='^LG\\.|^LN\\.', replacement='')
        # Site-defining color
        circos.rect(
            xleft=xlim[1], ybottom=0.66, xright=xlim[2], ytop=1,
            col=an.site.cols[tmp.site]
        )
        # Cluster-defining color
        circos.rect(
            xleft=xlim[1], ybottom=0, xright=xlim[2], ytop=0.66,
            col=clusts.cols[[tmp.set]][tmp.clust]
        )
        # Lines showing clonotype relative frequency.
        tmp.data <- segment.data[t.pop==sector & freq.rel>freq.rel.thold, .N]>0
        if(tmp.data){
            circos.lines(
                x=segment.data[t.pop==sector & freq.rel>freq.rel.thold, freq.cumm],
                y=rep(x=0.4, times=segment.data[t.pop==sector & freq.rel>freq.rel.thold, .N]),
                type='h', col='black'
            )
            # Rectangle distinguishing between singletons and clonotypes (expanded TCRs).
            tmp.lim <- segment.data[t.pop==sector & freq.rel>freq.rel.thold, max(freq.cumm)]
            circos.rect(
                xleft=xlim[1], ybottom=0, xright=tmp.lim, ytop=0.4,
                col=adjustcolor("#FFFFFF", alpha.f=0)
            )
            circos.rect(
                xleft=tmp.lim, ybottom=0, xright=xlim[2], ytop=0.4,
                col=adjustcolor("#FFFFFF", alpha.f=1)
            )
        }else{
            # Rectangle indicating singletons across segment.
            circos.rect(
                xleft=xlim[1], ybottom=0, xright=xlim[2], ytop=0.4,
                col=adjustcolor("#FFFFFF", alpha.f=1)
            )
        }
    }, track.height=0.22, bg.border=NA)
}
# To add links.
lvl.att.1.add.links <- function(link.data){
    for(idx in 1:nrow(link.data)){
        sector.1 <- link.data[idx, t.pop.lg]
        sector.2 <- link.data[idx, t.pop.ln]
        coords.1 <- c(link.data[idx, start.coord.lg], link.data[idx, end.coord.lg])
        coords.2 <- c(link.data[idx, start.coord.ln], link.data[idx, end.coord.ln])
        circos.link(
            sector.index1=sector.1, point1=coords.1,
            sector.index2=sector.2, point2=coords.2,
            col=adjustcolor("#FDFD96", alpha.f=0.3),
            border = "black"
        )
    }
}
# @ Process per lineage and antigen specificity
subset.rel.spcs <- list(
    # `CD4`=setdiff(x=rel.spcs, y=c('SARS-CoV-2', 'PIV', 'MPV')),
    `CD4`=c('CMV'),
    # `CD8`=c('CMV', 'EBV', 'IAV', 'SARS-CoV-2')
    `CD8`=c('CMV')
)
subset.rel.pops <- list(
    `c40.lg.cd4.d40`=c('0', '1', '2'),
    # `c40.lg.cd8.d40`=c('0', '1', '2', '3'),
    `c40.lg.cd8.d40`=c('0', '1', '2'),
    `c40.ln.cd4.d08`=c('0', '1', '2', '3', '4', '5', '6'),
    # `c40.ln.cd8.d08`=c('0', '1', '2', '3', '4', '5')
    `c40.ln.cd8.d08`=c('0', '1', '2')
)
tmp.lins <- sort(gen.tcr.data.1[[1]][, unique(gex.lin.tag)])
# tmp.tholds <- c(0, 0.05, 0.1, 0.2)
tmp.tholds <- c(0)
donor.set.1 <- c('L004', 'L027')
donor.set.2 <- setdiff(x=site.donor.ovlp, donor.set.1)
donor.sets <- list(
    `L004-L0027`=donor.set.1,
    `Rest`=donor.set.2
)
for(tmp.lin in tmp.lins[2]){
    cat(tmp.lin, '\n')
    tmp.spcs <- subset.rel.spcs[[tmp.lin]]
    for(tmp.spc in tmp.spcs){
        # ---> Preflights
        # Determine specificity-relevant donors.
        tmp.data <- meta.data.l[['c40.lg']][
            !is.na(clonotype.tag) & !is.na(donor.id.tag) &
            donor.id.tag %in% site.donor.ovlp &
            gex.lin.tag==tmp.lin &
            ag.spc.tag==tmp.spc,
            .(freq.abs=.N),
            by=.(donor.id.tag)
        ]
        for(donor.set.lab in names(donor.sets)){
            # spc.donor.set <- tmp.data[freq.abs>=limit.of.detect, donor.id.tag]
            spc.donor.set <- donor.sets[[donor.set.lab]]
            # ---> Segment-specific info
            site.orders <- c(
                `c40.lg`=-1,
                `c40.ln`=1
            )
            segment.data <- lapply(X=names(meta.data.l), FUN=function(tmp.aggr){
                tmp.set <- aggr.set.main[[tmp.aggr]]
                tmp.set <- tmp.set[gen.cell.types[tmp.set]==tmp.lin]
                tmp.data <- meta.data.l[[tmp.aggr]][
                    !is.na(clonotype.tag) & !is.na(donor.id.tag) &
                    donor.id.tag %in% spc.donor.set &
                    gex.lin.tag==tmp.lin &
                    get(clust.labs[tmp.set]) %in% subset.rel.pops[[tmp.set]] &
                    ag.spc.tag==tmp.spc
                ]
                tmp.data.1 <- tmp.data[,
                    .(freq.abs=uniqueN(barcode)),
                    by=.(
                        clonotype.tag, donor.id.tag,
                        t.pop=as.character(get(clust.labs[tmp.set]))
                    )
                ]
                tmp.data.2 <- tmp.data[,
                    .(freq.total=uniqueN(barcode)),
                    by=.(donor.id.tag)
                ]
                segment.data <- merge(x=tmp.data.1, y=tmp.data.2, by='donor.id.tag')
                segment.data[, freq.rel:=freq.abs/freq.total]
                setorderv(x=segment.data, cols=c('t.pop', 'freq.rel'), order=site.orders[tmp.aggr])
                # # Determine cummulative relative frequency.
                # segment.data$freq.cumm <- Reduce(x=segment.data$freq.rel, f=sum, accumulate=TRUE)
                return(segment.data)
            })
            names(segment.data) <- aggr.simp.labs[names(meta.data.l)]
            segment.data <- rbindlist(l=segment.data, use.names=TRUE, idcol='site')
            segment.data[, t.pop:=paste(site, t.pop, sep='.')]
            # Determine cummulative relative frequency.
            segment.data$freq.cumm <- Reduce(x=segment.data$freq.rel, f=sum, accumulate=TRUE)
            # ---> Link-specific info.
            # General links between anatomical sites.
            site.link.data <- lapply(X=gen.tcr.data.2, FUN=function(tmp.data){
                tmp.data <- tmp.data[
                    site.ovlp=='Overlapped' &
                    gex.lin.tag==tmp.lin, 
                    comb.ids
                ]
                tmp.data <- str_split(string=tmp.data, pattern='\\/')
                tmp.data <- unlist(tmp.data)
                tmp.data <- tmp.data[str_detect(string=tmp.data, pattern=';')]
                return(tmp.data)
            })
            tmp.check <- all(site.link.data[[1]] %in% site.link.data[[2]]) & all(site.link.data[[2]] %in% site.link.data[[1]])
            if(!tmp.check) stop('Unexpected error!\n')
            site.link.data <- site.link.data[[1]]
            site.link.data <- data.table(link=site.link.data)
            site.link.data <- as.data.table(separate(data=site.link.data, col=link, into=c('lg.id', 'ln.id'), sep=';'))
            # Population-specific links (ones to actually be depicted on p;ot)
            lg.pops <- segment.data[site=='LG', unique(t.pop)]
            ln.pops <- segment.data[site=='LN', unique(t.pop)]
            link.data <- lapply(X=lg.pops, FUN=function(lg.pop){
                tmp.data <- lapply(X=ln.pops, FUN=function(ln.pop){
                    tmp.cols <- c('donor.id.tag', 'clonotype.tag', 't.pop', 'freq.rel', 'freq.cumm')
                    tmp.data.1 <- segment.data[
                        site=='LG' & t.pop==lg.pop &
                        clonotype.tag %in% site.link.data[, lg.id],
                        ..tmp.cols
                    ]
                    tmp.data.1[, `:=`(start.coord=freq.cumm-freq.rel, end.coord=freq.cumm)]
                    tmp.data.2 <- segment.data[
                        site=='LN' & t.pop==ln.pop &
                        clonotype.tag %in% site.link.data[, ln.id],
                        ..tmp.cols
                    ]
                    tmp.data.2[, `:=`(start.coord=freq.cumm-freq.rel, end.coord=freq.cumm)]
                    tmp.data <- merge(
                        x=site.link.data, y=tmp.data.1,
                        by.x='lg.id', by.y='clonotype.tag',
                        all=FALSE
                    )
                    tmp.data <- merge(
                        x=tmp.data, y=tmp.data.2,
                        by.x='ln.id', by.y='clonotype.tag',
                        all=FALSE, suffixes=c('.lg', '.ln')
                    )
                    return(tmp.data)
                })
                tmp.data <- rbindlist(l=tmp.data, use.names=TRUE)
                return(tmp.data)
            })
            link.data <- rbindlist(l=link.data, use.names=TRUE)
            # ---> Process per relative frequency threshold.
            for(freq.rel.thold in tmp.tholds){
                # ---> Circos plot, option 1
                tmp.file.name <- paste0(
                    tmp.reports.path, '/',
                    'SimBetSites_Circos-1', '_Lin-', tmp.lin,
                    '_Spc-', tmp.spc,
                    '_Donors-', donor.set.lab,
                    '_PctThold-', freq.rel.thold*100,
                    '.B.pdf'
                )
                pdf(file=tmp.file.name, width=14, height=14)
                # Initialize the Circos plot
                lvl.att.1.init.circos(segment.data=segment.data, freq.rel.thold=freq.rel.thold)
                # Create links between sectors (representing shared TCRs)
                to.plot <- link.data[freq.rel.lg>freq.rel.thold]
                lvl.att.1.add.links(link.data=to.plot)
                # Clear the plot
                dev.off()
                circos.clear()
            }
        }
    }
}

# ---> Percentage of T cell clonally overlapped between anatomic sites per subset.
for(data.set in names(aggr.set.labs)){
    tmp.vals <- meta.data.l[[data.set]][, unique(data.set)]
    for(t.subset in tmp.vals){
        anti.set <- if(data.set=='c40.lg') 'c40.ln' else 'c40.lg'
        anti.subset <- aggr.set.main[[anti.set]]
        anti.subset <- names(gen.cell.types[anti.subset])[gen.cell.types[anti.subset]==gen.cell.types[t.subset]]
        # Collect details on anti dataset, including bias info.
        tmp.data.1 <- get.multi.site.info(
            data.set=data.set, t.subset=t.subset,
            anti.set=anti.set, anti.subset=anti.subset,
            rel.cols=NULL
        )
        # Add specificity info.
        tmp.data.2 <- meta.data.l[[data.set]][
            data.set==t.subset,
            .(barcode, donor.id.tag, clonotype.tag, meta.tag=TRUE, ag.spc.tag)
        ]
        tmp.data <- merge(
            x=tmp.data.1, y=tmp.data.2,
            by=c('barcode', 'donor.id.tag', 'clonotype.tag'),
            all.x=TRUE, all.y=FALSE
        )
        tmp.check <- tmp.data[, all(meta.tag)]
        if(!tmp.check) stop('Unexpected error!\n')
        # @ Set cluster order
        tmp.lvls <- rev(names(clusts.cols[[t.subset]]))
        tmp.vals <- factor(x=as.character(tmp.data$cluster), levels=tmp.lvls)
        set(x=tmp.data, j='cluster', value=tmp.vals)
        # Define file width according to donor amount.
        tmp.width <- (tmp.data[, uniqueN(cluster)] * 6.68/40) + 0.32
        tmp.width <- (tmp.width * 5)/1.8
        # ---> Donor-wise information.
        tmp.vals <- c('Cell', 'Clone')
        for(tmp.item in tmp.vals){
            # Overlap fraction per subset.
            if(tmp.item=='Cell'){
                to.plot <- tmp.data[
                    !is.na(ag.spc.tag),
                    .(ovlp.freq.rel=ifelse(
                            test=.N>limit.of.detect,
                            yes=.SD[site.ovlp=='Overlapped', .N]/.N,
                            no=Inf
                        )
                    ),
                    by=.(donor.id.tag, ag.spc.tag, cluster)
                ]
                to.plot <- to.plot[!is.infinite(ovlp.freq.rel)]
            }else{
                to.plot <- tmp.data[
                    !is.na(ag.spc.tag),
                    .(ovlp.freq.rel=ifelse(
                            test=.N>limit.of.detect,
                            yes=.SD[site.ovlp=='Overlapped', uniqueN(clonotype.tag)]/uniqueN(clonotype.tag),
                            no=Inf
                        )
                    ),
                    by=.(donor.id.tag, ag.spc.tag, cluster)
                ]
                to.plot <- to.plot[!is.infinite(ovlp.freq.rel)]
            }
            # # Set factors
            # tmp.lvls <- names(clusts.defs[[t.subset]])
            # tmp.vals <- factor(x=as.character(to.plot$cluster), levels=tmp.lvls)
            # set(x=to.plot, j='cluster', value=tmp.vals)
            # Plot.
            tmp.ggplot <- ggplot(data=to.plot, aes(x=cluster, y=ovlp.freq.rel)) +
                geom_boxplot(color='black', width=0.7, linewidth=4, outlier.shape=NA) +
                geom_jitter(aes(color=ag.spc.tag), shape=1, stroke=5, size=10, height=0, width=0.1) +
                scale_color_manual(values=pp.cols) +
                scale_x_discrete(drop=FALSE) +
                scale_y_continuous(expand=c(0, 0), limits=c(0, 1), breaks=scales::pretty_breaks(n=3)) +
                labs(x='Cluster', y=paste0(tmp.item, ' fraction'), fill='Cluster')
            tmp.lab <- paste0('/', obj.extended.names[t.subset], '_Ovlped', tmp.item, 'Fract_Cluster-Spc-Donor')
            publish.plot(
                tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
                blank.comp=blank.complement.3, do.legend=FALSE, do.rotate=TRUE, width=tmp.width*3, height=14
            )
        }
        # ---> Summary across donors.
        # Would need to chhange to show one bar per pathogen, if needed.
        # tmp.vals <- c('Cell', 'Clone')
        # for(tmp.item in tmp.vals){
        #     # Overlap fraction per subset.
        #     if(tmp.item=='Cell'){
        #         to.plot <- tmp.data[
        #             !is.na(ag.spc.tag),
        #             .(ovlp.freq.rel=.SD[site.ovlp=='Overlapped', .N]/.N),
        #             by=.(cluster, ag.spc.tag)
        #         ]
        #     }else{
        #         to.plot <- tmp.data[
        #             !is.na(ag.spc.tag),
        #             .(ovlp.freq.rel=.SD[site.ovlp=='Overlapped', uniqueN(clonotype.tag)]/uniqueN(clonotype.tag)),
        #             by=.(cluster, ag.spc.tag)
        #         ]
        #     }
        #     # Plot.
        #     tmp.ggplot <- ggplot(data=to.plot, aes(x=cluster, y=ovlp.freq.rel, fill=cluster)) +
        #         # geom_bar(stat='identity', width=0.8, linewidth=1.2, alpha=1, color='black') +
        #         geom_bar(stat='identity', width=0.8, linewidth=0, alpha=1, color='black') +
        #         scale_fill_manual(values=clusts.cols[[t.subset]]) +
        #         scale_y_continuous(expand=c(0, 0), limits=c(0, 1), breaks=scales::pretty_breaks(n=3)) +
        #         labs(x='Cluster', y=paste0(tmp.item, ' fraction'), fill='Cluster')
        #     tmp.lab <- paste0('/', obj.extended.names[t.subset], '_Ovlped', tmp.item, 'Fract_Cluster-Spc')
        #     publish.plot(
        #         tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
        #         blank.comp=blank.complement.7, do.legend=FALSE, do.rotate=TRUE, width=tmp.width, height=5
        #     )
        # }
    }
}


############    -----------------------------------------    ############
### --------------------------- Figure on --------------------------- ###
### ------------ Pathogen relatedness by immune response ------------ ###
############    -----------------------------------------    ############

# ----> Output directory.
fig.pat.rels.path <- paste0(reports.path, '/figure_on_pathogen_rels')
if(!dir.exists(fig.pat.rels.path)) dir.create(fig.pat.rels.path)
# ---> Comments: Main dataset for this section, CD8 T cell set.
main.obj.name <- NULL


### -------------------------- Main Figure -------------------------- ###

### ----------------- Associations among pathogens ------------------ ###

# ---> Function to define stats columns and specificities under consideration.
get.stat.cols <- function(these.data){
    p.cols <- colnames(these.data)[colnames(these.data) %like% '^p.adj[\\|\\.]{1}']
    tmp.rm <- c('p.adj..max.gen', 'p.adj..max.sig')
    p.cols <- p.cols[!p.cols %in% tmp.rm]
    names(p.cols) <- str_replace(string=p.cols, pattern='^p.adj[\\|\\.]{1}', replacement='')
    lfc.cols <- colnames(these.data)[colnames(these.data) %like% '^lfc[\\|\\.]{1}']
    tmp.rm <- c('lfc..min.gen', 'lfc..min.sig')
    lfc.cols <- lfc.cols[!lfc.cols %in% tmp.rm]
    names(lfc.cols) <- str_replace(string=lfc.cols, pattern='^lfc[\\|\\.]', replacement='')
    to.check <- length(c(
        setdiff(x=names(p.cols), y=names(lfc.cols)),
        setdiff(x=names(lfc.cols), y=names(p.cols)))
    ) == 0
    if(!to.check) stop('Unexpected error while defining unique specficities to consider.\n')
    ag.spcs <- names(p.cols)
    return(list(
        p.cols,
        lfc.cols,
        ag.spcs
    ))
}

# ---> Relatedness according to Jaccard index.

# @ Distance matrix calculation and corr plot.
clade.no <- c(
    `c40.lg.cd4.d40`=5,
    `c40.lg.cd8.d40`=2
)
for(t.subset in names(main.extended.names)[1]){
    # Retrieve DEA data.
    ind.ivb.res <- ivb.res[[t.subset]]
    # Define stat columns.
    tmp.data <- get.stat.cols(ind.ivb.res=ind.ivb.res)
    p.cols <- tmp.data[[1]]; lfc.cols <- tmp.data[[2]]; ag.spcs <- tmp.data[[3]]
    # Retrieve DEG set for each specificity.
    deg.sets <- lapply(X=ag.spcs, FUN=function(ag.spc){
        p.col <- p.cols[ag.spc]
        lfc.col <- lfc.cols[ag.spc]
        tmp.data <- ind.ivb.res[
            get(p.col)<=0.05 & get(lfc.col)>=0.25,
            gene.id
        ]
        return(tmp.data)
    })
    names(deg.sets) <- ag.spcs
    # Pairwise distance calculation
    dist.mat <- sapply(X=ag.spcs, FUN=function(ag.spc.1){
        dist.v <- unlist(lapply(X=ag.spcs, FUN=function(ag.spc.2){
            tmp.data <- get.jaccard(
                x=deg.sets[[ag.spc.1]],
                y=deg.sets[[ag.spc.2]]
            )
            return(tmp.data)
        }))
        names(dist.v) <- ag.spcs
        return(dist.v)
    })
    # Overlap significance estimation.
    p.mat <- sapply(X=ag.spcs, FUN=function(ag.spc.1){
        # cat(ag.spc.1, '\n')
        dist.v <- unlist(lapply(X=ag.spcs, FUN=function(ag.spc.2){
            # DEG sub-universe.
            deg.union <- unique(c(deg.sets[[ag.spc.1]], deg.sets[[ag.spc.2]]))
            # Overlapped DEGs.
            deg.ovlp <- intersect(x=deg.sets[[ag.spc.1]], y=deg.sets[[ag.spc.2]])
            # DEGs unique in G1.
            deg.u1 <- setdiff(x=deg.sets[[ag.spc.1]], y=deg.sets[[ag.spc.2]])
            # DEGs unique in G2.
            deg.u2 <- setdiff(x=deg.sets[[ag.spc.2]], y=deg.sets[[ag.spc.1]])
            # Non-DEG subuniverse.
            deg.non <- setdiff(x=ind.ivb.res[, unique(gene.id)], y=deg.union)
            # Define test input.
            tmp.data <- matrix(
                data=c(
                    length(deg.ovlp), length(deg.u2),
                    length(deg.u1), length(deg.non)
                ),
                nrow=2, byrow=TRUE
            )
            tmp.check <- sum(tmp.data) == ind.ivb.res[, uniqueN(gene.id)]
            if(!tmp.check) stop('Unexpected error while estimating overlap significance.\n')
            # Perform test
            tmp.data <- fisher.test(x=tmp.data)
            tmp.data <- tmp.data$p.value
            return(tmp.data)
        }))
        names(dist.v) <- ag.spcs
        return(dist.v)
    })
    # Scale interval.
    # tmp.int <- copy(dist.mat)
    diag(dist.mat) <- 0
    tmp.int <- max(dist.mat, na.rm=TRUE)
    tmp.int <- c(0, tmp.int)
    # Correlation plot-like plot
    tmp.file.name <- paste0(
        fig.pat.rels.path, '/',
        obj.extended.names[t.subset], '_CorrPlot_Dist-Jaccard',
        '.pdf'
    )
    pdf(file=tmp.file.name, width=15, height=15)
    corrplot(
        corr=dist.mat,
        col.lim=tmp.int,
        is.corr=FALSE,
        diag=FALSE,
        p.mat=p.mat, sig.level=0.01,
        order='hclust', hclust.method='average', addrect=clade.no[t.subset]
    )
    dev.off()
}


# ---> Relatedness according to DEG metric correlation.

# @ Distance matrix calculation and corr plot.
clade.no <- c(
    `c40.lg.cd4.d40`=4,
    `c40.lg.cd8.d40`=2
)
for(t.subset in names(main.extended.names)[1]){
    # Retrieve DEA data.
    ind.ivb.res <- ivb.res[[t.subset]]
    # Define stat columns.
    tmp.data <- get.stat.cols(ind.ivb.res=ind.ivb.res)
    p.cols <- tmp.data[[1]]; lfc.cols <- tmp.data[[2]]; ag.spcs <- tmp.data[[3]]
    # Define DEG universe.
    #       Genes that came up as a DEG in at least one of the comparisons.
    deg.sets <- lapply(X=ag.spcs, FUN=function(ag.spc){
        p.col <- p.cols[ag.spc]
        lfc.col <- lfc.cols[ag.spc]
        tmp.data <- ind.ivb.res[
            get(p.col)<=0.05 &
            (get(lfc.col)>=0.25 | get(lfc.col)<=(-0.25)),
            gene.id
        ]
        return(tmp.data)
    })
    deg.uni <- unique(unlist(deg.sets))
    # Pairwise distance and significance calculation
    tmp.data <- ind.ivb.res[
        gene.id %in% deg.uni,
        ..lfc.cols
    ]
    corr.data <- rcorr(as.matrix(tmp.data), type='spearman')
    # Correlation plot-like plot
    tmp.file.name <- paste0(
        fig.pat.rels.path, '/',
        obj.extended.names[t.subset], '_CorrPlot_Dist-Spearman',
        '.pdf'
    )
    pdf(file=tmp.file.name, width=15, height=15)
    corrplot(
        corr=corr.data$r,
        col=colorRampPalette(c("blue","white","red"))(200),
        diag=FALSE,
        p.mat=corr.data$P, sig.level=0.01,
        order='hclust', hclust.method='average', addrect=clade.no[t.subset]
    )
    dev.off()
}

# ---> Relatedness according to linear dimentionality reduction based on gene expression or DEA stats.
# @ Function to perform PCA and output reports.
perform.pca <- function(
    pca.data, app.lab,
    tmp.reports.path
){
    # Define donor groups aesthetics.
    donor.group.cols <- c(
        '#A6CEE3', # pastel sky blue
        '#B2DF8A', # pastel green
        '#FDBF6F', # pastel orange
        '#CAB2D6', # pastel lavender
        '#FFD967', # pastel yellow
        '#FB9A99', # pastel coral
        '#80B1D3', # pastel steel blue
        '#B3DE69'  # pastel lime green
    )
    names(donor.group.cols) <- paste0('group.', 1:length(donor.group.cols))
    donor.group.shapes <- c(1, 2, 6, 4, 5)
    names(donor.group.shapes) <- paste0('group.', 1:length(donor.group.shapes))
    donor.group.aes <- donor.meta[, .(donor.id=donor.id.tag)]
    setorderv(x=donor.group.aes, cols='donor.id')
    donor.group.aes[, color.group:=rep(x=names(donor.group.cols), length.out=.N)]
    donor.group.aes[, shape.group:=rep(x=names(donor.group.shapes), each=length(donor.group.cols))]
    # Perform PCA
    pca.res <- prcomp(x=t(pca.data), scale.=TRUE)
    tmp.file.name <- paste0(
        tmp.reports.path, '/',
        obj.extended.names[t.subset], '_PCA_App-', app.lab, '_Var',
        '.pdf'
    )
    pdf(file=tmp.file.name)
    print(factoextra::fviz_eig(X=pca.res, addlabels=TRUE))
    dev.off()
    # Retrieve data to plot.
    pc.no <- 10
    if(ncol(pca.data)<pc.no) pc.no <- ncol(pca.data)
    pc.cols <- paste0('PC', 1:pc.no)
    exp.var <- head((pca.res$sdev**2)*100 / sum(pca.res$sdev**2), pc.no)
    names(exp.var) <- pc.cols
    exp.var <- round(x=exp.var, digits=2)
    to.plot <- as.data.table(pca.res$x[, pc.cols])
    to.plot$obs <- row.names(pca.res$x)
    if(app.lab=='DEA-Mean' | app.lab=='HVG-Mean'){
        to.plot <- separate(data=to.plot, col='obs', into=c('ag', 'donor.id'), sep='\\.')
        tmp.vals <- c(`Antigen`='ag', `Donor ID`='donor.id')
        # Add donor aesthetic groups.
        to.plot <- merge(x=to.plot, y=donor.group.aes, by='donor.id', all.x=TRUE, all.y=FALSE)
    }else{
        to.plot <- separate(data=to.plot, col='obs', into=c('lfc', 'ag'), sep='\\|')
        tmp.vals <- c(`Antigen`='ag')
    }
    # General PC plots
    tmp.file.name <- paste0(
        tmp.reports.path, '/',
        obj.extended.names[t.subset], '_PCA_App-', app.lab, '_PCs',
        '.pdf'
    )
    pdf(file=tmp.file.name, width=8)
    for(tmp.pc in pc.cols[2:length(pc.cols)]){
        for(tmp.val in names(tmp.vals)){
            tmp.var <- tmp.vals[tmp.val]
            tmp.aes <- if(tmp.var=='donor.id') aes(x=PC1, y=get(tmp.pc), color=color.group, shape=shape.group) else aes(x=PC1, y=get(tmp.pc), color=get(tmp.var), shape='mock')
            tmp.ggplot <- ggplot(data=to.plot, tmp.aes) +
                geom_vline(xintercept=0, linewidth=2, color='black') +
                geom_hline(yintercept=0, linewidth=2, color='black') +
                geom_point(size=3, stroke=2) +
                scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
                labs(
                    title=paste0('PC1 vs ', tmp.pc, '; ', tmp.val),
                    x=paste0('PC1 (', exp.var['PC1'], '%)'),
                    y=paste0(tmp.pc, ' (', exp.var[tmp.pc], '%)'),
                    color=tmp.val
                ) +
                theme_bw()
            tmp.ggplot <- if(tmp.var=='ag') tmp.ggplot + scale_color_manual(values=pp.cols) else tmp.ggplot + scale_color_manual(values=donor.group.cols)
            tmp.ggplot <- if(tmp.var=='donor.id') tmp.ggplot + scale_shape_manual(values=donor.group.shapes) else tmp.ggplot + scale_shape_manual(values=c(`mock`=1))
            print(tmp.ggplot)
        }
    }
    dev.off()
    # Publication-quality PC1 vs PC2 plots.
    for(tmp.val in names(tmp.vals)){
        tmp.var <- tmp.vals[tmp.val]
        tmp.aes <- if(tmp.var=='donor.id') aes(x=PC1, y=PC2, color=color.group, shape=shape.group) else aes(x=PC1, y=PC2, color=get(tmp.var), shape='mock')
        tmp.ggplot <- ggplot(data=to.plot, tmp.aes) +
            geom_vline(xintercept=0, linewidth=1.5, color='black') +
            geom_hline(yintercept=0, linewidth=1.5, color='black') +
            geom_point(size=6, stroke=4) +
            scale_y_continuous(breaks=scales::pretty_breaks(n=3)) +
            labs(
                x=paste0('PC1 (', exp.var['PC1'], '%)'),
                y=paste0('PC2 (', exp.var['PC2'], '%)'),
                color=tmp.val
            )
        tmp.ggplot <- if(tmp.var=='ag') tmp.ggplot + scale_color_manual(values=pp.cols) else tmp.ggplot + scale_color_manual(values=donor.group.cols)
        tmp.ggplot <- if(tmp.var=='donor.id') tmp.ggplot + scale_shape_manual(values=donor.group.shapes) else tmp.ggplot + scale_shape_manual(values=c(`mock`=1))
        tmp.lab <- paste0(obj.extended.names[t.subset], '_PCA_App-', app.lab, '_PCs-1-vs-2_Var-', tmp.val)
        publish.plot(
            tmp.ggplot=tmp.ggplot, output.path=tmp.reports.path, file.name=tmp.lab, type='pdf',
            blank.comp=blank.complement.4, do.legend=FALSE,
            width=10, height=10
        )
    }
    return(NULL)
}
# @ Define relevant specificities.
subset.rel.spcs <- list(
    `c40.lg.cd4.d40`=rel.spcs,
    `c40.lg.cd8.d40`=c('CMV', 'EBV', 'IAV', 'SARS-CoV-2')
)
subset.rel.pops <- list(
    `c40.lg.cd4.d40`=c('0', '2'),
    `c40.lg.cd8.d40`=c('0', '1')
)
# @ Process per strategy
for(t.subset in names(main.extended.names)){
    for(tmp.rel.pop in subset.rel.pops[[t.subset]]){
        # Retrieve DEA data.
        ind.ivb.res <- ivb.res[[t.subset]]
        # Define stat columns.
        tmp.data <- get.stat.cols(these.data=ind.ivb.res)
        p.cols <- tmp.data[[1]]; lfc.cols <- tmp.data[[2]]; ag.spcs <- tmp.data[[3]]
        # Define DEG universe.
        #       Genes that came up as a DEG in at least one of the comparisons.
        deg.sets <- lapply(X=ag.spcs, FUN=function(ag.spc){
            p.col <- p.cols[ag.spc]
            lfc.col <- lfc.cols[ag.spc]
            tmp.data <- ind.ivb.res[
                cluster==tmp.rel.pop &
                get(p.col)<=0.05 &
                # (get(lfc.col)>=0.25 | get(lfc.col)<=(-0.25)),
                (get(lfc.col)>=0 | get(lfc.col)<=(-0)),
                gene.id
            ]
            return(tmp.data)
        })
        deg.uni <- unique(unlist(deg.sets))
        # ---> Process based on DEA stats.
        # Pairwise distance and significance calculation
        tmp.data <- ind.ivb.res[
            cluster==tmp.rel.pop &
            gene.id %in% deg.uni,
            ..lfc.cols
        ]
        tmp.data <- as.matrix(tmp.data)
        tmp.check <- sum(is.na(tmp.data))/(dim(tmp.data)[1]*dim(tmp.data)[2])
        if(tmp.check>0.1) stop('Unexpected error!\n')
        tmp.data[is.na(tmp.data)] <- 0
        tmp.data <- tmp.data[rowSums(tmp.data)!=0, ]
        tmp.data <- tmp.data[, colSums(tmp.data)!=0]
        # PCA
        tmp.reports.path <- paste0(fig.pat.rels.path, '/dea_', obj.extended.names[t.subset], '_', tmp.rel.pop)
        if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
        perform.pca(pca.data=tmp.data, app.lab='DEA-LFC', tmp.reports.path=tmp.reports.path)
        # ---> Process based on gene expression stats of DEGs.
        # Define specificity/donor tag in seurat object.
        tmp.vals <- as.character(srt.objs.list[[t.subset]]@meta.data[, 'ag.spc.tag'])
        tmp.vals[!is.na(tmp.vals)] <- paste(
            tmp.vals[!is.na(tmp.vals)],
            srt.objs.list[[t.subset]]@meta.data[
                !is.na(tmp.vals), 'donor.id.tag'
            ],
            sep='.'
        )
        srt.objs.list[[t.subset]]@meta.data[, 'tmp.tag'] <- tmp.vals
        # Retrieve expression stats
        tmp.data <- get.exp.stats(
            seurat.obj=srt.objs.list[[t.subset]], norm.method='cpms',
            group.tag='tmp.tag', na.rm=TRUE,
            filter.tags=c(
                clust.labs[t.subset],
                'ag.spc.tag'
            ),
            groups.to.filter=list(
                c(tmp.rel.pop),
                subset.rel.spcs[[t.subset]]
            ),
            to.keep=TRUE,
            # min.cell.count=30, do.pos.only=TRUE
            min.cell.count=10, do.pos.only=TRUE
        )
        exp.stats <- tmp.data[['exp.stats']]
        cell.counts <- tmp.data[['group.counts']]
        srt.objs.list[[t.subset]]@meta.data[, 'tmp.tag'] <- NULL
        # PCA data.
        tmp.cols <- c('group', 'feature', 'mean')
        tmp.data <- exp.stats[
            feature %in% deg.uni,
            ..tmp.cols
        ]
        tmp.data <- spread(data=tmp.data, key=group, value=mean)
        tmp.row.names <- tmp.data$feature; tmp.data$feature <- NULL
        row.names(tmp.data) <- tmp.row.names
        tmp.data <- as.matrix(tmp.data)
        tmp.check <- sum(is.na(tmp.data))/(dim(tmp.data)[1]*dim(tmp.data)[2])
        if(tmp.check>0.01) stop('Unexpected error!\n')
        tmp.data <- tmp.data[rowSums(tmp.data)!=0, ]
        tmp.data <- tmp.data[, colSums(tmp.data)!=0]
        # PCA
        tmp.reports.path <- paste0(fig.pat.rels.path, '/dea_', obj.extended.names[t.subset], '_', tmp.rel.pop)
        if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
        perform.pca(pca.data=tmp.data, app.lab='DEA-Mean', tmp.reports.path=tmp.reports.path)
        # ---> Process based on gene expression stats of HVGs.
        # Try for different % of variance.
        tmp.vars <- c(10, 15, 20, 25, 30, 35, 40)
        # tmp.vars <- c(15)
        for(tmp.var in tmp.vars){
            # HVG definition.
            tmp.reports.path <- paste0(fig.pat.rels.path, '/hvg-', tmp.var, '_', obj.extended.names[t.subset], '_', tmp.rel.pop)
            if(!dir.exists(tmp.reports.path)) dir.create(tmp.reports.path)
            tmp.cells <- as.character(srt.objs.list[[t.subset]]@meta.data[, clust.labs[t.subset]])=='0'
            tmp.cells <- Cells(srt.objs.list[[t.subset]])[tmp.cells]
            tmp.obj <- subset(x=srt.objs.list[[t.subset]], cells=tmp.cells)
            tmp.obj <- do.hvg.analysis(
                seurat.obj=tmp.obj,
                FVFs.method='vst', feats.for.dsa=tmp.var,
                mean.cutoff=0.01, prop.cutoff=0.001,
                is.five.prime=FALSE,
                feat.select.path=tmp.reports.path
            )
            tmp.hvg <- VariableFeatures(tmp.obj)
            # PCA data.
            tmp.cols <- c('group', 'feature', 'mean')
            tmp.data <- exp.stats[
                feature %in% tmp.hvg,
                ..tmp.cols
            ]
            tmp.data <- spread(data=tmp.data, key=group, value=mean)
            tmp.row.names <- tmp.data$feature; tmp.data$feature <- NULL
            row.names(tmp.data) <- tmp.row.names
            tmp.data <- as.matrix(tmp.data)
            tmp.check <- sum(is.na(tmp.data))/(dim(tmp.data)[1]*dim(tmp.data)[2])
            if(tmp.check>0.01) stop('Unexpected error!\n')
            tmp.data <- tmp.data[rowSums(tmp.data)!=0, ]
            tmp.data <- tmp.data[, colSums(tmp.data)!=0]
            # PCA
            perform.pca(pca.data=tmp.data, app.lab='HVG-Mean', tmp.reports.path=tmp.reports.path)
        }
    }
}


### ------------------------ Systematic GSEA ------------------------ ###

# ---> General function for GSEA based on general DEA results.

plot.gsea.hmap.main <- function(
    dea.res, data.lab,
    module.feat.sets, set.lab,
    gen.module.set, spc.module.set,
    clust.rows=TRUE, clust.cols=TRUE,
    plot.blank=FALSE
){
    # Define stat columns.
    tmp.data <- get.stat.cols(these.data=dea.res)
    p.cols <- tmp.data[[1]]; lfc.cols <- tmp.data[[2]]; ag.spcs <- tmp.data[[3]]
    # Collect gene signatures
    tmp.module.set <- c(
        gen.module.set,
        spc.module.set
    )
    tmp.module.set <- module.feat.sets[tmp.module.set]
    names(tmp.module.set) <- c(
        names(gen.module.set),
        names(spc.module.set)
    )
    # GSEA stat fetch
    fgsea.results <- lapply(X=ag.spcs, FUN=function(ag.spc){
        cat(ag.spc, '\n')
        tmp.data <- dea.res[, get(lfc.cols[ag.spc])]
        names(tmp.data) <- dea.res[, gene.id]
        tmp.data <- tmp.data[!is.na(tmp.data)]
        fgsea.results <- fgsea(
            pathways=tmp.module.set,
            stats=tmp.data,
            minSize=10, maxSize=500,
            nproc=total.workers
            #, nperm=10000
        )
        return(fgsea.results)
    })
    names(fgsea.results) <- ag.spcs
    fgsea.results <- rbindlist(l=fgsea.results, use.names=TRUE, idcol='ag.spc')
    # Remove significance-wise uninformative signatures
    tmp.data <- fgsea.results[, .(to.check=sum(padj<=0.05)), by=pathway]
    tmp.data <- tmp.data[to.check==0, pathway]
    fgsea.results <- fgsea.results[!pathway %in% tmp.data]
    # Tidy data.
    tmp.data <- fgsea.results[, .(ref.tag=ag.spc, signature=pathway, p.adj=padj, NES)]
    tmp.data[p.adj>0.05, NES:=NA]
    tmp.data[, p.adj:=NULL]
    tmp.data <- spread(data=tmp.data, key=ref.tag, value=NES)
    tmp.data <- as.data.frame(tmp.data, stringsAsFactors=FALSE)
    row.names(tmp.data) <- tmp.data$signature; tmp.data$signature <- NULL
    # Set color scale and breaks for heatmap.
    col.breaks <- seq(from=min(range(tmp.data, na.rm=TRUE)), to=max(range(tmp.data, na.rm=TRUE)), length.out=100)
    mid.point <- which.min(abs(col.breaks - 0))
    # hmap.col.scale.1 <- colorRampPalette(c('#233773', '#3553ae', '#72bcd4', '#ffffff'))(mid.point)
    # hmap.col.scale.2 <- colorRampPalette(c('#ffffff', '#ffff00', '#ffa500', '#ff664d', '#ff2500', '#b31a00'))(100-(mid.point+1))
    hmap.col.scale.1 <- colorRampPalette(c('#8FD4C8', '#9FDAD0', '#ABE0D6', '#BFE8E0', '#CEEBE5', '#DEF2EE', '#EAF8F5', '#FFFFFF'))(mid.point)
    hmap.col.scale.2 <- colorRampPalette(c('#FFFFFF', '#F6EDF4', '#ECDEEC', '#E2C7E2', '#DBBADC', '#CFA5CF', '#C58FC6', '#BC7FBE'))(100-(mid.point+1))
    hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
    # Perform clustering.
    if(clust.rows | clust.cols){
        clust.data <- fgsea.results[, .(ref.tag=ag.spc, signature=pathway, NES)]
        clust.data <- spread(data=clust.data, key=ref.tag, value=NES)
        clust.data <- as.data.frame(clust.data, stringsAsFactors=FALSE)
        row.names(clust.data) <- clust.data$signature; clust.data$signature <- NULL
        if(clust.rows){
            clust.row.obj <- dist(x=clust.data, method='euclidean')
            clust.row.obj <- hclust(d=clust.row.obj, method='complete')
        }else{
            clust.row.obj <- clust.rows
            tmp.module.set <- c(
                gen.module.set,
                spc.module.set
            )
            tmp.module.set <- names(tmp.module.set)[names(tmp.module.set) %in% row.names(tmp.data)]
            tmp.data <- tmp.data[tmp.module.set, ]
        }
        if(clust.cols){
            clust.col.obj <- dist(x=t(clust.data), method='euclidean')
            clust.col.obj <- hclust(d=clust.col.obj, method='complete')
        }else{
            clust.col.obj <- clust.cols
        }
    }
    # Set columns metadata.
    col.meta <- data.frame(
        row.names=colnames(tmp.data),
        `Specificity`=colnames(tmp.data),
        stringsAsFactors=FALSE
    )
    meta.cols <- list(
        `Specificity`=pp.cols[names(pp.cols) %in% colnames(tmp.data)]
    )
    # @ Heatmap.
    # Plot 'complete' version
    tmp.width <- (ncol(tmp.data)*1) + 0.6
    tmp.height <- (nrow(tmp.data)*1) + 0.8
    tmp.file.name <- paste0(
        fig.pat.rels.path, '/',
        data.lab, '_SystematicGSEA_Set-', set.lab,
        '.C.pdf'
    )
    pheatmap(
        mat=tmp.data,
        color=hmap.col.scale, breaks=col.breaks, scale='none',
        na_col='#b3b3b3', border_color='black',
        cluster_rows=clust.row.obj, cluster_cols=clust.col.obj,
        annotation_col=col.meta, annotation_names_col=FALSE, annotation_colors=meta.cols,
        show_rownames=TRUE, show_colnames=TRUE, legend=TRUE, annotation_legend=TRUE,
        filename=tmp.file.name, height=tmp.height, width=tmp.width+1.5
    )
    if(plot.blank){
        # Plot 'complete' version
        tmp.file.name <- paste0(
            fig.pat.rels.path, '/',
            data.lab, '_SystematicGSEA_Set-', set.lab,
            '.B.pdf'
        )
        pheatmap(
            mat=tmp.data,
            color=hmap.col.scale, breaks=col.breaks, scale='none',
            na_col='#b3b3b3', border_color='black',
            cluster_rows=clust.row.obj, cluster_cols=clust.col.obj,
            annotation_col=col.meta, annotation_names_col=FALSE, annotation_colors=meta.cols,
            show_rownames=FALSE, show_colnames=FALSE, legend=FALSE, annotation_legend=FALSE,
            filename=tmp.file.name, height=tmp.height, width=tmp.width
        )
    }
    # @ Color legend only.
    tmp.data <- as.data.frame(tmp.data)
    tmp.data$signature <- row.names(tmp.data)
    tmp.data <- gather(data=tmp.data, key='pop', value='nes', -`signature`)
    tmp.data <- tmp.data[!is.na(tmp.data$nes), ]
    # W/ labels
    tmp.ggplot <- ggplot(data=tmp.data, aes(x=nes, y=nes, col=nes)) + 
        geom_point() +
        scale_color_gradientn(colors=hmap.col.scale, breaks=scales::pretty_breaks(n=5)) +
        theme(
            legend.ticks=element_line(color='black', linewidth=0.6),
            legend.ticks.length=unit(0.22, "cm"),
            legend.frame=element_rect(color='black', linewidth=0.6)
        )
    tmp.ggplot <- get_legend(p=tmp.ggplot)
    tmp.file.name <- paste0(
        fig.pat.rels.path, '/',
        data.lab, '_SystematicGSEA_Set-', set.lab,
        '.L1.pdf'
    )
    pdf(file=tmp.file.name, height=2, width=2)
    print(as_ggplot(tmp.ggplot))
    dev.off()
    # W/out labels
    tmp.ggplot <- ggplot(data=as.data.frame(tmp.data), aes(x=nes, y=nes, col=nes)) + 
        geom_point() +
        scale_color_gradientn(colors=hmap.col.scale, breaks=scales::pretty_breaks(n=5), name=NULL, labels=NULL) +
        theme(
            legend.ticks=element_line(color='black', linewidth=0.6),
            legend.ticks.length=unit(0.22, "cm"),
            legend.frame=element_rect(color='black', linewidth=0.6)
        )
    tmp.ggplot <- get_legend(p=tmp.ggplot)
    tmp.file.name <- paste0(
        fig.pat.rels.path, '/',
        data.lab, '_SystematicGSEA_Set-', set.lab,
        '.L2.pdf'
    )
    pdf(file=tmp.file.name, height=1.25, width=0.3)
    print(as_ggplot(tmp.ggplot))
    dev.off()
    return(NULL)
}


# ----> Define sets of signatures.
# MSIGDB set.
msigdb.set <- msigdb.data[, unique(name)]
names(msigdb.set) <- msigdb.set
msigdb.module.feat <- lapply(X=msigdb.set, FUN=function(x) msigdb.data[name==x, unique(gene_symbol)])
dice.set <- c(
    names(modules.feats)[names(modules.feats) %like% 'RESTING'],
    names(modules.feats)[names(modules.feats) %like% 'Act-']
)
names(dice.set) <- dice.set
# Lists of sets
gen.module.set.l <- list(
    `PubComp`=c(
        `Interferon`='type1n2.interferon.broadreact',
        `Glycolysis`='hallmark.glycolysis',
        `Cell cycle`='cell.cycling.best',
        `TRM (Clarke)`='trm.signature.clarke',
        `TRM (Scott)`='trm.up.scott',
        `Non-TRM (Scott)`='trm.down.scott',
        `TCM (Pauken)`='TCM.Pauken.Tumor'
    ),
    `MSigDB`=msigdb.set,
    `DICE`=dice.set
)
spc.module.set.l <- list(
    `PubComp`=list(
        `c40.lg.cd4.d40`=c(
            `Cytotoxicity (Patil)`='cell.cytotoxicity.patil',
            `TFH (Locci)`='tfh.signature.locci',
            `TH1 (Arlehamn)`='th1.signature1.arlehamn',
            `TH17 (Arlehamn)`='th17.signature2.arlehamndan',
            `TREG (Schmiedel)`='treg.signature.schmiedle',
            `Exhaustion`='dixhaust.consensus2.vj',
            `Unhelped`='unhelped.cullen'
        ),
        `c40.lg.cd8.d40`=c(
            `Exhaustion`='dixhaust.consensus2.vj',
            `Unhelped`='unhelped.cullen',
            `Cytotoxicity (Guo)`='tcell.cytotoxicity.guo'
        )
    ),
    `MSigDB`=list(),
    `DICE`=list()
)



### ------------------------ Systematic GSEA ------------------------ ###
### --------------- for background-vs-index approach  --------------- ###

# ---> Function to create GSEA heatmap for a set of signatures.
plot.gsea.hmap.1 <- function(
    module.feat.sets, set.lab,
    gen.module.set, spc.module.sets,
    clust.rows=TRUE, clust.cols=TRUE
){
    for(t.subset in names(main.extended.names)[1]){
        # Retrieve DEA data.
        tmp.res <- ivb.res[[t.subset]]
        # Run GSEA
        plot.gsea.hmap.main(
            dea.res=tmp.res,
            data.lab=paste0(
                'IvBApp_',
                obj.extended.names[t.subset]
            ),
            module.feat.sets=module.feat.sets, set.lab=set.lab,
            gen.module.set=gen.module.set, spc.module.set=spc.module.sets[[t.subset]],
            clust.rows=clust.rows, clust.cols=clust.cols
        )
    }
    return(NULL)
}

# ---> Systematic GSEA, process for each signature collection.
for(tmp.set.lab in names(gen.module.set.l)){
    cat('Process for', tmp.set.lab, '\n')
    plot.gsea.hmap.1(
        module.feat.sets=c(modules.feats, msigdb.module.feat), set.lab=tmp.set.lab,
        gen.module.set=gen.module.set.l[[tmp.set.lab]], spc.module.sets=spc.module.set.l[[tmp.set.lab]]
    )
}

# ---> Consensus set 1 for GSEA.
gen.module.set <- c(
)
spc.module.sets <- list(
    `c40.lg.cd4.d40`=c(
        `Hypoxia (Hall.)`='HYPOXIA',
        `Apoptosis (Hall.)`='APOPTOSIS',
        `IL2/STAT5 signaling (Hall.)`='IL2_STAT5_SIGNALING',
        `Inflammatory response (Hall.)`='INFLAMMATORY_RESPONSE',
        `TNFA/NFKB signaling (Hall.)`='TNFA_SIGNALING_VIA_NFKB',
        `OXPHOS (Hall.)`='OXIDATIVE_PHOSPHORYLATION',
        `Epithelial/Mesenchimal transition (Hall.)`='EPITHELIAL_MESENCHYMAL_TRANSITION',
        `mTORC1 signaling (Hall.)`='MTORC1_SIGNALING',
        # `Apical junction (Hall.)`='APICAL_JUNCTION',
        `UPR (Hall.)`='UNFOLDED_PROTEIN_RESPONSE',
        `IFN-I response (Hall.)`='INTERFERON_ALPHA_RESPONSE',
        `IFN-II response (Hall.)`='INTERFERON_GAMMA_RESPONSE',
        `TH1 (Arlehamn)`='th1.signature1.arlehamn',
        `TH17 (Arlehamn)`='th17.signature2.arlehamndan',
        `Resting TH17 (DICE)`='R24.TH17.RESTING',
        `Resting TH1 (DICE)`='R24.TH1.RESTING',
        `Resting TH2 (DICE)`='R24.TH2.RESTING',
        `Activated TH17 (DICE)`='Act−3−TH17',
        `Activated TH1 (DICE)`='Act−2−TH1',
        `Activated TH2 (DICE)`='Act−5−TH2',
        `Activated THIFNR (DICE)`='Act−17−THIFNR'
    )#,
    # `c40.lg.cd8.d40`=c(
    #     `Non-TRM (Scott)`='trm.down.scott',
    #     `Cytotoxicity (Guo)`='tcell.cytotoxicity.guo',
    #     `Interferon`='type1n2.interferon.broadreact',
    #     `Unhelped`='unhelped.cullen',
    #     `Hypoxia (Hall.)`='HYPOXIA',
    #     `IL2/JAK/STAT5 signaling (Hall.)`='IL2_STAT5_SIGNALING',
    #     `Exhaustion`='dixhaust.consensus2.vj',
    #     `TRM (Scott)`='trm.up.scott',
    #     `TRM (Clarke)`='trm.signature.clarke'
    # )
)
tmp.data.1 <- msigdb.data[, unique(name)]
tmp.data <- lapply(X=tmp.data.1, FUN=function(x) msigdb.data[name==x, unique(gene_symbol)])
names(tmp.data) <- tmp.data.1
tmp.data <- c(
    modules.feats,
    tmp.data
)
plot.gsea.hmap.1(
    module.feat.sets=tmp.data, set.lab='Consens-1',
    gen.module.set=gen.module.set, spc.module.sets=spc.module.sets,
    clust.rows=TRUE
)

# ---> Consensus set 2 for GSEA.
gen.module.set <- c(
)
spc.module.sets <- list(
    `c40.lg.cd4.d40`=c(
        `Hypoxia (Hall.)`='HYPOXIA',
        `Apoptosis (Hall.)`='APOPTOSIS',
        `IL2/STAT5 signaling (Hall.)`='IL2_STAT5_SIGNALING',
        `Inflammatory response (Hall.)`='INFLAMMATORY_RESPONSE',
        `TNFA/NFKB signaling (Hall.)`='TNFA_SIGNALING_VIA_NFKB',
        `Epithelial/Mesenchimal transition (Hall.)`='EPITHELIAL_MESENCHYMAL_TRANSITION',
        `mTORC1 signaling (Hall.)`='MTORC1_SIGNALING',
        `UPR (Hall.)`='UNFOLDED_PROTEIN_RESPONSE',
        `IFN-I response (Hall.)`='INTERFERON_ALPHA_RESPONSE',
        `IFN-II response (Hall.)`='INTERFERON_GAMMA_RESPONSE'
    )#,
    # `c40.lg.cd8.d40`=c(
    #     `Non-TRM (Scott)`='trm.down.scott',
    #     `Cytotoxicity (Guo)`='tcell.cytotoxicity.guo',
    #     `Interferon`='type1n2.interferon.broadreact',
    #     `Unhelped`='unhelped.cullen',
    #     `Hypoxia (Hall.)`='HYPOXIA',
    #     `IL2/JAK/STAT5 signaling (Hall.)`='IL2_STAT5_SIGNALING',
    #     `Exhaustion`='dixhaust.consensus2.vj',
    #     `TRM (Scott)`='trm.up.scott',
    #     `TRM (Clarke)`='trm.signature.clarke'
    # )
)
tmp.data.1 <- msigdb.data[, unique(name)]
tmp.data <- lapply(X=tmp.data.1, FUN=function(x) msigdb.data[name==x, unique(gene_symbol)])
names(tmp.data) <- tmp.data.1
tmp.data <- c(
    modules.feats,
    msigdb.module.feat
)
plot.gsea.hmap.1(
    module.feat.sets=tmp.data, set.lab='Consens-2',
    gen.module.set=gen.module.set, spc.module.sets=spc.module.sets,
    clust.rows=TRUE
)


### ------------------------ Systematic GSEA ------------------------ ###
### --------------- for pairwise comparison approach  --------------- ###

# ---> Function to create GSEA heatmap for a set of signatures.

plot.gsea.hmap.2 <- function(
    dea.res, spc.of.int, data.subset.lab,
    filter.tags, groups.to.filter, to.keep=TRUE,
    module.feat.sets, set.lab,
    gen.module.set, spc.module.sets,
    clust.rows=TRUE, clust.cols=TRUE,
    plot.blank=FALSE
){
    # Preflights.
    if(!is.list(groups.to.filter)) groups.to.filter <- list(groups.to.filter)
    tmp.check <- length(filter.tags) == length(groups.to.filter)
    if(!tmp.check) stop('Faulty filtering inputs. `filter.tags` must be a character vector listing tags defined in metadata and `groups.to.filter` must be a list with the same number of entries as the length of `filter.tags`.\n')
    tmp.check <- filter.tags %in% colnames(dea.res)
    if(!all(tmp.check)){
        tmp.check <- filter.tags[!tmp.check]
        tmp.check <- paste0('Following tag(s) not properly defined in seurat object\'s metadata:\n', paste0(tmp.check, collapse='\n'), '\n')
        stop(tmp.check)
    }
    if(!'consensus.pred' %in% filter.tags){
        filter.tags <- c(filter.tags, 'consensus.pred')
        groups.to.filter <- c(
            groups.to.filter,
            list(spc.of.int)
        )
    }
    # Define rows to be maintained for each group.
    idxs.to.keep <- lapply(X=1:length(filter.tags), FUN=function(tag.idx){
        this.tmp.tag <- filter.tags[tag.idx]
        tmp.groups <- as.character(groups.to.filter[[tag.idx]])
        tmp.data <- dea.res[, .(idx=1:.N, this.tmp.tag=as.character(get(this.tmp.tag)))]
        if(to.keep){
            tmp.data <- tmp.data[this.tmp.tag %in% tmp.groups, idx]
        }else{
            tmp.data <- tmp.data[!this.tmp.tag %in% tmp.groups, idx]
        }
        return(tmp.data)
    })
    idxs.to.keep <- Reduce(x=idxs.to.keep, f=intersect)
    # Filter DEA data.
    dea.res <- dea.res[idxs.to.keep]
    tmp.cols <- colnames(dea.res)[!colnames(dea.res) %like% paste0(spc.of.int, '$')]
    dea.res <- dea.res[, ..tmp.cols]
    # Run GSEA
    plot.gsea.hmap.main(
        dea.res=dea.res,
        data.lab=paste0(
            'PWApp_', data.subset.lab,
            obj.extended.names[t.subset]
        ),
        module.feat.sets=module.feat.sets, set.lab=set.lab,
        gen.module.set=gen.module.set, spc.module.set=spc.module.sets[[t.subset]],
        clust.rows=clust.rows, clust.cols=clust.cols,
        plot.blank=plot.blank
    )
    return(NULL)
}

# ------------> Process for CD8 TEM, CMV
main.obj.name <- 'c40.lg.cd8.d40'
spc.of.int <- 'CMV'
pop.of.int <- '1-TEM'
gen.filter.tags <- c('pop.tag')
gen.groups.to.filter <- list(
    pop.of.int
)
gen.data.subset.lab <- paste0(spc.of.int, '_', pop.of.int)

# ---> Systematic GSEA, process for each signature collection.
t.subset <- main.obj.name
for(tmp.set.lab in names(gen.module.set.l)){
    cat('Process for', tmp.set.lab, '\n')
    plot.gsea.hmap.2(
        dea.res=apc.res[[main.obj.name]], spc.of.int=spc.of.int, data.subset.lab=gen.data.subset.lab,
        filter.tags=gen.filter.tags, groups.to.filter=gen.groups.to.filter,
        module.feat.sets=c(modules.feats, msigdb.module.feat), set.lab=tmp.set.lab,
        gen.module.set=gen.module.set.l[[tmp.set.lab]], spc.module.sets=spc.module.set.l[[tmp.set.lab]]
    )
}

# ---> Consensus set.
tmp.signatures <- c(
    'IL2_STAT5_SIGNALING', 'APOPTOSIS',
    'UV_RESPONSE_UP', 'INFLAMMATORY_RESPONSE', 'KRAS_SIGNALING_DN', 'TNFA_SIGNALING_VIA_NFKB', 'INTERFERON_ALPHA_RESPONSE', 'INTERFERON_GAMMA_RESPONSE',
    'OXIDATIVE_PHOSPHORYLATION'
)
names(tmp.signatures) <- tmp.signatures
plot.gsea.hmap.2(
    dea.res=apc.res[[main.obj.name]], spc.of.int=spc.of.int, data.subset.lab=gen.data.subset.lab,
    filter.tags=gen.filter.tags, groups.to.filter=gen.groups.to.filter,
    module.feat.sets=c(modules.feats, msigdb.module.feat), set.lab='Consens-1',
    gen.module.set=tmp.signatures, spc.module.sets=list(),
    plot.blank=TRUE
)


# ------------> Process for CD4-CTL, CMV
main.obj.name <- 'c40.lg.cd4.d40'
spc.of.int <- 'CMV'
pop.of.int <- '2-CTL'
gen.filter.tags <- c('pop.tag')
gen.groups.to.filter <- list(
    pop.of.int
)
gen.data.subset.lab <- paste0(spc.of.int, '_', pop.of.int)

# ---> Systematic GSEA, process for each signature collection.
t.subset <- main.obj.name
for(tmp.set.lab in names(gen.module.set.l)){
    cat('Process for', tmp.set.lab, '\n')
    plot.gsea.hmap.2(
        dea.res=apc.res[[main.obj.name]], spc.of.int=spc.of.int, data.subset.lab=gen.data.subset.lab,
        filter.tags=gen.filter.tags, groups.to.filter=gen.groups.to.filter,
        module.feat.sets=c(modules.feats, msigdb.module.feat), set.lab=tmp.set.lab,
        gen.module.set=gen.module.set.l[[tmp.set.lab]], spc.module.sets=spc.module.set.l[[tmp.set.lab]]
    )
}

# ---> Consensus set.
tmp.signatures <- c(
    'INTERFERON_ALPHA_RESPONSE', 'INTERFERON_GAMMA_RESPONSE', 'PEROXISOME', 'OXIDATIVE_PHOSPHORYLATION',
    'HYPOXIA', 'TNFA_SIGNALING_VIA_NFKB', 'APOPTOSIS'
    
)
names(tmp.signatures) <- tmp.signatures
plot.gsea.hmap.2(
    dea.res=apc.res[[main.obj.name]], spc.of.int=spc.of.int, data.subset.lab=gen.data.subset.lab,
    filter.tags=gen.filter.tags, groups.to.filter=gen.groups.to.filter,
    module.feat.sets=c(modules.feats, msigdb.module.feat), set.lab='Consens-1',
    gen.module.set=tmp.signatures, spc.module.sets=list(),
    plot.blank=TRUE
)


### --------------------- Relevant gene examples -------------------- ###

tmp.data <- meta.data[
    gex.lin.tag=='CD4' &
    `RNA_snn_res.0.2`=='2' &
    !is.na(ag.spc.tag),
    .(cell.count=.N),
    by=.(ag.spc.tag, donor.id.tag)
]
tmp.data <- spread(data=tmp.data, key=ag.spc.tag, value=cell.count, fill=0)
tmp.file.name <- paste0(fig.pat.rels.path, '/CD4-CTL_CellCounts.csv')
fwrite(file=tmp.file.name, x=tmp.data)

tmp.data <- meta.data[
    gex.lin.tag=='CD4' &
    `RNA_snn_res.0.2`=='0' &
    !is.na(ag.spc.tag),
    .(cell.count=.N),
    by=.(ag.spc.tag, donor.id.tag)
]
tmp.data <- spread(data=tmp.data, key=ag.spc.tag, value=cell.count, fill=0)
tmp.file.name <- paste0(fig.pat.rels.path, '/CD4-TRM_CellCounts.csv')
fwrite(file=tmp.file.name, x=tmp.data)

t.subset <- 'c40.lg.cd4.d40'
this.cluster <- '2'
these.markers <- c(
    'GNLY', 'GZMB',
    'KLRD1', 'BHLHE40', 'IFNGR1',
    'CD7', 'SELL'
    # 'GZMK',
)
# these.groups <- rel.spcs
these.groups <- c('Aspergillus', 'CMV', 'EBV', 'IAV', 'SARS-CoV-2', 'MPV')
# Subset seurat object.
tmp.obj <- subset.srt.obj(
    seurat.obj=srt.objs.list[[t.subset]],
    filter.tags=c(
        'ag.spc.tag',
        clust.labs[t.subset]
    ),
    groups.to.filter=list(
        these.groups,
        this.cluster
    ),
    to.keep=TRUE#,
    # na.rm=TRUE
)
# Plot
tmp.ggplot <- dot.plot(
    seurat.obj=tmp.obj, features=these.markers, slot='data', ensembl=FALSE,
    do.norm=TRUE,
    groups.tag='ag.spc.tag', groups.order=these.groups, groups.of.int=NULL,
    na.rm=TRUE, feature.thold=NULL,
    filter.tags=NULL, groups.to.filter=NULL, keep=TRUE,
    # this.color.scale=signatures.col.scale,
    this.color.scale='hot.and.cold',
    # col.min=col.min, col.max=col.max,
    scale.by='radius', dot.scale=12,
    size.min=0, #size.max=size.max,
    file.name=NULL
)
tmp.ggplot <- tmp.ggplot + theme(legend.position='bottom')
# Output.
width.add <- 2
tmp.width <- (length(these.groups) * 0.5) + width.add
tmp.height <- (length(these.markers) * 0.4) + 0.3
tmp.lab <- paste0(obj.extended.names[t.subset], '_Markers_DotPlot_Pop-', this.cluster)
publish.plot(
    tmp.ggplot=tmp.ggplot, output.path=fig.pat.rels.path, file.name=tmp.lab, type='pdf', blank.comp=blank.complement.2, do.legend=TRUE,
    width=tmp.width, height=tmp.height, legend.height=2.6, legend.width=0.5,
    do.rotate=TRUE
)


############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############
############    ----------    End of script    ----------    ############
############    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ############