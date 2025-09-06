cat('\n\n')
############    ------   Prolif. assay TCR data    -------    ############
############    -----------   QC and prep.    ------------    ############
cat('############    ------   Prolif. assay TCR data    -------    ############\n')
cat('############    -----------   QC and prep.    ------------    ############\n')

# By: Vicente Fajardo
# Best module: R/4.2.2


### -------------------------- Description -------------------------- ###
# Script to apply post-preprocessing QC to the bulk TCR data from the proliferation assay experiments.


cat('\n\n')
### -------------------------- Dependencies ------------------------- ###
cat('### -------------------------- Dependencies ------------------------- ###\n')
library(optparse)
library(parallel)
library(data.table)
library(tidyr)
library(stringr)
library(ggplot2)
library(pheatmap)
library(UpSetR)
source('/path/to/user/R24/paper_developments/R24_Cancer/paper_items/jobs_scripts/functions_collection.0.4.R')


cat('\n\n')
### --------------------------- Arguments --------------------------- ###
cat('### --------------------------- Arguments --------------------------- ###\n')

# ---> Arguments from command line
option.list <- list(
    # For output.
    make_option(opt_str="--ArtifactThold", type="numeric", default=0.5, dest="art.call.thold", help="Numeric, sample fraction threshold to call artifacts within a plate.\n")
)
# Getting arguments from command line and setting their values to their respective variable names.
opt.parser <- OptionParser(option_list=option.list)
opt <- parse_args(opt.parser)
# Moving options to their own variables
art.call.thold <- opt$art.call.thold

# ---> Hardcoded arguments
# @ Paper.
main.prj <- 'HIPC'
this.prj <- 'HIPC-Lung'
this.prj.lab <- 'prol_assay_reps'
batches.lab <- 'Batches-1-to-4'
gex.batches.lab <- 'Batches-1-to-9'
aggr.sffx <- batches.lab
aggr.id <- paste0(this.prj, '_', aggr.sffx)
donor.vers <- '1.0'
# Analysis-specific arguments
seq.dates <- c('06-07-2024', '04-08-2024', '05-14-2024')
seq.run.ids <- c('NV137', 'NV132', 'NV135')
names(seq.dates) <- seq.run.ids
sample.prjs <- c('BULK_TCR23', 'BULK_TCR24a', 'BULK_TCR24b', 'CTY2_BULKTCR25', 'CTY2_BULKTCR26')
# @ Reports date.
reports.date <- NULL
# @ CPU number
tmp.cmmd <- "echo $SLURM_JOB_CPUS_PER_NODE"
cpu.no <- as.integer(system(command=tmp.cmmd, intern=TRUE))
# @ Seed.
set.seed(seed=1)
# ---> Color defintions.
# @ Peptide pool colors.
chain.cols <- c(
  TCRA='#800080',
  TCRB='#008000'
)
pp.cols <- c(
    CMV4='#B22222',
    CMV8='#B22222',
    EBV='#FF8C00',
    IAV='#0095B3',
    `SARS-CoV-2`='#9932CC',
    MPV='#FF5EA2',
    PIV='#FFD700',
    RSV='#59B300',
    RV='#000000',
    HIV='#720000',
    `B-pertussis`='#754600',
    `B-pertussis-Vax`='#754600',
    Aspergillus='#005746',
    Alternaria='#319272'
)
ag.cols <- c(
    `CMV`='#B22222',
    `EBV`='#FF8C00',
    `IAV`='#0095B3',
    `SARS-CoV-2`='#9932CC',
    `MPV`='#FF5EA2',
    `PIV`='#FFD700',
    `RSV`='#59B300',
    `B-pertussis-Vax`='#754600',
    `B-pertussis-Rest`='#C37600',
    `Aspergillus`='#005746',
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
cell.type.cols <- c(
  CD4='#00BFFF',
  CD8='#EE82EE'
)
rep.cols <- c(
  `1`='#000000',
  `2`='#808080',
  `3`='#d3d3d3',
  `4`='#ffffff'
)
plate.cols <- c(
  `NV132-P1`='#00606E',
  `NV132-P2`='#01A3BA',
  `NV135-P1`='#070926',
  `NV135-P2`='#281259',
  `NV137-P1`='#6E3101',
  `NV137-P2`='#FF7F1B',
  `NV137-P3`='#BA5201'
)
# ---> Others for plotting
# blank.complement.1 <- theme(line=element_blank(), text=element_blank(), panel.background=element_blank(), panel.border=element_blank(), legend.position='none', axis.ticks=element_blank(), axis.line=element_blank()) #
# blank.complement.1.1 <- theme(line=element_blank(), text=element_blank(), panel.background=element_blank(), panel.border=element_blank(), legend.position='none', axis.ticks=element_line(size=1), axis.line=element_line(size=1), axis.ticks.length=unit(0.4, "cm"), axis.ticks.x=element_blank())
# # Alternative to 1 with no ticks in the x axis.
# blank.complement.1.2 <- theme(line=element_blank(), text=element_blank(), panel.background=element_blank(), panel.border=element_blank(), legend.position='none', axis.ticks=element_line(size=1), axis.line=element_line(size=1), axis.ticks.length=unit(0.4, "cm"))
# ---> Path definitions.
prj.gen.path <- paste0('/path/to/user/', main.prj, '/paper_developments/', this.prj)
donor.data.path <- paste0(prj.gen.path, '/donors_metadata')
data.path <- paste0(prj.gen.path, '/exploratory_analyses/', this.prj.lab, '/data_to_start')
reports.path <- paste0(prj.gen.path, '/exploratory_analyses/', this.prj.lab)
if(is.null(reports.date)){
  reports.path <- paste0(reports.path, '/prol_assay_reps_art-thold-', art.call.thold, '_', Sys.Date())
  if(!dir.exists(reports.path)) dir.create(reports.path)
}else{
  reports.path <- paste0(reports.path, '/prol_assay_reps_art-thold-', art.call.thold, '_', reports.date)
  if(!dir.exists(reports.path)) stop('A reports date was provided, but we failed to define an already existing reports path.\n')
}
gen.seq.path <- '/path/to/user/sequencing_data'
exp.layout.paths <- paste0(gen.seq.path, '/', seq.dates, '/bcl2fastq/data/experiment_layout')
migec.paths <- paste0(gen.seq.path, '/', seq.dates, '/migec')
cdr.filter.paths <- paste0(migec.paths, '/cdrfinal')
names(cdr.filter.paths) <- names(seq.dates)
# ---> File definitions.
# Barcodes file.
bc.info.files <- paste0(exp.layout.paths, '/', names(seq.dates), '_barcode_info_for_migec.csv')
names(bc.info.files) <- names(seq.dates)
# Mislabel fix file.
name.fix.file <- paste0(data.path, '/NameFixInfo.csv')
# Ex vivo lung data.
lung.data.file <- '/path/to/user/HIPC/paper_developments/HIPC-Lung/preliminary_process/Batches-1-to-8_2024-11-06/ClonotypeBasedTCRData.csv'
# ---> Check directories and files.
if(!dir.exists(reports.path)) dir.create(reports.path)
essential.files <- list(
  `Barcode info files`=bc.info.files
)
essential.paths <- c(
  `CDR filter paths`=cdr.filter.paths
)
essential.files <- essential.files[!unlist(lapply(X=essential.files, FUN=function(x) all(file.exists(x))))]
essential.paths <- essential.paths[!unlist(lapply(X=essential.paths, FUN=function(x) all(dir.exists(x))))]
if(length(essential.files) > 0) stop(paste0('Next essential files are not appropriately defined -make sure you\'ve got their names right-:\n', paste0(names(essential.files), collapse='\n'), '\n'))
if(length(essential.paths) > 0) stop(paste0('Next essential paths are not appropriately defined -make sure you\'ve got their names right-:\n', paste0(names(essential.paths), collapse='\n'), '\n'))


### --------------------------- Functions --------------------------- ###

# ---> Function to obtain jaccard index.

get.jaccard <- function(x, y){
  num.val <- length(intersect(x=x, y=y))
  den.val <- length(union(x=x, y=y))
  return(num.val/den.val)
}


cat('\n\n')
### --------------------------- Load data --------------------------- ###
### ---------------------- & data preprocessing --------------------- ###
cat('### --------------------------- Load data --------------------------- ###\n')
cat('### ---------------------- & data preprocessing --------------------- ###\n')

# ---> Mislabel fix info
name.fix.info <- fread(file=name.fix.file)

# ---> Barcode info files.
bc.info <- lapply(X=bc.info.files, FUN=fread)
bc.info <- rbindlist(l=bc.info, use.names=TRUE, idcol='seq.run')

# ---> Lung TCR data.
lung.data <- fread(file=lung.data.file)


cat('\n\n')
### ---------------------- Data preprocessing ----------------------- ###
### ----------------------- of bulk TCR data ------------------------ ###
cat('### ---------------------- Data preprocessing ----------------------- ###\n')
cat('### ----------------------- of bulk TCR data ------------------------ ###\n')

prepro.path <- paste0(reports.path, '/for_manual_check')
if(!dir.exists(prepro.path)) dir.create(prepro.path)

# ---> Barcode info files.
# Keep info only from relevant projects.
bc.info <- bc.info[sample.prj %in% sample.prjs, ]

# ---> TCR chain type.
bc.info[, tcr.chain:=str_extract(string=seq.lib.name, pattern='ACR|BCR')]
tmp.vals <- c(`ACR`='TCRA', `BCR`='TCRB')
tmp.check <- bc.info[, all(tcr.chain%in%names(tmp.vals))]
if(!tmp.check) stop('Failed to define type of chain.\n')
bc.info[, tcr.chain:=tmp.vals[tcr.chain]]

# ---> Mislabel fix.
# Add fix info to full barcode info table.
merge.cols <- c(
    'seq.run', 'sample.prj', 'tcr.lib.id', 'tcr.lib.name', 'tcr.chain',
    'master.id', 'slave.id', 'master.seq', 'slave.seq'
)
tmp.data <- copy(name.fix.info)
tmp.data[, tcr.lib.id:=as.character(tcr.lib.id)]
tmp.data[, tcr.lib.name:=orig.tcr.lib.name]; tmp.data[, orig.tcr.lib.name:=NULL]
tmp.data[, final.lib.name:=as.character(fix.tcr.lib.name)]; tmp.data[, fix.tcr.lib.name:=NULL]
bc.info <- merge(
    x=bc.info, y=tmp.data,
    by=merge.cols,
    all.x=TRUE, all.y=TRUE
)
tmp.check <- bc.info[!is.na(accordance.record) & is.na(seq.lib.name), .N==0]
if(!tmp.check) stop('Info for mislabel fix does not match up original barcode information.\n')
# Fix name as necessary.
bc.info[is.na(final.lib.name), final.lib.name:=tcr.lib.name]

# ---> Plate information (custom process).
# @ Extract pool number per sequencing library.
bc.info[, pool.no:=str_extract(string=seq.lib.name, pattern='\\d+$')]
# @ For seq. run 132
# Pools 1 to 24 correspond to plate 1.
tmp.vals <- 1:24
bc.info[seq.run=='NV132' & pool.no%in%tmp.vals, plate.id:='P1']
# Pools 25 to 30 correspond to plate 2.
tmp.vals <- 25:30
bc.info[seq.run=='NV132' & pool.no%in%tmp.vals, plate.id:='P2']
# @ For seq. run 135
# Pools 1 to 18 correspond to plate 1.
tmp.vals <- 1:18
bc.info[seq.run=='NV135' & pool.no%in%tmp.vals, plate.id:='P1']
# Pools 19 to 36 correspond to plate 2.
tmp.vals <- 19:36
bc.info[seq.run=='NV135' & pool.no%in%tmp.vals, plate.id:='P2']
# @ For seq. run 137
# Sample project TCR023 corresponds to a single plate.
bc.info[
    seq.run=='NV137' &
    str_detect(string=seq.lib.name, pattern='TCR023'),
    plate.id:='P1'
]
# Sample project TCR024 plate info is specified in the seq. library name
bc.info[
    seq.run=='NV137' &
    str_detect(string=seq.lib.name, pattern='TCR024'),
    # unique(plate.id)
    plate.id:=str_replace_all(
        string=str_extract(
            string=seq.lib.name,
            pattern='_P[12]+_'
        ),
        pattern='_P*', replacement=''
    )
]
# Adjust for the naming of plate from projetc 23 as P1.
bc.info[
    seq.run=='NV137' &
    str_detect(string=seq.lib.name, pattern='TCR024'),
    plate.id:=paste0('P', as.integer(plate.id)+1)
]
# @ Final check (manual).
# bc.info[, uniqueN(seq.lib.name), by=.(seq.run, plate.id)]
bc.info[, pool.no:=NULL]
# @ Full plate ID.
bc.info[, full.plate.id:=paste(seq.run, plate.id, sep='-')]

# ---> Extract pieces of information from TCR library name.
bc.info[, tcr.lib.info:=final.lib.name]
# Identify samples that represent replicates of the same experiment (those for donors "148" and "141") and extract their replicate ID if available.
bc.info[
    str_detect(string=tcr.lib.info, pattern='^141|^148'),
    exp.sample.rep:=str_extract(
        string=str_extract(string=tcr.lib.info, pattern='^141[-_]\\d+|^148[-_]\\d+'),
        pattern='\\d+$'
    )
]
# Fix barcode library name pattern accordingly.
bc.info[
    str_detect(string=tcr.lib.info, pattern='^141|^148') & !is.na(exp.sample.rep),
    tmp.donor.id:=str_extract(
        string=tcr.lib.info,
        pattern='^\\d+'
    )
]
bc.info[
    str_detect(string=tcr.lib.info, pattern='^141|^148') & !is.na(exp.sample.rep),
    tcr.lib.info:=str_replace(
        string=tcr.lib.info,
        pattern=paste0('^\\d+[_-]', exp.sample.rep),
        replacement=tmp.donor.id
    )
]
bc.info[, tmp.donor.id:=NULL]
# Remove "Redo" suffix for samples that need it, while keeping track of such information as part of a separate column.
bc.info[, is.redo:=str_detect(string=tcr.lib.info, pattern='_Redo$')]
bc.info[, tcr.lib.info:=str_replace(string=tcr.lib.info, pattern='_Redo$', replacement='')]
# For CMV libraries from run NV135, fix library name.
bc.info[
  seq.run=='NV135' & tcr.lib.info%like%'_CMV_',
  tmp.col:=str_extract(string=tcr.lib.info, pattern='[84]$')
]
bc.info[
  seq.run=='NV135' & final.lib.name%like%'_CMV_',
  tcr.lib.info:=str_replace(string=tcr.lib.info, pattern='_[84]$', replacement='')
]
bc.info[
  !is.na(tmp.col),
  tcr.lib.info:=str_replace(string=tcr.lib.info, pattern='_CMV_', replacement=paste0('_CMV', tmp.col, '_'))
]
bc.info[, tmp.col:=NULL]
# For run NV135 (terribly labeled), extract cell subset from rest of info when necessary.
bc.info[
  seq.run=='NV135' & tcr.lib.info%like%'[^_(CD)][48]$',
  tmp.col:=str_extract(string=tcr.lib.info, pattern='[48]$')
]
bc.info[
  !is.na(tmp.col),
  tcr.lib.info:=str_replace(string=tcr.lib.info, pattern=paste0(tmp.col, '$'), replacement=paste0('_CD', tmp.col))
]
# Check for any other samples whose library name do not meet the pattern criterion.
tmp.check <- str_split(string=bc.info[, tcr.lib.info], pattern='_|-', simplify=FALSE)
tmp.check <- unlist(lapply(X=tmp.check, FUN=length))!=3
if(any(tmp.check)){
  tmp.warn <- bc.info[tmp.check, paste0(unique(tcr.lib.info), collapse='\n')]
  tmp.warn <- paste0('Not all of the TCR library names follow the expected pattern (listed below).\n', tmp.warn, '\nWe will attempt to correct them by assuming that the first dash/underscore was dropped.')
  warning(tmp.warn)
  bc.info[
    tmp.check,
    tcr.lib.info:=paste(
      str_extract(string=tcr.lib.info, pattern='^\\d+'),
      str_extract(string=tcr.lib.info, pattern='[^\\d]+(_|-)[^_-]+$'),
      sep='-'
    )
  ]
  tmp.check <- str_split(string=bc.info[, tcr.lib.info], pattern='_|-', simplify=FALSE)
  tmp.check <- unlist(lapply(X=tmp.check, FUN=length))!=3
  if(any(tmp.check)) stop('Failed to correct the TCR library name. Please manually correct them.')
}
# Extract info from TCR library name.
bc.info <- as.data.table(
  separate(
    data=bc.info, 
    col='tcr.lib.info', into=c('donor.id', 'peptide.pool', 'cell.subset'),
    sep='-|_'
  )
)

# ---> Harmonize library metadata.

# @ Set T-cell subset for which the peptide pool was developed.
bc.info[, pp.cell.subset:=str_extract(string=peptide.pool, pattern='[48]$')]
bc.info[is.na(pp.cell.subset), pp.cell.subset:='4']
tmp.vals <- c('4'='CD4', '8'='CD8')
tmp.check <- bc.info[, all(pp.cell.subset %in% names(tmp.vals))]
if(!tmp.check) stop('Failed to define peptide pool T-cell subset for all libraries.\n')
bc.info[, pp.cell.subset:=tmp.vals[pp.cell.subset]]

# @ Harmonize peptide pool.
bc.info[, peptide.pool:=toupper(peptide.pool)]
bc.info[, .N, by=peptide.pool]
pp.vals <- c(
  `CMV4`='CMV4',
  `CMV8`='CMV8',
  `EBV`='EBV',
  `HA`='IAV',
  `FLU`='IAV',
  `SPIKE`='SARS-CoV-2',
  `RSV`='RSV',
  `MPV`='MPV',
  `PERTUSIS`='B-pertussis-Vax',
  `PERT`='B-pertussis-Vax',
  `PR`='B-pertussis-Vax',
  `PT`='B-pertussis-Vax',
  `ASP`='Aspergillus',
  `PI`='PIV',
  `RH`='RV',
  `RHINO`='RV'
)
tmp.check <- bc.info[, all(peptide.pool %in% names(pp.vals))]
if(!tmp.check) stop('Failed to define peptide pool identity all libraries.\n')
bc.info[, peptide.pool:=pp.vals[peptide.pool]]

# @ Set antigen according to peptide pool.
bc.info[, .N, by=peptide.pool]
ag.vals <- c(
    `CMV4`='CMV',
    `CMV8`='CMV',
    `EBV`='EBV',
    `RSV`='RSV',
    `IAV`='IAV',
    `PIV`='PIV',
    `MPV`='MPV',
    `RV`='RV',
    `SARS-CoV-2`='SARS-CoV-2',
    `B-pertussis-Vax`='B-pertussis-Vax',
    `Aspergillus`='Aspergillus'
)
tmp.check <- bc.info[, all(peptide.pool %in% names(ag.vals))]
if(!tmp.check) stop('Failed to define peptide pool identity all libraries.\n')
bc.info[, pool.ag:=ag.vals[peptide.pool]]

# @ Set T-cell subset.
bc.info[, cell.subset:=toupper(cell.subset)]
bc.info[cell.subset=='8', cell.subset:='CD8']
bc.info[cell.subset=='4', cell.subset:='CD4']
tmp.vals <- c('CD4', 'CD8')
tmp.check <- bc.info[, all(cell.subset %in% tmp.vals)]
if(!tmp.check) stop('Unexpected error when defining T-cell subset.\n')

# @ Set donor ID.
#   For SDBB donors.
tmp.vals <- c('141', '148')
bc.info[
    donor.id %in% tmp.vals,
    donor.id:=paste0('SDBB', donor.id)
]
#   For HIPC lung donors.
tmp.vals <- paste0('SDBB', tmp.vals)
bc.info[
    !donor.id %in% tmp.vals,
    tmp.col:=3-str_length(donor.id)
]
bc.info$tmp.col <- unlist(lapply(X=bc.info[, tmp.col], FUN=function(x){
    if(is.na(x)) return(NA) else return(paste0(rep(x='0', times=x), collapse=''))
}))
bc.info[!donor.id %in% tmp.vals, donor.id:=paste0('L', tmp.col, donor.id)]
bc.info[, tmp.col:=NULL]
# bc.info[, gtools::mixedsort(unique(donor.id))]

# @ Set sample status regarding control experiments.
bc.info[donor.id %like% '^SDBB', sample.status:='Control']
bc.info[is.na(sample.status), sample.status:='Index']

# @ Set final condition labels.
bc.info[,
    exp.condition:=paste(tcr.chain, donor.id, peptide.pool, cell.subset, sep='.')
]

# @ Define experimental replicates.
bc.info[,
    sample.rep:=1:.SD[, .N],
    by=exp.condition
]
tmp.conds <- bc.info[
    sample.rep>1,
    exp.condition
]
bc.info[,
    has.rep:=exp.condition %in% tmp.conds
]
# Specific sanity check: all sample replicates should have been performed with the same type of peptide pool.
tmp.check <- bc.info[,
    .(pp.count=uniqueN(pp.cell.subset)),
    by=exp.condition
][pp.count>1, .N==0]
if(!tmp.check) stop('Unexpected error.\n')

# ---> Output preflight checks.
# Output file for manual inspection of variable value assignment.
tmp.data <- unique(bc.info[, .(exp.condition, seq.run, sample.prj, tcr.lib.id, tcr.chain, tcr.lib.name=final.lib.name, donor.id, peptide.pool, cell.subset, pp.cell.subset, exp.sample.rep, has.rep, sample.rep, is.redo)])
setorderv(x=tmp.data, cols=c('exp.condition', 'seq.run', 'sample.prj', 'tcr.chain', 'donor.id', 'peptide.pool', 'pp.cell.subset', 'sample.rep'))
tmp.file.name <- paste0(prepro.path, '/PreflightChecks_1.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)
# Alternative version.
tmp.data <- unique(bc.info[, .(seq.run, sample.prj, tcr.chain, tcr.lib.name=final.lib.name, donor.id, peptide.pool, cell.subset, pp.cell.subset, exp.sample.rep, has.rep, sample.rep)])
tmp.data.1 <- tmp.data[tcr.chain=='TCRA']; tmp.data.1[, tcr.chain:=NULL]
tmp.data.2 <- tmp.data[tcr.chain=='TCRB']; tmp.data.2[, tcr.chain:=NULL]
tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by=c('seq.run', 'sample.prj', 'tcr.lib.name', 'donor.id', 'peptide.pool', 'cell.subset', 'pp.cell.subset', 'exp.sample.rep', 'sample.rep'), suffixes=c('.alpha', '.beta'), all=TRUE)
setorderv(x=tmp.data, cols=c('seq.run', 'sample.prj', 'donor.id', 'peptide.pool', 'pp.cell.subset', 'sample.rep'))
tmp.file.name <- paste0(prepro.path, '/PreflightChecks_2.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)
# Alternative version.
tmp.vals <- c('SDBB141', 'SDBB148')
tmp.data <- bc.info[
  !donor.id %in% tmp.vals,
  .(sample.count=uniqueN(tcr.lib.id)),
  by=.(donor.id, peptide.pool)
]
tmp.data <- spread(data=tmp.data, key='peptide.pool', value='sample.count', fill=0)
tmp.file.name <- paste0(prepro.path, '/PreflightChecks_3.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)


cat('\n\n')
### --------------------- Load sample-based data -------------------- ###
cat('### --------------------- Load sample-based data -------------------- ###\n')

# ---> Load MIGEC outputs
# Define final sample name prefix.
bc.info[, tcr.lib.file.prfx:=paste(tcr.lib.id, tcr.lib.name, tcr.chain, sep='.')]
# File definition.
tcr.clone.files <- lapply(X=bc.info[, unique(seq.run)], FUN=function(tmp.run){
    tcr.clone.files <- paste0(
        cdr.filter.paths[tmp.run], '/',
        bc.info[seq.run==tmp.run, tcr.lib.file.prfx],
        '.filtered.cdrblast.txt'
    )
    names(tcr.clone.files) <- bc.info[seq.run==tmp.run, tcr.lib.file.prfx]
    return(tcr.clone.files)
})
tcr.clone.files <- unlist(tcr.clone.files)
# Checkpoint.
tmp.check <- all(file.exists(tcr.clone.files))
if(!tmp.check){
    tmp.check <- tcr.clone.files[!file.exists(tcr.clone.files)]
    tmp.check <- paste0('Failed to define the MIGEC output files. Next are the failed attempted definitions:\n', paste0(tmp.check, collapse='\n'))
    stop(tmp.check)
}
# Load data.
tcr.data <- lapply(X=tcr.clone.files, FUN=fread, na.strings='.')
tcr.data <- rbindlist(l=tcr.data, use.names=TRUE, idcol='tcr.lib.file.prfx')
tmp.cols <- c(
    `tcr.lib.file.prfx`='tcr.lib.file.prfx',
    `cdr3.nt.seq`='CDR3 nucleotide sequence',
    `cdr3.aa.seq`='CDR3 amino acid sequence',
    `trv.allele`='V segments',
    `trj.allele`='J segments',
    `trd.allele`='D segments',
    `umi.count`='Good events',
    `read.count`='Good reads'
)
tcr.data <- tcr.data[, ..tmp.cols]; colnames(tcr.data) <- names(tmp.cols)
# Merge w/ relevant pieces from barcode info.
tmp.data <- bc.info[, .(
    tcr.lib.file.prfx,
    seq.run, sample.prj, plate.id, full.plate.id, seq.lib.name,
    tcr.lib.id, tcr.chain,
    sample.status, exp.condition,
    donor.id, cell.subset,
    peptide.pool, pp.cell.subset, pool.ag,
    sample.rep, has.rep
)]
tcr.data <- merge(
  x=tmp.data, y=tcr.data,
  by='tcr.lib.file.prfx', all=TRUE
)
tcr.data[, tcr.lib.file.prfx:=NULL]
# # TO REMOVE. JUST FOR GEN. TRIAL. -------
# tmp.rows <- 1:tcr.data[, .N]
# tmp.rows <- sample(x=tmp.rows, size=floor(length(tmp.rows)/10))
# tcr.data <- tcr.data[tmp.rows]
# # TO REMOVE. JUST FOR GEN. TRIAL. -------

# ---> Define base v and j genes (w/out allele info).
these.cols <- paste0('tr', c('v', 'j'), '.allele')
for(full.col in these.cols){
    tcr.data[, allele:=get(full.col)]
    allele.count <- tcr.data[, max(str_count(string=allele, pattern=','))+1]
    allele.cols <- paste0('allele.', 1:allele.count)
    tcr.data <- as.data.table(separate(
        data=tcr.data, sep=',',
        col=allele, into=allele.cols,
        remove=TRUE, fill='right'
    ))
    for(allele.col in allele.cols){
        tcr.data[[allele.col]] <- str_extract(string=tcr.data[[allele.col]], pattern='^[^*]+')
        tcr.data[[allele.col]][is.na(tcr.data[[allele.col]])] <- ''
    }
    simp.col <- str_extract(string=full.col, pattern='tr[vdj]')
    tcr.data[[simp.col]] <- apply(
            X=tcr.data[, ..allele.cols],
            MARGIN=1, FUN=paste, collapse=';'
        )
    tcr.data[[simp.col]] <- str_replace(string=tcr.data[[simp.col]], pattern=';$', replacement='')
    # Remove irrelevant columns.
    for(allele.col in allele.cols) tcr.data[[allele.col]] <- NULL
}
# Manual inspection of corresponding values.
tmp.data <- lapply(X=these.cols, FUN=function(full.col){
    simp.col <- str_extract(string=full.col, pattern='tr[vdj]')
    tcr.data[, .(entry.count=.N), by=
        .(
            allele=get(full.col),
            gene=get(simp.col)
        )
    ]
})
names(tmp.data) <- str_extract(string=these.cols, pattern='^tr[vj]')
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='gene.type')
setorderv(x=tmp.data, cols=c('gene.type', 'entry.count'), order=-1)
tmp.file.name <- paste0(prepro.path, '/PreflightChecks_4.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)

# ---> Set clonotype ID.
# Define sets of IDs, one for each set of T-cell subset and TCR chain type.
chain.types <- tcr.data[, unique(tcr.chain)]
tmp.cols <- c('tcr.chain', 'cdr3.nt.seq', 'trv', 'trj')
tmp.data <- lapply(X=chain.types, FUN=function(chain.type){
    tmp.data <- tcr.data[
        tcr.chain==chain.type,
        .(total.count=sum(umi.count)),
        by=tmp.cols
    ]
    setorderv(x=tmp.data, cols='total.count', order=-1)
    tmp.data[, total.count:=NULL]
    tmp.data[, tcr.id:=paste0(
        paste0('000000', collapse=''),
        1:.N
    )]
    tmp.data[, tcr.id:=str_extract(string=tcr.id, pattern='\\d{6}$')]
    tmp.data[,
        tcr.id:=paste0(
            'ATCR', # Preffix "ATCR" stands for "Assay TCR"
            ifelse(test=chain.type=='TCRA', yes='A', no='B'), '-',
            tcr.id
        )
    ]
})
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE)
tcr.data <- merge(x=tmp.data, y=tcr.data, by=tmp.cols, all=TRUE)
tmp.check <- tcr.data[, sum(is.na(tcr.id))==0]
if(!tmp.check) stop('Failed to assign unique clonotype IDs.')
# Make sure that there's unique entries for clonotype IDs within donors.
tmp.check <- tcr.data[has.rep==FALSE, .N, by=.(tcr.id, donor.id, peptide.pool)][, sum(N>2)==0] # There might be up to 2 entries for each clonoype given that the same clonotype can be captured in either or both T cell lineages.
if(!tmp.check) stop('Failed to assign unique IDs to TCRs. Some of them ended up being redundant.\n')


cat('\n\n')
### ---------------------- Artifact detection ----------------------- ###
cat('### ---------------------- Artifact detection ----------------------- ###\n')

art.hunt.path <- paste0(reports.path, '/artifact_detection')
if(!dir.exists(art.hunt.path)) dir.create(art.hunt.path)

# ---> Fetch data to work with
# Consider only index samples w/ properly labeled repeats # Note that at this point, as shown on the PreflightCheks tables, none of the index samples should have any technical repeats.
merge.cols <- c('tcr.chain', 'cell.subset', 'donor.id', 'peptide.pool')
tmp.data <- unique(bc.info[has.rep==FALSE & sample.status=='Index', ..merge.cols])
tmp.data <- merge(
  x=tcr.data, y=tmp.data,
  by=merge.cols, all=FALSE
)
# Condense into uniquely defined clonotypes.
clone.cols <- c('cdr3.nt.seq', 'trv', 'trj')
tmp.cols <- c('seq.run', 'plate.id', 'full.plate.id', 'exp.condition', merge.cols, clone.cols)
tmp.data <- tmp.data[,
  .(
    umi.abs.freq=sum(umi.count),
    read.abs.freq=sum(read.count)
  ),
  by=tmp.cols
]
# Add relative frequencies
tmp.data.1 <- tmp.data[,
  .(
    umi.total.freq=sum(umi.abs.freq),
    read.total.freq=sum(read.abs.freq)
  ),
  by=exp.condition
]
tmp.data <- merge(x=tmp.data, y=tmp.data.1, by='exp.condition', all=TRUE)
tmp.data[, umi.rel.freq:=umi.abs.freq/umi.total.freq]
tmp.data[, read.rel.freq:=read.abs.freq/read.total.freq]
# Add putative publicity degree to each TCR.
tmp.data.1 <- tmp.data[,
  .(public.degree=uniqueN(donor.id)),
  by=clone.cols
]
tmp.data.1 <- merge(x=tmp.data, y=tmp.data.1, by=clone.cols, all=TRUE)
tmp.check <- tmp.data.1[, .N]==tmp.data[, .N]
if(!tmp.check) stop('Something went wrong.\n')
tmp.data <- tmp.data.1

# ---> Plate-based info
tmp.ggplot.1 <- ggplot(data=tmp.data, aes(x=plate.id, y=public.degree)) +
  geom_boxplot(outlier.shape=NA, linewidth=2, alpha=0) +
  facet_wrap(facets=~seq.run, nrow=1, scales='free_x') +
  labs(title='Degree distribution on a clonotype basis', x='Plate ID', y='Publicity degree') +
  coord_cartesian(ylim=c(0, 28)) +
  scale_y_continuous(expand=c(0, 0)) +
  theme_bw()
to.plot <- tmp.data[,
  .(sample.count=uniqueN(exp.condition)),
  by=.(seq.run, plate.id)
]
tmp.ggplot.2 <- ggplot(data=to.plot, aes(x=plate.id, y=sample.count)) +
  geom_bar(stat='identity', linewidth=2, color='black', alpha=0) +
  facet_wrap(facets=~seq.run, nrow=1, scales='free_x') +
  labs(title='Sample count per plate', x='Plate ID', y='Sample count') +
  scale_y_continuous(expand=c(0, 1)) +
  theme_bw()
to.plot <- tmp.data[,
  .(clonotype.count=.SD[, .N, by=clone.cols][, .N]),
  by=.(exp.condition, seq.run, plate.id)
]
tmp.ggplot.3 <- ggplot(data=to.plot, aes(x=plate.id, y=clonotype.count)) +
  geom_boxplot(outlier.shape=NA, linewidth=2, alpha=0) +
  geom_jitter(shape=21, size=2, width=0.3) +
  facet_wrap(facets=~seq.run, nrow=1, scales='free_x') +
  labs(title='Clonotype count per sample when segregated by plate', x='Plate ID', y='No. of clonotypes', caption='Every dot represents a different sample') +
  scale_y_continuous(expand=c(0, 1)) +
  theme_bw()
to.plot <- tmp.data[,
  .(umi.total.freq=sum(umi.abs.freq)),
  by=.(exp.condition, seq.run, plate.id)
]
tmp.ggplot.4 <- ggplot(data=to.plot, aes(x=plate.id, y=umi.total.freq)) +
  geom_boxplot(outlier.shape=NA, linewidth=2, alpha=0) +
  geom_jitter(shape=21, size=2, width=0.3) +
  facet_wrap(facets=~seq.run, nrow=1, scales='free_x') +
  labs(title='Total UMI frequency per sample when segregated by plate', x='Plate ID', y='Total UMI frequency') +
  scale_y_continuous(expand=c(0, 1)) +
  theme_bw()
tmp.file.name <- paste0(art.hunt.path, '/PlateBasedGenInfo.pdf')
pdf(file=tmp.file.name, height=5)
print(tmp.ggplot.1)
print(tmp.ggplot.2)
print(tmp.ggplot.3)
print(tmp.ggplot.4)
dev.off()

# ---> Artifact TCR abundance
# Process as a general function dependent on TCR degree and chain under evaluation
show.artifact.abund <- function(
  degree.thold, these.chains
){
  # Retrieve general heatmap format. 
  tmp.cols <- c(clone.cols, 'exp.condition', 'umi.rel.freq')
  gen.to.plot <- tmp.data[
    public.degree>=degree.thold & tcr.chain%in%these.chains,
    ..tmp.cols
  ]
  tmp.cols <- gen.to.plot[, unique(exp.condition)]
  gen.to.plot <- spread(data=gen.to.plot, key=exp.condition, value=umi.rel.freq, fill=0)
  gen.to.plot <- as.matrix(gen.to.plot[, tmp.cols])
  # Process for different upper bounds.
  # tmp.bounds <- c(`None`=NA, `0.05`=0.05, `0.01`=0.01)
  tmp.bounds <- c(`0.05`=0.05, `0.01`=0.01)
  for(bound.lab in names(tmp.bounds)){
    # To prevent outliers from drowning the signal.
    upper.bound <- tmp.bounds[bound.lab]
    to.plot <- copy(gen.to.plot)
    if(!is.na(upper.bound)) to.plot[to.plot>upper.bound] <- upper.bound
    # Set color scale and breaks for heatmap.
    col.breaks <- seq(from=min(range(to.plot)), to=max(range(to.plot)), length.out=100)
    mid.point <- which.min(abs(col.breaks < 0.02))
    hmap.col.scale.1 <- colorRampPalette(c('#ffffff', '#ffff9a', '#ffff4d', '#ffd34d', '#ffc000', '#ff4d4d'))(mid.point)
    hmap.col.scale.2 <- colorRampPalette(c('#ff4d4d', '#ff0000', '#890000', '#140000'))(100-(mid.point+1))
    hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
    # Set metadata.
    tmp.cols <- c(
      `sample`='exp.condition',
      `Run`='seq.run', `Plate`='plate.id',
      `Chain`='tcr.chain', `Lineage`='cell.subset'
    )
    col.meta <- as.data.frame(unique(tmp.data[, ..tmp.cols]))
    colnames(col.meta) <- names(tmp.cols)
    row.names(col.meta) <- col.meta$sample; col.meta$sample <- NULL
    # Define colors for metadata tracks
    meta.colors <- list(
        Lineage=c(
            CD4='#00BFFF', CD8='#EE82EE'
        ),
        Plate=c(
            `P1`='#ffe500', `P2`='#009aff', `P3`='#ff001b'
        ),
        Run=c(
            NV132='#b300b3', NV135='#59b300', NV137='#0100b3'
        ),
        Chain=c(
            TCRA='#0095B3', TCRB='#808080'
        )
    )
    # Plot.
    tmp.file.name <- paste0(art.hunt.path, '/_ClonotypeUMIRelFreq_DegreeThold-', degree.thold, '_Chains-', paste0(these.chains, collapse='-'), '_Bound-', bound.lab, '.pdf')
    pheatmap(
        mat=to.plot,
        color=hmap.col.scale, breaks=col.breaks, scale='none',
        cluster_rows=TRUE, cluster_cols=TRUE,
        annotation_col=col.meta, annotation_colors=meta.colors,
        show_rownames=FALSE, show_colnames=FALSE,
        filename=tmp.file.name
    )
  }
}

# Process for both chains.
show.artifact.abund(degree.thold=10, these.chains=c('TCRA', 'TCRB'))
# Process for alpha chain only
show.artifact.abund(degree.thold=5, these.chains=c('TCRA'))
show.artifact.abund(degree.thold=10, these.chains=c('TCRA'))
# Process for beta chain only
show.artifact.abund(degree.thold=5, these.chains=c('TCRB'))
show.artifact.abund(degree.thold=10, these.chains=c('TCRB'))


# ---> Artifact TCR sharing
# @ Preflights
tmp.data[,
  `:=`(
    sample.id=paste(tcr.chain, full.plate.id, sep='-'),
    clone.id=paste(cdr3.nt.seq, trv, trj, sep='_')
  )
]

# @ TCR sharing between samples according to TCR chains.
tmp.vals <- tmp.data[, unique(tcr.chain)]
tmp.file.name <- paste0(art.hunt.path, '/SharingAccordingToChains.pdf')
pdf(file=tmp.file.name, height=10, width=15)
size.vals <- 1:2
for(size.thold in size.vals){
  to.plot <- lapply(X=tmp.vals, FUN=function(tmp.val){
    tmp.data <- tmp.data[
      tcr.chain==tmp.val & umi.abs.freq>=size.thold,
      clone.id
    ]
    return(tmp.data)
  })
  names(to.plot) <- tmp.vals
  print(upset(
      fromList(to.plot),
      nsets=length(to.plot),
      nintersects=NA, order.by='freq',
      main.bar.color='lightgray', att.color='black',
      show.numbers='no',
      # scale.intersections='log10',
      mainbar.y.label=paste0('Intersection size for clonotypes w/ UMI size equal to or larger than ', size.thold)
  ))
} 
dev.off()

# @ TCR sharing between samples according to plate and on a chain basis
for(tmp.chain in tmp.data[, unique(tcr.chain)]){
  tmp.vals <- tmp.data[tcr.chain==tmp.chain, unique(sample.id)]
  tmp.file.name <- paste0(art.hunt.path, '/SharingAccordingToPlates_Chain-', tmp.chain, '.pdf')
  pdf(file=tmp.file.name, height=10, width=15)
  size.vals <- 2:10
  for(size.thold in size.vals){
    to.plot <- lapply(X=tmp.vals, FUN=function(tmp.val){
      tmp.data <- tmp.data[
        sample.id==tmp.val & umi.abs.freq>=size.thold,
        clone.id
      ]
      return(tmp.data)
    })
    names(to.plot) <- tmp.vals
    print(upset(
        fromList(to.plot),
        nsets=length(to.plot),
        nintersects=NA, order.by='freq',
        main.bar.color='lightgray', att.color='black',
        show.numbers='no',
        # scale.intersections='log10',
        mainbar.y.label=paste0('Intersection size for clonotypes w/ UMI size equal to or larger than ', size.thold)
    ))
  } 
  dev.off()
}


# ----> Artifact calling on a plate basis.
# Criterion: TCRs found in 50% of all samples in a plate only when index samples are considered.
# art.call.thold <- 0.5
tmp.vals <- tcr.data[, sort(unique(full.plate.id))]
tmp.data <- lapply(X=tmp.vals, FUN=function(tmp.plate){
  tmp.data <- tcr.data[full.plate.id==tmp.plate & sample.status=='Index']
  chain.types <- tmp.data[, unique(tcr.chain)]
  tmp.data <- lapply(X=chain.types, FUN=function(chain.type){
    tmp.data <- tmp.data[tcr.chain==chain.type]
    sample.thold <- floor(tmp.data[, uniqueN(exp.condition)]*art.call.thold)
    tmp.data <- tmp.data[,
      .(is.artifact=uniqueN(exp.condition)>=sample.thold),
      by=.(tcr.id)
    ]
    return(tmp.data)
  })
  tmp.data <- rbindlist(l=tmp.data, use.names=TRUE)
  return(tmp.data)
})
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE)
# Harmonize calls across plates. A clonotype will be called an artifact if it was identified as such for any of the plates.
tmp.data <- tmp.data[, .(is.artifact=any(is.artifact)), by=tcr.id]
tcr.data <- merge(x=tcr.data, y=tmp.data, by='tcr.id')

# Fraction of artifactual clonotypes and UMIs per sample.
to.plot <- tcr.data[
  sample.status=='Index',
  .(
    art.clone.rel.freq=.SD[is.artifact==TRUE, uniqueN(tcr.id)]/.SD[, uniqueN(tcr.id)],
    art.umi.rel.freq=.SD[is.artifact==TRUE, sum(umi.count)]/.SD[, sum(umi.count)]
  ),
  by=.(exp.condition, seq.run, plate.id, tcr.chain)
]
tmp.ggplot.1 <- ggplot(data=to.plot, aes(x=plate.id, y=art.clone.rel.freq)) +
  geom_boxplot(outlier.shape=NA, linewidth=2, alpha=0) +
  geom_jitter(aes(color=tcr.chain), shape=21, size=2, width=0.3) +
  facet_wrap(facets=~seq.run, nrow=1, scales='free_x') +
  labs(title='Artifactual clonotype relative frequency per sample', x='Plate ID', y='Artifactual clonotype relative frequency', color='Chain', caption='Every dot represents a different sample') +
  # scale_y_continuous(expand=c(0, 1)) +
  theme_bw()
tmp.ggplot.2 <- ggplot(data=to.plot, aes(x=plate.id, y=art.umi.rel.freq)) +
  geom_boxplot(outlier.shape=NA, linewidth=2, alpha=0) +
  geom_jitter(aes(color=tcr.chain), shape=21, size=2, width=0.3) +
  facet_wrap(facets=~seq.run, nrow=1, scales='free_x') +
  labs(title='Artifactual UMI relative frequency per sample', x='Plate ID', y='Artifactual UMI relative frequency', color='Chain', caption='Every dot represents a different sample') +
  # scale_y_continuous(expand=c(0, 1)) +
  theme_bw()
tmp.file.name <- paste0(art.hunt.path, '/PlateBasedArtifactInfo.pdf')
pdf(file=tmp.file.name, height=5)
print(tmp.ggplot.1)
print(tmp.ggplot.2)
dev.off()
# Report.
tmp.file.name <- paste0(art.hunt.path, '/PlateBasedArtifactInfo.csv')
fwrite(file=tmp.file.name, x=to.plot, quote=TRUE)

# ---> Artifact removal
tcr.data <- tcr.data[is.artifact==FALSE]
tcr.data[, is.artifact:=NULL]


### ------------------------- Main program -------------------------- ###


cat('\n\n')
### ----------------- Antigen specificity definition ---------------- ###
cat('### ----------------- Antigen specificity definition ---------------- ###\n')

# Here, we have the advantage that specificity is given as a binary prediction (reactive vs non-reactive) for every single peptide, allowing us to better better discern between monoreactive and  crossreactive TCRs.

ag.spc.path <- paste0(reports.path, '/ag_specificity_determination')
if(!dir.exists(ag.spc.path)) dir.create(ag.spc.path)

# ---> Define data for which the process is necessary.
# Subset for TCRs from samples that meet the following criteria:
# The sample does not have a replicate (rationale: there were no replicates for index samples).
# The sample comes from an index donor.
# The peptide pool used for the sample matches its corresponding T-cell subset in the case of CMV.
# Data from the RV peptide pool are fully disregarded.
this.tcr.data <- tcr.data[
    has.rep==FALSE & sample.status=='Index' &
    !(pool.ag=='CMV' & cell.subset!=pp.cell.subset) &
    pool.ag!='RV'
]
this.tcr.data <- this.tcr.data[, .(tcr.id, cell.subset, tcr.chain, donor.id.tag=donor.id, peptide.pool.tag=pool.ag, fad=umi.count)]
this.tcr.data[, long.id.tag:=paste(tcr.id, cell.subset, donor.id.tag, sep=';')]
# Identify peptide pools under consideration.
peptide.pools <- this.tcr.data[, unique(peptide.pool.tag)]

# ----> General clonotype status.
# Define clonotype general status among: 
#   - Singleton: Clonotype ID is observed w/ a single UMI count within a given donor.
#   - Ambiguous non-singleton: Clonotype ID is observed w/ a UMI count greater than 1 within a given donor, but with these coming up for two or more peptide pool antigens.
#   - Unambiguous non-singleton: Clonotype ID is observed w/ a UMI count greater than 1 for a single peptide pool antigen (and 0 for any other antigen) within a given donor.
tmp.data <- this.tcr.data[,
    .(
        gen.clonotype.status=ifelse(
            test=sum(fad)==1,
            yes='Singleton',
            no=ifelse(
                test=uniqueN(peptide.pool.tag)==1,
                yes='Unambiguous',
                no='Ambiguous'
            )
        )
    ),
    by=.(long.id.tag)
]
tmp.data$gen.clonotype.status <- factor(x=tmp.data$gen.clonotype.status, levels=c('Singleton', 'Ambiguous', 'Unambiguous'))
this.tcr.data <- merge(x=this.tcr.data, y=tmp.data, by='long.id.tag', all=TRUE)
# Explore frequency of different statuses in the whole dataset and for each peptide pool antigen per cell type.
tmp.data <- this.tcr.data[, .(gen.clonotype.status, cell.subset, peptide.pool.tag)]
tmp.ggplot.1 <- ggplot(data=tmp.data, aes(x=gen.clonotype.status)) +
  geom_bar(width=0.7, fill='lightblue', color='black', linewidth=1.5, alpha=0.8) +
  scale_y_continuous(expand=c(0, NA)) +
  labs(x='General clonotype status', y='Clonotype count') +
  theme_bw()
tmp.ggplot.2 <- ggplot(data=tmp.data, aes(x=gen.clonotype.status)) +
  geom_bar(width=0.7, fill='lightblue', color='black', linewidth=1.5, alpha=0.8) +
  facet_wrap(facets=~cell.subset, scales='free_y', ncol=1) +
  scale_y_continuous(expand=c(0, NA)) +
  labs(x='General clonotype status', y='Clonotype count') +
  theme_bw()
tmp.ggplot.3 <- ggplot(data=tmp.data, aes(x=gen.clonotype.status, fill=peptide.pool.tag)) +
  geom_bar(width=0.7, color='black', linewidth=1.5, alpha=0.8) +
  facet_wrap(facets=~cell.subset+peptide.pool.tag, scales='free_y') +
  scale_y_continuous(expand=c(0, NA)) +
  scale_fill_manual(values=ag.cols) +
  labs(x='General clonotype status', y='Clonotype count') +
  theme_bw() + theme(legend.position='none', axis.text.x=element_text(angle=-45))
tmp.file.name <- paste0(ag.spc.path, '/GeneralClonotypeStatus_Freqs.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot.1)
print(tmp.ggplot.2)
print(tmp.ggplot.3)
dev.off()

# ---> Within-donor relative frequencies (FRD).
# Remove singletons at this point.
this.tcr.data <- this.tcr.data[gen.clonotype.status!='Singleton']
# Define within-donor absolute and relative frequencies (FADs and FRDs, respectively).
tmp.data <- this.tcr.data[, .(faw=sum(fad)), by=.(tcr.chain, cell.subset, donor.id.tag, peptide.pool.tag)]
this.tcr.data <- merge(x=this.tcr.data, y=tmp.data, by=c('tcr.chain', 'cell.subset', 'donor.id.tag', 'peptide.pool.tag'))
this.tcr.data[, frd:=fad/faw]
# Per (long) clonotype ID, define top first and second peptide pool antigens when they're ranked according to the FRD they show. Keep track of their FRDs.
setorderv(x=this.tcr.data, cols='frd', order=-1)
tmp.data.1 <- this.tcr.data[,
  .(
    pp.1=.SD[1, peptide.pool.tag],
    frd.1=.SD[1, frd],
    fad.1=.SD[1, fad]
  ),
  by=.(long.id.tag)
]
tmp.data.1[, is.top.1:=TRUE]
tmp.data.2 <- this.tcr.data[, .(long.id.tag, peptide.pool.tag, frd, fad)]
tmp.data <- merge(
  x=tmp.data.2, y=tmp.data.1,
  by.x=c('long.id.tag', 'peptide.pool.tag'),
  by.y=c('long.id.tag', 'pp.1'),
  all=TRUE
)
tmp.check <- tmp.data[, .N]==tmp.data.2[, .N]
if(!tmp.check) stop('Something went wrong while defining antigen specificities (1).\n')
tmp.data[is.na(is.top.1), is.top.1:=FALSE]
tmp.data <- tmp.data[is.top.1==FALSE]
setorderv(x=tmp.data, cols='frd', order=-1)
tmp.data.2 <- tmp.data[,
  .(
    pp.2=.SD[1, peptide.pool.tag],
    frd.2=.SD[1, frd],
    fad.2=.SD[1, fad]
  ),
  by=.(long.id.tag)
]
tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by='long.id.tag', all=TRUE)
tmp.data[, is.top.1:=NULL]
tmp.data.1 <- unique(this.tcr.data[, .(long.id.tag, gen.clonotype.status)])
tmp.data <- merge(x=tmp.data.1, y=tmp.data, by='long.id.tag')

# ---> Specificity fold change (FC)
# Calculate the fold change between the top first and second FRDs; we will refer to this FC as the "specificity" FC.
tmp.data[, spc.fc:=frd.1/frd.2]
spc.data <- copy(tmp.data)
# Explore FRD distribution for non-singletons.
tmp.data <- spc.data[gen.clonotype.status!='Singleton']
frd.tholds <- tmp.data[, quantile(x=frd.1, probs=c(0.5, 0.75, 0.9))]
tmp.ggplot.1 <- ggplot(data=tmp.data, aes(x=frd.1)) +
  geom_density(linewidth=1.2, color='black', alpha=0.2, fill='pink') +
    scale_y_continuous(expand=c(0, 0)) + scale_x_continuous(expand=c(0, 0)) +
    labs(x='Maximum within-donor relative frequency', y='Density') +
    theme_bw()
tmp.ggplot.2 <- ggplot(data=tmp.data, aes(x=frd.1)) +
    geom_density(linewidth=1.2, color='black', alpha=0.2, fill='pink') +
    geom_vline(xintercept=frd.tholds[1], linewidth=1, linetype='dashed', color='red') +
    geom_vline(xintercept=frd.tholds[2], linewidth=1, linetype='dashed', color='red') +
    scale_y_continuous(expand=c(0, 0)) + scale_x_continuous(expand=c(0, 0)) +
    coord_cartesian(xlim=c(NA, frd.tholds[3])) +
    labs(
        x='Maximum within-donor relative frequency',
        y='Density',
        caption=paste0('Zoomed in on x axis to avoid showing the long tail of the distribution.\nThe ', names(frd.tholds)[1], ' and ', names(frd.tholds)[2], ' percentiles are ', round(x=frd.tholds[1], digits=8), ' and ', round(x=frd.tholds[2], digits=8), ', respectively.')
    ) +
    theme_bw()
tmp.ggplot.3 <- ggplot(data=tmp.data, aes(x=frd.1)) +
    geom_density(linewidth=1.2, color='black', alpha=0.2, fill='pink') +
    geom_vline(xintercept=frd.tholds[1], linewidth=1, linetype='dashed', color='red') +
    geom_vline(xintercept=frd.tholds[2], linewidth=1, linetype='dashed', color='red') +
    facet_wrap(facets=~gen.clonotype.status, ncol=1) +
    scale_y_continuous(expand=c(0, 0)) + scale_x_continuous(expand=c(0, 0)) +
    coord_cartesian(xlim=c(NA, frd.tholds[3])) +
    labs(
        x='Maximum within-donor relative frequency',
        y='Density',
        caption=paste0('Zoomed in on x axis to avoid showing the long tail of the distribution.\nThe ', names(frd.tholds)[1], ' and ', names(frd.tholds)[2], ' percentiles are ', round(x=frd.tholds[1], digits=8), ' and ', round(x=frd.tholds[2], digits=8), ', respectively.')
    ) +
    theme_bw()
tmp.file.name <- paste0(ag.spc.path, '/MaxWithinDonorFreq_Dists.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot.1)
print(tmp.ggplot.2)
print(tmp.ggplot.3)
dev.off()
# Explore specificity FC distribution for ambiguous clonotypes.
tmp.data <- spc.data[gen.clonotype.status=='Ambiguous']
fc.tholds <- tmp.data[, quantile(x=spc.fc, probs=c(0.20, 0.60, 0.9))]
tmp.ggplot.1 <- ggplot(data=tmp.data, aes(x=spc.fc)) +
    geom_density(linewidth=1.2, color='black', alpha=0.2, fill='pink') +
    scale_y_continuous(expand=c(0, 0)) + scale_x_continuous(expand=c(0, 0)) +
    labs(x='Relative frequency fold change (top vs second)', y='Density') +
    theme_bw()
tmp.ggplot.2 <- ggplot(data=tmp.data, aes(x=spc.fc)) +
    geom_density(linewidth=1.2, color='black', alpha=0.2, fill='pink') +
    geom_vline(xintercept=fc.tholds[1], linewidth=1, linetype='dashed', color='red') +
    geom_vline(xintercept=fc.tholds[2], linewidth=1, linetype='dashed', color='red') +
    scale_y_continuous(expand=c(0, 0)) + scale_x_continuous(expand=c(0, 0)) +
    coord_cartesian(xlim=c(NA, fc.tholds[3])) +
    labs(
        x='Relative frequency fold change (top vs second)',
        y='Density',
        caption=paste0('Zoomed in on x axis to avoid showing the long tail of the distribution.\nThe ', names(fc.tholds)[1], ' and ', names(fc.tholds)[2], ' percentiles are ', round(x=fc.tholds[1], digits=2), ' and ', round(x=fc.tholds[2], digits=2), ', respectively.')
    ) +
    theme_bw()
tmp.file.name <- paste0(ag.spc.path, '/SpecificityFoldChange_Dists.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot.1)
print(tmp.ggplot.2)
dev.off()

# ---> Assign antigen specificities.
# Define final specificity.
spc.data[
  gen.clonotype.status=='Unambiguous',
  specificity.class.tag:=pp.1
]
spc.data[
  gen.clonotype.status=='Ambiguous' &
  spc.fc>=2 &
  fad.1>1,
  specificity.class.tag:=pp.1
]
# Define list of all antigens each clonotype had UMI counts for.
tmp.data <- this.tcr.data[,
  .(specificity.ags.tag=paste(sort(unique(peptide.pool.tag)), collapse=';')),
  by=long.id.tag
]
spc.data <- merge(x=spc.data, y=tmp.data, by='long.id.tag')
# Define a raw specificity tag that remains the same as the one defined above for clonotypes w/ a final value different than NA, but with the value borrowed from the top peptide pool for the rest of the clonotypes.
spc.data[, specificity.raw.tag:=specificity.class.tag]
spc.data[is.na(specificity.raw.tag), specificity.raw.tag:=pp.1]

# ---> Frequency-based confidence level.
tmp.data <- spc.data[!is.na(specificity.class.tag)]
frd.tholds <- tmp.data[, quantile(x=frd.1, probs=c(0.2, 0.4, 0.6, 0.8))]
spc.data[, freq.conf.lvl:=1]
spc.data[frd.1>frd.tholds['20%'], freq.conf.lvl:=2]
spc.data[frd.1>frd.tholds['40%'], freq.conf.lvl:=3]
spc.data[frd.1>frd.tholds['60%'], freq.conf.lvl:=4]
spc.data[frd.1>frd.tholds['80%'], freq.conf.lvl:=5]
spc.data[is.na(specificity.class.tag), freq.conf.lvl:=0]

# ---> Info in metadata.
spc.data <- separate(data=spc.data, col='long.id.tag', into=c('tcr.id', 'cell.subset', 'donor.id'), sep=';')
tcr.data <- merge(x=tcr.data, y=spc.data, by=c('tcr.id', 'cell.subset', 'donor.id'), all=TRUE)
tmp.lvls <- names(ag.cols)[names(ag.cols)%in%tcr.data$specificity.class.tag]
tcr.data$specificity.class.tag <- factor(x=tcr.data$specificity.class.tag, levels=tmp.lvls)
tcr.data$specificity.raw.tag <- factor(x=tcr.data$specificity.raw.tag, levels=tmp.lvls)

# ---> Exploration of final specificity assignments.
tmp.caption <- 'Each facet indicates confidence level based on frequency.'
# Summary across donors.
tmp.data <- unique(tcr.data[
  !is.na(freq.conf.lvl) & freq.conf.lvl!=0,
  .(tcr.id, donor.id, cell.subset, tcr.chain, specificity.class.tag, freq.conf.lvl)
])
tmp.ggplot.1 <- ggplot(data=tmp.data[cell.subset=='CD4'], aes(x=specificity.class.tag, fill=specificity.class.tag)) +
  geom_bar(width=0.7, color='black', linewidth=1.5, alpha=0.8) +
  facet_wrap(facets=~freq.conf.lvl+tcr.chain, scales='free_y', nrow=5) +
  scale_y_continuous(expand=c(0, NA)) +
  scale_fill_manual(values=ag.cols) +
  labs(
    title='For CD4 T cells, summary across donors',
    x='Antigen specificity', y='Clonotype count',
    caption=tmp.caption
  ) +
  theme_bw() + theme(legend.position='none', axis.text.x=element_text(angle=90))
tmp.ggplot.2 <- ggplot(data=tmp.data[cell.subset=='CD8'], aes(x=specificity.class.tag, fill=specificity.class.tag)) +
  geom_bar(width=0.7, color='black', linewidth=1.5, alpha=0.8) +
  facet_wrap(facets=~freq.conf.lvl+tcr.chain, scales='free_y', nrow=5) +
  scale_y_continuous(expand=c(0, NA)) +
  scale_fill_manual(values=ag.cols) +
  labs(
    title='For CD8 T cells, summary across donors',
    x='Antigen specificity', y='Clonotype count',
    caption=tmp.caption
  ) +
  theme_bw() + theme(legend.position='none', axis.text.x=element_text(angle=90))
# Summary considering donor-specific counts.
tmp.data <- tmp.data[, .(clonotype.count=.N), by=.(freq.conf.lvl, specificity.class.tag, cell.subset, tcr.chain, donor.id)]
tmp.ggplot.3 <- ggplot(data=tmp.data[cell.subset=='CD4'], aes(x=specificity.class.tag, y=clonotype.count)) +
  geom_boxplot(aes(fill=specificity.class.tag), width=0.7, color='black', linewidth=1.5, alpha=0.6, outlier.shape=NA) +
  geom_jitter(width=0.3, size=1.2) +
  facet_wrap(facets=~freq.conf.lvl+tcr.chain, scales='free_y', nrow=5) +
  scale_y_continuous(expand=c(0, NA)) +
  scale_fill_manual(values=ag.cols) +
  labs(
    title='For CD4 T cells, summary considering donor-specific counts',
    x='Antigen specificity', y='Clonotype count',
    caption=tmp.caption
  ) +
  theme_bw() + theme(legend.position='none', axis.text.x=element_text(angle=90))
tmp.ggplot.4 <- ggplot(data=tmp.data[cell.subset=='CD8'], aes(x=specificity.class.tag, y=clonotype.count)) +
  geom_boxplot(aes(fill=specificity.class.tag), width=0.7, color='black', linewidth=1.5, alpha=0.6, outlier.shape=NA) +
  geom_jitter(width=0.3, size=1.2) +
  facet_wrap(facets=~freq.conf.lvl+tcr.chain, scales='free_y', nrow=5) +
  scale_y_continuous(expand=c(0, NA)) +
  scale_fill_manual(values=ag.cols) +
  labs(
    title='For CD8 T cells, summary considering donor-specific counts',
    x='Antigen specificity', y='Clonotype count',
    caption=tmp.caption
  ) +
  theme_bw() + theme(legend.position='none', axis.text.x=element_text(angle=90))
tmp.file.name <- paste0(ag.spc.path, '/FinalAssignments_Summ.pdf')
pdf(file=tmp.file.name, height=14, width=8)
print(tmp.ggplot.1)
print(tmp.ggplot.2)
print(tmp.ggplot.3)
print(tmp.ggplot.4)
dev.off()


cat('\n\n')
### ------------------ Processed sample-based data ------------------ ###
cat('### ------------------ Processed sample-based data ------------------ ###\n')

# ---> Output
# Final column order.
tmp.cols <- c(
    'seq.run', 'sample.prj', 'plate.id', 'full.plate.id', 'seq.lib.name', 'tcr.lib.id',
    'sample.status', 'exp.condition', 'has.rep', 'sample.rep',
    'donor.id', 'cell.subset',
    'peptide.pool', 'pp.cell.subset', 'pool.ag',
    'tcr.id', 'tcr.chain',
    'cdr3.nt.seq', 'cdr3.aa.seq',
    'trv', 'trj',
    'trv.allele', 'trd.allele', 'trj.allele',
    'umi.count', 'read.count'
)
tmp.cols <- c(
  tmp.cols,
  setdiff(x=colnames(spc.data), y=tmp.cols)
)
tmp.check <- c(
  length(x=setdiff((tmp.cols), y=colnames(tcr.data))),
  length(x=setdiff(colnames(tcr.data), y=tmp.cols))
)
tmp.check <- all(tmp.check==0)
if(!tmp.check) stop('Failed to put columns in order for general sample-based table.\n')
tcr.data <- tcr.data[, ..tmp.cols]
# Set order according to TCR ID and output.
setorderv(x=tcr.data, cols='tcr.id')
tmp.file.name <- paste0(reports.path, '/SampleBasedTCRData.csv')
fwrite(file=tmp.file.name, x=tcr.data)


cat('\n\n')
### --------------------- Clonotype-based info ---------------------- ###
cat('### --------------------- Clonotype-based info ---------------------- ###\n')

# ---> Define data for which the process is necessary.
# Subset for TCRs from samples that meet the following criteria:
# The sample does not have a replicate (rationale: there were no replicates for index samples).
# The sample comes from an index donor.
# The peptide pool used for the sample matches its corresponding T-cell subset in the case of CMV.
# Data from the RV peptide pool are fully disregarded.
this.tcr.data <- tcr.data[
    has.rep==FALSE & sample.status=='Index' &
    !(pool.ag=='CMV' & cell.subset!=pp.cell.subset) &
    pool.ag!='RV'
]

# ---> General information.
tmp.data.1 <- unique(this.tcr.data[,
    .(
        tcr.chain,
        cdr3.nt.seq, cdr3.aa.seq,
        trv, trj,
        trv.allele, trj.allele, trd.allele
    ),
    by=tcr.id
])
tmp.check <- tmp.data.1[, .N]==tmp.data.1[, uniqueN(tcr.id)]
if(!tmp.check) stop('Unexpected error while compiling clonotype-based info.\n')

# ---> Identify potentially public TCRs.
tmp.data.2 <- this.tcr.data[,
    .(
        public.status=ifelse(
            test=uniqueN(donor.id)>1,
            yes='likely public',
            no='private'
        )
    ),
    by=.(tcr.id)
]

# ---> Total frequencies.
tmp.data.3 <- this.tcr.data[,
    .(
      total.abs.freq=sum(umi.count)
    ),
    by=.(tcr.id, cell.subset, donor.id)
]
tmp.data.4 <- tmp.data.3[, .(total=sum(total.abs.freq)), by=.(cell.subset, donor.id)]
tmp.data.3 <- merge(
  x=tmp.data.3, y=tmp.data.4,
  by=c('cell.subset', 'donor.id')
)
tmp.data.3[, total.rel.freq:=total.abs.freq/total]
tmp.data.3[, total:=NULL]

# ---> Specificity per clonotype and donor.
tmp.data.4 <- unique(this.tcr.data[,
    .(
        gen.clonotype.status,
        frd.1, fad.1, frd.2, fad.2, spc.fc,
        specificity.class.tag, specificity.raw.tag, specificity.ags.tag,
        freq.conf.lvl
    ),
    by=.(tcr.id, cell.subset, donor.id)
])
tmp.check <- tmp.data.4[, .N]==tmp.data.4[, unique(length(paste0(tcr.id, donor.id, sep=';')))]
if(!tmp.check) stop('Unexpected error\n')

# ---> Absolute counts per peptide pool.
tmp.data.5 <- this.tcr.data[,
    .(size=sum(umi.count)),
    by=.(tcr.id, cell.subset, donor.id, pool.ag)
]
tmp.data.5[, pool.ag:=paste0('abs.freq.', pool.ag)]
tmp.data.5 <- as.data.table(spread(data=tmp.data.5, key=pool.ag, value=size, fill=0))

# ---> Relative counts per peptide pool.
tmp.data.6.1 <- as.data.table(gather(data=tmp.data.5, key='pool.ag', value='size', -`donor.id`, -`cell.subset`, -`tcr.id`))
tmp.data.6.2 <- tmp.data.6.1[, .(total.freq=sum(size)), by=.(donor.id, cell.subset, pool.ag)]
tmp.data.6 <- merge(x=tmp.data.6.1, y=tmp.data.6.2, by=c('donor.id', 'cell.subset', 'pool.ag'), all.x=TRUE)
tmp.data.6[, r.val:=size/total.freq]
tmp.data.6[, `:=`(total.freq=NULL, size=NULL)]
tmp.data.6 <- as.data.table(spread(data=tmp.data.6, key=pool.ag, value=r.val, fill=0))
colnames(tmp.data.6) <- str_replace(string=colnames(tmp.data.6), pattern='^abs.freq', replacement='rel.freq')

# Merge all data.
clone.meta <- merge(x=tmp.data.1, y=tmp.data.2, by='tcr.id', all=TRUE)
tmp.data <- merge(x=tmp.data.3, y=tmp.data.4, by=c('tcr.id', 'cell.subset', 'donor.id'), all=TRUE)
clone.meta <- merge(x=clone.meta, y=tmp.data, by='tcr.id')
tmp.data <- merge(x=tmp.data.6, y=tmp.data.5, by=c('tcr.id', 'cell.subset', 'donor.id'), all=TRUE)
clone.meta <- merge(x=clone.meta, y=tmp.data, by=c('tcr.id', 'cell.subset', 'donor.id'), all=TRUE)
# Sanity check
to.check <- clone.meta[, .(tcr.id, cell.subset, donor.id)]
to.check.1 <- to.check[, .N]
to.check.2 <- unique(to.check); to.check.2 <- to.check.2[, .N]
tmp.check <- to.check.1==to.check.2
if(!tmp.check) stop('Unexpected error. There should be as many lines as there are unique pairs of clonotypes and donor IDs.\n')

# Final format.
setorderv(x=clone.meta, col=c('cell.subset', 'tcr.chain', 'total.rel.freq'), order=c(1, 1, -1))
# @ Whole dataset.
tmp.file.name <- paste0(reports.path, '/ClonotypeBasedTCRData.csv')
fwrite(file=tmp.file.name, x=clone.meta, quote=TRUE, na=NA)


cat('\n\n')
### ---------------------- General exploration ---------------------- ###
cat('### ---------------------- General exploration ---------------------- ###\n')

exp.path <- paste0(reports.path, '/gen_exploration')
if(!dir.exists(exp.path)) dir.create(exp.path)

# ---> Control donor definition.
ctrl.donors <- tcr.data[sample.status=='Control', unique(donor.id)]

# ---> On lung samples, basic assessments.
# @ Number of donors covered per peptide pool per chain.
tmp.data <- tcr.data[
  !donor.id %in% ctrl.donors,
  .(donor.count=uniqueN(donor.id)),
  by=.(tcr.chain, peptide.pool)
]
tmp.ggplot.1 <- ggplot(data=tmp.data, aes(x=peptide.pool, y=donor.count, fill=peptide.pool)) +
  geom_bar(stat='identity', width=0.6, color='black', linewidth=1) +
  facet_wrap(facets=~tcr.chain, ncol=1, scales='free_y') +
  scale_fill_manual(values=pp.cols) +
  scale_y_continuous(expand=c(0, 0)) +
  labs(title='Donors w/ available samples per peptide pool', x='Peptide pool', y='Count of unique donors with available samples') +
  theme_bw() +
  theme(legend.position='none', axis.text.x=element_text(angle=-45))
# @ Distribution of number of CDR3 NUCLEOTIDE and AMINOACID sequences per donor per peptide pool per TCR chain.
tmp.data <- tcr.data[
  !donor.id %in% ctrl.donors,
  .(
    cdr3.nt.count=uniqueN(cdr3.nt.seq),
    cdr3.aa.count=uniqueN(cdr3.aa.seq)
  ),
  by=.(tcr.chain, peptide.pool, donor.id, sample.rep)
]
tmp.ggplot.2 <- ggplot(data=tmp.data, aes(x=peptide.pool, y=cdr3.nt.count, fill=peptide.pool)) +
  geom_boxplot(alpha=0.6, width=0.6, outlier.shape=NA) +
  geom_jitter(size=1, width=0.3) +
  facet_wrap(facets=~tcr.chain, ncol=1, scales='free_y') +
  scale_fill_manual(values=pp.cols) +
  scale_y_continuous(expand=c(0, 0)) +
  labs(title='Number of CDR3 NUCLEOTIDE sequences per donor per peptide pool', x='Peptide pool', y='Count of unique CDR3 NUCLEOTIDE sequences', caption='Every dot represents a different experimental sample.') +
  theme_bw() +
  theme(legend.position='none', axis.text.x=element_text(angle=-45))
tmp.ggplot.3 <- ggplot(data=tmp.data, aes(x=peptide.pool, y=cdr3.aa.count, fill=peptide.pool)) +
  geom_boxplot(alpha=0.6, width=0.6, outlier.shape=NA) +
  geom_jitter(size=1, width=0.3) +
  facet_wrap(facets=~tcr.chain, ncol=1, scales='free_y') +
  scale_fill_manual(values=pp.cols) +
  scale_y_continuous(expand=c(0, 0)) +
  labs(title='Number of CDR3 AMINOACID sequences per donor per peptide pool', x='Peptide pool', y='Count of unique CDR3 AMINOACID sequences', caption='Every dot represents a different experimental sample.') +
  theme_bw() +
  theme(legend.position='none', axis.text.x=element_text(angle=-45))
# @ Output
tmp.file.name <- paste0(exp.path, '/GenReports.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot.1)
print(tmp.ggplot.2)
print(tmp.ggplot.3)
dev.off()


# ---> On lung samples, on clonal expansion.
tmp.data <- tcr.data[
  !donor.id %in% ctrl.donors
]
# General UMI count distribution (across all samples)
tmp.ggplot.1 <- ggplot(data=tmp.data, aes(x=umi.count)) +
  geom_histogram(bins=100, alpha=0, linewidth=0.5, color='black') +
  facet_wrap(facets=~cell.subset, ncol=1) +
  scale_y_continuous(expand=c(0, 0)) +
  labs(title='General UMI count distribution (across all samples)', subtitle='No count limit.', x='UMI count', y='Within-sample CDR3 frequency') +
  theme_bw()
tmp.data[, tmp.val:=umi.count]
tmp.pctl <- 0.75
tmp.data.2 <- tmp.data[, .(val.limit=quantile(x=tmp.val, probs=tmp.pctl)), by=.(cell.subset)]
tmp.data <- merge(x=tmp.data, y=tmp.data.2, by='cell.subset')
tmp.data[tmp.val>val.limit, tmp.val:=val.limit]
tmp.caption <- paste0(tmp.data.2[, paste(cell.subset, val.limit, sep=': ')], collapse=', ')
tmp.ggplot.2 <- ggplot(data=tmp.data, aes(x=tmp.val)) +
  geom_histogram(bins=10, alpha=0, linewidth=0.5, color='black') +
  facet_wrap(facets=~cell.subset, ncol=1) +
  scale_y_continuous(expand=c(0, 0)) +
  labs(title='General UMI count distribution (across all samples)', subtitle='With percentile limit.', x='UMI count', y='Within-sample CDR3 frequency', caption=paste0('Bars beyond the ', tmp.pctl*100, ' percentile (', tmp.caption, ') were stacked into the last bar.')) +
  theme_bw()
# Output all UMI count distributions.
tmp.file.name <- paste0(exp.path, '/UMICountDists.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot.1)
print(tmp.ggplot.2)
tmp.pps <- tmp.data[, gtools::mixedsort(unique(peptide.pool))]
for(tmp.pp in tmp.pps){
    to.plot <- tmp.data[peptide.pool==tmp.pp]
    to.plot[, `:=`(tmp.val=umi.count, val.limit=NULL)]
    tmp.pctl <- 0.75
    tmp.data.2 <- to.plot[, .(val.limit=quantile(x=tmp.val, probs=tmp.pctl)), by=.(cell.subset)]
    to.plot <- merge(x=to.plot, y=tmp.data.2, by='cell.subset')
    to.plot[tmp.val>val.limit, tmp.val:=val.limit]
    tmp.caption <- paste0(tmp.data.2[, paste(cell.subset, val.limit, sep=': ')], collapse=', ')
    tmp.ggplot.3 <- ggplot(data=to.plot, aes(x=tmp.val)) +
        geom_histogram(bins=10, alpha=0.5, linewidth=0.5, color='black', fill=pp.cols[tmp.pp]) +
        facet_wrap(facets=~cell.subset, ncol=1) +
        scale_y_continuous(expand=c(0, 0)) +
        labs(title=paste0('UMI count distribution across ', tmp.pp, ' samples'), subtitle='With percentile limit.', x='UMI count', y='Within-sample CDR3 frequency', caption=paste0('Bars beyond the ', tmp.pctl*100, ' percentile (', tmp.caption, ') were stacked into the last bar.')) +
        theme_bw()
    print(tmp.ggplot.3)
}
dev.off()


cat('\n\n')
### -------------- Assay validation based on controls --------------- ###
cat('### -------------- Assay validation based on controls --------------- ###\n')

assay.val.path <- paste0(reports.path, '/assay_validation')
if(!dir.exists(assay.val.path)) dir.create(assay.val.path)

# Define sample ID
tmp.data <- tcr.data[
  sample.status=='Control'
]
tmp.data[, sample.id:=paste(donor.id, cell.subset, peptide.pool, sample.rep, sep='|')]
# Calculate jaccard indexes in a pairwise manner either or not accounting for clone size, per TCR chain and per type of CDR3 sequence molecule.
# molecule.vals <- c(`NT`='cdr3.nt.seq', `AA`='cdr3.aa.seq')
molecule.vals <- c(`NT`='cdr3.nt.seq')
chain.vals <- c('TCRA', 'TCRB')
to.plot <- lapply(X=molecule.vals, FUN=function(mol.col){
  cat(paste0(mol.col), '\n')
  tmp.data <- lapply(X=chain.vals, FUN=function(chain.type){
    cat(paste0(chain.type), '\n')
    tmp.data <- tmp.data[tcr.chain==chain.type]
    tmp.vals <- tmp.data[, unique(sample.id)]
    tmp.data <- lapply(X=tmp.vals, FUN=function(sample.id.1){
      # cat(paste0(sample.id.1), '\n')
      tmp.data <- lapply(X=tmp.vals, FUN=function(sample.id.2){
        sample.1.data <- tmp.data[sample.id==sample.id.1]
        sample.2.data <- tmp.data[sample.id==sample.id.2]
        # # Jaccard index when not accounting for clone size.
        # sample.1.seqs.1 <- sample.1.data[, get(mol.col)]
        # sample.2.seqs.1 <- sample.2.data[, get(mol.col)]
        # # Jaccard index when accounting for clone size.
        # sample.1.seqs.2 <- unlist(lapply(X=1:length(sample.1.seqs.1), FUN=function(idx){
        #   paste(sample.1.seqs.1[idx], 1:sample.1.data[idx, umi.count], sep='.')
        # }))
        # sample.2.seqs.2 <- unlist(lapply(X=1:length(sample.2.seqs.1), FUN=function(idx){
        #   paste(sample.2.seqs.1[idx], 1:sample.2.data[idx, umi.count], sep='.')
        # }))
        # Cosine similarity
        tmp.cols <- c(mol.col, 'umi.count')
        sample.1.seqs.3 <- sample.1.data[, ..tmp.cols]
        colnames(sample.1.seqs.3)[1] <- c('mol.col')
        sample.1.seqs.3 <- sample.1.seqs.3[, .(umi.count=sum(umi.count)), by=mol.col]
        sample.2.seqs.3 <- sample.2.data[, ..tmp.cols]
        colnames(sample.2.seqs.3)[1] <- c('mol.col')
        sample.2.seqs.3 <- sample.2.seqs.3[, .(umi.count=sum(umi.count)), by=mol.col]
        sample.seqs.3 <- merge(x=sample.1.seqs.3, y=sample.2.seqs.3, by='mol.col', all=TRUE)
        sample.1.seqs.3 <- sample.seqs.3[, umi.count.x]
        sample.1.seqs.3[is.na(sample.1.seqs.3)] <- 0
        sample.2.seqs.3 <- sample.seqs.3[, umi.count.y]
        sample.2.seqs.3[is.na(sample.2.seqs.3)] <- 0
        # Calculate indexes and return.
        cosine.idx <- lsa::cosine(x=sample.1.seqs.3, y=sample.2.seqs.3)[1, 1]
        if(length(cosine.idx)>1){
          warning('Issue found.\n')
          cosine.idx <- cosine.idx[1]
        }
        tmp.data <- data.table(
          # raw.jaccard.idx=get.jaccard(x=sample.1.seqs.1, y=sample.2.seqs.1),
          # nor.jaccard.idx=get.jaccard(x=sample.1.seqs.2, y=sample.2.seqs.2),
          cosine.idx=cosine.idx
        )
        return(tmp.data)
      })
      names(tmp.data) <- tmp.vals
      tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='sample.id.2')
      return(tmp.data)
    })
    names(tmp.data) <- tmp.vals
    tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='sample.id.1')
    return(tmp.data)
  })
  names(tmp.data) <- chain.vals
  tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='tcr.chain')
  cat('End of process for this chain.\n')
  cat('\n\n')
  return(tmp.data)
})
names(to.plot) <- names(molecule.vals)
to.plot <- rbindlist(l=to.plot, use.names=TRUE, idcol='mol.type')
# Output heatmap per type of similarity index, molecule and chain.
gen.tmp.cols <- c('sample.id.1', 'sample.id.2')
idx.types <- c(
  # 'UnweightJaccard'='raw.jaccard.idx',
  # 'WeightedJaccard'='nor.jaccard.idx',
  'Cosine'='cosine.idx'
)
for(idx.type in names(idx.types)){
  tmp.cols <- c(gen.tmp.cols, idx.types[idx.type])
  for(molecule.type in names(molecule.vals)){
    for(chain.type in chain.vals){
      # Matrix to plot.
      to.plot.tmp <- to.plot[
        mol.type==molecule.type & tcr.chain==chain.type,
        ..tmp.cols
      ]
      to.plot.tmp <- spread(data=to.plot.tmp, key=sample.id.2, value=idx.types[idx.type])
      tmp.row.names <- to.plot.tmp$sample.id.1
      to.plot.tmp$sample.id.1 <- NULL
      to.plot.tmp <- as.matrix(to.plot.tmp)
      row.names(to.plot.tmp) <- tmp.row.names
      to.plot.tmp <- to.plot.tmp[tmp.row.names, tmp.row.names]
      diag(to.plot.tmp) <- NA
      # Set metadata.
      col.meta <- data.frame(
        row.names=colnames(to.plot.tmp),
        others=colnames(to.plot.tmp)
      )
      col.meta <- separate(data=col.meta, col='others', into=c('Donor', 'Cell type', 'PP', 'Replicate'), sep='\\|')
      # Define colors for metadata tracks.
      ann.colors <- list(
        `Cell type`=cell.type.cols,
        `PP`=pp.cols[names(pp.cols)%in%col.meta[, 'PP']],
        `Replicate`=rep.cols,
        `Donor`=c(`SDBB141`='#ffc0cb', `SDBB148`='#ffff00')
      )
      # Set color scale and breaks for heatmap.
      col.breaks <- seq(from=min(range(to.plot.tmp, na.rm=TRUE)), to=max(range(to.plot.tmp, na.rm=TRUE)), length.out=100)
      mid.point <- which.min(abs(col.breaks - 0))
      hmap.col.scale.1 <- colorRampPalette(c('#ffffff', '#fffff6', '#ffffed', '#ffffb2'))(mid.point)
      hmap.col.scale.2 <- colorRampPalette(c('#ffff77', '#ffff00', '#ffffb2', '#ffff77', '#ffd280', '#ffd280', '#ffbc40', '#ffa500', '#ff5300', '#ff0000', '#8b0000', '#3f0000'))(100-(mid.point+1))
      hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
      # Plot.
      tmp.file.name <- paste0(assay.val.path, '/ControlSamples_Index', idx.type, '_CDR3-', molecule.type, '_', chain.type, '.pdf')
      pheatmap(
        mat=to.plot.tmp, 
        color=hmap.col.scale, breaks=col.breaks, scale='none',
        cluster_rows=FALSE, cluster_cols=FALSE,
        annotation_col=col.meta, annotation_row=col.meta, 
        annotation_colors=ann.colors,
        show_colnames=FALSE, show_rownames=FALSE,
        filename=tmp.file.name, width=25, height=22
      )
      # Raw data.
      tmp.file.name <- paste0(assay.val.path, '/ControlSamples_Index', idx.type, '_CDR3-', molecule.type, '_', chain.type, '.csv')
      write.csv(file=tmp.file.name, x=to.plot.tmp, row.names=TRUE)
    }
  }
}
tmp.file.name <- paste0(assay.val.path, '/__GenIdxInfo.csv')
fwrite(x=to.plot, file=tmp.file.name)


cat('\n\n')
### ---------------- Sample misassignment correction ---------------- ###
### ----------------------- External strategy ----------------------- ###
cat('### ---------------- Sample misassignment correction ---------------- ###\n')
cat('### ----------------------- External strategy ----------------------- ###\n')

correction.1.path <- paste0(reports.path, '/misassignment_correction_1')
if(!dir.exists(correction.1.path)) dir.create(correction.1.path)

# Identify donors for whom some samples are indicated to have replicates and discard from analysis.
discard.donors <- bc.info[has.rep==TRUE, unique(donor.id)]
# Define sample ID
tmp.data <- tcr.data[
  !donor.id %in% discard.donors
]
tmp.data[, sample.id:=paste(donor.id, peptide.pool, sep='|')]
# Calculate similarity indexes in a pairwise manner either or not accounting for clone size, per TCR chain and per type of CDR3 sequence molecule.
mol.1.vals <- c(`NT`='cdr3.nt.seq', `AA`='cdr3.aa.seq')
mol.2.vals <- list(
  'TCRA'=c(`cdr3.nt.seq`='cdr3a.nt.seq', `cdr3.aa.seq`='cdr3a.aa.seq'),
  'TCRB'=c(`cdr3.nt.seq`='cdr3b.nt.seq', `cdr3.aa.seq`='cdr3b.aa.seq')
)
cell.subsets <- tmp.data[, unique(cell.subset)]
chain.vals <- c('TCRA', 'TCRB')
tmp.data <- lapply(X=mol.1.vals, FUN=function(mol.1.col){
  cat(paste0(mol.1.col), '\n')
  tmp.data <- lapply(X=cell.subsets, FUN=function(subset.id){
    tmp.data <- tmp.data[cell.subset==subset.id]
    cat(paste0(subset.id), '\n')
    tmp.data <- lapply(X=chain.vals, FUN=function(chain.type){
        cat(paste0(chain.type), '\n')
        mol.2.col <- mol.2.vals[[chain.type]][mol.1.col]
        tmp.data <- tmp.data[tcr.chain==chain.type]
        tmp.donors <- tmp.data[, unique(donor.id)]
        tmp.data <- lapply(X=tmp.donors, FUN=function(donor.id.1){
          cat(paste0(donor.id.1), '\n')
          tmp.data <- tmp.data[donor.id==donor.id.1]
          tmp.vals <- tmp.data[, unique(sample.id)]
          tmp.data <- lapply(X=tmp.vals, FUN=function(tmp.sample.id){
            tmp.data <- mclapply(X=tmp.donors, mc.cores=cpu.no, FUN=function(donor.id.2){
              sample.1.data <- tmp.data[sample.id==tmp.sample.id]
              sample.2.data <- lung.data[!is.na(get(mol.2.col)) & donor.id.tag==donor.id.2]
              # # Jaccard index when not accounting for clone size.
              # sample.1.seqs.1 <- sample.1.data[, get(mol.1.col)]
              # sample.2.seqs.1 <- sample.2.data[, get(mol.2.col)]
              # # Jaccard index when accounting for clone size.
              # sample.1.seqs.2 <- unlist(lapply(X=1:length(sample.1.seqs.1), FUN=function(idx){
              #   paste(sample.1.seqs.1[idx], 1:sample.1.data[idx, umi.count], sep='.')
              # }))
              # sample.2.seqs.2 <- unlist(lapply(X=1:nrow(sample.2.data), FUN=function(idx){
              #   paste(sample.2.seqs.1[idx], 1:sample.2.data[idx, abs.freq], sep='.')
              # }))
              # Cosine similarity
              tmp.cols <- c(mol.1.col, 'umi.count')
              sample.1.seqs.3 <- sample.1.data[, ..tmp.cols]
              colnames(sample.1.seqs.3) <- c('mol.col', 'freq')
              sample.1.seqs.3 <- sample.1.seqs.3[, .(freq=sum(freq)), by=mol.col]
              tmp.cols <- c(mol.2.col, 'abs.freq')
              sample.2.seqs.3 <- sample.2.data[, ..tmp.cols]
              colnames(sample.2.seqs.3) <- c('mol.col', 'freq')
              sample.2.seqs.3 <- sample.2.seqs.3[, .(freq=sum(freq)), by=mol.col]
              sample.seqs.3 <- merge(x=sample.1.seqs.3, y=sample.2.seqs.3, by='mol.col', all=TRUE)
              sample.1.seqs.3 <- sample.seqs.3[, freq.x]
              sample.1.seqs.3[is.na(sample.1.seqs.3)] <- 0
              sample.2.seqs.3 <- sample.seqs.3[, freq.y]
              sample.2.seqs.3[is.na(sample.2.seqs.3)] <- 0
              # Intersection size.
              sample.1.seqs.4 <- sample.1.data[, get(mol.1.col)]
              sample.2.seqs.4 <- sample.2.data[, get(mol.2.col)]
              # Calculate indexes and return.
              cosine.idx <- lsa::cosine(x=sample.1.seqs.3, y=sample.2.seqs.3)
              if(length(cosine.idx)>1){
                warning('Issue found.\n')
                cosine.idx <- cosine.idx[1]
              }
              cosine.idx <- cosine.idx[1, 1]
              tmp.data <- data.table(
                # raw.jaccard.idx=get.jaccard(x=sample.1.seqs.1, y=sample.2.seqs.1),
                # nor.jaccard.idx=get.jaccard(x=sample.1.seqs.2, y=sample.2.seqs.2),
                raw.jaccard.idx=NA,
                nor.jaccard.idx=NA,
                cosine.idx=cosine.idx,
                int.size=length(intersect(x=sample.1.seqs.4, y=sample.2.seqs.4))
              )
              return(tmp.data)
            })
            names(tmp.data) <- tmp.donors
            tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='donor.id.2')
            return(tmp.data)
          })
          names(tmp.data) <- tmp.vals
          tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='sample.id')
          return(tmp.data)
        })
        names(tmp.data) <- tmp.donors
        tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='donor.id.1')
        return(tmp.data)
    })
    names(tmp.data) <- chain.vals
    tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='tcr.chain')
    return(tmp.data)
  })
  names(tmp.data) <- cell.subsets
  tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='cell.subset')
  return(tmp.data)
})
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='mol.type')
# Define peptide pool condition.
tmp.data[, peptide.pool:=str_replace(
  string=sample.id,
  pattern='^[^|]+\\|',
  replacement=''
)]
# Include plate info.
tmp.data.1 <- unique(bc.info[, .(
  sample.id=paste(donor.id, peptide.pool, sep='|'),
  cell.subset, tcr.chain,
  full.plate.id=paste(seq.run, plate.id, sep='-')
)])
tmp.data.1 <- tmp.data.1[,
  .(full.plate.id=paste0(full.plate.id, collapse='|')),
  by=.(cell.subset, tcr.chain, sample.id)
]
tmp.data <- merge(x=tmp.data, y=tmp.data.1, by=c('sample.id', 'cell.subset', 'tcr.chain'))
# Output heatmap per type of similarity index, molecule and chain.
summ.idxs <- c(
  `CosineIdx`='cosine.idx',
  `IntSize`='int.size'
)
for(this.summ.idx in names(summ.idxs)){
  # Define similarity index of choice.
  tmp.data[, summ.idx:=get(summ.idxs[this.summ.idx])]
  donor.ids <- tmp.data[, unique(donor.id.1)]
  for(donor.id in donor.ids){
    to.plot <- tmp.data[
      mol.type=='NT' & donor.id.1==donor.id,
      .(full.plate.id, cell.subset, tcr.chain, peptide.pool, donor.id=donor.id.2, summ.idx)
    ]
    to.plot[, tmp.id:=paste(cell.subset, tcr.chain, peptide.pool, sep='-')]
    to.plot.tmp <- spread(data=to.plot[, .(tmp.id, donor.id, summ.idx)], key=donor.id, value=summ.idx)
    tmp.row.names <- to.plot.tmp$tmp.id
    to.plot.tmp$tmp.id <- NULL
    to.plot.tmp <- as.matrix(to.plot.tmp)
    row.names(to.plot.tmp) <- tmp.row.names
    if(all(to.plot.tmp==0) | is.null(dim(to.plot.tmp))) next
    # Set color scale and breaks for heatmap.
    col.breaks <- seq(from=min(range(to.plot.tmp, na.rm=TRUE)), to=max(range(to.plot.tmp, na.rm=TRUE)), length.out=100)
    mid.point <- which.min(abs(col.breaks - 0))
    hmap.col.scale.1 <- colorRampPalette(c('#ffffff', '#fffff6', '#ffffed', '#ffffb2'))(mid.point)
    hmap.col.scale.2 <- colorRampPalette(c('#ffff77', '#ffff00', '#ffffb2', '#ffff77', '#ffd280', '#ffd280', '#ffbc40', '#ffa500', '#ff5300', '#ff0000', '#8b0000', '#3f0000'))(100-(mid.point+1))
    hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
    # Row metadata.
    row.meta <- as.data.frame(unique(to.plot[, .(row.names=tmp.id, Lineage=cell.subset, Chain=tcr.chain, PP=peptide.pool, Plate=full.plate.id)]))
    row.names(row.meta) <- row.meta$row.names; row.meta$row.names <- NULL
    tmp.cols <- plate.cols[names(plate.cols)%in%row.meta[, 'Plate']]
    row.cols <- list(
      `PP`=pp.cols,
      `Lineage`=cell.type.cols,
      `Chain`=chain.cols,
      `Plate`=tmp.cols
    )
    # Plot.
    show.nos <- this.summ.idx=='IntSize'
    tmp.file.name <- paste0(correction.1.path, '/', this.summ.idx, '_', donor.id, '.pdf')
    pheatmap(
        mat=to.plot.tmp, 
        color=hmap.col.scale, breaks=col.breaks, scale='none',
        cluster_rows=FALSE, cluster_cols=FALSE,
        annotation_col=NULL, annotation_row=row.meta, 
        annotation_colors=row.cols,
        display_numbers=show.nos, number_format = "%.0f",
        show_colnames=TRUE, show_rownames=FALSE,
        filename=tmp.file.name, width=17, height=6
    )
    if(this.summ.idx=='IntSize'){
      tmp.file.name <- paste0(correction.1.path, '/', this.summ.idx, '_', donor.id, '.csv')
      write.csv(x=to.plot.tmp, file=tmp.file.name, quote=TRUE)
    }
  }
}
tmp.data[, summ.idx:=NULL]
tmp.file.name <- paste0(correction.1.path, '/__GenIdxInfo.csv')
fwrite(x=tmp.data, file=tmp.file.name)


cat('\n\n')
### ---------------- Sample misassignment correction ---------------- ###
### ----------------------- Internal strategy ----------------------- ###
cat('### ---------------- Sample misassignment correction ---------------- ###\n')
cat('### ----------------------- Internal strategy ----------------------- ###\n')

correction.2.path <- paste0(reports.path, '/misassignment_correction_2')
if(!dir.exists(correction.2.path)) dir.create(correction.2.path)

# Identify donors for whom some samples are indicated to have replicates and discard from analysis.
discard.donors <- bc.info[has.rep==TRUE, unique(donor.id)]
# Define sample ID
tmp.data <- tcr.data[
  !donor.id %in% discard.donors
]
tmp.data[, sample.id:=paste(donor.id, peptide.pool, sep='|')]
# Calculate jaccard indexes in a pairwise manner either or not accounting for clone size, per TCR chain and per type of CDR3 sequence molecule.
seq.runs <- tmp.data[, unique(seq.run)]
# molecule.vals <- c(`NT`='cdr3.nt.seq', `AA`='cdr3.aa.seq')
molecule.vals <- c(`NT`='cdr3.nt.seq')
cell.subsets <- tmp.data[, unique(cell.subset)]
chain.vals <- c('TCRA', 'TCRB')
to.plot <- lapply(X=seq.runs, FUN=function(run.id){
  cat(paste0(run.id), '\n')
  tmp.data <- tmp.data[seq.run==run.id]
  tmp.data <- lapply(X=molecule.vals, FUN=function(mol.col){
    cat(paste0(mol.col), '\n')
    tmp.data <- lapply(X=cell.subsets, FUN=function(subset.id){
      tmp.data <- tmp.data[cell.subset==subset.id]
      cat(paste0(subset.id), '\n')
      tmp.data <- lapply(X=chain.vals, FUN=function(chain.type){
        cat(paste0(chain.type), '\n')
        tmp.data <- tmp.data[tcr.chain==chain.type]
        tmp.vals <- tmp.data[, unique(sample.id)]
        tmp.data <- lapply(X=tmp.vals, FUN=function(sample.id.1){
          # cat(paste0(sample.id.1), '\n')
          tmp.data <- lapply(X=tmp.vals, FUN=function(sample.id.2){
            sample.1.data <- tmp.data[sample.id==sample.id.1]
            sample.2.data <- tmp.data[sample.id==sample.id.2]
            # Jaccard index when not accounting for clone size.
            sample.1.seqs.1 <- sample.1.data[, get(mol.col)]
            sample.2.seqs.1 <- sample.2.data[, get(mol.col)]
            # Jaccard index when accounting for clone size.
            sample.1.seqs.2 <- unlist(lapply(X=1:length(sample.1.seqs.1), FUN=function(idx){
              paste(sample.1.seqs.1[idx], 1:sample.1.data[idx, umi.count], sep='.')
            }))
            sample.2.seqs.2 <- unlist(lapply(X=1:length(sample.2.seqs.1), FUN=function(idx){
              paste(sample.2.seqs.1[idx], 1:sample.2.data[idx, umi.count], sep='.')
            }))
            # Cosine similarity
            tmp.cols <- c(mol.col, 'umi.count')
            sample.1.seqs.3 <- sample.1.data[, ..tmp.cols]
            colnames(sample.1.seqs.3)[1] <- c('mol.col')
            sample.1.seqs.3 <- sample.1.seqs.3[, .(umi.count=sum(umi.count)), by=mol.col]
            sample.2.seqs.3 <- sample.2.data[, ..tmp.cols]
            colnames(sample.2.seqs.3)[1] <- c('mol.col')
            sample.2.seqs.3 <- sample.2.seqs.3[, .(umi.count=sum(umi.count)), by=mol.col]
            sample.seqs.3 <- merge(x=sample.1.seqs.3, y=sample.2.seqs.3, by='mol.col', all=TRUE)
            sample.1.seqs.3 <- sample.seqs.3[, umi.count.x]
            sample.1.seqs.3[is.na(sample.1.seqs.3)] <- 0
            sample.2.seqs.3 <- sample.seqs.3[, umi.count.y]
            sample.2.seqs.3[is.na(sample.2.seqs.3)] <- 0
            # Calculate indexes and return.
            cosine.idx <- lsa::cosine(x=sample.1.seqs.3, y=sample.2.seqs.3)
            if(length(cosine.idx)>1){
              warning('Issue found.\n')
              cosine.idx <- cosine.idx[1]
            }
            cosine.idx <- cosine.idx[1, 1]
            tmp.data <- data.table(
              raw.jaccard.idx=get.jaccard(x=sample.1.seqs.1, y=sample.2.seqs.1),
              nor.jaccard.idx=get.jaccard(x=sample.1.seqs.2, y=sample.2.seqs.2),
              cosine.idx=cosine.idx
            )
            return(tmp.data)
          })
          names(tmp.data) <- tmp.vals
          tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='sample.id.2')
          return(tmp.data)
        })
        names(tmp.data) <- tmp.vals
        tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='sample.id.1')
        return(tmp.data)
      })
      names(tmp.data) <- chain.vals
      tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='tcr.chain')
      return(tmp.data)
    })
    names(tmp.data) <- cell.subsets
    tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='cell.subset')
    return(tmp.data)
  })
  names(tmp.data) <- names(molecule.vals)
  tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='mol.type')
  cat('End of process for this run.\n')
  cat('\n\n')
  return(tmp.data)
})
names(to.plot) <- seq.runs
to.plot <- rbindlist(l=to.plot, use.names=TRUE, idcol='seq.run')
# Output heatmap per type of similarity index, molecule and chain.
gen.tmp.cols <- c('sample.id.1', 'sample.id.2')
idx.types <- c(
  'UnweightJaccard'='raw.jaccard.idx',
  'WeightedJaccard'='nor.jaccard.idx',
  'Cosine'='cosine.idx'
)
for(idx.type in names(idx.types)){
  tmp.cols <- c(gen.tmp.cols, idx.types[idx.type])
  for(run.id in seq.runs){
    for(molecule.type in names(molecule.vals)){
      for(subset.id in cell.subsets){
        for(chain.type in chain.vals){
          # Matrix to plot.
          to.plot.tmp <- to.plot[
            seq.run==run.id & mol.type==molecule.type & cell.subset==subset.id & tcr.chain==chain.type,
            ..tmp.cols
          ]
          to.plot.tmp <- spread(data=to.plot.tmp, key=sample.id.2, value=idx.types[idx.type])
          tmp.row.names <- to.plot.tmp$sample.id.1
          to.plot.tmp$sample.id.1 <- NULL
          to.plot.tmp <- as.matrix(to.plot.tmp)
          row.names(to.plot.tmp) <- tmp.row.names
          to.plot.tmp <- to.plot.tmp[tmp.row.names, tmp.row.names]
          diag(to.plot.tmp) <- NA
          # Set metadata.
          col.meta <- data.frame(
            row.names=colnames(to.plot.tmp),
            others=colnames(to.plot.tmp)
          )
          col.meta <- separate(data=col.meta, col='others', into=c('Donor', 'PP'), sep='\\|')
          # Define colors for metadata tracks.
          donor.count <- length(unique(col.meta$Donor))
          donor.cols.1 <- RColorBrewer::brewer.pal(n=floor(donor.count/2), name='Paired')
          donor.cols.2 <- RColorBrewer::brewer.pal(n=ceiling(donor.count/2), name='Set3')
          donor.cols <- c(donor.cols.1, donor.cols.2)
          donor.cols <- sample(x=donor.cols, size=donor.count, replace=FALSE)
          names(donor.cols) <- gtools::mixedsort(unique(col.meta$Donor))
          ann.colors <- list(
            `PP`=pp.cols[names(pp.cols)%in%col.meta[, 'PP']],
            `Donor`=donor.cols
          )
          # Set color scale and breaks for heatmap.
          col.breaks <- seq(from=min(range(to.plot.tmp, na.rm=TRUE)), to=max(range(to.plot.tmp, na.rm=TRUE)), length.out=100)
          mid.point <- which.min(abs(col.breaks - 0))
          hmap.col.scale.1 <- colorRampPalette(c('#ffffff', '#fffff6', '#ffffed', '#ffffb2'))(mid.point)
          hmap.col.scale.2 <- colorRampPalette(c('#ffff77', '#ffff00', '#ffffb2', '#ffff77', '#ffd280', '#ffd280', '#ffbc40', '#ffa500', '#ff5300', '#ff0000', '#8b0000', '#3f0000'))(100-(mid.point+1))
          hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
          # Plot.
          tmp.file.name <- paste0(correction.2.path, '/Index', idx.type, '_', run.id, '_CDR3-', molecule.type, '_', subset.id, '_', chain.type, '.pdf')
          pheatmap(
            mat=to.plot.tmp, 
            color=hmap.col.scale, breaks=col.breaks, scale='none',
            cluster_rows=FALSE, cluster_cols=FALSE,
            annotation_col=col.meta, annotation_row=col.meta, 
            annotation_colors=ann.colors,
            show_colnames=FALSE, show_rownames=FALSE,
            filename=tmp.file.name, width=25, height=22
          )
        }
      }
    }
  }
}
tmp.file.name <- paste0(correction.2.path, '/__GenIdxInfo.csv')
fwrite(x=to.plot, file=tmp.file.name)


cat('\n\n')
### --------------------------------- END -------------------------------- ###
cat('### --------------------------------- END -------------------------------- ###\n')
# ---> Print session info.
sessionInfo()

cat('PROGRAM FINISHED!\n\n')
