############    ---   Preliminary process to prepare    ---    ############
############    --------   the DICE TCR database    --------    ############

# By: Vicente Fajardo

### -------------------------- Description -------------------------- ###
# Perform all preliminary processes required to obtain clean tables of the DICE TCR database data.

### -------------------------- Dependencies ------------------------- ###
library(optparse)
library(parallel)
library(data.table)
library(tidyr)
library(stringr)
library(ggplot2)
library(UpSetR)
library(VennDiagram)
library(pheatmap)
# library(rstatix)
# library(Hmisc)
# library(corrplot)
# source('/home/vfajardo/scripts/functions/R_handy_functions.0.4.R')
source('/path/to/user/R24/paper_developments/R24_Cancer/paper_items/jobs_scripts/functions_collection.0.4.R')


### --------------------------- Arguments --------------------------- ###

# ---> Arguments from command line
option.list <- list(
    # For output.
    make_option(opt_str="--BatchesLab", type="character", default=NULL, dest="batches.lab", help="Character, prefix label for aggregation.\n"),
    make_option(opt_str="--SeqDate", type="character", default=NULL, dest="seq.date", help="Character, sequencing date (ID) for folder where aggr results were stored.\n"),
    make_option(opt_str="--AggrVers", type="character", default=NULL, dest="aggr.vers", help="Character, aggr version.\n"),
    make_option(opt_str="--AggrDate", type="character", default=NULL, dest="aggr.date", help="Character, aggr version.\n")
)
# Getting arguments from command line and setting their values to their respective variable names.
opt.parser <- OptionParser(option_list=option.list)
opt <- parse_args(opt.parser)
# Moving options to their own variables
batches.lab <- opt$batches.lab
seq.date <- opt$seq.date
aggr.vers <- opt$aggr.vers
aggr.date <- opt$aggr.date

# ---> Hardcoded arguments
# @ Paper.
main.prj <- 'R24'
this.prj <- 'DICE-VirusSpec'
this.prj.lab <- 'preliminary_process'
batches.lab <- batches.lab
aggr.sffx <- paste0(batches.lab, '_CD4')
aggr.id <- paste0(this.prj, '_', aggr.sffx)
seq.date <- seq.date
aggr.date <- aggr.date
aggr.vers <- aggr.vers
donor.vers <- '1.0'
core.no <- as.integer(system(command='echo $SLURM_JOB_CPUS_PER_NODE', intern=TRUE))
# @ Reports date.
reports.date <- NULL
# @ Seed.
set.seed(seed=1)
# ---> Color defintions.
# @ Color per peptide pool
peptide.pool.cols <- c(
  'CMV'='#B22222',
  'EBV'='#FF8C00',
  'SARS-CoV-2'='#9932CC',
  'IAV'='#0095B3',
  'RSV'='#59B300',
  'MPV'='#FF5EA2',
  'PIV'='#FFD700',
  'B-pertussis-Vax'='#754600',
  'B-pertussis-Rest'='#C37600',
  'Aspergillus'='#005746',
  'Alternaria'='#319272'
)
# @ Other terms' colors.
signatures.col.scale <- c('#ffffff', '#ffffe0', '#ffffad', '#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
# ---> Path definitions.
prj.gen.path <- paste0('/path/to/user/', main.prj, '/paper_developments/', this.prj)
donor.data.path <- paste0(prj.gen.path, '/donors_metadata')
reports.path <- paste0(prj.gen.path, '/exploratory_analyses/', this.prj.lab)
corr.data.path <- paste0(reports.path, '/data')
if(is.null(reports.date)){
    reports.path <- paste0(reports.path, '/', batches.lab, '_v-', aggr.vers, '_', Sys.Date())
    if(!dir.exists(reports.path)) dir.create(reports.path)
}else{
    reports.path <- paste0(reports.path, '/', batches.lab, '_', reports.date)
    if(!dir.exists(reports.path)) stop('A reports date was provided, but we failed to define an already existing reports path.\n')
}
gen.seq.path <- '/path/to/user/sequencing_data'
aggr.data.path <- paste0(gen.seq.path, '/', seq.date, '/aggr_vdj/data')
aggr.res.path <- paste0(gen.seq.path, '/', seq.date, '/aggr_vdj/', aggr.id)
aggr.res.path <- list.dirs(path=aggr.res.path, recursive=FALSE)
aggr.res.path <- aggr.res.path[str_detect(string=aggr.res.path, pattern=aggr.date)]
tmp.check <- length(aggr.res.path)==1 # Sanity check.
if(!tmp.check) stop('Unexpected error!\n')
# ---> File definitions.
# TCR data.
tcr.clones.file <- paste0(aggr.res.path, '/clonotypes_aggr.csv')
tcr.anns.file <- paste0(aggr.res.path, '/filtered_contig_annotations_aggr.csv')
# Metadata
aggr.table.file <- paste0(aggr.data.path, '/', aggr.id, '_aggr_table_annotated.', aggr.vers, '.csv')
donor.meta.file <- paste0(donor.data.path, '/GeneralMetaDataProcessedTable_For', this.prj, '_', batches.lab, '.', donor.vers, '.csv')
# For sample swap correction.
#       For peptide pool-specific donor sample swap correction.
pd.swap.corr.file <- paste0(corr.data.path, '/DonorPPSwapFix.0.2.csv')
#       For repeat sample swap correction.
rep.swap.corr.file <- paste0(corr.data.path, '/RepeatSwapFix.0.2.csv')
# Previous aggr
prev.meta.file <- '/home/fcastaneda/fcastaneda-temp/R24/GLIPH_per_virus/results/TCR_analysis/4viruses/hto_tcr_cd3norm_info.csv'
# ---> Check directories and files.
if(!dir.exists(reports.path)) dir.create(reports.path)
essential.files <- c(
  tcr.clones.file, tcr.anns.file,
  aggr.table.file, donor.meta.file,
  pd.swap.corr.file, rep.swap.corr.file
)
essential.files <- essential.files[!unlist(lapply(X=essential.files, FUN=file.exists))]
if(length(essential.files) > 0) stop(paste0('Next essential files are not appropriately defined -make sure you\'ve got their names right-:\n', paste0(essential.files, collapse='\n'), '\n'))


### --------------------------- Functions --------------------------- ###

# ---> Function to obtain jaccard index.
get.jaccard <- function(x, y){
  num.val <- length(intersect(x=x, y=y))
  den.val <- length(union(x=x, y=y))
  return(num.val/den.val)
}

# ---> Function to obtain sample pairwise similarity indexes.
get.sim.idxs <- function(sample.freqs, extra.verbose=TRUE){
    # Identity unique sample IDs.
    sample.ids <- sample.freqs[, unique(sample.id)]
    # Iterate over samples in a pairwise fashion.
    tmp.data <- lapply(X=sample.ids, FUN=function(sample.id.1){
        cat(paste0(sample.id.1, '\n'))
        tmp.data <- mclapply(X=sample.ids, mc.cores=core.no, FUN=function(sample.id.2){
            if(extra.verbose) cat(paste0('\t\t', sample.id.2, '\n'))
            s1.data <- sample.freqs[sample.id==sample.id.1]
            s2.data <- sample.freqs[sample.id==sample.id.2]
            # Jaccard index when not accounting for clone size.
            s1.seqs.1 <- s1.data[, clonotype.tag]
            s2.seqs.1 <- s2.data[, clonotype.tag]
            # Cosine similarity
            tmp.cols <- c('clonotype.tag', 'rel.freq')
            s1.seqs.3 <- s1.data[, ..tmp.cols]
            s2.seqs.3 <- s2.data[, ..tmp.cols]
            ss.seqs.3 <- merge(x=s1.seqs.3, y=s2.seqs.3, by='clonotype.tag', all=TRUE)
            s1.seqs.3 <- ss.seqs.3[, rel.freq.x]
            s1.seqs.3[is.na(s1.seqs.3)] <- 0
            s2.seqs.3 <- ss.seqs.3[, rel.freq.y]
            s2.seqs.3[is.na(s2.seqs.3)] <- 0
            # Calculate indexes and return.
            cosine.idx <- lsa::cosine(x=s1.seqs.3, y=s2.seqs.3)
            cosine.idx <- cosine.idx[1, 1]
            tmp.data <- data.table(
                raw.jaccard.idx=get.jaccard(x=s1.seqs.1, y=s2.seqs.1),
                cosine.idx=cosine.idx
            )
            return(tmp.data)
        })
        names(tmp.data) <- sample.ids
        tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='sample.id.2')
        return(tmp.data)
    })
    names(tmp.data) <- sample.ids
    tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='sample.id.1')
    return(tmp.data)
}


### -------------------------- Load data & -------------------------- ###
### ---------------------- Data preprocessing ----------------------- ###

# ---> TCR data.
tcr.clones <- fread(file=tcr.clones.file)
tcr.anns <- fread(file=tcr.anns.file)

# ---> Metadata
aggr.table <- fread(file=aggr.table.file)
donor.meta <- fread(file=donor.meta.file)

# For sample swap correction.
#       For peptide pool-specific donor sample swap correction.
pd.swap.corr <- fread(file=pd.swap.corr.file)
#       For repeat sample swap correction.
rep.swap.corr <- fread(file=rep.swap.corr.file)

# ---> Previous results (for comparison against the new results).
prev.meta <- fread(file=prev.meta.file)


### ------------------------- Main program -------------------------- ###

### ------------------------- Preprocessing ------------------------- ###

# Find all suffixes
all.sffxs <- tcr.anns[, str_extract(string=barcode, pattern='-.+$')]
all.sffxs <- unique(all.sffxs)
all.sffxs <- gtools::mixedsort(str_replace(string=all.sffxs, pattern='^-', replacement=''))

# ---> Barcode-clone associations.
cell.clone.info <- tcr.anns[is_cell & high_confidence & full_length & productive]
cols.to.keep <- c(
    `barcode`='barcode',
    `chain`='chain', `v.gene`='v_gene', `j.gene`='j_gene',
    `cdr3.nt.seq`='cdr3_nt',
    `cdr3.aa.seq`='cdr3',
    `umi.count`='umis',
    `clonotype.tag`='raw_clonotype_id'
)
cell.clone.info <- cell.clone.info[, ..cols.to.keep]
colnames(cell.clone.info) <- names(cols.to.keep)
cell.clone.info <- unique(cell.clone.info)

# Perform only if necessary.
if(is.null(reports.date)){
    # ---> vdj gene usage info.
    # General vdj gene usage info.
    # We obtain unique v and j genes per clonotype. When cellranger has defined multiple v and j genes for a given clonotype, we take the one that is most supported. The most supported is the gene that has the largest amount of entries (contigs) in the table below.
    #       We first split the dataset into clonotypes w/ a single V gene and a single J gene definition from the rest, for which we need to pick one per the criterion established right above.
    tmp.data <- cell.clone.info[,
        .(
            v.gene.count=uniqueN(v.gene),
            j.gene.count=uniqueN(j.gene)
        ),
        by=.(
            clonotype.tag,
            chain, cdr3.nt.seq
        )
    ]
    tmp.clonotypes <- tmp.data[v.gene.count>1 | j.gene.count>1, unique(clonotype.tag)]
    #          For well-defined clonotypes, collect genes.
    vdj.gene.info.1 <- cell.clone.info[
        !clonotype.tag %in% tmp.clonotypes,
        .(
            v.gene=unique(v.gene),
            j.gene=unique(j.gene)
        ),
        by=.(
            clonotype.tag,
            chain, cdr3.nt.seq
        )
    ]
    #          For the rest of the clonotypes, select best gene according to the established criterion.
    vdj.gene.info.2 <- cell.clone.info[
        clonotype.tag %in% tmp.clonotypes,
        .(
            v.gene=.SD[,
                .(support=sum(umi.count)),
                by=v.gene
            ][support==max(support), paste0(unique(v.gene), collapse='|')],
            j.gene=.SD[,
                .(support=sum(umi.count)),
                by=j.gene
            ][support==max(support), paste0(unique(j.gene), collapse='|')]
        ),
        by=.(
            clonotype.tag,
            chain, cdr3.nt.seq
        )
    ]
    #           Out of this second group, for those clonotypes with multiple v and j genes that were equally supported according to our criteorion, make a guess and keep only one gene per set of clonotype, chain and CDR3 nt seq.
    vdj.gene.info.2[str_detect(string=v.gene, pattern='\\|'), v.gene:=str_replace(string=v.gene, pattern='\\|.+$', replacement='')]
    vdj.gene.info.2[str_detect(string=j.gene, pattern='\\|'), j.gene:=str_replace(string=j.gene, pattern='\\|.+$', replacement='')]
    #           Concatenate info from individual approaches.
    vdj.gene.info <- rbindlist(l=list(vdj.gene.info.1, vdj.gene.info.2), use.names=TRUE)
    # Sanity check. Confirm we have unique entries for each combination of clonotype ID, chain and CDR3 sequence.
    tmp.check <- vdj.gene.info[, .N, by=.(clonotype.tag, chain, cdr3.nt.seq)][N>1, .N==0]
    if(!tmp.check) stop('Unexpected error post vdj data retrieval.\n')

    # ---> Unique sequences per clonotype.
    tmp.data <- unique(cell.clone.info[, .(clonotype.tag, chain, cdr3.nt.seq, cdr3.aa.seq)])
    # Add vdj gene info.
    tmp.data <- merge(
        x=tmp.data, y=vdj.gene.info,
        by=c('clonotype.tag', 'chain', 'cdr3.nt.seq'),
        all=TRUE
    )
    tmp.check <- tmp.data[is.na(v.gene) | is.na(j.gene), .N==0]
    if(!tmp.check) stop('Error while defining unique sequences per clonotype (1).\n')
    # Collapse chain-specific info per clonotype. That is, if several, collapse info from multiple CDR3 and V/J gene sequences into single entries while keeping track of the correspondance between CDR3 and gene sequences.
    tmp.data <- tmp.data[,
        .(
            cdr3.nt.seq=paste0(cdr3.nt.seq, collapse=';'),
            cdr3.aa.seq=paste0(cdr3.aa.seq, collapse=';'),
            v.gene=paste0(v.gene, collapse=';'),
            j.gene=paste0(j.gene, collapse=';')
        ),
        by=.(clonotype.tag, chain)
    ]
    # For clonotypes w/ multiple sequences of either chain, confirm that we have complete info for all four relevant pieces: CDR3 nt and aa sequences and V and J genes.
    tmp.check <- tmp.data[,
        .(
            cdr3.nt.seq=str_count(string=cdr3.nt.seq, pattern=';'),
            cdr3.aa.seq=str_count(string=cdr3.aa.seq, pattern=';'),
            v.gene=str_count(string=v.gene, pattern=';'),
            j.gene=str_count(string=j.gene, pattern=';')
        )
    ]
    tmp.check <- apply(X=tmp.check, MARGIN=1, FUN=function(x) length(unique(x)))
    tmp.check <- all(tmp.check==1)
    if(!tmp.check) stop('Error while defining unique sequences per clonotype (2).\n')
    # Spread chain-specific info into multiple columns.
    tmp.cols <- c('cdr3.nt.seq', 'cdr3.aa.seq', 'v.gene', 'j.gene')
    these.col.names <- list(
        `cdr3.nt.seq`=c('TRA'='cdr3a.nt.seq', 'TRB'='cdr3b.nt.seq'),
        `cdr3.aa.seq`=c('TRA'='cdr3a.aa.seq', 'TRB'='cdr3b.aa.seq'),
        `v.gene`=c('TRA'='tra.v', 'TRB'='trb.v'),
        `j.gene`=c('TRA'='tra.j', 'TRB'='trb.j')
    )
    these.col.names <- lapply(X=these.col.names, FUN=function(x) c(c(`clonotype.tag`='clonotype.tag'), x))
    tmp.data <- lapply(X=tmp.cols, FUN=function(tmp.col){
        tmp.cols <- c('clonotype.tag', 'chain', tmp.col)
        tmp.data <- tmp.data[, ..tmp.cols]
        tmp.data <- as.data.table(spread(data=tmp.data, key=chain, value=tmp.col, fill=NA))
        colnames(tmp.data) <- these.col.names[[tmp.col]][colnames(tmp.data)]
        return(tmp.data)
    })
    clone.info <- Reduce(x=tmp.data, f=function(x, y) merge(x=x, y=y, by='clonotype.tag', all=TRUE)); rm(tmp.data)
    # Final sanity check.
    tmp.check <- clone.info[, uniqueN(clonotype.tag)]==cell.clone.info[, uniqueN(clonotype.tag)]
    if(!tmp.check) stop('Error while defining unique sequences per clonotype (3).\n')
}


### ------------------------ Cell-based info ------------------------ ###

# Perform only if necessary.
if(is.null(reports.date)){

    # ---> Metadata.
    # @ Lane tags and rest of info.
    lane.info <- aggr.table[, .(sample.id=sample_id, chrom.batch.tag, seq.batch.tag, lane.tag)]
    lane.info <- as.data.table(separate(data=lane.info, col=lane.tag, into=c('cell.type.tag', 'peptide.pool.tag'), sep=';'))
    
    # @ Hashtag information.
    hto.info <- lapply(X=aggr.table[, hto.tag], FUN=function(hto.file){
        if(is.na(hto.file) | hto.file=='') return(NA)
        tmp.data <- fread(file=hto.file, header=TRUE)
        return(tmp.data)
    })
    names(hto.info) <- aggr.table[, sample_id]
    hto.info <- hto.info[!is.na(hto.info)]
    hto.info <- rbindlist(l=hto.info, use.names=TRUE, idcol='sample.id')
    colnames(hto.info)[2:3] <- c('barcode', 'hto.tag')
    # Replace '-' in hto classes.
    hto.info[, hto.tag:=str_replace_all(string=hto.tag, pattern='-', replacement='_')]
    # Update barcode suffix according to cellranger's style during aggr. See documentation: https://support.10xgenomics.com/single-cell-vdj/software/pipelines/latest/using/aggr
    tmp.data <- aggr.table[, .(sample.id=sample_id, aggr.sffx=as.character(1:.N))]
    hto.info <- merge(x=hto.info, y=tmp.data, by='sample.id')
    hto.info[, barcode:=str_replace(string=barcode, pattern='\\d+$', replacement=aggr.sffx)]
    hto.info[, aggr.sffx:=NULL]
    # Further collect lane info.
    hto.info <- merge(x=hto.info, y=lane.info[, .(sample.id, chrom.batch.tag)], by='sample.id', all=TRUE)
    hto.info[, sample.id:=NULL]
    # Create column to collect with donor info.
    hto.info[hto.tag!='Doublet' & hto.tag!='Negative', donor.tag:=paste(hto.tag, chrom.batch.tag, sep='-')]
    hto.info[, chrom.batch.tag:=NULL]

    # @ Donor metadata.
    # Sanity check. Donor tags from both sources should be a perfect match after joining.
    tmp.check.1 <- hto.info[!is.na(donor.tag), unique(donor.tag)]
    tmp.check.2 <- donor.meta[, unique(donor.tag)]
    # all(tmp.check.1 %in% tmp.check.2)
    # tmp.check.2[!tmp.check.2 %in% tmp.check.1] # NOTE: THIS IS INFORMATION ABOUT MISSING DONORS.
    hto.info <- merge(x=hto.info, y=donor.meta, by='donor.tag', all.x=TRUE, all.y=FALSE)

    # @ CD3 surface expression information.
    cite.info <- lapply(X=1:aggr.table[, .N], FUN=function(sffx){
        cat(sffx, '\n')
        # Load CITE-seq data and subset for CD3E counts.
        cite.file <- paste0(aggr.table[sffx, cite.tag], '/raw_feature_bc_matrix')
        tmp.data <- read.10X.data(file=cite.file)
        tmp.data <- tmp.data['CD3E', ]
        tmp.data <- data.table(
            barcode=names(tmp.data),
            cd3e.umi=tmp.data
        )
        # Keep only cells found in the TCR dataset and update barcode suffix.
        tmp.bcs <- cell.clone.info[
            str_detect(string=barcode, pattern=paste0('-', sffx, '$')),
            str_replace(string=barcode, pattern='-\\d+', replacement='-1')
        ]
        tmp.data <- tmp.data[barcode %in% tmp.bcs]
        tmp.data[,
            barcode:=str_replace(string=barcode, pattern='-\\d$', replacement=paste0('-', sffx))
        ]
        return(tmp.data)
    })
    cite.info <- rbindlist(l=cite.info, use.names=TRUE)

    # ---> General merge and collection of further information.
    # @ General merge
    meta.data <- unique(
        cell.clone.info[, .(
            barcode, clonotype.tag, sffx=str_extract(string=barcode, pattern='\\d+$')
        )]
    )
    lane.info[, sffx:=as.character(1:.N)]
    meta.data <- merge(x=meta.data, y=lane.info, by='sffx'); meta.data[, sffx:=NULL]
    meta.data <- merge(x=meta.data, y=hto.info, by='barcode', all.x=TRUE, all.y=FALSE)
    meta.data <- merge(x=meta.data, y=cite.info, by='barcode', all.x=TRUE, all.y=FALSE)
    meta.data <- merge(x=meta.data, y=clone.info, by='clonotype.tag')
}

# ---> General sanity checks.
# @1
tmp.data.1 <- meta.data[
    !is.na(donor.id.tag),
    .(
        pp.no=uniqueN(peptide.pool.tag),
        sample.no=uniqueN(sample.id),
        cell.count=.N,
        pps=paste0(sort(unique(peptide.pool.tag)), collapse=';')
    ),
    by=.(donor.id.tag)
]
setorderv(x=tmp.data.1, cols='pp.no', order=-1)
tmp.file.name <- paste0(reports.path, '/ManualInspection01.csv')
fwrite(file=tmp.file.name, x=tmp.data.1, quote=TRUE, na=NA)
# @2
tmp.data.1 <- meta.data[
    !is.na(donor.id.tag),
    .(
        sample.no=uniqueN(sample.id),
        sample.ids=paste0(sort(unique(sample.id)), collapse=';'),
        cell.count=.N
    ),
    by=.(donor.id.tag, peptide.pool.tag)
]
setorderv(x=tmp.data.1, cols='sample.no', order=-1)
tmp.file.name <- paste0(reports.path, '/ManualInspection02.csv')
fwrite(file=tmp.file.name, x=tmp.data.1, quote=TRUE, na=NA)
# @3
tmp.data.1 <- meta.data[
    !is.na(donor.id.tag),
    .(
        donor.no=uniqueN(donor.id.tag),
        batch.count=uniqueN(chrom.batch.tag),
        batches=paste0(sort(unique(chrom.batch.tag)), collapse=';'),
        cell.count=.N
    ),
    by=.(peptide.pool.tag)
]
setorderv(x=tmp.data.1, cols='peptide.pool.tag', order=-1)
tmp.file.name <- paste0(reports.path, '/ManualInspection03.csv')
fwrite(file=tmp.file.name, x=tmp.data.1, quote=TRUE, na=NA)


### ------------- Experimental misassignment correction ------------- ###

exp.corr.path <- paste0(reports.path, '/exp_misassignment_correction')
if(!dir.exists(exp.corr.path)) dir.create(exp.corr.path)

# ---> Sample pairwise similarity indexes for index peptide pools BEFORE correction.
# -------- 
# @ Rationale
# By exploring the results from the previous section (the one right above), we found that there might be issues for the CMV repeats. CMV repeats correspond to FACS batches DICE-TCR018 and DICE-TCR028, with them being repeats of samples originally processed in DICE-TCR001. Therefore, we will explore these dissimilarities in more detail.
# --------
exp.corr.rep.bf.path <- paste0(exp.corr.path, '/repeat_idxs_bf_corr')
if(!dir.exists(exp.corr.rep.bf.path)) dir.create(exp.corr.rep.bf.path)
# ------        Process defined as a function to be repeated before and after correction.
assess.rep.sim <- function(meta.data, exp.corr.rep.path){
    rep.facs.batches <- list(
        CMV=c('DICE_TCR018', 'DICE_TCR028'),
        PIV='DICE_TCR028'
    )
    for(peptide.pool in names(rep.facs.batches)){
        # ------        Retrieve data
        # Identify donors w/ repeats.
        donors.w.rep <- meta.data[
            !is.na(donor.id.tag) & chrom.batch.tag%in%rep.facs.batches[[peptide.pool]] & peptide.pool.tag==peptide.pool,
            unique(donor.id.tag)
        ]
        # Identify samples that correspond to repeats (either for the original or the repeated library) and relabel.
        tmp.data <- meta.data[
            peptide.pool.tag==peptide.pool & !is.na(donor.id.tag) & donor.id.tag%in%donors.w.rep,
            .(donor.id.tag, clonotype.tag, lib.long.id=sample.id)
        ]
        tmp.names <- paste(peptide.pool, 1:tmp.data[, uniqueN(lib.long.id)], sep='.')
        names(tmp.names) <- tmp.data[, gtools::mixedsort(unique(lib.long.id))]
        tmp.data[, lib.id:=tmp.names[lib.long.id]]
        # Keep track of library labels.
        tmp.names <- data.table(
            lib.id.new=tmp.names,
            lib.id.original=names(tmp.names),
            is.repeat=names(tmp.names) %in% meta.data[chrom.batch.tag%in%rep.facs.batches[[peptide.pool]], unique(sample.id)]
        )
        tmp.file.name <- paste0(exp.corr.rep.path, '/LibraryIDDict_', peptide.pool, '.csv')
        fwrite(file=tmp.file.name, x=tmp.names, quote=TRUE, na=NA)
        # Calculate relative frequencies per sample (defined as a set of donor and peptide pool).
        tmp.data.1 <- tmp.data[,
            .(abs.freq=.N),
            by=.(lib.id, donor.id.tag, clonotype.tag)
        ]
        tmp.data.2 <- tmp.data[,
            .(sample.freq=.N),
            by=.(lib.id, donor.id.tag)
        ]
        sample.freqs <- merge(x=tmp.data.1, y=tmp.data.2, by=c('lib.id', 'donor.id.tag'))
        sample.freqs[, rel.freq:=abs.freq/sample.freq]
        tmp.check <- sample.freqs[, .(check=round(sum(rel.freq), digits=10)==1), by=.(lib.id, donor.id.tag)][, all(check)]
        if(!tmp.check) stop('Unexpected error during experimental misassignment correction (2).\n')
        sample.freqs[, sample.id:=paste(donor.id.tag, lib.id, sep=';')]
        # Calculate jaccard indexes in a pairwise manner either or not accounting for clone size, per TCR chain and per type of CDR3 sequence molecule.
        tmp.data <- get.sim.idxs(sample.freqs=sample.freqs)
        tmp.data <- as.data.table(separate(data=tmp.data, col='sample.id.1', into=c('donor.id.1', 'lib.1'), sep=';', remove=FALSE))
        tmp.data <- as.data.table(separate(data=tmp.data, col='sample.id.2', into=c('donor.id.2', 'lib.2'), sep=';', remove=FALSE))
        # ------        Viz results.
        # Depict summary index as a heatmap.
        to.plot <- tmp.data[, .(sample.id.1, sample.id.2, cosine.idx)]
        to.plot <- spread(data=to.plot, key=sample.id.2, value=cosine.idx)
        tmp.row.names <- to.plot$sample.id.1
        to.plot$sample.id.1 <- NULL
        to.plot <- as.matrix(to.plot)
        row.names(to.plot) <- tmp.row.names
        diag(to.plot) <- NA
        # Set color scale and breaks for heatmap.
        col.breaks <- seq(from=min(range(to.plot, na.rm=TRUE)), to=max(range(to.plot, na.rm=TRUE)), length.out=100)
        mid.point <- which.min(abs(col.breaks - 0))
        hmap.col.scale.1 <- colorRampPalette(c('#ffffff', '#fffff6', '#ffffed', '#ffffb2'))(mid.point)
        hmap.col.scale.2 <- colorRampPalette(c('#ffff77', '#ffff00', '#ffffb2', '#ffff77', '#ffd280', '#ffd280', '#ffbc40', '#ffa500', '#ff5300', '#ff0000', '#8b0000', '#3f0000'))(100-(mid.point+1))
        hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
        # Row metadata.
        row.col.meta <- as.data.frame(unique(tmp.data[, .(sample.id=sample.id.1, Donor=donor.id.1, Library=lib.1)]))
        row.col.meta <- merge(
            x=row.col.meta,
            y=as.data.frame(tmp.names[, .(
                Library=lib.id.new,
                Repeat=ifelse(test=is.repeat, yes='Repeat', no='Original')
            )]),
            by='Library'
        )
        row.names(row.col.meta) <- row.col.meta$sample.id; row.col.meta$sample.id <- NULL
        meta.cols <- list(
            Library=viridis::magma(n=length(unique(row.col.meta$Library))),
            Donor=viridis::viridis(n=length(unique(row.col.meta$Donor))),
            Repeat=viridis::turbo(n=length(unique(row.col.meta$Repeat)))
        )
        names(meta.cols$Library) <- unique(row.col.meta$Library)
        names(meta.cols$Donor) <- unique(row.col.meta$Donor)
        names(meta.cols$Repeat) <- unique(row.col.meta$Repeat)
        # Plot.
        tmp.file.name <- paste0(exp.corr.rep.path, '/CosineIdx_RepeatExps_', peptide.pool, '.pdf')
        pheatmap(
            mat=to.plot, 
            color=hmap.col.scale, breaks=col.breaks, scale='none',
            cluster_rows=FALSE, cluster_cols=FALSE,
            annotation_col=row.col.meta, annotation_row=row.col.meta, 
            annotation_colors=meta.cols,
            show_colnames=TRUE, show_rownames=TRUE,
            filename=tmp.file.name, width=15, height=12.5
        )
        # @ In the case of the CMV peptide pool, subset for donors in libraries CMV.1 and CMV.5.
        if(peptide.pool=='CMV'){
            tmp.samples <- rev(gtools::mixedsort(row.names(row.col.meta)[row.col.meta$Library=='CMV.1' | row.col.meta$Library=='CMV.5']))
            row.col.meta <- row.col.meta[tmp.samples, ]
            meta.cols <- list(
                Library=viridis::magma(n=length(unique(row.col.meta$Library))),
                Donor=viridis::viridis(n=length(unique(row.col.meta$Donor))),
                Repeat=viridis::turbo(n=length(unique(row.col.meta$Repeat)))
            )
            names(meta.cols$Library) <- unique(row.col.meta$Library)
            names(meta.cols$Donor) <- unique(row.col.meta$Donor)
            names(meta.cols$Repeat) <- unique(row.col.meta$Repeat)
            to.plot <- to.plot[tmp.samples, tmp.samples]
            # Re-set color scale and breaks for heatmap.
            col.breaks <- seq(from=min(range(to.plot, na.rm=TRUE)), to=max(range(to.plot, na.rm=TRUE)), length.out=100)
            mid.point <- which.min(abs(col.breaks - 0))
            hmap.col.scale.1 <- colorRampPalette(c('#ffffff', '#fffff6', '#ffffed', '#ffffb2'))(mid.point)
            hmap.col.scale.2 <- colorRampPalette(c('#ffff77', '#ffff00', '#ffffb2', '#ffff77', '#ffd280', '#ffd280', '#ffbc40', '#ffa500', '#ff5300', '#ff0000', '#8b0000', '#3f0000'))(100-(mid.point+1))
            hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
            # Plot.
            tmp.file.name <- paste0(exp.corr.rep.path, '/CosineIdx_RepeatExps_', peptide.pool, '_Subset.pdf')
            pheatmap(
                mat=to.plot, 
                color=hmap.col.scale, breaks=col.breaks, scale='none',
                cluster_rows=FALSE, cluster_cols=FALSE,
                annotation_col=row.col.meta, annotation_row=row.col.meta, 
                annotation_colors=meta.cols,
                show_colnames=TRUE, show_rownames=TRUE,
                filename=tmp.file.name, width=12, height=10
            )
        }
    }
}
assess.rep.sim(meta.data=meta.data, exp.corr.rep.path=exp.corr.rep.bf.path)

# ---> Apply necessary modifications for repeat sample swap correction.
tmp.data <- rep.swap.corr[,
    .(
        sample.id=repeat.sample.id,
        donor.id.tag=repeat.donor.id,  # This might seem counterintuitive because of the labeling, but it must be correct.
        new.donor.id.tag=orig.donor.id
    )
]
tmp.data <- unique(tmp.data) # Not an error. Donor SDBB-120 was repeated twice, but we mean to make the correction only for the last repeat (i.e., for sample 'DICE_TCR028_Hu_CMV_2_17D_TCR)
meta.data <- merge(
    x=meta.data, y=tmp.data,
    by=c('sample.id', 'donor.id.tag'),
    all.x=TRUE, all.y=FALSE
)
meta.data[!is.na(new.donor.id.tag), donor.id.tag:=new.donor.id.tag]
meta.data[, new.donor.id.tag:=NULL]

# ---> Sample pairwise similarity indexes for index peptide pools AFTER correction.
exp.corr.rep.afr.path <- paste0(exp.corr.path, '/repeat_idxs_afr_corr')
if(!dir.exists(exp.corr.rep.afr.path)) dir.create(exp.corr.rep.afr.path)
assess.rep.sim(meta.data=meta.data, exp.corr.rep.path=exp.corr.rep.afr.path)

# ---> Peptide pool-specific donor similarity indexes BEFORE correction.
exp.corr.ppd.bf.path <- paste0(exp.corr.path, '/donor_pairwise_idxs_bf_corr')
if(!dir.exists(exp.corr.ppd.bf.path)) dir.create(exp.corr.ppd.bf.path)
# ------        Process defined as a function to be repeated before and after correction.
assess.ppd.sim <- function(meta.data, exp.corr.ppd.path){
    # ------        Retrieve data
    # Calculate relative frequencies per sample (defined as a set of donor and peptide pool).
    tmp.data.1 <- meta.data[
        !is.na(donor.id.tag),
        .(abs.freq=.N),
        by=.(donor.id.tag, peptide.pool.tag, clonotype.tag)
    ]
    tmp.data.2 <- meta.data[
        !is.na(donor.id.tag),
        .(sample.freq=.N),
        by=.(donor.id.tag, peptide.pool.tag)
    ]
    sample.freqs <- merge(x=tmp.data.1, y=tmp.data.2, by=c('donor.id.tag', 'peptide.pool.tag'))
    sample.freqs[, rel.freq:=abs.freq/sample.freq]
    tmp.check <- sample.freqs[, .(check=round(sum(rel.freq), digits=10)==1), by=.(donor.id.tag, peptide.pool.tag)][, all(check)]
    if(!tmp.check) stop('Unexpected error during experimental misassignment correction (1).\n')
    sample.freqs[, sample.id:=paste(peptide.pool.tag, donor.id.tag, sep=';')]
    # Calculate jaccard indexes in a pairwise manner either or not accounting for clone size, per TCR chain and per type of CDR3 sequence molecule.
    tmp.data <- get.sim.idxs(sample.freqs=sample.freqs, extra.verbose=FALSE)
    tmp.data <- as.data.table(separate(data=tmp.data, col='sample.id.1', into=c('peptide.pool.1', 'donor.id.1'), sep=';', remove=FALSE))
    tmp.data <- as.data.table(separate(data=tmp.data, col='sample.id.2', into=c('peptide.pool.2', 'donor.id.2'), sep=';', remove=FALSE))
    # ------        Viz results.
    # For each donor and peptide pool ("sample"), we will determine a single similarity index per donor, calculated as the median of the cosine index across all peptide pools for a given donor.
    to.plot <- tmp.data[,
        .(
            summ.idx=median(cosine.idx)
        ),
        by=.(peptide.pool=peptide.pool.1, donor.id.1, donor.id.2)
    ]
    # Depict summary index as a heatmap per donor.
    donor.ids <- to.plot[, unique(donor.id.1)]
    for(donor.id in donor.ids){
        to.plot.tmp <- to.plot[
            donor.id.1==donor.id,
            .(peptide.pool, donor.id=donor.id.2, summ.idx)
        ]
        to.plot.tmp <- spread(data=to.plot.tmp, key=donor.id, value=summ.idx)
        tmp.row.names <- to.plot.tmp$peptide.pool
        to.plot.tmp$peptide.pool <- NULL
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
        row.meta <- data.frame(
            row.names=row.names(to.plot.tmp),
            `PP`=row.names(to.plot.tmp)
        )
        row.cols <- list(
            `PP`=peptide.pool.cols
        )
        # Plot.
        tmp.file.name <- paste0(exp.corr.ppd.path, '/SummIdx_', donor.id, '.pdf')
        pheatmap(
            mat=to.plot.tmp, 
            color=hmap.col.scale, breaks=col.breaks, scale='none',
            cluster_rows=FALSE, cluster_cols=FALSE,
            annotation_col=NULL, annotation_row=row.meta, 
            annotation_colors=row.cols,
            show_colnames=TRUE, show_rownames=FALSE,
            filename=tmp.file.name, width=12, height=4
        )
    }
}
assess.ppd.sim(meta.data=meta.data, exp.corr.ppd.path=exp.corr.ppd.bf.path)

# ---> Apply necessary modifications for peptide pool-specific donor sample swap correction.
tmp.data <- pd.swap.corr[,
    .(
        peptide.pool.tag=peptide.pool,
        donor.id.tag=donor.id.1,  # This might seem counterintuitive because of the labeling, but it must be correct.
        new.donor.id.tag=donor.id.2
    )
]
# tmp.data <- unique(tmp.data) # Not an error. Donor SDBB-120 was repeated twice, but we mean to make the correction only for the last repeat (i.e., for sample 'DICE_TCR028_Hu_CMV_2_17D_TCR)
meta.data <- merge(
    x=meta.data, y=tmp.data,
    by=c('peptide.pool.tag', 'donor.id.tag'),
    all.x=TRUE, all.y=FALSE
)
meta.data[!is.na(new.donor.id.tag), donor.id.tag:=new.donor.id.tag]
meta.data[, new.donor.id.tag:=NULL]

# ---> Peptide pool-specific donor similarity indexes BEFORE correction.
exp.corr.ppd.afr.path <- paste0(exp.corr.path, '/donor_pairwise_idxs_afr_corr')
if(!dir.exists(exp.corr.ppd.afr.path)) dir.create(exp.corr.ppd.afr.path)
assess.ppd.sim(meta.data=meta.data, exp.corr.ppd.path=exp.corr.ppd.afr.path)


### ----------------- TCR-based donor deconvolution ----------------- ###
### ------- HTO-based deconvolution-misassignment correction -------- ###

correction.path <- paste0(reports.path, '/deconvolution_correction')
if(!dir.exists(correction.path)) dir.create(correction.path)

# Perform only if necessary.
if(is.null(reports.date)){

    # ---> Identify potentially public TCRs.
    # A clonotype is said to be potentially public when it was found to be expressed by cells from 2 or more donors, either for one or more specificities, and w/ clone sie >1 in at least one of the donors.
    # Identify all clonotypes found in more than one donor.
    tmp.data <- meta.data[
        !is.na(donor.id.tag),
        .(donor.count=uniqueN(donor.id.tag)),
        by=.(clonotype.tag)
    ]
    # Save for report.
    tmp.report <- tmp.data[, .SD[donor.count>1, .N]*100/.N]
    # Exclude all clonotypes that have clone size=1 in all donors they were found in.
    tmp.data.2 <- meta.data[
        !is.na(donor.id.tag) & clonotype.tag %in% tmp.data[donor.count>1, clonotype.tag],
        .(clone.size=.N),
        by=.(donor.id.tag, clonotype.tag)
    ]
    tmp.data.2 <- tmp.data.2[,
        .(check=sum(clone.size>1)),
        by=clonotype.tag
    ][check>0, clonotype.tag]
    tmp.data <- tmp.data[clonotype.tag %in% tmp.data.2]

    # ---> Determine major and minor sizes for a potentially public TCRs.
    # There exists only one major size, the one for the donor(s) with size==max(size).
    # There might be multiple donors with the major size. If that's the case, we pick randomly one of them to be the major donor. Otherwise, the single donor w/ size==max(size) is taken as the major donor.
    # All of the donors that were not determined to be the major donor will be taken as the minor donors. Their corresponding sizes (even if they are the same as the one for the major donor -tie, but it cannot be larger-) will be taken as the minor sizes.
    tmp.data.1 <- meta.data[
        !is.na(donor.id.tag) & clonotype.tag %in% tmp.data[, clonotype.tag],
        .(
            size=.N
        ),
        by=.(clonotype.tag, donor.id.tag)
    ]
    tmp.data.2 <- tmp.data.1[,
        .SD[size==max(size), ][1,],
        by=clonotype.tag
    ]
    tmp.check <- tmp.data.2[, uniqueN(clonotype.tag)==.N]
    if(!tmp.check) stop('We must have kept the same number of entries as we have unique potentially public clonotypes.\n')
    pp.clone.info <- merge(x=tmp.data.2, y=tmp.data.1, by='clonotype.tag', suffixes=c('.major', '.minor'))
    pp.clone.info <- pp.clone.info[donor.id.tag.major!=donor.id.tag.minor]

    # ---> Determine difference in size between major and minor donors and difference fraction
    pp.clone.info[, size.diff:=size.major-size.minor]
    pp.clone.info[, diff.prop:=size.diff/size.major]

    # ---> Visual exploration BEFORE correction
    # NOTE: In an earlier version, I used this part to determine what would be the best steps for a well-performer correction algorithm (i.e., to develope the 2-step correction process).
    tmp.legend <- paste0('n=', pp.clone.info[, .N], ' major vs minor donor relations')
    tmp.ggplot.1 <- ggplot(data=pp.clone.info, aes(x=size.major, y=size.minor)) +
        geom_point(size=0.6) +
        labs(x='Clone size in major donor', y='Clone size in minor donor', caption=tmp.legend) +
        theme_bw()
    tmp.ggplot.2 <- ggplot(data=pp.clone.info, aes(x=size.major, y=size.minor)) +
        geom_point(size=0.6) +
        scale_y_continuous(limits=c(0, pp.clone.info[, max(size.major)])) +
        labs(x='Clone size in major donor', y='Clone size in minor donor', caption=tmp.legend) +
        theme_bw()
    tmp.ggplot.3 <- ggplot(data=pp.clone.info, aes(x=size.major, y=size.minor)) +
        geom_point(size=0.6) +
        scale_x_log10() +
        scale_y_log10() +
        labs(x='Clone size in major donor', y='Clone size in minor donor', caption=tmp.legend) +
        theme_bw()
    tmp.ggplot.4 <- ggplot(data=pp.clone.info, aes(x=size.major, y=size.minor)) +
        geom_point(size=0.6) +
        scale_x_log10() +
        scale_y_log10(limits=c(1, pp.clone.info[, max(size.major)])) +
        labs(x='Clone size in major donor', y='Clone size in minor donor', caption=tmp.legend) +
        theme_bw()
    tmp.ggplot.5 <- ggplot(data=pp.clone.info, aes(x=size.major, y=size.minor)) +
        geom_point(size=0.6) +
        geom_abline(linewidth=1.5, alpha=0.5, color='red', linetype='dashed') +
        scale_x_log10() +
        scale_y_log10(limits=c(1, pp.clone.info[, max(size.major)])) +
        labs(x='Clone size in major donor', y='Clone size in minor donor', caption=tmp.legend) +
        theme_bw()
    tmp.ggplot.6 <- ggplot(data=pp.clone.info[size.minor>1], aes(x=diff.prop)) +
        geom_density(fill='lightblue', alpha=0.6) +
        geom_vline(
            xintercept=quantile(x=pp.clone.info[size.minor>1, diff.prop], probs=0.25),
            linewidth=1.5, alpha=0.5, color='red', linetype='dashed'
        ) +
        geom_vline(
            xintercept=quantile(x=pp.clone.info[size.minor>1, diff.prop], probs=0.50),
            linewidth=1.5, alpha=0.5, color='black'
        ) +
        geom_vline(
            xintercept=quantile(x=pp.clone.info[size.minor>1, diff.prop], probs=0.75),
            linewidth=1.5, alpha=0.5, color='red', linetype='dashed'
        ) +
        scale_x_continuous(limits=c(0, NA), expand=c(0, 0)) +
        scale_y_continuous(limits=c(0, NA), expand=c(0, 0)) +
        labs(x='Difference in size between major and minor donors', y='Density', caption='Here we removed entries taken care of by step 1.') +
        theme_bw()
    tmp.file.name <- paste0(correction.path, '/PotPublicTCRsBeforeCorrection.pdf')
    pdf(file=tmp.file.name)
    print(tmp.ggplot.1)
    print(tmp.ggplot.2)
    print(tmp.ggplot.3)
    print(tmp.ggplot.4)
    print(tmp.ggplot.5)
    print(tmp.ggplot.6)
    dev.off()

    # ---> Apply correction.
    # See the two-step correction process documentation for correction criteria.

    # Determine difference fraction threshold (75% quantile of the farction threshold distribution)
    diff.frac.thold <- pp.clone.info[size.minor>1, quantile(x=diff.prop, probs=0.6)]
    # tmp.data.1[size.minor>1][
    #   , .SD[diff.prop>=diff.frac.thold, .N]/.N
    # ]

    # Do correction process per potentially public clonotype. 
    # --------| NOTE: See decision chart in the supplementary information.
    pp.clones <- pp.clone.info[, unique(clonotype.tag)]
    meta.data[, final.donor.id.tag:=donor.id.tag]
    for(tmp.clonotype in pp.clones){
        # cat(tmp.clonotype, '\n')
        tmp.data <- pp.clone.info[clonotype.tag==tmp.clonotype]
        major.donor <- tmp.data[, unique(donor.id.tag.major)]
        # Determine set B (TCR libraries from which we retrieved major donor cells that express the pp TCR).
        set.b <- meta.data[
            donor.id.tag==major.donor & clonotype.tag==tmp.clonotype,
            unique(sample.id)
        ]
        # Apply step 1 correction.
        minor.donors <- tmp.data[size.minor==1, donor.id.tag.minor]
        meta.data[
            donor.id.tag%in%minor.donors & clonotype.tag==tmp.clonotype & sample.id%in%set.b,
            # donor.id.tag%in%minor.donors & clonotype.tag==tmp.clonotype,
            final.donor.id.tag:=major.donor
        ]
        # Apply step 2 correction.
        minor.donors <- tmp.data[diff.prop>=diff.frac.thold, donor.id.tag.minor]
        meta.data[
            donor.id.tag%in%minor.donors & clonotype.tag==tmp.clonotype & sample.id%in%set.b,
            final.donor.id.tag:=major.donor
        ]
    }
    meta.data[, uncorrected.donor.id.tag:=donor.id.tag]
    meta.data[, donor.id.tag:=final.donor.id.tag]
    meta.data[, final.donor.id.tag:=NULL]

    # ---> Update list of potentially public TCRs.
    # Same definition as above.
    # Identify all clonotypes found in more than one donor.
    tmp.data <- meta.data[
        !is.na(donor.id.tag),
        .(donor.count=uniqueN(donor.id.tag)),
        by=.(clonotype.tag)
    ]
    # Save for report.
    tmp.report <- c(tmp.report, tmp.data[, .SD[donor.count>1, .N]*100/.N])
    # Exclude all clonotypes that have clone size=1 in all donors they were found in.
    tmp.data.2 <- meta.data[
        !is.na(donor.id.tag) & clonotype.tag %in% tmp.data[donor.count>1, clonotype.tag],
        .(clone.size=.N),
        by=.(donor.id.tag, clonotype.tag)
    ]
    tmp.data.2 <- tmp.data.2[,
        .(check=sum(clone.size>1)),
        by=clonotype.tag
    ][check>0, clonotype.tag]
    tmp.data <- tmp.data[clonotype.tag %in% tmp.data.2]

    # ---> Determine major and minor sizes for updated list of potentially public TCRs.
    # Do as above.
    tmp.data.1 <- meta.data[
        !is.na(donor.id.tag) & clonotype.tag %in% tmp.data[, clonotype.tag],
        .(
            size=.N
        ),
        by=.(clonotype.tag, donor.id.tag)
    ]
    tmp.data.2 <- tmp.data.1[,
        .SD[size==max(size), ][1,],
        by=clonotype.tag
    ]
    tmp.check <- tmp.data.2[, uniqueN(clonotype.tag)==.N]
    if(!tmp.check) stop('We must have kept the same number of entries as we have unique potentially public clonotypes.\n')
    tmp.data.1 <- merge(x=tmp.data.2, y=tmp.data.1, by='clonotype.tag', suffixes=c('.major', '.minor'))
    tmp.data.1 <- tmp.data.1[donor.id.tag.major!=donor.id.tag.minor]

    # ---> Merge original vs updated info of potentially public TCRs (i.e., before and after correction)
    tmp.data.2 <- merge(x=pp.clone.info, y=tmp.data.1, by=c('clonotype.tag', 'donor.id.tag.major', 'donor.id.tag.minor'), suffixes=c('.before', '.after'), all=TRUE)
    tmp.check <- pp.clone.info[, .N]==tmp.data.2[, .N] # Interestingly enough, there aren't any new entries of potentially public TCRs where a major-minor donor relationship would be uncovered (hard to explain why this might be possible, yet I think it would). Setting a condition just in case this case comes up at some point when adding new data.
    if(!tmp.check) stop('Unexpected error.\nPlease check why new pairs of minor and minor donors are being uncovered after the two-step correction process.\n')
    pp.clone.info <- tmp.data.2
    # Tag corrected pairs.
    # pp.clone.info[is.na(size.major.after), all(is.na(size.minor.after))]
    pp.clone.info[, was.corrected:=FALSE]
    pp.clone.info[is.na(size.major.after), was.corrected:=TRUE]

    # ---> Visual exploration after correction
    tmp.cols <- c('FALSE'='#0495AD', 'TRUE'='#AD4E00')
    tmp.legend <- paste0('n=', pp.clone.info[, .N], ' major vs minor donor relations')
    tmp.ggplot.5 <- ggplot(data=pp.clone.info, aes(x=size.major.before, y=size.minor.before, col=was.corrected)) +
        geom_point(size=0.35, alpha=0.6) +
        geom_abline(linewidth=1.5, alpha=0.5, color='red', linetype='dashed') +
        scale_x_log10() +
        scale_y_log10(limits=c(1, pp.clone.info[, max(size.major.before)])) +
        scale_color_manual(values=tmp.cols) +
        labs(x='Clone size in major donor', y='Clone size in minor donor', color='Was\ncorrected?', caption=tmp.legend) +
        theme_bw()
    tmp.file.name <- paste0(correction.path, '/PotPublicTCRsAfterCorrection.pdf')
    pdf(file=tmp.file.name)
    print(tmp.ggplot.5)
    dev.off()

    # ---> Output report.
    tmp.report <- data.table(
        `Percent of potentially public clonotypes before correction`=tmp.report[1],
        `Percent of potentially public clonotypes after correction`=tmp.report[2],
        `Corrected clonotypes`=pp.clone.info[was.corrected==TRUE, uniqueN(clonotype.tag)],
        `Corrected clonotype and donor pairs`=pp.clone.info[was.corrected==TRUE, .N],
        `Corrected cells`=pp.clone.info[was.corrected==TRUE, sum(size.minor.before)]
    )
    tmp.report <- gather(data=tmp.report, key='Description', value='Number Or Fraction')
    tmp.file.name <- paste0(correction.path, '/CorrectionProcessReport.csv')
    fwrite(file=tmp.file.name, x=tmp.report, na=NA)
}


### ----------------- Antigen specificity definition ---------------- ###

# Here, we have the advantage that specificity is given as a binary prediction (reactive vs non-reactive) for every single peptide, potentially allowing us to better discern the cross-reactive TCRs.

ag.spc.path <- paste0(reports.path, '/ag_specificity_determination')
if(!dir.exists(ag.spc.path)) dir.create(ag.spc.path)

if(is.null(reports.date)){
    # ---> Clonotype ID.
    this.tcr.data <- meta.data[!is.na(donor.id.tag), .(clonotype.tag, donor.id.tag, peptide.pool.tag)]
    this.tcr.data <- this.tcr.data[, .(fad=.N), by=.(clonotype.tag, donor.id.tag, peptide.pool.tag)]
    this.tcr.data[, long.id.tag:=paste(clonotype.tag, donor.id.tag, sep=';')]

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
    tmp.data <- this.tcr.data[, .(gen.clonotype.status, peptide.pool.tag)]
    tmp.ggplot.1 <- ggplot(data=tmp.data, aes(x=gen.clonotype.status)) +
        geom_bar(width=0.7, fill='lightblue', color='black', linewidth=1.5, alpha=0.8) +
        scale_y_continuous(expand=c(0, NA)) +
        labs(x='General clonotype status', y='Clonotype count') +
        theme_bw()
    tmp.ggplot.2 <- ggplot(data=tmp.data, aes(x=gen.clonotype.status, fill=peptide.pool.tag)) +
        geom_bar(width=0.7, color='black', linewidth=1.5, alpha=0.8) +
        facet_wrap(facets=~peptide.pool.tag, scales='free_y') +
        scale_y_continuous(expand=c(0, NA)) +
        scale_fill_manual(values=peptide.pool.cols) +
        labs(x='General clonotype status', y='Clonotype count') +
        theme_bw() + theme(legend.position='none', axis.text.x=element_text(angle=-45))
    tmp.file.name <- paste0(ag.spc.path, '/GeneralClonotypeStatus_Freqs.pdf')
    pdf(file=tmp.file.name)
    print(tmp.ggplot.1)
    print(tmp.ggplot.2)
    dev.off()

    # ---> Within-donor relative frequencies (FRD).
    # Remove singletons at this point.
    this.tcr.data <- this.tcr.data[gen.clonotype.status!='Singleton']
    # Define within-donor absolute and relative frequencies (FADs and FRDs, respectively).
    tmp.data <- this.tcr.data[, .(faw=sum(fad)), by=.(donor.id.tag, peptide.pool.tag)]
    this.tcr.data <- merge(x=this.tcr.data, y=tmp.data, by=c('donor.id.tag', 'peptide.pool.tag'))
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
    # Sanity check.
    tmp.check <- tmp.data[,
        .(
            check.1=.SD[
                gen.clonotype.status=='Unambiguous',
                all(is.na(pp.2) & is.na(frd.2) & is.na(fad.2))
            ],
            check.2=.SD[
                gen.clonotype.status=='Ambiguous',
                all(!is.na(pp.2) & !is.na(frd.2) & !is.na(fad.2))
            ]
        )
    ][, check.1 & check.2]
    if(!tmp.check) stop('Something went wrong while defining antigen specificities (2).\n')

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
    fc.tholds <- tmp.data[, quantile(x=spc.fc, probs=c(0.20, 0.6, 0.9))]
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
    spc.data <- separate(data=spc.data, col='long.id.tag', into=c('clonotype.tag', 'donor.id.tag'), sep=';')
    meta.data <- merge(x=meta.data, y=spc.data, by=c('clonotype.tag', 'donor.id.tag'), all=TRUE)
    tmp.lvls <- names(peptide.pool.cols)[names(peptide.pool.cols)%in%meta.data$specificity.class.tag]
    meta.data$specificity.class.tag <- factor(x=meta.data$specificity.class.tag, levels=tmp.lvls)
    meta.data$specificity.raw.tag <- factor(x=meta.data$specificity.raw.tag, levels=tmp.lvls)
}


### ------------------------ Cell-based info ------------------------ ###

# Save if necessary. Otherwise, read.
if(is.null(reports.date)){
    # Sort according to barcode.
    meta.data[, to.order:=as.numeric(str_extract(string=barcode, pattern='\\d+$'))]
    setorderv(x=meta.data, col=c('to.order', 'barcode'), order=1)
    meta.data[, to.order:=NULL]
    # @ Whole dataset.
    tmp.file.name <- paste0(reports.path, '/CellBasedTCRData.csv')
    fwrite(file=tmp.file.name, x=meta.data, quote=TRUE, na=NA)
    # @ Cells with donor IDs.
    tmp.file.name <- paste0(reports.path, '/CellBasedTCRData_Clean.csv')
    tmp.data <- meta.data[!is.na(donor.id.tag)]
    fwrite(file=tmp.file.name, x=tmp.data, quote=TRUE, na=NA)
}else{
    tmp.file.name <- paste0(reports.path, '/CellBasedTCRData.csv')
    if(!file.exists(tmp.file.name)) stop('Failure at attempt to define cell-based info file.\n')
    meta.data <- fread(file=tmp.file.name)
    rm(tmp.file.name)
}


### --------------------- Clonotype-based info ---------------------- ###

# Perform only if necessary.
if(is.null(reports.date)){

    # ---> General information.
    tmp.data.1 <- unique(meta.data[
        !is.na(donor.id.tag),
        .(
            cdr3b.nt.seq, cdr3a.nt.seq,
            cdr3b.aa.seq, cdr3a.aa.seq,
            trb.v, trb.j,
            tra.v, tra.j
        ),
        by=clonotype.tag
    ])

    # ---> Identify potentially public TCRs.
    tmp.data.2 <- meta.data[
        !is.na(donor.id.tag),
        .(
            public.status=ifelse(
                test=uniqueN(donor.id.tag)>1,
                yes='likely public',
                no='private'
            )
        ),
        by=.(clonotype.tag)
    ]

    # ---> Total frequencies.
    tmp.data.3 <- meta.data[
        !is.na(donor.id.tag),
        .(total.abs.freq=.N),
        by=.(clonotype.tag, donor.id.tag)
    ]
    tmp.data.4 <- tmp.data.3[, .(total=sum(total.abs.freq)), by=.(donor.id.tag)]
    tmp.data.3 <- merge(
        x=tmp.data.3, y=tmp.data.4,
        by='donor.id.tag'
    )
    tmp.data.3[, total.rel.freq:=total.abs.freq/total]
    tmp.data.3[, total:=NULL]

    # ---> Specificity per clonotype and donor.
    tmp.data.4 <- unique(meta.data[
        !is.na(donor.id.tag) & !is.na(clonotype.tag),
        .(
            gen.clonotype.status,
            frd.1, frd.2, pp.1, pp.2, spc.fc,
            specificity.class.tag, specificity.raw.tag, specificity.ags.tag,
            freq.conf.lvl
        ),
        by=.(clonotype.tag, donor.id.tag)
    ])
    tmp.check <- tmp.data.4[, .N]==tmp.data.4[, unique(length(paste0(clonotype.tag, donor.id.tag, sep=';')))]
    if(!tmp.check) stop('Unexpected error\n')

    # ---> Absolute counts per peptide pool.
    tmp.data.5 <- meta.data[
        !is.na(donor.id.tag),
        .(size=.N),
        by=.(clonotype.tag, donor.id.tag, peptide.pool.tag)
    ]
    tmp.data.5[, peptide.pool.tag:=paste0('abs.freq.', peptide.pool.tag)]
    tmp.data.5 <- as.data.table(spread(data=tmp.data.5, key=peptide.pool.tag, value=size, fill=0))

    # ---> Relative counts per peptide pool.
    tmp.data.6.1 <- as.data.table(gather(data=tmp.data.5, key='peptide.pool.tag', value='size', -`donor.id.tag`, -`clonotype.tag`))
    tmp.data.6.2 <- tmp.data.6.1[, .(total.freq=sum(size)), by=.(donor.id.tag, peptide.pool.tag)]
    tmp.data.6 <- merge(x=tmp.data.6.1, y=tmp.data.6.2, by=c('donor.id.tag', 'peptide.pool.tag'), all.x=TRUE)
    tmp.data.6[, r.val:=size/total.freq]
    tmp.data.6[, `:=`(total.freq=NULL, size=NULL)]
    tmp.data.6 <- as.data.table(spread(data=tmp.data.6, key=peptide.pool.tag, value=r.val, fill=0))
    colnames(tmp.data.6) <- str_replace(string=colnames(tmp.data.6), pattern='^abs.freq', replacement='rel.freq')

    # Scaled CD3E surface expression rank per clonotype and donor on a peptide pool basis.
    tmp.data.7 <- meta.data[
        !is.na(donor.id.tag) & !is.na(clonotype.tag) & !is.na(cd3e.umi) &
        !is.na(specificity.class.tag),
        .SD[,
            .(
                abs.freq=.N,
                cd3e.mean=mean(cd3e.umi),
                cd3e.median=median(cd3e.umi)
            ),
            by=.(clonotype.tag, donor.id.tag)
        ],
        by=.(sample.id)
    ]
    tmp.data.7 <- tmp.data.7[abs.freq>2] # Consider only clonotypes w/ absolute freq>2
    tmp.cols <- colnames(tmp.data.7)
    # Ranking according to either stat.
    tmp.data.7 <- lapply(X=c('mean', 'median'), FUN=function(tmp.stat){
        setorderv(x=tmp.data.7, cols=paste0('cd3e.', tmp.stat), order=-1)
        tmp.data <- tmp.data.7[,
            .SD[, .(
                rank=1:.N,
                clonotype.tag, donor.id.tag
            )],
            by=sample.id
        ]
        tmp.data <- merge(
            x=tmp.data,
            y=tmp.data.7[, .(scale.factor=.SD[, 100/.N]), by=sample.id],
            by='sample.id', sort=FALSE
        )
        tmp.data[, scaled.rank:=rank*scale.factor]
        tmp.data <- tmp.data[, .(sample.id, clonotype.tag, donor.id.tag, scaled.rank)]
        colnames(tmp.data)[colnames(tmp.data)=='scaled.rank'] <- paste0(tmp.stat, '.scaled.rank')
        tmp.data <- merge(x=tmp.data.7, y=tmp.data, by=c('sample.id', 'clonotype.tag', 'donor.id.tag'))
        return(tmp.data)
    })
    tmp.data.7 <- merge(x=tmp.data.7[[1]], y=tmp.data.7[[2]], by=tmp.cols)
    scaled.rank.data <- copy(tmp.data.7)
    # Take consensus ranks for each peptide pool according to clonotype/donor ID pairs.
    tmp.data.7 <- merge(
        x=tmp.data.7,
        y=unique(meta.data[, .(sample.id, peptide.pool.tag)]),
        by='sample.id',
        all.x=TRUE, all.y=FALSE
    )
    tmp.data.7 <- tmp.data.7[,
        .(
            cd3e.mean=median(cd3e.mean),
            cd3e.median=median(cd3e.median),
            mean.scaled.rank=median(mean.scaled.rank),
            median.scaled.rank=median(median.scaled.rank)
        ),
        by=.(clonotype.tag, donor.id.tag, peptide.pool.tag)
    ]
    # Tidy format according to clonotype/donor ID pairs and peptide pools.
    tmp.cols <- c('cd3e.mean', 'cd3e.median', 'mean.scaled.rank', 'median.scaled.rank')
    tmp.data.7 <- lapply(X=tmp.cols, FUN=function(tmp.col){
        tmp.cols <- c('peptide.pool.tag', 'clonotype.tag', 'donor.id.tag', tmp.col)
        tmp.data <- tmp.data.7[, ..tmp.cols]
        tmp.data[, peptide.pool.tag:=paste0(tmp.col, '.', peptide.pool.tag)]
        tmp.data <- as.data.table(spread(data=tmp.data, key=peptide.pool.tag, value=tmp.col, fill=NA))
        return(tmp.data)
    })
    tmp.data.7 <- Reduce(x=tmp.data.7, f=function(x, y){
        merge(x=x, y=y, by=c('clonotype.tag', 'donor.id.tag'))
    })

    # Merge all data.
    clone.meta <- merge(x=tmp.data.1, y=tmp.data.2, by='clonotype.tag', all=TRUE)
    tmp.data <- merge(x=tmp.data.3, y=tmp.data.4, by=c('clonotype.tag', 'donor.id.tag'), all=TRUE)
    clone.meta <- merge(x=clone.meta, y=tmp.data, by='clonotype.tag')
    tmp.data <- merge(x=tmp.data.6, y=tmp.data.5, by=c('clonotype.tag', 'donor.id.tag'), all=TRUE)
    clone.meta <- merge(x=clone.meta, y=tmp.data, by=c('clonotype.tag', 'donor.id.tag'), all=TRUE)
    clone.meta <- merge(x=clone.meta, y=tmp.data.7, by=c('clonotype.tag', 'donor.id.tag'), all=TRUE)
    # Sanity check
    to.check <- clone.meta[, .(clonotype.tag, donor.id.tag)]
    to.check.1 <- to.check[, .N]
    to.check.2 <- unique(to.check); to.check.2 <- to.check.2[, .N]
    tmp.check <- to.check.1==to.check.2
    if(!tmp.check) stop('Unexpected error. There should be as many lines as there are unique pairs of clonotypes and donor IDs.\n')

    # Final format.
    clone.meta[, to.order:=as.integer(str_extract(string=clonotype.tag, pattern='\\d+$'))]
    setorderv(x=clone.meta, col=c('to.order', 'total.rel.freq'), order=c(1, -1))
    clone.meta[, to.order:=NULL]

    # @ Whole dataset.
    tmp.file.name <- paste0(reports.path, '/ClonotypeBasedTCRData.csv')
    fwrite(file=tmp.file.name, x=clone.meta, quote=TRUE, na=NA)
    # @ For expanded clonotypes only
    tmp.file.name <- paste0(reports.path, '/ClonotypeBasedTCRData_ExpandedOnly.csv')
    fwrite(file=tmp.file.name, x=clone.meta[freq.conf.lvl>0], quote=TRUE, na=NA)
}


### -------------------- Preliminary exploration -------------------- ###
### ------------------------ of hashtag data ------------------------ ###

hto.reports.path <- paste0(reports.path, '/hto_assessment')
if(!dir.exists(hto.reports.path)) dir.create(hto.reports.path)

# Perform only if necessary.
if(is.null(reports.date)){

    # ---> General overlap between modalities.
    tmp.data <- meta.data[, .(
        mod.ovlp.tag=ifelse(
            test=is.na(hto.tag),
            yes='Missed',
            no='Overlapped'
        ),
        sample.id=str_replace_all(
            string=sample.id, pattern='^DICE_TCR\\d+_Hu_|_\\d+D_TCR$', replacement=''
        ),
        chrom.batch.tag
    )]
    tmp.ggplot <- ggplot(data=tmp.data, aes(x=sample.id, fill=mod.ovlp.tag)) +
        geom_bar(position='fill') +
        facet_wrap(~chrom.batch.tag, ncol=9, scales='free_x') +
        labs(x='Sample ID', y='Proportion', fill='Overlap\nclass') +
        theme_bw() + theme(legend.position='bottom', axis.text.x=element_text(angle=30))
    tmp.file.name <- paste0(hto.reports.path, '/OverlapBetweenModalities.pdf')
    pdf(file=tmp.file.name, width=10)
    print(tmp.ggplot)
    dev.off()

    # ---> Rates. For singlets and for general classes.
    tmp.data <- meta.data[, .(
        hto.tag, donor.id.tag, 
        class.tag=ifelse(
            test=is.na(hto.tag) | hto.tag=='Negative',
            yes='Negative',
            no=ifelse(
            hto.tag=='Doublet',
            yes='Doublet',
            no='Singlet'
            )
        ),
        sample.id=str_replace_all(string=sample.id, pattern='^DICE_TCR\\d+_Hu_|_\\d+D_TCR$', replacement=''),
        chrom.batch.tag
    )]
    tmp.ggplot.1 <- ggplot(data=tmp.data, aes(x=sample.id, fill=class.tag)) +
        geom_bar(position='fill') +
        facet_wrap(~chrom.batch.tag, ncol=9, scales='free_x') +
        labs(title='Fraction of general HTO-based deconvolution classes for all barcodes', x='Sample minimal ID', y='Proportion', caption='Barcodes without HTO information entered the class Negative.') +
        theme_bw() + theme(legend.position='bottom', axis.text.x=element_text(angle=30))
    tmp.ggplot.2 <- ggplot(data=tmp.data[!is.na(hto.tag)], aes(x=sample.id, fill=class.tag)) +
        geom_bar(position='fill') +
        facet_wrap(~chrom.batch.tag, ncol=9, scales='free_x') +
        labs(title='Fraction of general HTO-based deconvolution classes for barcodes w/ HTO information', x='Sample minimal ID', y='Proportion', caption='Barcodes without HTO information were excluded here.') +
        theme_bw() + theme(legend.position='bottom', axis.text.x=element_text(angle=30))
    tmp.ggplot.3 <- ggplot(data=tmp.data[!is.na(hto.tag) & hto.tag!='Doublet' & hto.tag!='Negative'], aes(x=sample.id, fill=hto.tag)) +
        geom_bar(position='fill') +
        facet_wrap(~chrom.batch.tag, ncol=9, scales='free_x') +
        labs(title='For singlets, classes\' fractions', x='Sample minimal ID', y='Proportion') +
        theme_bw() + theme(legend.position='bottom', axis.text.x=element_text(angle=30))
    tmp.file.name <- paste0(hto.reports.path, '/ClassRates.pdf')
    pdf(tmp.file.name, width=10)
    print(tmp.ggplot.1)
    print(tmp.ggplot.2)
    print(tmp.ggplot.3)
    dev.off()
}

### -------------------- Preliminary exploration -------------------- ###
### ---------------- of CD3E surface expression data ---------------- ###

cite.reports.path <- paste0(reports.path, '/cite_assessment')
if(!dir.exists(cite.reports.path)) dir.create(cite.reports.path)

# Perform only if necessary.
if(is.null(reports.date)){

    # ---> General overlap between modalities.
    tmp.data <- meta.data[, .(
        mod.ovlp.tag=ifelse(
            test=is.na(cd3e.umi),
            yes='Missed',
            no='Overlapped'
        ),
        sample.id=str_replace_all(
            string=sample.id, pattern='^DICE_TCR\\d+_Hu_|_\\d+D_TCR$', replacement=''
        ),
        chrom.batch.tag
    )]
    tmp.ggplot <- ggplot(data=tmp.data, aes(x=sample.id, fill=mod.ovlp.tag)) +
        geom_bar(position='fill') +
        facet_wrap(~chrom.batch.tag, ncol=9, scales='free_x') +
        labs(x='Sample ID', y='Proportion', fill='Overlap\nclass') +
        theme_bw() + theme(legend.position='bottom', axis.text.x=element_text(angle=30))
    tmp.file.name <- paste0(cite.reports.path, '/OverlapBetweenModalities.pdf')
    pdf(file=tmp.file.name, width=10)
    print(tmp.ggplot)
    dev.off()

    # ---> CD3E UMI counts, comparison between batches.
    tmp.get.plot <- function(batch.var, batch.lab, col.no){
        tmp.ggplot <- ggplot(data=meta.data, aes(x=cd3e.umi+1)) +
            geom_density(linewidth=2, color='black') +
            geom_vline(xintercept=meta.data[, median(cd3e.umi, na.rm=TRUE)], linewidth=2, color='red', linetype='dashed') +
            facet_wrap(facets=paste0('~', batch.var), ncol=col.no) +
            scale_x_log10() +
            labs(title=paste0('Barcode-specific CD3E UMI expression distribution\nSplit according to ', batch.lab), x='log10(CD3E UMI counts + 1)', y='Density') +
            theme_bw()
        return(tmp.ggplot)
    }
    tmp.ggplot.1 <- tmp.get.plot(batch.var='seq.batch.tag', batch.lab='sequencing run', col.no=5)
    tmp.ggplot.2 <- tmp.get.plot(batch.var='chrom.batch.tag', batch.lab='FACS-sorting batch', col.no=9)
    tmp.ggplot.3 <- tmp.get.plot(batch.var='sample.id', batch.lab='sequencing library', col.no=11)
    tmp.file.name <- paste0(cite.reports.path, '/CD3EExpression_CompBetweenBatches.pdf')
    pdf(file=tmp.file.name, width=15)
    print(tmp.ggplot.1)
    print(tmp.ggplot.2)
    print(tmp.ggplot.3)
    dev.off()

    # ---> Exploration of rank scaling approaches.
    # Comparison between stats for rank scaling.
    tmp.ggplot.1 <- ggplot(data=scaled.rank.data, aes(x=mean.scaled.rank, y=median.scaled.rank)) +
        geom_point(size=0.3) +
        geom_smooth(method=lm, se=FALSE) +
        ggpubr::stat_cor(method='pearson') +
        labs(title='Comparison of approaches for rank scaling', x='Mean-based approach', y='Median-based approach') +
        theme_bw()
    # For clonotypes, comparison between scaled ranks and absolute frequencies.
    tmp.ggplot.2 <- ggplot(data=scaled.rank.data, aes(x=mean.scaled.rank, y=abs.freq)) +
        geom_point(size=0.3) +
        geom_smooth(method=lm, se=FALSE) +
        ggpubr::stat_cor(method='pearson') +
        scale_y_log10() +
        labs(title='Comparison between mean-based scaled ranks and absolute freq.', x='Mean-based scaled ranks', y='Clone absolute frequency (log10)') +
        theme_bw()
    tmp.ggplot.3 <- ggplot(data=scaled.rank.data, aes(x=median.scaled.rank, y=abs.freq)) +
        geom_point(size=0.3) +
        geom_smooth(method=lm, se=FALSE) +
        ggpubr::stat_cor(method='pearson') +
        scale_y_log10() +
        labs(title='Comparison between median-based scaled ranks and absolute freq.', x='Median-based scaled ranks', y='Clone absolute frequency (log10)') +
        theme_bw()
    # Output
    tmp.file.name <- paste0(cite.reports.path, '/RankScalingMethodExploration.pdf')
    pdf(file=tmp.file.name)
    print(tmp.ggplot.1)
    print(tmp.ggplot.2)
    print(tmp.ggplot.3)
    dev.off()

    # ---> CD3E UMI counts, contrast between cells expressing singleton clonotypes and cells expressing rest of the clonotypes.
    # Multiple tests.
    # tmp.vals <- meta.data[, sort(unique(sample.id))]
    # test.data <- lapply(X=tmp.vals, FUN=function(tmp.val){
    #     tmp.data <- meta.data[
    #         sample.id==tmp.val & !is.na(gen.clonotype.status),
    #         .(
    #             cd3e.umi=cd3e.umi,
    #             pseudo.count=log2(cd3e.umi+1),
    #             singleton=ifelse(test=gen.clonotype.status=='Singleton', yes='Singleton', no='Rest'),
    #             peptide.pool=unique(peptide.pool.tag)
    #         )
    #     ]
    #     if(tmp.data[, .N]<10) return(NA)
    #     tmp.data <- wilcox.test(
    #         x=tmp.data[singleton=='Singleton', cd3e.umi],
    #         y=tmp.data[singleton=='Rest', cd3e.umi],
    #         paired=FALSE
    #     )
    #     tmp.data <- data.table(p.val=tmp.data$p.value)
    #     return(tmp.data)
    # })
    # names(test.data) <- tmp.vals
    # test.data <- test.data[!is.na(test.data)]
    # test.data <- rbindlist(l=test.data, use.names=TRUE, idcol='sample.id')
    # test.data[, adj.p:=p.adjust(p=p.val, method='BH')]
    # # Plots.
    # tmp.cols <- c('Singleton'='#8b0000', 'Rest'='#ffbf00')
    # tmp.file.name <- paste0(cite.reports.path, '/UMIExp_Singleton-vs-Rest.pdf')
    # pdf(file=tmp.file.name)
    # for(tmp.val in tmp.vals){
    #     tmp.data.1 <- meta.data[
    #         sample.id==tmp.val & !is.na(gen.clonotype.status),
    #         .(
    #             cd3e.umi=cd3e.umi,
    #             pseudo.count=log2(cd3e.umi+1),
    #             singleton=ifelse(test=gen.clonotype.status=='Singleton', yes='Singleton', no='Rest'),
    #             peptide.pool=unique(peptide.pool.tag)
    #         )
    #     ]
    #     tmp.data.2 <- tmp.data.1[, .(pseudo.count=median(pseudo.count)), by=singleton]
    #     tmp.data.3 <- paste0('Adj. P-value: ', round(x=test.data[sample.id==tmp.val, adj.p], digits=5))
    #     tmp.ggplot <- ggplot(data=tmp.data.1, aes(color=singleton)) +
    #         geom_density(linewidth=1.5, aes(x=pseudo.count, color=singleton), alpha=0.5) +
    #         scale_y_continuous(expand=c(0, NA)) +
    #         scale_color_manual(values=tmp.cols) +
    #         labs(
    #             title=paste0('Sequencing library: ', tmp.val, '\nPeptide pool: ', tmp.data.1[, unique(peptide.pool)]),
    #             x='log2(CD3E UMI counts + 1)', y='Density', color='', caption=tmp.data.3
    #         ) +
    #         theme_bw()
    #     for(tmp.group in names(tmp.cols)){
    #         tmp.ggplot <- tmp.ggplot + geom_vline(xintercept=tmp.data.2[singleton==tmp.group, pseudo.count], linewidth=2.5, linetype='dashed', color=tmp.cols[tmp.group])
    #     }
    #     print(tmp.ggplot)
    # }
    # dev.off()
}


### -------------------- Preliminary exploration -------------------- ###
### -------------------------- of TCR data -------------------------- ###

tcr.reports.path <- paste0(reports.path, '/tcr_assessment')
if(!dir.exists(tcr.reports.path)) dir.create(tcr.reports.path)

# # ---> Potentially public TCRs that are likely to be mislabeling errors.
# tmp.data <- meta.data[
#     !is.na(donor.id.tag),
#     .(donor.count=uniqueN(donor.id.tag)),
#     clonotype.tag
# ]

# ---> Define sets of peptide pools to apply process to.
pp.sets <- list(
    `peptide_pools-ALL`=meta.data[, unique(specificity.class.tag)]
    # `peptide_pools-SARS-CMV-FLU`=c('SARS-CoV-2', 'CMV', 'FLU')
)

for(pp.set in names(pp.sets)){
    # ---> Define metadata and sub-reports path.
    pp.tcr.reports.path <- paste0(tcr.reports.path, '/', pp.set)
    if(!dir.exists(pp.tcr.reports.path)) dir.create(pp.tcr.reports.path)
    this.pp.set <- pp.sets[[pp.set]]
    sub.meta.data <- meta.data[specificity.class.tag %in% this.pp.set]

    # ---> Count summary.
    # Collect data to plot.
    tmp.data.1 <- sub.meta.data[
        !is.na(donor.id.tag),
        .(
            Cells=.N,
            Clonotypes=uniqueN(clonotype.tag)
        ),
        by=.(donor.id.tag, specificity.class.tag)
    ]
    tmp.data.2 <- as.data.table(gather(data=tmp.data.1, key='key', value='count', `Cells`, `Clonotypes`))
    # General summary stats (to be shown on plot, too)
    tmp.stats <- tmp.data.2[,
        .(
            quartile.25=quantile(x=count, probs=0.25),
            quartile.50=median(count),
            quartile.75=quantile(x=count, probs=0.75)
        ),
        by=key
    ]
    # @ Plots showing stats on a peptide pool basis.
    tmp.ggplot.1 <- ggplot(data=tmp.data.2, aes(x=specificity.class.tag, y=count, fill=specificity.class.tag)) +
        geom_boxplot(width=0.5, alpha=0.4) +
        geom_jitter(size=0.12, width=0.25) +
        geom_hline(data=tmp.stats, linetype='dashed', color='red', aes(yintercept=quartile.25)) +
        geom_hline(data=tmp.stats, linetype='dashed', color='red', aes(yintercept=quartile.50)) +
        geom_hline(data=tmp.stats, linetype='dashed', color='red', aes(yintercept=quartile.75)) +
        facet_wrap(facets=~key, nrow=2, scales='free_y') +
        scale_y_log10() +
        coord_cartesian(ylim=c(10, NA)) +
        scale_fill_manual(values=peptide.pool.cols) +
        labs(title='Cell and clonotype count per donor', x='ARTE peptide pool', y='Count (log10)', fill='', caption='Zoomed on y axis.\nDashed lines represent the 25, 50 and 75% quartiles of the general distribution (i.e., across peptide pools).') +
        theme_bw() + theme(legend.position='none', axis.text.x=element_text(angle=90))
    tmp.ggplot.2 <- ggplot(data=tmp.data.2, aes(x=specificity.class.tag, y=count, fill=specificity.class.tag)) +
        geom_boxplot(width=0.5, alpha=0.4) +
        geom_jitter(size=0.12, width=0.25) +
        geom_hline(data=tmp.stats, linetype='dashed', color='red', aes(yintercept=quartile.25)) +
        geom_hline(data=tmp.stats, linetype='dashed', color='red', aes(yintercept=quartile.50)) +
        geom_hline(data=tmp.stats, linetype='dashed', color='red', aes(yintercept=quartile.75)) +
        facet_wrap(facets=~key, nrow=2, scales='free_y') +
        coord_cartesian(ylim=c(0, 5000)) +
        scale_fill_manual(values=peptide.pool.cols) +
        labs(title='Cell and clonotype count per donor', x='ARTE peptide pool', y='Count', fill='', caption='Zoomed on y axis.\nDashed lines represent the 25, 50 and 75% quartiles of the general distribution (i.e., across peptide pools).') +
        theme_bw() + theme(legend.position='none', axis.text.x=element_text(angle=90))
    # @ Plots comparing number of clonotypes obtained when compared with the number of cells obatined.
    tmp.ggplot.3 <- ggplot(data=tmp.data.1, aes(x=Cells, y=Clonotypes)) +
        geom_point(aes(color=specificity.class.tag)) +
        geom_smooth(formula=y~x, se=FALSE, method='lm', color='black', linewidth=0.6) +
        ggpubr::stat_cor(method='pearson') +
        scale_color_manual(values=peptide.pool.cols) +
        scale_x_continuous(expand=c(0, 0)) +
        scale_y_continuous(expand=c(0, 0)) +
        labs(x='Cell count', y='Clonotype count', color='Peptide pool') +
        theme_bw()
    # Output plot.
    tmp.file.name <- paste0(pp.tcr.reports.path, '/CountSummary.pdf')
    pdf(tmp.file.name)
    print(tmp.ggplot.1)
    print(tmp.ggplot.2)
    print(tmp.ggplot.3)
    dev.off()

    # ---> Clonal expansion per peptide pool.
    # Retrieve data.
    tmp.data.1 <- sub.meta.data[
        !is.na(donor.id.tag),
        .(clone.size=.N),
        by=.(clonotype.tag, specificity.class.tag, donor.id.tag)
    ]
    up.thold <- quantile(x=tmp.data.1[, clone.size], probs=0.99)
    tmp.data.1[clone.size>up.thold, clone.size:=up.thold+1]
    # Colors and captions
    tmp.caption.1 <- paste0('The bar for the maximum value in the x axis represents the frequency of clonotypes\nwith size equal to or greater than ', up.thold+1, ' (which happens to be the 99% quartile).\nNotice that the y axis scale was trimmed for better viz, but the number of clonotypes\n w/ size equal to 1 is ', tmp.data.1[clone.size==1, .N], '. See next page for a look without trimming.')
    tmp.caption.2 <- paste0('The bar for the maximum value in the x axis represents the frequency of clonotypes\nwith size equal to or greaterthan ', up.thold+1, ' (which happens to be the 99% quartile).')
    # Plot.
    tmp.ggplot.1 <- ggplot(data=tmp.data.1, aes(x=clone.size, fill=specificity.class.tag)) +
        geom_histogram(alpha=0.7, bins=up.thold+1, linewidth=1, color='black') +
        facet_wrap(facets=~specificity.class.tag, ncol=1) +
        scale_y_continuous(expand=c(0, 10), n.breaks=3) +
        coord_cartesian(ylim=c(0, 1e4)) +
        scale_fill_manual(values=peptide.pool.cols) +
        labs(x='Clone size', y='Absolute frequency', fill='', caption=tmp.caption.1) +
        theme_bw() + theme(legend.position='none')
    tmp.ggplot.2 <- ggplot(data=tmp.data.1, aes(x=clone.size, fill=specificity.class.tag)) +
        geom_histogram(alpha=0.7, bins=up.thold+1, linewidth=1, color='black') +
        facet_wrap(facets=~specificity.class.tag, ncol=1) +
        scale_y_continuous(expand=c(0, 10), n.breaks=3) +
        scale_fill_manual(values=peptide.pool.cols) +
        labs(x='Clone size', y='Absolute frequency', fill='', caption=tmp.caption.2) +
        theme_bw() + theme(legend.position='none')
    tmp.file.name <- paste0(pp.tcr.reports.path, '/ExpansionSummary_PerPeptidePool.pdf')
    pdf(file=tmp.file.name, height=12, width=5)
    print(tmp.ggplot.1)
    print(tmp.ggplot.2)
    dev.off()

    # ---> Clonotype sharing between antigens.
    # For expanded clonotypes only.
    tmp.data <- sub.meta.data[
        !is.na(donor.id.tag),
        .(
            clonotype.tag=paste(clonotype.tag, donor.id.tag, sep='_'),
            peptide.pool.tag
        )
    ]
    tmp.data <- tmp.data[, .(clone.size=.N), by=.(clonotype.tag, peptide.pool.tag)]
    uniq.viruses <- tmp.data[, unique(peptide.pool.tag)]
    # Define plot.
    tmp.file.name <- paste0(pp.tcr.reports.path, '/CloneSharingAmongPeptidePools.pdf')
    pdf(tmp.file.name, height=10, width=10)
    # Plot according to size threshold.
    size.vals <- 1:10
    these.cols <- colorRampPalette(colors=c('#00ffff', '#0000ff'))(10)
    for(size.thold in size.vals){
        data.to.plot <- lapply(X=uniq.viruses, FUN=function(tmp.virus){
            tmp.data[peptide.pool.tag==tmp.virus & clone.size>=size.thold, unique(clonotype.tag)]
        })
        names(data.to.plot) <- uniq.viruses
        print(upset(
            fromList(data.to.plot),
            nsets=length(data.to.plot),
            main.bar.color=these.cols[size.thold], att.color='black',
            show.numbers='no',
            mainbar.y.label=paste0('Intersection size for clonotypes w/ size equal to or larger than ', size.thold)
        ))
    }
    dev.off()
    # Define plot for all intersections.
    tmp.file.name <- paste0(pp.tcr.reports.path, '/CloneSharingAmongPeptidePools_AllInts.pdf')
    pdf(tmp.file.name, height=10, width=10)
    # Plot according to size threshold.
    size.vals <- 1:10
    these.cols <- colorRampPalette(colors=c('#00ffff', '#0000ff'))(10)
    for(size.thold in size.vals){
        data.to.plot <- lapply(X=uniq.viruses, FUN=function(tmp.virus){
            tmp.data[peptide.pool.tag==tmp.virus & clone.size>=size.thold, unique(clonotype.tag)]
        })
        names(data.to.plot) <- uniq.viruses
        print(upset(
            fromList(data.to.plot),
            nsets=length(data.to.plot),
            nintersects=NA,
            main.bar.color=these.cols[size.thold], att.color='black',
            show.numbers='no',
            mainbar.y.label=paste0('Intersection size for clonotypes w/ size equal to or larger than ', size.thold)
        ))
    }
    dev.off()

    # ---> Fraction accounted for by each reactivity for clonotypes overall.
    tmp.data <- sub.meta.data[, .(clone.count=uniqueN(clonotype.tag)), by=specificity.class.tag]
    # tmp.data[specificity.class.tag=='Cross-reactive', specificity.class.tag:='pR']
    tmp.ggplot <- ggplot(data=tmp.data[!is.na(specificity.class.tag)], aes(x='', y=clone.count, fill=specificity.class.tag)) +
        geom_bar(stat='identity', color='black', linewidth=0.8, alpha=0.8) +
        scale_fill_manual(values=peptide.pool.cols) +
        coord_polar("y", start=0) +
        labs(x='', y='', fill='Specificity') +
        theme_minimal() + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
    tmp.file.name <- paste0(pp.tcr.reports.path, '/ReactivityFractionForClones.pdf')
    pdf(file=tmp.file.name)
    print(tmp.ggplot)
    dev.off()

    # ---> Clonotype count summary for clonotypes captured w/ confidence of certain clonality.
    # Collect data to plot.
    tmp.data.1 <- sub.meta.data[
        !is.na(donor.id.tag),
        .(
            size=uniqueN(barcode)
        ),
        by=.(donor.id.tag, specificity.class.tag, clonotype.tag)
    ]
    size.conf.vals <- 1:10
    tmp.data.2 <- lapply(X=size.conf.vals, FUN=function(conf.val){
        tmp.data <- tmp.data.1[,
            .(
            clone.frac=.SD[
                size>=conf.val, .N
            ] / .N
            ),
            by=.(donor.id.tag, specificity.class.tag)
        ]
        return(tmp.data)
    })
    names(tmp.data.2) <- as.character(size.conf.vals)
    tmp.data.2 <- rbindlist(l=tmp.data.2, use.names=TRUE, idcol='conf.level')
    tmp.data.2$conf.level <- factor(x=tmp.data.2$conf.level, levels=as.character(size.conf.vals))
    # General summary stats per confidence level
    tmp.stats <- tmp.data.2[
        clone.frac!=0,
        .(
            quartile.50=median(clone.frac),
            condition.no=.N,
            group='median'
        ),
        by=conf.level
    ]
    tmp.caption.1 <- paste0('Each dot represents a condition where a condition is defined by the combination of a donor ID and a peptide pool.\nDots are colored according to condition\'s peptide pool.\nBoxplot is defined only for conditions w/ fractions greater than 0.\nNumber of conditions w/ fractions greater than 0 is depicted on top of the boxes.')
    tmp.caption.2 <- paste0('Same plot as above, but cutoff equals 1 was excluded and y scale is linear for better viz.\nEach dot represents a condition where a condition is defined by the combination of a donor ID and a peptide pool.\nDots are colored according to condition\'s peptide pool.\nBoxplot is defined only for conditions w/ fractions greater than 0.\nNumber of conditions w/ fractions greater than 0 is depicted on top of the boxes.')
    # Plot
    tmp.ggplot.1 <- ggplot(data=tmp.data.2[clone.frac>0], aes(x=conf.level, y=clone.frac)) +
        geom_boxplot(width=0.5, alpha=0.4, outlier.shape=NA) +
        geom_jitter(size=0.12, width=0.25, aes(color=specificity.class.tag)) +
        geom_line(data=tmp.stats, linewidth=0.8, color='#ffff00', aes(y=quartile.50, group=group)) +
        geom_text(data=tmp.stats, color='black', aes(y=0.9, label=condition.no)) +
        scale_y_log10() +
        scale_color_manual(values=peptide.pool.cols) +
        labs(title='Fraction of total clonotypes maintained after setting a clone size confidence\nthreshold to call antigen-reactive cells', x='Clone size confidence threshold', y='Fraction of total clonotypes (log10)', color='Peptide pool', caption=tmp.caption.1) +
        theme_bw() + theme(legend.position='bottom')
    tmp.ggplot.2 <- ggplot(data=tmp.data.2[clone.frac>0 & conf.level!='1'], aes(x=conf.level, y=clone.frac)) +
        geom_boxplot(width=0.5, alpha=0.4, outlier.shape=NA) +
        geom_jitter(size=0.12, width=0.25, aes(color=specificity.class.tag)) +
        geom_line(data=tmp.stats[conf.level!='1'], linewidth=0.8, color='#ffff00', aes(y=quartile.50, group=group)) +
        geom_text(data=tmp.stats[conf.level!='1'], color='black', aes(y=0.2, label=condition.no)) +
        coord_cartesian(ylim=c(NA, 0.2)) +
        scale_color_manual(values=peptide.pool.cols) +
        labs(title='Fraction of total clonotypes maintained after setting a clone size confidence\nthreshold to call antigen-reactive cells', x='Clone size confidence threshold', y='Fraction of total clonotypes (linear)', color='Peptide pool', caption=tmp.caption.2) +
        theme_bw() + theme(legend.position='bottom')
    tmp.file.name <- paste0(pp.tcr.reports.path, '/PickingSizeConfCutoff.pdf')
    pdf(tmp.file.name)
    print(tmp.ggplot.1)
    print(tmp.ggplot.2)
    dev.off()
}

cat('\n\n')
### --------------------------------- END -------------------------------- ###
# ---> Print session info.
sessionInfo()

cat('PROGRAM FINISHED!\n\n')
