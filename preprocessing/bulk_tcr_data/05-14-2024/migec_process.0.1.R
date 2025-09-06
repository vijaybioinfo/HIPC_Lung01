############    -----   Sample sheet sanity checks    -----    ############
############    -------------   NV127 run    --------------    ############

### -------------------------- Description -------------------------- ###
# In general, we must cross-compare the content from the sample sheet and the TCR barcode info files. The following are the checks to perform:
#       Sequence library names match. There can't be library names in the TCR barcode info that left undefined in the sample sheet.
#       There has to be a unique set of illumina barcodes defined per sequencing library across rows in the TCR barcode info file.
#       There have to be unique TCR barcode IDs across TCR samples within each sequencing library. This applies to both master and slave barcoeds.
#       There has to be a 1:1 correspondence between TCR barcode ID and sequence across samples for both barcode types, master and slave.
#       Illumina seqs. must match between the sample sheet and the barcode info.


### -------------------------- Dependencies ------------------------- ###
library(data.table)
library(stringr)
library(tidyr)
library(rjson)
library(ggplot2)


### --------------------------- Arguments --------------------------- ###

# ---> General definitions.
seq.date <- '05-14-2024'
nv.run.id <- 'NV135'
# @ Master and slave barcode constructs.
bc.seqs <- list(
    'TCRA'=c(
        master='cagtggtatcaacgcagagtNNNNtNNNNtNNNNtct',
        slave='NNNN<BARCODE>gggtcagggttctggatat'
    ),
    'TCRB'=c(
        master='cagtggtatcaacgcagagtNNNNtNNNNtNNNNtct',
        slave='NNNN<BARCODE>acacsttkttcaggtcctc'
    )
)
# @ About job definition.
job.head.pttn <- '#!/bin/sh
#SBATCH --job-name=MIGEC_{ROUTINE}
#SBATCH --output={OUT_FILE}
#SBATCH --error={ERR_FILE}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem={MEM}g
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00'
# ---> Path definitions.
gen.seq.path <- paste0('/path/to/user/sequencing_data/', seq.date)
bcl2f.path <- paste0(gen.seq.path, '/bcl2fastq')
fastq.path <- paste0(bcl2f.path, '/', nv.run.id, '/outs/fastq_path')
if(!dir.exists(fastq.path)) stop('Output directory from bcl2fastq does not exist.')
stats.path <- paste0(fastq.path, '/Stats')
data.path <- paste0(bcl2f.path, '/data/experiment_layout')
if(!dir.exists(data.path)) stop('Sequencing directory does not exist.\n')
migec.path <- paste0(gen.seq.path, '/migec')
if(!dir.exists(migec.path)) dir.create(migec.path)
jobs.path <- paste0(gen.seq.path, '/jobs_scripts/migec')
if(!dir.exists(jobs.path)) dir.create(jobs.path)
# ---> File and sepcific definitions.
bc.info.file <- paste0(data.path, '/', nv.run.id, '_barcode_info_for_migec.csv')
dec.stats.file <- paste0(stats.path, '/Stats.json')
migec.tool <- '/home/vfajardo/bin/MIGEC/migec-1.2.9.jar' # Not a file tho.
essential.files <- c(
  bc.info.file,
  dec.stats.file
)
essential.files <- essential.files[!unlist(lapply(X=essential.files, FUN=file.exists))]
if(length(essential.files) > 0) stop(paste0('Next essential files are not appropriately defined -make sure you\'ve got their names right-:\n', paste0(essential.files, collapse='\n'), '\n'))



### -------------------------- Load data & -------------------------- ###
### ---------------------- Data preprocessing ----------------------- ###

# ---> Barcode info file.
bc.info <- fread(file=bc.info.file)
# bc.info <- bc.info[sample.prj=='BulkTCR021_restim']

# ---> Determine cahin type for each sequencing library.
bc.info[, chain.type:=str_extract(string=seq.lib.name, pattern='ACR|BCR')]
tmp.vals <- c(`ACR`='TCRA', `BCR`='TCRB')
bc.info[, chain.type:=tmp.vals[chain.type]]

# ---> Retrieve names of the fastq files.
# List valid sample projects
fastq.info <- basename(list.dirs(path=fastq.path, full.names=TRUE, recursive=FALSE))
fastq.info <- fastq.info[fastq.info %in% bc.info[, sample.prj]]
# Sanity check. All sample projects previously defined under the TCR barcode info file should have fastq files available.
tmp.check <- bc.info[, all(sample.prj %in% fastq.info)]
if(!tmp.check) stop('Unexpected error!\n')
# Capture fastq files per sample project
tmp.data <- lapply(X=fastq.info, FUN=function(x){
    tmp.data <- data.table(
        fastq.file=list.files(path=paste0(fastq.path, '/', x), pattern='fastq', all.files=FALSE, full.names=TRUE, recursive=FALSE)
    )
    return(tmp.data)
})
names(tmp.data) <- fastq.info
fastq.info <- rbindlist(l=tmp.data, use.names=TRUE, idcol='sample.prj')
# Read type.
fastq.info[,
    read.type:=str_replace_all(
        string=str_extract(string=fastq.file, pattern='_R[12]_'),
        pattern='_', replacement=''
    )
]
# Fastq file basename.
fastq.info[,
    fastq.base:=str_replace(
        string=basename(fastq.file),
        pattern=paste0('_', read.type, '_001.fastq(.gz)*$'), replacement=''
    )
]
# Lane info.
fastq.info[,
    lane.no:=str_extract(string=fastq.base, pattern='L\\d{3}')
]
# Extract sequencing sample name.
fastq.info[,
    seq.lib.name:=str_replace(
        string=fastq.base,
        pattern=paste0('_S\\d+_', lane.no), replacement=''
    )
]
# Determine R1 to R2 relationship for each sample name.
fastq.info <- as.data.table(spread(
    data=fastq.info, key=read.type, value=fastq.file
))
# Sanity check. All sequencing libraries previously defined under the TCR barcode info file should have fastq files available.
tmp.check <- bc.info[, all(seq.lib.name %in% fastq.info[, seq.lib.name])]
if(!tmp.check) stop('Unexpected error!\n')

# ---> Deconvolution stats file.
dec.stats <- fromJSON(file=dec.stats.file)
dec.stats <- dec.stats[['ConversionResults']]
names(dec.stats) <- sapply(X=dec.stats, FUN=function(x) return(x$LaneNumber))
tmp.data <- lapply(X=dec.stats, FUN=function(lane.info){
    lane.info <- lane.info$DemuxResults
    tmp.data <- lapply(X=lane.info, FUN=function(sample.info){
        tmp.data <- data.table(
            seq.lib.name=sample.info$SampleId,
            total.reads=sample.info$NumberReads
        )
        return(tmp.data)
    })
    tmp.data <- rbindlist(l=tmp.data, use.names=TRUE)
    return(tmp.data)
})
names(tmp.data) <- names(dec.stats)
dec.stats <- rbindlist(l=tmp.data, use.names=TRUE, idcol='lane.no')
# Add read count to general barcode info table.
# dec.stats[!sample.id%in%bc.info[, seq.lib.name]]
tmp.data <- dec.stats[, .(total.reads=sum(total.reads)), by=seq.lib.name]
tmp.check <- bc.info[, all(seq.lib.name%in%tmp.data[, seq.lib.name])]
if(!tmp.check) stop('Could not retrieve deconvolution stats for certain libraries listed in the barcode info table.\n')
bc.info <- merge(x=bc.info, y=tmp.data, by='seq.lib.name', all.x=TRUE, all.y=FALSE)


### ------------------------- Main program -------------------------- ###
### ------------------------- MIGEC process ------------------------- ###


### ---------------------- Basic job definition --------------------- ###

# ---> Define the basic structure of a job given a job output directory and a module name. Then, run the job from command line.

run.job <- function(module.name, module.path, line.rest, memory=100){
    # ---> Define and run job.
    # Job files.
    job.file <- paste0(tmp.path, '/', module.name, '_', nv.run.id, '.job.sh')
    out.file <- paste0(tmp.path, '/', module.name, '_', nv.run.id, '.out.txt')
    err.file <- paste0(tmp.path, '/', module.name, '_', nv.run.id, '.err.txt')
    # Custom job jeader
    job.head <- job.head.pttn
    job.head <- str_replace(string=job.head, pattern='\\{ROUTINE\\}', replacement=module.name)
    job.head <- str_replace(string=job.head, pattern='\\{OUT_FILE\\}', replacement=out.file)
    job.head <- str_replace(string=job.head, pattern='\\{ERR_FILE\\}', replacement=err.file)
    job.head <- str_replace(string=job.head, pattern='\\{ERR_FILE\\}', replacement=err.file)
    job.head <- str_replace(string=job.head, pattern='\\{MEM\\}', replacement='100')
    # Job content.
    migec.cmd <- paste(
        'source activate MIGEC\n',
        'java -jar ', migec.tool, module.name,
        line.rest, sep=' '
    )
    # Create and save job file.
    tmp.cmd <- paste0(job.head, '\n\n', migec.cmd)
    write(x=tmp.cmd, file=job.file)
    tmp.cmd <- paste0('sbatch ', job.file)
    system(command=tmp.cmd)
}


### ------------------- TCR sample demultiplexing ------------------- ###

# ---> Predefinitions.
tmp.path <- paste0(jobs.path, '/checkout')
if(!dir.exists(tmp.path)) dir.create(tmp.path)

# ---> Create barcode file.
tmp.data.1 <- bc.info[, .(sample.prj, seq.lib.name, tcr.lib.id, tcr.lib.name, chain.type, master.seq, slave.seq)]
tmp.data.2 <- fastq.info[, .(sample.prj, seq.lib.name, R1, R2)]
tmp.data <- merge(x=tmp.data.1, y=tmp.data.2, by=c('sample.prj', 'seq.lib.name'), all=TRUE, allow.cartesian=TRUE)
tmp.check <- any(is.na(tmp.data))
if(tmp.check) stop('Unexpected error!\n')
tmp.data <- tmp.data[, .(
    sample.id=paste(tcr.lib.id, tcr.lib.name, chain.type, sep='.'),
    master.seq, slave.seq,
    R1, R2,
    chain.type
)]
# Sanity check: There should be unique entries in the table.
tmp.check <- tmp.data[, .N]
tmp.data <- unique(tmp.data)
tmp.check <- tmp.check == tmp.data[, .N]
if(!tmp.check) stop('Unexpected error!\n')
# Final definition of 5' and 3' fragment sequence patterns.
for(tmp.chain in c('TCRA', 'TCRB')){
    tmp.data[
        chain.type==tmp.chain,
        master.seq:=paste0(master.seq, bc.seqs[[tmp.chain]]['master'])
    ]
    tmp.data[
        chain.type==tmp.chain,
        slave.seq:=str_replace(
            string=bc.seqs[[tmp.chain]]['slave'],
            pattern='<BARCODE>',
            replacement=slave.seq
        )
    ]
}
# Final format and save.
# colnames(tmp.data) <- c(
#     'Sample ID', 'Master barcode sequence', 'Slave barcode sequence', 'Read#1 FASTQ', 'Read#2 FASTQ'
# )
barcodes.file <- paste0(migec.path, '/BarcodesFile_Checkout.tsv')
fwrite(file=barcodes.file, x=tmp.data, sep='\t', quote=FALSE, col.names=FALSE)

# ---> Define and run job.
run.job(
    module.name='CheckoutBatch',
    module.path=tmp.path,
    line.rest=paste0(' -cute ', barcodes.file, ' ', migec.path, '/checkout/'),
    memory=250
)

# ---> Results exploration.
# Read results and preprocess.
tmp.file.name <- paste0(migec.path, '/checkout/checkout.log.txt')
file.exists(tmp.file.name)
tmp.data <- fread(file=tmp.file.name, header=TRUE)
colnames(tmp.data) <- c('fastq.r1', 'fastq.r2', 'tcr.lib.name', 'reads.w.master.bc', 'reads.w.both.bcs', 'ovlped.reads')
tmp.data[, ovlped.reads:=NULL]
tmp.data[,
    seq.lib.name:=str_replace(
        string=basename(fastq.r1),
        pattern='_[^_]+_L\\d+_R1_\\d+.fastq.gz',
        replacement=''
    )
]
tmp.data[,
    lane.no:=str_extract(
        string=str_extract(string=fastq.r1, pattern='_L\\d+'),
        pattern='\\d$'
    )
]
tmp.data <- tmp.data[, .(lane.no, seq.lib.name, tcr.lib.name, reads.w.master.bc, reads.w.both.bcs)]
# Sanity check. Total number of reads should greatly approximate the total number of reads from the deconvolution stats.
tmp.data.1 <- tmp.data[, .(migec.reads=sum(reads.w.master.bc)), by=.(seq.lib.name, lane.no)]
tmp.data.1 <- merge(x=tmp.data.1, y=dec.stats, by=c('lane.no', 'seq.lib.name'), all.x=TRUE)
tmp.check <- tmp.data.1[, all(migec.reads==total.reads)]
if(!tmp.check) stop('Inconsistencies between MIGEC read counts and reads provided by bcl2fastq might have been encountered.\n')
# Final reports.
#       Per TCR library.
tmp.file.name <- paste0(migec.path, '/checkout/ReportPerTCRLibrary.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA)
#       Per sequencing library.
tmp.data <- merge(x=tmp.data, y=dec.stats, by=c('lane.no', 'seq.lib.name'), all.x=TRUE, all.y=FALSE)
colnames(tmp.data)[colnames(tmp.data)=='total.reads'] <- 'tcr.lib.total.read'
tmp.data[, `:=`(
    reads.w.master.bc.percent=reads.w.master.bc/tcr.lib.total.read,
    reads.w.both.bcs.percent=reads.w.both.bcs/tcr.lib.total.read
)]
tmp.data.1 <- tmp.data[, 
    .(
        undefined.reads=.SD[
            tcr.lib.name %like% '^undef', sum(reads.w.master.bc.percent)
        ],
        reads.assigned.to.master.bc=.SD[
            !tcr.lib.name %like% '^undef', sum(reads.w.master.bc.percent)
        ],
        reads.assigned.to.both.bcs=.SD[
            !tcr.lib.name %like% '^undef', sum(reads.w.both.bcs.percent)
        ]
    ),
    by=.(lane.no, seq.lib.name)
]
tmp.file.name <- paste0(migec.path, '/checkout/ReportPerSeqLibrary.csv')
fwrite(file=tmp.file.name, x=tmp.data.1, na=NA)
# Viz.
tmp.data.1 <- as.data.table(gather(data=tmp.data.1, key='read.group', value='read.fraction', -`lane.no`, -`seq.lib.name`))
to.plot <- tmp.data.1[read.group=='undefined.reads' | read.group=='reads.assigned.to.master.bc']
tmp.vals <- c('undefined.reads'='Undefined', 'reads.assigned.to.master.bc'='Assigned to master barcodes')
to.plot$read.group <- factor(x=tmp.vals[to.plot$read.group], levels=tmp.vals)
tmp.ggplot.1 <- ggplot(data=to.plot, aes(x=read.group, y=read.fraction)) +
    geom_boxplot(fill='lightblue', linewidth=1.5, color='black', outlier.shape=NA) +
    geom_jitter(size=1, color='black', width=0.2) +
    labs(title='Reads assigned to master barcodes versus the ones left unassigned', x='', y='Fraction of reads') +
    theme_bw()
to.plot <- tmp.data.1[read.group=='reads.assigned.to.master.bc' | read.group=='reads.assigned.to.both.bcs']
tmp.vals <- c('reads.assigned.to.master.bc'='Assigned to master barcodes', 'reads.assigned.to.both.bcs'='Assigned to both barcodes')
to.plot$read.group <- factor(x=tmp.vals[to.plot$read.group], levels=tmp.vals)
to.plot[, tmp.group:=paste(seq.lib.name, lane.no, sep='.')]
tmp.ggplot.2 <- ggplot(data=to.plot, aes(x=read.group, y=read.fraction)) +
    geom_boxplot(fill='lightblue', linewidth=1.5, color='black', outlier.shape=NA) +
    geom_jitter(size=1, color='black', width=0.05) +
    geom_line(aes(group=tmp.group), color='black', linewidth=0.5) +
    labs(title='Reads assigned either only to master barcodes or to both types of barcodes', x='', y='Fraction of reads') +
    theme_bw()
tmp.file.name <- paste0(migec.path, '/checkout/ReportSummary.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot.1)
print(tmp.ggplot.2)
dev.off()


### --------------------- Demultiplexing stats ---------------------- ###

# ---> Predefinitions.
tmp.path <- paste0(jobs.path, '/histogram')
if(!dir.exists(tmp.path)) dir.create(tmp.path)

# ---> Define and run job.
run.job(
    module.name='Histogram',
    module.path=tmp.path,
    line.rest=paste0(migec.path, '/checkout/', ' ', migec.path, '/histogram/'),
    memory=100
)

# ---> Results exploration.
tmp.file.name <- paste0(migec.path, '/histogram/overseq.txt')
file.exists(tmp.file.name)
tmp.data <- fread(file=tmp.file.name, header=TRUE)
colnames(tmp.data)[1:2] <- c('tcr.lib.name', 'sample.type')
tmp.data <- as.data.table(gather(data=tmp.data, key='read.count', value='umi.count', -`tcr.lib.name`, -`sample.type`))
to.keep <- tmp.data[, .(to.check=all(umi.count==0)), by=read.count][to.check==FALSE, read.count]
tmp.data <- tmp.data[read.count %in% to.keep]
tmp.data$read.count <- factor(x=tmp.data$read.count, levels=gtools::mixedsort(unique(tmp.data$read.count)))
# Identify bad-quality TCR samples, where one of them is defined as a TCR sample w/ less than 1000 UMI molecules regardless of read support.
tmp.thold <- 10000
tmp.data.2 <- tmp.data[, .(umi.count=sum(umi.count)), by=tcr.lib.name]
tmp.data.2[, sample.qual:=ifelse(test=umi.count<tmp.thold, yes='bad', no='rest')]
tmp.data <- merge(x=tmp.data, y=tmp.data.2[, .(tcr.lib.name, sample.qual)], by='tcr.lib.name', all=TRUE, sort=FALSE)
# tmp.data$read.count <- as.integer(tmp.data$read.count)
tmp.capt <- paste0('Threshold of total UMI count to call a sample as bad is ', tmp.thold)
tmp.ggplot.1 <- ggplot(data=tmp.data.2, aes(x=umi.count)) +
    geom_histogram(bins=30, alpha=0, color='black', linewidth=1) +
    geom_vline(xintercept=tmp.thold, linewidth=1.5, color='darkred', linetype='dashed') +
    scale_y_continuous(expand=expansion(add=0), n.breaks=6) + 
    scale_x_continuous(expand=expansion(add=10), n.breaks=6) +
    labs(title='Distribution of total UMI count per TCR sample', x='Sample total UMI count', y='Frequency', caption=tmp.capt) +
    theme_bw()
tmp.ggplot.2 <- ggplot(data=tmp.data, aes(x=read.count, y=umi.count+1)) +
    geom_boxplot(fill='lightblue', linewidth=1.5, color='black') +
    geom_jitter(size=1, color='black', width=0.6) +
    scale_y_log10(expand=expansion(add=0)) +
    labs(title='UMI count distribution when UMIs are split accordind to read support', x='Number of reads', y='Number of UMIs (log10)') +
    theme_bw()
tmp.ggplot.3 <- ggplot(data=tmp.data[sample.qual!='bad'], aes(x=read.count, y=umi.count+1)) +
    geom_boxplot(fill='lightblue', linewidth=1.5, color='black') +
    geom_jitter(size=1, color='black', width=0.6) +
    scale_y_log10(expand=expansion(add=0)) +
    labs(title='UMI count distribution when UMIs are split accordind to read support\nSamples previously called bad were excluded here', x='Number of reads', y='Number of UMIs (log10)') +
    theme_bw()
tmp.file.name <- paste0(migec.path, '/histogram/ReadSupportOfUMIs.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot.1)
print(tmp.ggplot.2)
print(tmp.ggplot.3)
dev.off()



### ---------------------- Consensus assembly ----------------------- ###

# ---> Predefinitions.
tmp.path <- paste0(jobs.path, '/assemble')
if(!dir.exists(tmp.path)) dir.create(tmp.path)

# ---> Define and run job.
run.job(
    module.name='AssembleBatch',
    module.path=tmp.path,
    line.rest=paste0(
        '--force-collision-filter --force-overseq 2 ',
        migec.path, '/checkout/', ' ',
        migec.path, '/histogram/', ' ',
        migec.path, '/assemble/'),
    memory=250
)
# Particularly when samples were shallowly sequenced, force the MIG UMI count threshold to be X with the following line:
# java -jar migec.jar AssembleBatch --force-collision-filter --force-overseq X <checkout_output_folder> <histogram_output_folder> <assemble_output_folder>


### -------------------------- VDJ mapping -------------------------- ###

# ---> Predefinitions.
tmp.path <- paste0(jobs.path, '/cdrblast')
if(!dir.exists(tmp.path)) dir.create(tmp.path)

# ---> Assemble results to filter out libraries w/out enough MIGs to continue with process.
# Define MIG count threshold.
mig.thold <- 100
# Filter out low-quality TCR libraries.
tmp.file.name <- paste0(migec.path, '/assemble/assemble.log.txt')
file.exists(tmp.file.name)
for.qc <- fread(file=tmp.file.name, header=TRUE)
for.qc <- for.qc[, .(sample.id=`#SAMPLE_ID`, mig.total=MIGS_GOOD_TOTAL)]
for.qc <- for.qc[mig.total>mig.thold]

# ---> Create sample meta file.
tmp.vals <- c(`TCRB`='TRB', `TCRA`='TRA')
tmp.data <- bc.info[,
    .(
        sample.id=paste(tcr.lib.id, tcr.lib.name, chain.type, sep='.'),
        species='HomoSapiens',
        gene=tmp.vals[chain.type],
        file.type='paired',
        mask='1:1',
        qual='20,25'
    )
]
# Filter out low-quality TCR libraries according to MIG number (see definition of "for.qc" above).
tmp.data <- tmp.data[sample.id %in% for.qc[, sample.id]]
# Sanity check: There should be unique entries in the table.
tmp.check <- tmp.data[, .N]
tmp.data <- unique(tmp.data)
tmp.check <- tmp.check == tmp.data[, .N]
if(!tmp.check) stop('Unexpected error!\n')
# Final format and save.
meta.file <- paste0(migec.path, '/SampleMeta_CdrBlast.tsv')
fwrite(file=meta.file, x=tmp.data, sep='\t', quote=FALSE, col.names=FALSE)

# ---> Define and run job.
run.job(
    module.name='CdrBlastBatch',
    module.path=tmp.path,
    line.rest=paste0(
        ' --sample-metadata ', meta.file, ' ',
        migec.path, '/checkout/', ' ',
        migec.path, '/assemble/', ' ',
        migec.path, '/cdrblast/'
    ),
    memory=100
)


### ----------------------- Results filtering ----------------------- ###

# ---> Predefinitions.
tmp.path <- paste0(jobs.path, '/cdrfinal')
if(!dir.exists(tmp.path)) dir.create(tmp.path)

# ---> Define and run job.
run.job(
    module.name='FilterCdrBlastResultsBatch',
    module.path=tmp.path,
    line.rest=paste0(
        ' -s -c ', # To filter out "singleton" MIGs.
        migec.path, '/cdrblast/', ' ',
        migec.path, '/cdrfinal/'
    ),
    memory=100
)

# ---> Results exploration.
# Read results and preprocess.
tmp.file.name <- paste0(migec.path, '/cdrfinal/cdrblastfilter.log.txt')
file.exists(tmp.file.name)
tmp.data <- fread(file=tmp.file.name, header=TRUE)
colnames(tmp.data)[1] <- 'tcr.lib.name'
colnames(tmp.data) <- str_replace_all(string=tolower(colnames(tmp.data)), pattern='_', replacement='.')
# Viz.
#       Clonotype count before and after filtering.
to.plot <- tmp.data[, .(tcr.lib.name, clonotypes.total, clonotypes.filtered)]
to.plot <- as.data.table(gather(data=to.plot, key='filt.stage', value='clonotype.count', -`tcr.lib.name`))
tmp.vals <- c('clonotypes.total'='Before', 'clonotypes.filtered'='After')
to.plot$filt.stage <- factor(x=tmp.vals[to.plot$filt.stage], levels=tmp.vals)
tmp.ggplot.1 <- ggplot(data=to.plot, aes(x=filt.stage, y=clonotype.count)) +
    geom_boxplot(fill='lightblue', linewidth=1.5, color='black') +
    geom_jitter(size=1, color='black', width=0.05) +
    geom_line(aes(group=tcr.lib.name), color='black', linewidth=0.5) +
    scale_y_log10() +
    labs(title='Number of clonotypes before and after filtering', x='Filtering stage', y='Number of clonotypes') +
    theme_bw()
#       UMI count before and after filtering.
to.plot <- tmp.data[, .(tcr.lib.name, events.total, events.filtered)]
to.plot <- as.data.table(gather(data=to.plot, key='filt.stage', value='umi.count', -`tcr.lib.name`))
tmp.vals <- c('events.total'='Before', 'events.filtered'='After')
to.plot$filt.stage <- factor(x=tmp.vals[to.plot$filt.stage], levels=tmp.vals)
tmp.ggplot.2 <- ggplot(data=to.plot, aes(x=filt.stage, y=umi.count)) +
    geom_boxplot(fill='lightblue', linewidth=1.5, color='black') +
    geom_jitter(size=1, color='black', width=0.05) +
    geom_line(aes(group=tcr.lib.name), color='black', linewidth=0.5) +
    scale_y_log10() +
    labs(title='Number of UMIs before and after filtering', x='Filtering stage', y='Number of UMIs') +
    theme_bw()
tmp.file.name <- paste0(migec.path, '/cdrfinal/ReportSummary.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot.1)
print(tmp.ggplot.2)
dev.off()