############    --------------   aggr vdj    --------------    ############
cat('############    --------------   aggr vdj    --------------    ############\n')

# By: Vicente Fajardo
# Version: 2
# Subversion: 2
# Updates.
# ---> Version updates:
#   Version 2 should always be used as the stable one compared the lower ones. A horrible bug (DON'T YOU FUCK*NG EVER USE FACTORS IN apply's, @VICENTE) ruined the relationships between individual samples' IDs and aggregated clonotypes IDs.
# ---> Subversion updates:
#   Consider the aggregation table may provide the suffixes each sample ID cells must contain. If none, it may just assign them according to their row number order.
#   Get a copy of the aggregation table within the output directory. Small change, yet essential.

### -------------------------- Description -------------------------- ###
# Given vdj output data for different samples, this program will aggregate their clonotypes.csv and [filtered_]contig_annotations.csv files, producing a new table to allow the connection between barcodes and new clonotype IDs.

### -------------------------- Dependencies ------------------------- ###
cat('### -------------------------- Dependencies ------------------------- ###\n')
library(optparse)
library(ggplot2)
library(stringr)
source('/root/BioHome/vfajardo/scripts/functions/R_handy_functions.R')
cat('Dependencies loaded!\n')
cat('\n\n')

### --------------------------- Arguments --------------------------- ###
cat('### --------------------------- Arguments --------------------------- ###\n')
# Declaring arguments to parse from command line ----------------------->
option.list <- list(
  make_option(opt_str="--ReportsPath", type="character", default=NULL, dest="reports.path", help="Absolute path to reports directory."),
  make_option(opt_str="--GenInputPath", type="character", default=NULL, dest="gen.input.path", help="Absolute path to directory saving the output for all vdj samples."),
  make_option(opt_str="--AggrTable", type="character", default=NULL, dest="aggr.table.file", help="Absolute path to file describing the aggregation order."),
  make_option(opt_str="--FreqThold", type="integer", default=2, dest="freq.thold", help="For a chain-related multiplet clonotype (see code), frequency threshold for it to be called well-supported."),
  make_option(opt_str="--SampleCountThold", type="integer", default=2, dest="samples.count.thold", help="For a chain-related multiplet clonotype (see code), minimum amount of samples it has to be seen in for it to be called well-supported.")
)
# Getting arguments from command line and setting their values to their respective variable names.
opt.parser = OptionParser(option_list=option.list);
opt = parse_args(opt.parser);
# Moving options to their own variables
reports.path <- opt$reports.path
gen.input.path <- opt$gen.input.path
aggr.table.file <- opt$aggr.table.file
freq.thold <- opt$freq.thold
samples.count.thold <- opt$samples.count.thold

# Show parameter values.
cat('General reports path:', reports.path, '\n')
cat('General input path:', gen.input.path, '\n')
cat('Aboslute path to aggregation table file:', aggr.table.file, '\n')
cat('\n\n')

# Defaults ------------------------------------------------------------->
# Others --------------------------------------------------------------->
### --------------------------- Functions --------------------------- ###
# ---> 1
# Name: Merge CDR3s.
# Description:
# Given a vector of TRA and TRB CDR3 sequences, this function sorts them separately and merges them returning them alltogether.
# Arguments ------------------------->
# cdr3s - Character vector of TRA and TRB CDR3 sequences, either aa or nt ones.
# Function:
merge.cdr3s <- function(cdr3s){
  TRA.matches <- grep(x=cdr3s, pattern='TRA:', value=TRUE)
  TRB.matches <- grep(x=cdr3s, pattern='TRB:', value=TRUE)
  return(paste(c(sort(TRA.matches), sort(TRB.matches)), collapse=';'))
}

### --------------------------- Load data --------------------------- ###
cat('### --------------------------- Load data --------------------------- ###\n')
# ---> Load aggr table.
aggr.df <- read.csv(file=aggr.table.file, stringsAsFactor=FALSE)
# Check if we've got a wide aggr table.
if('clonotypes' %in% colnames(aggr.df) & 'annotations' %in% colnames(aggr.df)) total.aggr.df <- TRUE else total.aggr.df <- FALSE
# Next may be redundant, essential though.
if('sample.sffx' %in% colnames(aggr.df)) row.names(aggr.df) <- aggr.df$sample.sffx else row.names(aggr.df) <- as.character(1:nrow(aggr.df))
cat('Aggregation table loaded...\n')

# ---> Load clonotypes data.
# If files are defined in the aggr table, we take them as input and sample files. Else, we take the ones defined in the general input path.
if(total.aggr.df){
  vdj.samples <- aggr.df$library_id
}else{
  # We gotta make sure that we have data files for all samples described in the table.
  vdj.samples <- list.files(gen.input.path)
  if(!all(vdj.samples %in% aggr.df$library_id)) stop('Not all files listed in the general input directory are described in the aggregation table.')
}

# ------ Clonotypes data.
# Loading according to our aggregation data frame.
clons.list <- lapply(X=vdj.samples, FUN=function(sample){
  if(!total.aggr.df) tmp.file.name <- paste(gen.input.path, sample, 'outs/clonotypes.csv', sep='/') else tmp.file.name <- aggr.df$clonotypes[aggr.df$library_id==sample]
  return(read.csv(file=tmp.file.name, stringsAsFactor=FALSE))
})
names(clons.list) <- vdj.samples
cat('Clonptypes data loaded...\n')

# ------ Cells-clonotypes relationships data.
cells.clons.list <- lapply(X=vdj.samples, FUN=function(sample){
  if(!total.aggr.df) tmp.file.name <- paste(gen.input.path, sample, 'outs/filtered_contig_annotations.csv', sep='/') else tmp.file.name <- aggr.df$annotations[aggr.df$library_id==sample]
  return(read.csv(file=tmp.file.name, stringsAsFactor=FALSE))
})
names(cells.clons.list) <- vdj.samples
cat('Cells-clonotypes relations data loaded...\n')
cat('\n\n')

### ------------------------- Main program -------------------------- ###
cat('### ------------------------- Main program -------------------------- ###\n')
### ------------------------ Pre-aggregation ------------------------ ###
cat('### ------------------------ Pre-aggregation ------------------------ ###\n')
# ---> Add/change barcode and clonotype ID suffixes for the cells-clonotypes relationships table.
for(sample in vdj.samples){
  # Get suffix according to the order in the aggregation table.
  suffix <- row.names(aggr.df)[aggr.df$library_id==sample]
  # Add/Change barcode suffix.
  if(all(grepl(pattern='[ACTG]+-[1-9]+$', x=cells.clons.list[[sample]][, 'barcode'], perl=TRUE))){
    cells.clons.list[[sample]][, 'barcode'] <- gsub(pattern='[1-9]+$', x=cells.clons.list[[sample]][, 'barcode'], perl=TRUE, replacement=suffix)
  }else{
    cells.clons.list[[sample]][, 'barcode'] <- paste0(cells.clons.list[[sample]][, 'barcode'], '-', suffix)
  }
  # Add clonotype ID suffix in both tables.
  cells.clons.list[[sample]][, 'raw_clonotype_id'] <- paste0(cells.clons.list[[sample]][, 'raw_clonotype_id'], '-', suffix)
  cells.clons.list[[sample]][, 'raw_consensus_id'] <- paste0(cells.clons.list[[sample]][, 'raw_consensus_id'], '-', suffix)
  clons.list[[sample]][, 'clonotype_id'] <- paste0(clons.list[[sample]][, 'clonotype_id'], '-', suffix)
}

# ---> Sort TRA and TRB sequences so that chain-multiple clonotypes can be considered in the aggregation process.
for(sample in vdj.samples){
  # ---> Merge nucleotide CDR3 sequences.
  nt.cdr3s <- str_split(string=clons.list[[sample]]$cdr3s_nt, pattern=';')
  nt.cdr3s <- sapply(X=nt.cdr3s, FUN=merge.cdr3s)
  # ---> Merge aa CDR3 sequences.
  aa.cdr3s <- str_split(string=clons.list[[sample]]$cdr3s_aa, pattern=';')
  aa.cdr3s <- sapply(X=aa.cdr3s, FUN=merge.cdr3s)
  # ---> Modify the values in the original tables.
  clons.list[[sample]]$cdr3s_nt <- nt.cdr3s
  clons.list[[sample]]$cdr3s_aa <- aa.cdr3s
}
cat('Data modified for aggregation!\n')
cat('\n\n')

### -------------------------- Aggregation -------------------------- ###
cat('### -------------------------- Aggregation -------------------------- ###\n')
# Three different approaches to be considered.
# ---> Get the unique nt CDR3 sequences.
nt.cdr3s <- unique(unlist(lapply(X=vdj.samples, FUN=function(sample) return(clons.list[[sample]]$cdr3s_nt))))
names(nt.cdr3s) <- paste0('clonotype', 1:length(nt.cdr3s))

# ---> Raw aggregated clonotypes table.
raw.clons.aggr <- Reduce(x=clons.list, f=rbind)

# ---> Get aggregated clonotypes table.
clean.clons.aggr <- t(sapply(X=names(nt.cdr3s), FUN=function(clon.id){
  # CDR3 nt sequence as such.
  cdr3 <- nt.cdr3s[clon.id]
  # All info subset for this clonotype.
  to.keep <- raw.clons.aggr$cdr3s_nt==cdr3
  this.clon <- raw.clons.aggr[to.keep, ]
  # Get the aa sequence.
  cdr3.aa <- unique(this.clon$cdr3s_aa)
  # Get the clonotype frequency.
  clon.freq <- sum(this.clon$frequency)
  # Get to know from how many samples we're recovering this CDR3.
  # Get samples.
  orig.samples <- unique(str_match(string=this.clon$clonotype_id, pattern='clonotype[0-9]+-([0-9]+)')[, 2]) # Check str_match at https://stringr.tidyverse.org
  # Get the total count.
  samples.count <- length(orig.samples)
  # Collapse samples.
  orig.samples <- paste0(orig.samples, collapse=';')
  # Output.
  to.output <- c(cdr3, clon.id, clon.freq, cdr3.aa, orig.samples, samples.count)
  return(to.output)
}))
# Turn the object into a data frame and calculate propotions.
# This is done here rather than inside the apply in order to add the proportion column and keep the original order.
clean.clons.aggr <- data.frame(clonotype_id=clean.clons.aggr[,2], frequency=as.integer(clean.clons.aggr[,3]), proportion=as.integer(clean.clons.aggr[,3])/sum(as.integer(clean.clons.aggr[,3])), cdr3s_aa=clean.clons.aggr[,4], cdr3s_nt=clean.clons.aggr[,1], samples=clean.clons.aggr[,5], samples_count=as.integer(clean.clons.aggr[,6]), stringsAsFactors=FALSE)

# ---> Filter non-well-supported CMCs.
# CMC stand for chain-multiplet clonotype.
# Filter CMCs according to input criteria, minimum frequency (freq.thold) and minimum amount of original samples (samples.count.thold).
# Idenitify CMCs.
tra.counts <- str_count(string=clean.clons.aggr$cdr3s_nt, pattern='TRA:')
trb.counts <- str_count(string=clean.clons.aggr$cdr3s_nt, pattern='TRB:')
is.cmc <- tra.counts>1 | trb.counts>1
# Identify which CMCs should be filtered out.
to.delete <- is.cmc & !(clean.clons.aggr$frequency >= freq.thold | clean.clons.aggr$samples_count >= samples.count.thold)
# Before deleting them, report a CMC table along a column indicating if it passed or not.
cmcs.table <- clean.clons.aggr
cmcs.table$pass <- !to.delete
cmcs.table <- cmcs.table[is.cmc, ]
# Then, delete non-well-supported CMCs.
clean.clons.aggr <- clean.clons.aggr[!to.delete,]
# Report how many CMCs were deleted.
cmcs.number <- sum(is.cmc)
cmcs.deleted <- sum(to.delete)
cat('\n\n# ------> CMCs report.\n\tCMCs detected:', cmcs.number, '\n\tNon-well supported CMCs:', cmcs.deleted, '\n\n')

# ---> Get a table for last id-new id relationships.
# Here, we iterate over the clean aggregation table since we've filtered the non-supported CMCs and we don't want to include them in downstream analyses anymore.
clons.ids.relationships <- lapply(X=clean.clons.aggr$clonotype_id, FUN=function(clon.id){
  # CDR3 nt sequence as such.
  cdr3 <- nt.cdr3s[clon.id]
  # All info subset for this clonotype.
  to.keep <- raw.clons.aggr$cdr3s_nt==cdr3
  this.clon <- raw.clons.aggr[to.keep, ]
  # Identify the previous clonotypes associated to this unique clonoptype.
  previous.ids <- data.frame(new.id=rep(clon.id, times=length(this.clon$clonotype_id)), last.id=this.clon$clonotype_id, stringsAsFactors=FALSE)
  return(previous.ids)
})
clons.ids.relationships <- Reduce(x=clons.ids.relationships, f=rbind)

# ---> Get the aggregation for the clonotypes-cells relationships table.
raw.cells.clons.aggr <- Reduce(x=cells.clons.list, f=rbind)

# Clean the results.
# First, we get rid of the relationships that involve a clonotype ID that, somehow (mainly because they have None values), was eliminated from the clonotypes part.
to.keep <- raw.cells.clons.aggr$raw_clonotype_id %in% clons.ids.relationships$last.id
raw.cells.clons.aggr <- raw.cells.clons.aggr[to.keep, ]

# Then, exchange the clonotype ID for the new one.
clean.cells.clons.aggr <- as.data.frame(t(apply(X=raw.cells.clons.aggr, MARGIN=1, FUN=function(cell.clon.rel){
  last.clon.id <- cell.clon.rel['raw_clonotype_id']
  new.clon.id <- as.character(unique(clons.ids.relationships$new.id[clons.ids.relationships$last.id==last.clon.id]))
  cell.clon.rel['raw_clonotype_id'] <- new.clon.id
  return(cell.clon.rel)
})), stringsAsFactors=FALSE)

### ----------------------------- Output ---------------------------- ###
# ---> Cells-clonotypes relationships.
tmp.file.name <- paste0(reports.path, '/filtered_contig_annotations_aggr.csv')
write.csv(x=clean.cells.clons.aggr, file=tmp.file.name, quote=FALSE, row.names=FALSE)

# ---> Clonotypes.
tmp.file.name <- paste0(reports.path, '/clonotypes_aggr.csv')
write.csv(x=clean.clons.aggr, file=tmp.file.name, quote=FALSE, row.names=FALSE)
cat('Aggregation finished and both aggregated files output to general reports path!\n\n')

# ---> CMCs table.
tmp.file.name <- paste0(reports.path, '/chain_related_multiplet_clonotypes.csv')
write.csv(x=cmcs.table, file=tmp.file.name, quote=FALSE, row.names=FALSE)

# ---> Aggr table.
tmp.file.name <- paste0(reports.path, '/vdj_aggr_table.csv')
write.csv(x=aggr.df, file=tmp.file.name, quote=FALSE, row.names=FALSE)

cat('PROGRAM FINISHED!\n')
