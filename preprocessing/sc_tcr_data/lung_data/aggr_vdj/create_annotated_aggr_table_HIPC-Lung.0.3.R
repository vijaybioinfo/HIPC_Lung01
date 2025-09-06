############    --------   Create vdj aggr table    -------    ############

# By: Vicente Fajardo
# Version: 0
# Subversion: 1
# Sub-subversion: 1
# Updates.
# ---> Version updates:
#   None, first version.
# ---> Subversion updates:
#   None, first subversion.
# ---> Sub-subversion updates:
#   None, first sub-subversion.

### -------------------------- Description -------------------------- ###
# Given the definition of an aggregation (run ID (run date) and sample projects that should be included) the script will attempt to create an aggregation table for the set of vdj samples within the project. In the end, it will also consider the order of the GenEx samples aggregation table for them to be very similar regarding sample suffixes.

### -------------------------- Dependencies ------------------------- ###
library(stringr)
library(data.table)
library(lubridate)

### --------------------------- Arguments --------------------------- ###
# # ---> Commands from command line.
arguments = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if(length(arguments)!=1) {
  stop("One (and only one) argument must be supplied.", call.=FALSE)
}

# ---> Sequencing run definition.
input.run.dates <- c(
  '05-25-2023', '05-26-2023',
  '06-27-2023', '06-28-2023',
  '10-20-2023', '10-28-2023',
  '11-10-2023', '11-17-2023'
)
output.run.date <- '11-17-2023'
aggr.id <- arguments[1]
aggr.table.version <- '0.1'
# Does the aggregation table contain annotations?
is.annotated <- TRUE
if(is.annotated) aggr.table.basename <- paste0(aggr.id, '_aggr_table_annotated.', aggr.table.version, '.csv') else aggr.table.basename <- paste0(aggr.id, '_aggr_table.', aggr.table.version, '.csv')
output.table.version <- '.1.0'
# ---> Constants.
vdj.file.sffx <- '_prj_vdj_sample_ids.txt'
# ---> Path definitions.
# Regarding input.
input.gen.paths <- paste0('/path/to/user/sequencing_data/', input.run.dates)
vdj.paths <- paste0(input.gen.paths, '/vdj')
if(!all(dir.exists(vdj.paths))) cat('WARNING: No vdj info or no info at all for some IDs provided.\n')
# Regarding output.
out.gen.path <- paste0('/path/to/user/sequencing_data/', output.run.date)
aggr.path <- paste0(out.gen.path, '/aggr')
aggr.data.path <- paste0(aggr.path, '/data')
aggr.table.file <- paste(aggr.data.path, aggr.table.basename, sep='/')
if(!file.exists(aggr.table.file)) stop('Could not define input aggr table properly.\n')
output.path <- paste0(out.gen.path, '/aggr_vdj/', aggr.id, '/data')
if(!dir.exists(output.path)) dir.create(output.path, recursive=TRUE)
#count.path <- paste0(run.gen.paths, '/count')
#layout.path <- paste0(count.path, '/data/experiment_layout')

### --------------------------- Functions --------------------------- ###

### --------------------------- Load data --------------------------- ###
# ---> Aggr table.
aggr.table <- read.csv(file=aggr.table.file, stringsAsFactors=FALSE)

### ---------------------- Data preprocessing ----------------------- ###
aggr.table$library_id <- str_replace(string=aggr.table$library_id, pattern='_Gex$|_GenEx|_GEX$', replacement='')

### ------------------------- Main program -------------------------- ###

# ---> Samples projects and IDs within the important runs.
# Get sample IDs of interest.
sample.ids <- aggr.table$library_id
# Get rid of type of library suffix if necessary.
sample.ids <- str_replace(string=sample.ids, pattern='_Gex|_GenEx|_GEX', replacement='')

# Get paths to all sample projects where our data could be located.
sample.prjs.info <- list.dirs(path=vdj.paths, full.names=TRUE, recursive=FALSE)
sample.ids.info <- list.dirs(path=sample.prjs.info, full.names=TRUE, recursive=FALSE)
sample.ids.info <- data.frame(sample.id=basename(path=sample.ids.info), id.path=sample.ids.info, stringsAsFactors=FALSE)

# Preserve only bonafide TCR sample IDs. This may make some troubles, though, be cautios.
sample.ids.info <- sample.ids.info[grepl(x=sample.ids.info$sample.id, pattern='_TCR'), ]

# Get sample ID basename.
sample.ids.info$sample.id.name <- str_replace(string=sample.ids.info$sample.id, pattern='_TCR$', replacement='')

# Turn sample IDs info dataframe into a datatable
sample.ids.info <- as.data.table(sample.ids.info)

# Define sequencing date. NOTE: When a sample ID is defined for more than one sequencing run, we assume that it was re-sequenced and we will take the data from the path with the latest date.
sample.ids.info[, date:=mdy(str_extract(string=id.path, pattern='\\d{2}-\\d{2}-\\d{2}'))]


# ---> Create vdj aggr table.
ids.table <- lapply(X=sample.ids, FUN=function(tmp.library.id){
  # Find latest definition for sample ID and define files.
  tmp.data <- sample.ids.info[sample.id.name==tmp.library.id][date==max(date)]
  sample.id.path <- tmp.data[, id.path]
  sample.id.path <- paste0(sample.id.path, '/outs')
  clons.file <- paste(sample.id.path, 'clonotypes.csv', sep='/')
  contigs.file <- paste(sample.id.path, 'filtered_contig_annotations.csv', sep='/')
  to.output <- data.frame(library_id=tmp.library.id, clonotypes=clons.file, annotations=contigs.file)
  return(to.output)
})
ids.table <- Reduce(x=ids.table, f=rbind)
ids.table$library_id <- as.character(ids.table$library_id)

# ---> Sort and filter aggr vdj table according to GenEx aggr table.
# Merge tables.
ids.table <- merge(x=aggr.table, y=ids.table, by='library_id', all.x=TRUE, all.y=FALSE, sort=FALSE)
# Keep only the columns of interest.
ids.table <- ids.table[, c('library_id', 'clonotypes', 'annotations')]
# Determine if there were any undefined libraries
ids.table$sample.sffx <- 1:nrow(ids.table)
tmp.check <-  nrow(ids.table[(is.na(ids.table$annotations) | is.na(ids.table$clonotypes)), ])
if(tmp.check>0) stop('There were some samples with undefined TCR files.\n')
# Finally, check integrity of all defined files.
tmp.check <- all(file.exists(as.character(ids.table[, 'clonotypes']))) & all(file.exists(as.character(ids.table[, 'annotations'])))
if(!tmp.check) stop('There were some samples with undefined TCR files.\n')

# ---> Output file.
output.file <- paste0(output.path, '/', aggr.id, '_vdj_aggr_table', output.table.version, '.csv')
# file.exists(output.file)
if(!file.exists(output.file)) write.csv(file=output.file, x=ids.table, row.names=FALSE, quote=FALSE) else stop('It was attempted to overwite an existing file.\n')
