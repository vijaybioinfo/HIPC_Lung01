############    --------   Create aggr tables    --------    ############

# By: Vicente Fajardo

### -------------------------- Dependencies ------------------------- ###
library(stringr)
library(data.table)
library(lubridate)


### --------------------------- Arguments --------------------------- ###
# ---> General definitions.
# Output version.
output.version <- 0.3
this.user <- system(command='echo $USER', intern=TRUE)
# Indicates whether the output GEx files are to be taken from scratch.
from.scratch <- FALSE
# Project definition.
gen.prj <- 'R24'
prj.name <- 'DICE-VirusSpec'
sample.prj.pttn <- '^DICE_TCR'
batches.lab <- '1-to-15'
aggr.prj <- paste0(prj.name, '_Batches-', batches.lab)
seq.dates <- c(
  '04-19-2023', '04-20-2023', '04-21-2023', '04-22-2023',
  '04-23-2023', '04-24-2023', '04-25-2023', '05-05-2023',
  '05-25-2023', '05-26-2023',
  '06-27-2023', '06-28-2023',
  '07-03-2023', '07-24-2023',
  '08-21-2023'
)
names(seq.dates) <- as.character(1:length(seq.dates))
sample.prj.pttn <- 'DICE_TCR'
hto.dates <- c('05-02-2023', '05-12-2023', '05-29-2023', '06-28-2023', '07-05-2023', '07-26-2023') # There might be one or more potential dates to looks patterns from.
# Samples with no HTO data.
ids.no.hto <- c('DICE_TCR018_Hu_CMV1_17D_CITE')
prjs.no.hto <- NA
# ---> Parameters to split main aggr table.
cell.types <- list(CD4='CD4', CD8='CD8', ALL=c('CD4', 'CD8'))
# ---> Path definitions.
gen.seq.path <- '/path/to/user/sequencing_data'
gen.hto.path <- paste0('/path/to/user/', gen.prj, '/deconvolution/HTO_based')
if(!dir.exists(gen.hto.path)) stop('Unexpected absence of HTO path.\n')
hto.sffx <- '_100umis_thold' # Leave blank if necessary.
# For output
output.date <- tail(seq.dates, 1)
output.path <- paste0(gen.seq.path, '/', output.date, '/aggr_vdj/data/')
if(!dir.exists(output.path)) dir.create(output.path)


### ------------------------- Main program -------------------------- ###

# ---> General annotated aggr table.

# @ General info per vdj library.
vdj.paths <- paste0(gen.seq.path, '/', seq.dates, '/vdj')
if(!all(dir.exists(vdj.paths))) stop('Unexpected absence of vdj paths.\n')
vdj.paths <- list.dirs(path=vdj.paths, full.names=TRUE, recursive=FALSE)
tmp.data <- data.table(
  dir.path=list.dirs(path=vdj.paths, full.names=TRUE, recursive=FALSE)
)
tmp.data[, sample.id:=basename(dir.path)]
# Keep only project-relevant samples.
tmp.data <- tmp.data[str_detect(string=sample.id, pattern=sample.prj.pttn)]
tmp.data[, sample.prj.path:=str_replace(string=dir.path, pattern=sample.id, replacement='')]
tmp.data[, chrom.batch.tag:=basename(sample.prj.path)]
tmp.data[,
  seq.batch.date:=basename(
    str_replace(
      string=sample.prj.path,
      pattern=paste0('vdj/', chrom.batch.tag, '/'),
      replacement=''
    )
  )
]
tmp.vals <- names(seq.dates); names(tmp.vals) <- seq.dates
tmp.data[, seq.batch.tag:=tmp.vals[seq.batch.date]]

# @ Define paths.
# PB file to aggr vdj.
tmp.data[, vdj_contig_info:=paste0(dir.path, '/outs/vdj_contig_info.pb')]
tmp.check <- tmp.data[, all(file.exists(vdj_contig_info))]
if(!tmp.check) stop('Not all PB files were properly defined.\n')
# Clonotype files.
tmp.data[, clonotypes:=paste0(dir.path, '/outs/clonotypes.csv')]
tmp.check <- tmp.data[, all(file.exists(clonotypes))]
if(!tmp.check) stop('Not all clonotype files were properly defined.\n')
# Annotation files.
tmp.data[, annotations:=paste0(dir.path, '/outs/filtered_contig_annotations.csv')]
tmp.check <- tmp.data[, all(file.exists(annotations))]
if(!tmp.check) stop('Not all clonotype annotation files were properly defined.\n')
# Keep only latest results per sample ID.
tmp.data[,
  tmp.batch.date:=mdy(seq.batch.date)
]
tmp.data.1 <- tmp.data[,
  .SD[
    tmp.batch.date==max(tmp.batch.date),
    .(
      dir.path, chrom.batch.tag, seq.batch.date, vdj_contig_info, clonotypes, annotations
    )
  ],
  by=sample.id
]
tmp.data.2 <- tmp.data[,
  .(seq.batch.tag=paste0(seq.batch.tag, collapse=';')),
  by=sample.id
]
tmp.data <- merge(x=tmp.data.2, y=tmp.data.1, by='sample.id')

# @ Donor deconvolution paths.
hto.paths <- paste(gen.hto.path, seq.dates, sep='/')
if(!all(dir.exists(hto.paths))) warning('There is no defined HTO-based deconvolution directory for all provided paths.\n')
hto.paths <- hto.paths[dir.exists(hto.paths)]
hto.paths <- list.dirs(path=hto.paths, full.names=TRUE, recursive=FALSE)
hto.paths <- list.dirs(path=hto.paths, full.names=TRUE, recursive=FALSE)
hto.paths <- list.dirs(path=hto.paths, full.names=TRUE, recursive=FALSE)
hto.paths <- hto.paths[str_detect(string=hto.paths, pattern=hto.sffx)]
hto.paths <- data.table(
  base.path=hto.paths,
  deconv.date=mdy(str_replace(
    string=str_extract(
      string=hto.paths, pattern='\\d{2}-\\d{2}-\\d{4}[^-]+$'
    ),
    pattern=hto.sffx,
    replacement=''
  ))
)
hto.paths[,
  sample.id:=basename(str_replace(
    string=base.path,
    pattern=basename(base.path),
    replacement=''
  ))
]
hto.paths[, hto.tag:=paste0(base.path, '/general_reports/ClassificationPerCell.csv')]
tmp.check <- hto.paths[, all(file.exists(hto.tag))]
if(!tmp.check) warning('Process might have not been completed for all samples. Make sure this is expected before continuing.\n')
hto.paths <- hto.paths[file.exists(hto.tag)]
# Identify latest seq. run for a given sample ID.
hto.paths[,
  seq.date:=mdy(
    str_replace(
      string=str_extract(string=hto.tag, pattern='HTO_based/\\d{2}-\\d{2}-\\d{4}'),
      pattern='HTO_based/', replacement=''
    )
  )
]
hto.paths <- hto.paths[,
  .SD[seq.date==max(seq.date),],
  by=sample.id
]
# Identify latest deconvolution run for a given sample ID.
hto.paths <- hto.paths[,
  .SD[deconv.date==max(deconv.date),],
  by=sample.id
]
# Merge accordingly.
tmp.data[,
  hto.sample.id:=str_replace(string=sample.id, pattern='_TCR$', replacement='_CITE')
] 
# tmp.data[!sample.id %in% ids.no.hto & !(hto.sample.id %in% hto.paths[, sample.id]), ]
tmp.check <- all(tmp.data[!(hto.sample.id %in% ids.no.hto) & !(chrom.batch.tag %in% prjs.no.hto), hto.sample.id] %in% hto.paths[, sample.id])
if(!tmp.check) stop('HTO data not defined for some samples.\n')
tmp.data <- merge(x=tmp.data, y=hto.paths[, .(sample.id, hto.tag)], all.x=TRUE, all.y=FALSE, by.x='hto.sample.id', by.y='sample.id')


# @ CD3 surface expression info.
count.paths <- paste0(gen.seq.path, '/', seq.dates, '/count')
if(!all(dir.exists(count.paths))) stop('Unexpected absence of count paths.\n')
count.paths <- list.dirs(path=count.paths, full.names=TRUE, recursive=FALSE)
# Only output directories labeled w/ suffix "_v2" contain CD3 surface expression information.
count.paths <- count.paths[str_detect(string=count.paths, pattern='_v2$')]
tmp.data.2 <- data.table(
  dir.path=list.dirs(path=count.paths, full.names=TRUE, recursive=FALSE)
)
tmp.data.2[, sample.id:=basename(dir.path)]
# Keep only project-relevant samples.
tmp.data.2 <- tmp.data.2[str_detect(string=sample.id, pattern=sample.prj.pttn)]
tmp.data.2[, sample.prj.path:=str_replace(string=dir.path, pattern=sample.id, replacement='')]
tmp.data.2[, chrom.batch.tag:=str_replace(string=basename(sample.prj.path), pattern='_v2$', replacement='')]
tmp.data.2[,
  seq.batch.date:=basename(
    str_replace(
      string=sample.prj.path,
      pattern=paste0('count/', chrom.batch.tag, '_v2/'),
      replacement=''
    )
  )
]
tmp.vals <- names(seq.dates); names(tmp.vals) <- seq.dates
tmp.data.2[, seq.batch.tag:=tmp.vals[seq.batch.date]]
# CITE-seq out directory.
tmp.data.2[, cite.path:=paste0(dir.path, '/outs')]
tmp.check <- tmp.data.2[, all(dir.exists(cite.path))]
if(!tmp.check) stop('Not all CITE-seq out directories were properly defined.\n')
# Keep only latest results per sample ID.
tmp.data.2[,
  tmp.batch.date:=mdy(seq.batch.date)
]
tmp.data.1 <- tmp.data.2[,
  .SD[
    tmp.batch.date==max(tmp.batch.date),
    .(
      dir.path, chrom.batch.tag, seq.batch.date, cite.path
    )
  ],
  by=sample.id
]
tmp.data.2 <- tmp.data.2[,
  .(seq.batch.tag=paste0(seq.batch.tag, collapse=';')),
  by=sample.id
]
tmp.data.2 <- merge(x=tmp.data.2, y=tmp.data.1, by='sample.id')
# Merge accordingly.
tmp.data.2[,
  tcr.sample.id:=str_replace(string=sample.id, pattern='_CITE$', replacement='_TCR')
] 
# tmp.data[!sample.id %in% ids.no.hto & !(hto.sample.id %in% hto.paths[, sample.id]), ]
tmp.check <- all(tmp.data[!(sample.id %in% ids.no.hto) & !(chrom.batch.tag %in% prjs.no.hto), hto.sample.id] %in% tmp.data.2[, sample.id])
if(!tmp.check) stop('HTO data not defined for some samples.\n')
tmp.data <- merge(x=tmp.data, y=tmp.data.2[, .(sample.id, cite.tag=cite.path)], all.x=TRUE, all.y=FALSE, by.x='hto.sample.id', by.y='sample.id')


# ---> Add lane tags.
lane.tags <- data.table(
  sample.id=tmp.data[, sample.id],
  cell.type=str_replace_all(
    string=str_extract(string=tmp.data[, sample.id], pattern='_CD4_|_CD8_'),
    pattern='_', replacement=''
  ),
  peptide.pool=str_extract(
    string=str_replace(
      string=tmp.data[, sample.id], pattern='DICE_TCR\\d+_Hu_', replacement=''
    ),
    pattern='[^_]+'
  )
)
# @ Fix cell type.
#       Minor fix
lane.tags[is.na(cell.type), cell.type:='CD4']
#       Labeling error during sequencing: All libraries from batch TCR008 are CD4 T cells.
tmp.samples <- c('DICE_TCR008_Hu_EBV_CD8_16D_TCR', 'DICE_TCR008_Hu_SC2_CD8_16D_TCR')
lane.tags[sample.id %in% tmp.samples, cell.type:='CD4']
# Fix peptide pool for pertussis pools.
lane.tags[sample.id=='DICE_TCR015_Hu_1_17D_TCR' & peptide.pool=='1', peptide.pool:='RSV']
lane.tags[sample.id=='DICE_TCR015_Hu_2_17D_TCR' & peptide.pool=='2', peptide.pool:='HIV']
# Fix peptide pool for pertussis pools.
lane.tags[peptide.pool=='1', peptide.pool:='R']
lane.tags[peptide.pool=='2', peptide.pool:='VAX']
# Fix peptide pool.
tmp.vals <- c(
  `CMV`='CMV',
  `FLU`='IAV',
  `EBV`='EBV',
  `SC2`='SARS-CoV-2',
  `R`='B-pertussis-Rest',
  `VAX`='B-pertussis-Vax',
  `RSV`='RSV',
  `HIV`='HIV',
  `CMV1`='CMV',
  `CMV2`='CMV',
  `hMPV`='MPV',
  `RSPV`='PIV',
  `MPV`='MPV',
  `ALT`='Alternaria',
  `ASP`='Aspergillus'
)
tmp.check <- lane.tags[, all(peptide.pool%in%names(tmp.vals))]
if(!tmp.check) stop('Not all peptide pool values were properly defined.\n')
lane.tags[, peptide.pool:=tmp.vals[peptide.pool]]
lane.tags[, sample.id:=NULL]
tmp.data$lane.tag <- apply(X=lane.tags, MARGIN=1, FUN=paste, collapse=';')

# ---> General subset.
# Remove HIV-specific  and alternaria-specific cell libraries.
tmp.vals <- c('HIV', 'Alternaria')
tmp.data <- tmp.data[lane.tags[, !peptide.pool%in%tmp.vals], ]
lane.tags <- lane.tags[!peptide.pool %in% tmp.vals]

# ---> Output main aggr table.
# Format to output.
tmp.data[, sample_id:=sample.id]
cols.order <- c(
  'sample_id',
  'vdj_contig_info',
  'clonotypes',
  'annotations',
  'hto.tag',
  'cite.tag',
  'chrom.batch.tag',
  'seq.batch.tag',
  'lane.tag'
)
tmp.data <- tmp.data[, ..cols.order]
tmp.file.name <- paste0(output.path, '/', aggr.prj, '_aggr_table_annotated.', output.version,'.csv')
if(!file.exists(tmp.file.name)) fwrite(file=tmp.file.name, quote=TRUE, x=tmp.data)

# ---> Split main aggr table/
# Define aggr IDs file.
aggr.ids.file <- paste0(output.path, '/aggr_ids.csv')
write(file=aggr.ids.file, x='id')
# Cell type.
for(cell.type in names(cell.types)){
  # Subset.
  tmp.table <- tmp.data[lane.tags$cell.type %in% cell.types[[cell.type]]]
  # Check if there are any entries for a given set.
  if(dim(tmp.table)[1]==0) next
  # Define ID and save.
  aggr.id <- paste0(aggr.prj, '_', cell.type)
  write(x=aggr.id, file=aggr.ids.file, append=TRUE)
  # Aggr table with annotations.
  tmp.file.name <- paste0(output.path, '/', aggr.id, '_aggr_table_annotated.', output.version, '.csv')
  fwrite(x=tmp.table, file=tmp.file.name)
  # Aggr table with no annotations.
  tmp.file.name <- str_replace(string=tmp.file.name, pattern='_annotated', replacement='')
  fwrite(x=tmp.table[, .(library_id=sample_id, clonotypes, annotations)], file=tmp.file.name)
}