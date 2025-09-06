cat('\n\n')
############    -----   Compile info from external    -----    ############
############    -------   TCR/epitope databases    --------    ############
cat('############    -----   Compile info from external    -----    ############\n')
cat('############    -------   TCR/epitope databases    --------    ############\n')

# By: Vicente Fajardo
# Best module: R/3.6.1

### -------------------------- Description -------------------------- ###
# This script was created to compile the information obtained in multiple files from external databases on TCR-epitope information. So far, we've included info from IEDB and VDJdb.
# For this process, we take unique clonotype/epitope pairs as unique entries for our final reference. Specifically, we aim to retrieve the following entries for each clonotype/epitope pair:
#       Beta and alpha CDR3 sequences. A beta sequence is sufficient to define a unique clonotype, but database entries defined only by an alpha chain without a beta chain will be discarded.
#       V and J genes for each chain. As above, beta alone is sufficient, but this does not apply to alpha.
#       Epitope sequence recognized by the clonotype.
#       Name of the epitope sequence recognized by the clonotype.
#       Antigen and organism that the epitope comes from.
#       Assay used to determine specificity.
#       HLA context if provided.


cat('\n\n')
### -------------------------- Dependencies ------------------------- ###
cat('### -------------------------- Dependencies ------------------------- ###\n')
library(data.table)
library(tidyr)
library(stringr)
library(ggplot2)
library(VennDiagram)
library(UpSetR)


cat('\n\n')
### --------------------------- Arguments --------------------------- ###
cat('### --------------------------- Arguments --------------------------- ###\n')

# @ Paper.
main.prj <- 'HIPC'
this.prj <- 'HIPC-Lung'
this.prj.lab <- 'ag-spc_tcr_reference'
reports.date <- NULL
# @ Seed.
set.seed(seed=1)
# ---> Path definitions.
prj.gen.path <- paste0('/path/to/user/', main.prj, '/paper_developments/', this.prj)
gen.reports.path <- paste0(prj.gen.path, '/', this.prj.lab)
reports.date <- if(!is.null(reports.date)) reports.date else Sys.Date()
reports.path <- paste0(gen.reports.path, '/compiled_reference_', reports.date)
if(!dir.exists(reports.path)) dir.create(reports.path)
ext.data.path <- paste0(gen.reports.path, '/external_info')
iedb.data.path <- paste0(ext.data.path, '/iedb')
vdjdb.data.path <- paste0(ext.data.path, '/vdjdb')
mcpas.data.path <- paste0(ext.data.path, '/mcpas')
imgt.data.path <- paste0(ext.data.path, '/imgt')
hla.allele.path <- paste0(imgt.data.path, '/IMGTHLA/allelelist')
# ---> File definition.
# IMGT files.
tr.gene.info.file <- paste0(imgt.data.path, '/tr_genes.csv')
# ---> Check directories and files.
if(!dir.exists(reports.path)) dir.create(reports.path)
essential.files <- c(
  tr.gene.info.file
)
essential.files <- essential.files[!unlist(lapply(X=essential.files, FUN=file.exists))]
if(length(essential.files) > 0) stop(paste0('Next essential files are not appropriately defined -make sure you\'ve got their names right-:\n', paste0(essential.files, collapse='\n'), '\n'))


cat('\n\n')
### --------------------------- Functions --------------------------- ###
cat('### --------------------------- Functions --------------------------- ###\n')

# ---> General vj gene info clean-up process
clean.vdj.info <- function(tcr.data, data.base){
    # @ Clean VDJ gene usage info.
    tcr.data[,
        `:=`(
            trav.base=str_replace(
                string=str_replace_all(string=trav.full, pattern='\\*\\d+(:)*|^\\s+|\\s+$', replacement=''),
                pattern='TCRAV', replacement='TRAV'
            ),
            traj.base=str_replace(
                string=str_replace_all(string=traj.full, pattern='\\*\\d+(:)*|^\\s+|\\s+$', replacement=''),
                pattern='TCRAJ', replacement='TRAJ'
            ),
            trbv.base=str_replace(
                string=str_replace_all(string=trbv.full, pattern='\\*\\d+(:)*|^\\s+|\\s+$', replacement=''),
                pattern='TCRBV', replacement='TRBV'
            ),
            trbj.base=str_replace(
                string=str_replace_all(string=trbj.full, pattern='\\*\\d+(:)*|^\\s+|\\s+$', replacement=''),
                pattern='TCRBJ', replacement='TRBJ'
            )
        )
    ]
    tcr.gene.cols <- c('trav.base', 'traj.base', 'trbv.base', 'trbj.base')
    # Let's give it a second chance to those genes w/ suffix "-1" that can't be found in the IMGT reference at this point, this by trimming such a suffix.
    for(gene.col in tcr.gene.cols){
        tcr.data[, tmp.col:=get(gene.col)]
        tcr.data[
            !is.na(tmp.col) & !(tmp.col %in% tr.gene.info[, IMGT.GENE.DB]),
            `:=`(
                tmp.col=str_replace(
                    string=tmp.col,
                    pattern='-0*1$', replacement=''
                )
            )
        ]
        cols.to.keep <- setdiff(x=colnames(tcr.data), y=gene.col)
        tcr.data <- tcr.data[, ..cols.to.keep]
        colnames(tcr.data)[colnames(tcr.data)=='tmp.col'] <- gene.col
    }
    # Let's give it a second chance to those genes w/ suffix "-0\d" that can't be found in the IMGT reference at this point, this by trimming the 0 before the last suffix digit.
    for(gene.col in tcr.gene.cols){
        tcr.data[, tmp.col:=get(gene.col)]
        tcr.data[
            !is.na(tmp.col) & !(tmp.col %in% tr.gene.info[, IMGT.GENE.DB]) & str_detect(string=tmp.col, pattern='-0\\d+$'),
            `:=`(
                tmp.col=str_replace(
                    string=tmp.col,
                    pattern='-0', replacement='-'
                )
            )
        ]
        cols.to.keep <- setdiff(x=colnames(tcr.data), y=gene.col)
        tcr.data <- tcr.data[, ..cols.to.keep]
        colnames(tcr.data)[colnames(tcr.data)=='tmp.col'] <- gene.col
    }
    # Let's give it a second chance to those genes w/ a colon for entries in the MCPas DB.
    if(data.base == 'MCPas'){
        for(gene.col in tcr.gene.cols){
            tcr.clean[, tmp.col:=get(gene.col)]
            tcr.clean[
                !is.na(tmp.col) & !(tmp.col %in% tr.gene.info[, IMGT.GENE.DB]) & tmp.col %like% ':',
                `:=`(
                    tmp.col=str_replace(
                        string=tmp.col,
                        pattern=':.+$', replacement=''
                    )
                )
            ]
            cols.to.keep <- setdiff(x=colnames(tcr.clean), y=gene.col)
            tcr.clean <- tcr.clean[, ..cols.to.keep]
            colnames(tcr.clean)[colnames(tcr.clean)=='tmp.col'] <- gene.col
        }
    }
    # Sanity check. Confirm validity of VDJ gene.
    tmp.check <- c(
        `trav`=tcr.data[!is.na(trav.base) & !trav.base %in% tr.gene.info[, IMGT.GENE.DB], .N]/tcr.data[!is.na(trav.base), .N],
        `traj`=tcr.data[!is.na(traj.base) & !traj.base %in% tr.gene.info[, IMGT.GENE.DB], .N]/tcr.data[!is.na(traj.base), .N],
        `trbv`=tcr.data[!is.na(trbv.base) & !trbv.base %in% tr.gene.info[, IMGT.GENE.DB], .N]/tcr.data[!is.na(trbv.base), .N],
        `trbj`=tcr.data[!is.na(trbj.base) & !trbj.base %in% tr.gene.info[, IMGT.GENE.DB], .N]/tcr.data[!is.na(trbj.base), .N]
    )
    if(any(tmp.check>0.25)){
        tmp.check <- tmp.check[tmp.check>0.25]
        stop(
            paste0(
                'Error for database ', data.base, '. Number of unrecognized TCR genes was too high for the following genes:\n',
                paste0(paste(names(tmp.check), tmp.check, sep=': '), collapse='\n'),
                '\n'
            )
        )
    }else{
        tcr.data[!is.na(trav.base) & !trav.base %in% tr.gene.info[, IMGT.GENE.DB], trav.base:=NA]
        tcr.data[!is.na(traj.base) & !traj.base %in% tr.gene.info[, IMGT.GENE.DB], traj.base:=NA]
        tcr.data[!is.na(trbv.base) & !trbv.base %in% tr.gene.info[, IMGT.GENE.DB], trbv.base:=NA]
        tcr.data[!is.na(trbj.base) & !trbj.base %in% tr.gene.info[, IMGT.GENE.DB], trbj.base:=NA]
    }
    return(tcr.data)
}


cat('\n\n')
### ------------------------- Loading data -------------------------- ###
cat('### ------------------------- Loading data -------------------------- ###\n')

# ---> IMGT files
# VDJ gene info.
tr.gene.info <- fread(file=tr.gene.info.file)
tr.gene.info <- unique(tr.gene.info) # Quite surprised that this is even necessary. It is, though.
colnames(tr.gene.info) <- str_replace_all(string=colnames(tr.gene.info), pattern='\\s+|/|-|\\*', replacement='.')
tr.gene.info[, V14:=NULL]
# HLA allele info.
# All of them have the same prefix.
hla.allele.info.files <- list.files(path=hla.allele.path, full.names=TRUE, pattern='^Allelelist')
hla.allele.info <- lapply(X=hla.allele.info.files, FUN=fread, skip='AlleleID,Allele')
hla.allele.info <- rbindlist(l=hla.allele.info, use.names=TRUE)
# Basic gene and allele info. Essentially, the first two pairs of digits.
hla.allele.info[,
    hla.base.gene:=paste0(
        'HLA-',
        str_replace(
            string=str_extract(string=Allele, pattern='^[A-Z0-9]+\\*'),
            pattern='\\*', replacement=''
        )
    )
]
hla.allele.info[,
    hla.base.allele:=paste0(
        'HLA-',
        str_replace(
            string=str_extract(string=Allele, pattern='^[A-Z0-9]+\\*(\\d{2}:*){2}'),
            pattern=':$', replacement=''
        )
    )
]
# hla.allele.info[is.na(hla.base.gene)]
# hla.allele.info[is.na(hla.base.allele)]
# hla.allele.info[, unique(hla.base.gene)] # All make sense based on manual inspection. :)
# hla.allele.info[, uniqueN(hla.base.allele)]


# ---> vdjdb files
# Full human dataset.


cat('\n\n')
### ------------------------- Main program -------------------------- ###
cat('### ------------------------- Main program -------------------------- ###\n')


### ------------------------- General stuff ------------------------- ###

# ---> Final reference info.
# List that'll store the final info from each DB.
db.info <- list()


cat('\n\n')
### ------------------------- IEDB process -------------------------- ###
cat('### ------------------------- IEDB process -------------------------- ###\n')

# ---> Files definition.
ref.info.file <- paste0(iedb.data.path, '/reference_full_v3.csv')
ag.info.file <- paste0(iedb.data.path, '/antigen_full_v3.csv')
assay.info.file <- paste0(iedb.data.path, '/assay_human_v3.csv')
epi.info.file <- paste0(iedb.data.path, '/epitope_full_v3.csv')
mhc.info.file <- paste0(iedb.data.path, '/mhc_ligand_full.csv')
tcell.info.file <- paste0(iedb.data.path, '/tcell_full_v3.csv')
tcr.info.file <- paste0(iedb.data.path, '/tcr_full_v3_modified.csv')

# ---> Loading data
ref.info <- fread(file=ref.info.file, na.strings=c('', 'NA', 'na'))
ag.info <- fread(file=ag.info.file, na.strings=c('', 'NA', 'na'))
assay.info <- fread(file=assay.info.file, na.strings=c('', 'NA', 'na'))
epi.info <- fread(file=epi.info.file, na.strings=c('', 'NA', 'na'))
mhc.info <- fread(file=mhc.info.file, na.strings=c('', 'NA', 'na'))
tcell.info <- fread(file=tcell.info.file, na.strings=c('', 'NA', 'na'))
tcr.info <- fread(file=tcr.info.file, na.strings=c('', 'NA', 'na'))

# ---> Data preprocessing
# @ Tidy columns per type of data.
# Function to turn double-line header into 
tidy.data.cols <- function(tmp.data){
    tmp.cols <- paste(
        str_replace_all(string=colnames(tmp.data), pattern='\\s+', replacement='.'),
        str_replace_all(string=tmp.data[1, ], pattern='\\s+', replacement='.'),
        sep='.')
    tmp.cols <- str_replace_all(string=tmp.cols, pattern='#', replacement='No')
    colnames(tmp.data) <- tmp.cols
    tmp.data <- tmp.data[2:.N]
    return(tmp.data)
}
# Tidy column names for all of the files.
ref.info <- tidy.data.cols(ref.info)
ag.info <- tidy.data.cols(ag.info)
assay.info <- tidy.data.cols(assay.info)
epi.info <- tidy.data.cols(epi.info)
mhc.info <- tidy.data.cols(mhc.info)
tcell.info <- tidy.data.cols(tcell.info)
tcr.info <- tidy.data.cols(tcr.info)

# @ Assay information.
# There's too many fields for this table so previous cleaning facilitates the rest of the proces.
cols.to.keep <- colnames(assay.info)[str_detect(string=colnames(assay.info), pattern='Assay.')]
cols.to.keep <- cols.to.keep[!str_detect(string=cols.to.keep, pattern='Assay.Antigen')]
assay.info <- assay.info[, ..cols.to.keep]
assay.info[, Assay.ID:=str_extract(string=Assay.ID.IEDB.IRI, pattern='\\d+$')]

# @ Exclude uninformative reference.
# Identify reference ID for uninformative reference based on its PMID.
tmp.data <- ref.info[Reference.PMID=='32793919', Reference.ID.IEDB.IRI]
tcr.info <- tcr.info[Reference.IEDB.IRI!=tmp.data]

# ---> Data exploration.

# tcr.info[, .N, by=.(Assay.Type)] # Useless entry for this particular table. All of them have the value "T cell"

# @ Number of DB entries for each type of chain pair
tmp.data <- tcr.info[, .N, by=.(Chain.1.Type, Chain.2.Type)]
tmp.data[, tcr.type:=paste(Chain.1.Type, Chain.2.Type, sep='/')]
tmp.data$tcr.type <- factor(x=tmp.data$tcr.type, levels=tmp.data$tcr.type)
setorderv(x=tmp.data, cols='tcr.type')
tmp.ggplot.1 <- ggplot(data=tmp.data, aes(x=tcr.type, y=N)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Number of DB entries for each type of chain pair', x='Chain pair type', y='Count') +
    theme_bw()

# @ Number of entries where FULL sequences are provided for each type of chain.
tmp.data <- tcr.info[
    Chain.1.Type=='alpha' | Chain.2.Type=='beta',
    .(
        alpha.full=!is.na(Chain.1.Protein.Sequence), beta.full=!is.na(Chain.2.Protein.Sequence)
    )
]
tmp.data <- tmp.data[, .N, by=.(alpha.full, beta.full)]
tmp.data[,
    full.seq.summ:=paste0(
        'alpha: ', alpha.full, ' - ',
        'beta: ', beta.full
    )
]
tmp.ggplot.2 <- ggplot(data=tmp.data, aes(x=full.seq.summ, y=N)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Number of DB entries with FULL chain sequence provided', x='Chain pair type', y='Count', caption='TCRs that were not alpha/beta were discarded here') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))

# @ Number of DB entries where CDR3 sequence is provided for each type of chain.
tmp.data <- tcr.info[
    Chain.1.Type=='alpha' | Chain.2.Type=='beta',
    .(
        alpha.determined.cdr3=!is.na(Chain.1.CDR3.Curated), alpha.calculated.cdr3=!is.na(Chain.1.CDR3.Calculated),
        beta.determined.cdr3=!is.na(Chain.2.CDR3.Curated), beta.calculated.cdr3=!is.na(Chain.2.CDR3.Calculated)
    )
]
tmp.data <- tmp.data[, .N, by=.(alpha.determined.cdr3, alpha.calculated.cdr3, beta.determined.cdr3, beta.calculated.cdr3)]
tmp.data[, beta.cdr3:=FALSE]
tmp.data[beta.determined.cdr3 | beta.calculated.cdr3, beta.cdr3:=TRUE]
tmp.data[, alpha.cdr3:=FALSE]
tmp.data[alpha.determined.cdr3 | alpha.calculated.cdr3, alpha.cdr3:=TRUE]
tmp.data <- tmp.data[, .(N=sum(N)), by=.(beta.cdr3, alpha.cdr3)]
tmp.data[,
    cdr3.seq.summ:=paste0(
        'alpha: ', alpha.cdr3, ' - ',
        'beta: ', beta.cdr3
    )
]
tmp.ggplot.3 <- ggplot(data=tmp.data, aes(x=cdr3.seq.summ, y=N)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Number of DB entries with CDR3 sequence provided considering processed CDR3 sequences as well', x='Chain pair type', y='Count', caption='TCRs that were not alpha/beta were discarded here\nProcessed CDR3 sequences were determined based on full chain sequence.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))


# @ Frequency of DB entries and unique epitopes recognized by unique alpha/beta clonotypes.
# A clonotype is defined as a set of alpha and beta CDR3 and TRAV, TRBV, TRAJ and TRBJ sequences.
tmp.data <- tcr.info[
    Chain.1.Type=='alpha' | Chain.2.Type=='beta',
    .(
        alpha.chain=ifelse(
            test=!is.na(Chain.1.CDR3.Curated), yes=Chain.1.CDR3.Curated, no=Chain.1.CDR3.Calculated
        ),
        beta.chain=ifelse(
            test=!is.na(Chain.2.CDR3.Curated), yes=Chain.2.CDR3.Curated, no=Chain.2.CDR3.Calculated
        ),
        TRAV=ifelse(
            test=!is.na(Chain.1.Curated.V.Gene), yes=Chain.1.Curated.V.Gene, no=Chain.1.Calculated.V.Gene
        ),
        TRAJ=ifelse(
            test=!is.na(Chain.1.Curated.J.Gene), yes=Chain.1.Curated.J.Gene, no=Chain.1.Calculated.J.Gene
        ),
        TRBV=ifelse(
            test=!is.na(Chain.2.Curated.V.Gene), yes=Chain.2.Curated.V.Gene, no=Chain.2.Calculated.V.Gene
        ),
        TRBJ=ifelse(
            test=!is.na(Chain.2.Curated.J.Gene), yes=Chain.2.Curated.J.Gene, no=Chain.2.Calculated.J.Gene
        ),
        Epitope.Name
    )
]
# Remove V and J gene allele information.
tmp.data[, TRAV:=str_replace(string=TRAV, pattern='\\*.+$', replacement='')]
tmp.data[, TRAJ:=str_replace(string=TRAJ, pattern='\\*.+$', replacement='')]
tmp.data[, TRBV:=str_replace(string=TRBV, pattern='\\*.+$', replacement='')]
tmp.data[, TRBJ:=str_replace(string=TRBV, pattern='\\*.+$', replacement='')]
# Identify total number of TCRs that react to each possible number of different epitopes.
tmp.data <- tmp.data[,
    .(entry.count=.N, epitope.count=uniqueN(Epitope.Name)),
    by=.(alpha.chain, beta.chain, TRAV, TRAJ, TRBV, TRBJ)
]
tmp.data[entry.count>10, entry.count:=10]
tmp.data[epitope.count>10, epitope.count:=10]
tmp.ggplot.4 <- ggplot(data=tmp.data, aes(x=entry.count)) +
    geom_histogram(bins=10, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Distribution of frequency of DB entries for unique clonotypes', x='Number of DB entries', y='Count', caption='TCRs that were not alpha/beta were discarded here\nNumber of entries was capped at 10.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))
tmp.ggplot.5 <- ggplot(data=tmp.data, aes(x=epitope.count)) +
    geom_histogram(bins=10, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Distribution of frequency of epitopes recognized by unique clonotypes', x='Number of recognized epitopes', y='Count', caption='TCRs that were not alpha/beta were discarded here\nNumber of epitopes was capped at 10.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))

# @ Distribution of frequency of epitopes per organism.
tmp.data <- tcr.info[
    Chain.1.Type=='alpha' | Chain.2.Type=='beta',
    .(Epitope.Name, Epitope.Source.Molecule, Epitope.Source.Organism)
]
tmp.data <- tmp.data[!is.na(Epitope.Name) & !is.na(Epitope.Source.Organism), .(N=uniqueN(Epitope.Name)), by=.(Epitope.Source.Organism)]
tmp.ggplot.6 <- ggplot(data=tmp.data, aes(x=N)) +
    geom_histogram(bins=10, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Distribution of frequency of epitopes per organism', x='Number of epitopes per organism', y='Count', caption='TCRs that were not alpha/beta were discarded here\nNumber of epitopes was capped at 10.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))

# @ Top 10 organisms w/ greatest number of epitopes assessed.
tmp.data <- tcr.info[
    Chain.1.Type=='alpha' | Chain.2.Type=='beta',
    .(Epitope.Name, Epitope.Source.Molecule, Epitope.Source.Organism)
]
tmp.data <- tmp.data[!is.na(Epitope.Name) & !is.na(Epitope.Source.Organism), .(N=uniqueN(Epitope.Name)), by=.(Epitope.Source.Organism)]
setorderv(x=tmp.data, cols='N', order=-1)
tmp.data <- tmp.data[1:10]
tmp.data$Epitope.Source.Organism <- factor(x=tmp.data$Epitope.Source.Organism, levels=tmp.data$Epitope.Source.Organism)
tmp.ggplot.7 <- ggplot(data=tmp.data, aes(x=Epitope.Source.Organism, y=N)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Top 10 organisms w/ greatest number of epitopes assessed', x='Organism', y='Count', caption='TCRs that were not alpha/beta were discarded here.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))

# @ Top 10 organisms w/ greatest number of DB entries.
tmp.data <- tcr.info[
    Chain.1.Type=='alpha' | Chain.2.Type=='beta',
    .(Epitope.Name, Epitope.Source.Molecule, Epitope.Source.Organism)
]
tmp.data <- tmp.data[!is.na(Epitope.Name) & !is.na(Epitope.Source.Organism), .(N=.N), by=.(Epitope.Source.Organism)]
setorderv(x=tmp.data, cols='N', order=-1)
tmp.data <- tmp.data[1:10]
tmp.data$Epitope.Source.Organism <- factor(x=tmp.data$Epitope.Source.Organism, levels=tmp.data$Epitope.Source.Organism)
tmp.ggplot.8 <- ggplot(data=tmp.data, aes(x=Epitope.Source.Organism, y=N)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Top 10 organisms w/ greatest number of DB entries', x='Organism', y='Count', caption='TCRs that were not alpha/beta were discarded here.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))

# @ Distribution of frequency of DB entries per epitope.
tmp.data <- tcr.info[
    Chain.1.Type=='alpha' | Chain.2.Type=='beta',
    .(Epitope.Name, Epitope.Source.Molecule, Epitope.Source.Organism)
]
tmp.data <- tmp.data[!is.na(Epitope.Name) & !is.na(Epitope.Source.Organism), .(N=.N), by=.(Epitope.Name, Epitope.Source.Organism)]
setorderv(x=tmp.data, cols='N', order=-1)
tmp.data[N>10, N:=10]
tmp.ggplot.9 <- ggplot(data=tmp.data, aes(x=N)) +
    geom_histogram(bins=10, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Distribution of frequency of DB entries per epitopes', x='Number of epitopes-specific DB entries', y='Count', caption='TCRs that were not alpha/beta were discarded here\nNumber of DB entries was capped at 10.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))

# @ Top 10 epitopes w/ greatest number of DB entries.
tmp.data <- tcr.info[
    Chain.1.Type=='alpha' | Chain.2.Type=='beta',
    .(Epitope.Name, Epitope.Source.Molecule, Epitope.Source.Organism)
]
tmp.data <- tmp.data[!is.na(Epitope.Name) & !is.na(Epitope.Source.Organism), .(N=.N), by=.(Epitope.Name, Epitope.Source.Organism)]
tmp.data[, full.epitope:=paste(Epitope.Source.Organism, Epitope.Name, sep=' / ')]
setorderv(x=tmp.data, cols='N', order=-1)
tmp.data <- tmp.data[1:10]
tmp.data$Epitope.Name <- factor(x=tmp.data$Epitope.Name, levels=tmp.data$Epitope.Name)
tmp.ggplot.10 <- ggplot(data=tmp.data, aes(x=full.epitope, y=N)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Top 10 epitopes w/ greatest number of DB entries', x='Epitope', y='DB entries', caption='TCRs that were not alpha/beta were discarded here\nNumber of epitopes was capped at 10.') +
    theme_bw() + theme(axis.text.x=element_text(angle=45))

# ---> Retrieve info required for compilation.
# @ Pairs of TCRs and epitopes, along with their corresponding relevant information.
tcr.clean <- tcr.info[
    Chain.1.Type=='alpha' | Chain.2.Type=='beta',
    .(
        alpha.chain=ifelse(
            test=!is.na(Chain.1.CDR3.Curated), yes=Chain.1.CDR3.Curated, no=Chain.1.CDR3.Calculated
        ),
        beta.chain=ifelse(
            test=!is.na(Chain.2.CDR3.Curated), yes=Chain.2.CDR3.Curated, no=Chain.2.CDR3.Calculated
        ),
        alpha.chain.full=Chain.1.Protein.Sequence,
        beta.chain.full=Chain.2.Protein.Sequence,
        trav.full=ifelse(
            test=!is.na(Chain.1.Curated.V.Gene), yes=Chain.1.Curated.V.Gene, no=Chain.1.Calculated.V.Gene
        ),
        traj.full=ifelse(
            test=!is.na(Chain.1.Curated.J.Gene), yes=Chain.1.Curated.J.Gene, no=Chain.1.Calculated.J.Gene
        ),
        trbv.full=ifelse(
            test=!is.na(Chain.2.Curated.V.Gene), yes=Chain.2.Curated.V.Gene, no=Chain.2.Calculated.V.Gene
        ),
        trbj.full=ifelse(
            test=!is.na(Chain.2.Curated.J.Gene), yes=Chain.2.Curated.J.Gene, no=Chain.2.Calculated.J.Gene
        ),
        epitope.seq=Epitope.Name,
        antigen=Epitope.Source.Molecule,
        epitope.organism=Epitope.Source.Organism,
        assay.id=Assay.IEDB.IDs,
        assay.mhcs=Assay.MHC.Allele.Names,
        entry.ref=Reference.IEDB.IRI
    )
]

# @ Relevant assay information.
assay.clean <- assay.info[,
    .(
        assay.id=Assay.ID,
        assay.method=Assay.Method
    )
]

# @ Create custom assay reference by considering "combined" assay IDs. Specifically, we'll combine the method description of individual assays only for the "combined" assay IDs that exists in the TCR table.
tmp.data <- tcr.clean[assay.id %like% ", ", unique(assay.id)]
tmp.data.1 <- str_split(string=tmp.data, pattern=', ')
tmp.data.2 <- unlist(lapply(X=tmp.data.1, FUN=function(x){
    tmp.data <- assay.clean[assay.id %in% x, unique(assay.method)]
    is.null
    tmp.data <- paste0(sort(tmp.data), collapse=';')
    return(tmp.data)
}))
tmp.data <- data.table(
    assay.id=tmp.data,
    assay.method=tmp.data.2
)
assay.clean <- rbind(assay.clean, tmp.data)
assay.clean <- assay.clean[assay.method!='']

# @ Merge TCR and assay info.
tcr.clean <- merge(x=tcr.clean, y=assay.clean, by='assay.id', all.x=TRUE)
# @ Visualization of the different number of DB entries we have per assay method.
tmp.data <- tcr.clean[, .(assay.method)]
tmp.data[is.na(assay.method), assay.method:='None']
tmp.data[str_detect(string=assay.method, pattern=';'), assay.method:='Multiple']
tmp.data.2 <- tmp.data[, .N, by=assay.method]
setorderv(x=tmp.data.2, cols='N', order=-1)
tmp.data$assay.method <- factor(x=tmp.data$assay.method, levels=tmp.data.2$assay.method)
tmp.ggplot.11 <- ggplot(data=tmp.data, aes(x=assay.method)) +
    geom_bar(width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Number of DB entries per assay method', x='Assay method', y='DB entries', caption='TCRs that were not alpha/beta were discarded here\nNumber of epitopes was capped at 10.') +
    theme_bw() + theme(axis.text.x=element_text(angle=45))

# @ Reference information.
# tcr.clean[!entry.ref %in% ref.info[, Reference.ID.IEDB.IRI]]
tmp.data <- ref.info[, .(entry.ref=Reference.ID.IEDB.IRI, pmid=Reference.PMID)]
tcr.clean <- merge(x=tcr.clean, y=tmp.data, by='entry.ref', all.x=TRUE, all.y=FALSE)
# Visualization.
tmp.data <- tcr.clean[, .(pmid)]
tmp.data[is.na(pmid), pmid:='None']
tmp.data.2 <- tmp.data[, .N, by=pmid]
tmp.data.2[N>20, N:=20]
tmp.ggplot.12 <- ggplot(data=tmp.data.2, aes(x=N)) +
    geom_histogram(bins=10, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Distribution of the number of DB entries per PMID', x='Number of DB entries', y='Frequency', caption='TCRs that were not alpha/beta were discarded here\nNumber of DB entries was capped at 20.') +
    theme_bw() + theme(axis.text.x=element_text(angle=45))

# @ Clean MHC allele information.
tcr.clean[assay.mhcs %like% "^HLA-[A-Z0-9]+(\\*\\d+:\\d+)*$", mhc:=assay.mhcs]
tcr.clean[is.na(mhc) & !is.na(assay.mhcs), unique(assay.mhcs)]
# Basic gene information if provided.
tcr.clean[,
    mhc.gene:=str_replace(
        string=str_extract(string=mhc, pattern='HLA-[A-Z0-9]+\\**'),
        pattern='\\*', replacement=''
    )
]
tmp.check <- tcr.clean[!is.na(mhc.gene) & !mhc.gene %in% hla.allele.info[, hla.base.gene], .N]/tcr.clean[!is.na(mhc), .N]
if(tmp.check>0.05) stop('We failed to validate the HLA gene name for a fairly great amount of DB entries with HLA information.\n')
# Remove MHC base gene info if MHC gene not reported on IMGT.
tcr.clean[!is.na(mhc.gene) & !mhc.gene %in% hla.allele.info[, hla.base.gene], mhc.gene:=NA]
# Basic allele information if provided.
tcr.clean[
    !is.na(mhc.gene),
    mhc.allele:=str_replace(
        string=str_extract(string=mhc, pattern='HLA-[A-Z0-9]+\\*(\\d{2}:*){2}'),
        pattern=':$', replacement=''
    )
]
# Sanity check. Confirm validity of HLA alleles.
tmp.check <- tcr.clean[!is.na(mhc.allele) & !mhc.allele %in% hla.allele.info[, hla.base.allele], .N==0]
if(!tmp.check) stop('Some of the identified alleles in the DB cannot be confirmed based on IMGT.\n')

# @ Visualization of the different number of DB entries we have per HLA gene.
tmp.data <- tcr.clean[, .(mhc.gene)]
tmp.data[is.na(mhc.gene), mhc.gene:='None']
tmp.data.2 <- tmp.data[, .N, by=mhc.gene]
setorderv(x=tmp.data.2, cols='N', order=-1)
tmp.data$mhc.gene <- factor(x=tmp.data$mhc.gene, levels=tmp.data.2$mhc.gene)
tmp.ggplot.13 <- ggplot(data=tmp.data, aes(x=mhc.gene)) +
    geom_bar(width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Number of DB entries per HLA gene', x='HLA gene', y='DB entries', caption='TCRs that were not alpha/beta were discarded here\nNumber of epitopes was capped at 10.') +
    theme_bw() + theme(axis.text.x=element_text(angle=45))

# @ Clean VDJ gene usage info.
tcr.clean <- clean.vdj.info(tcr.data=tcr.clean, data.base='IEDB')

# ---> Output all plots.
tmp.file.name <- paste0(reports.path, '/BasicExploration_IEDB.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot.1)
print(tmp.ggplot.2)
print(tmp.ggplot.3)
print(tmp.ggplot.4)
print(tmp.ggplot.5)
print(tmp.ggplot.6)
print(tmp.ggplot.7)
print(tmp.ggplot.8)
print(tmp.ggplot.9)
print(tmp.ggplot.10)
print(tmp.ggplot.11)
print(tmp.ggplot.12)
print(tmp.ggplot.13)
dev.off()

# ---> Final format of the table.
cols.order <- c(
    'alpha.chain', 'beta.chain',
    'trav.full', 'traj.full', 'trbv.full', 'trbj.full',
    'trav.base', 'traj.base', 'trbv.base', 'trbj.base',
    'epitope.organism', 'antigen', 'epitope.seq',
    'assay.method',
    'mhc', 'mhc.gene', 'mhc.allele',
    'pmid',
    'assay.mhcs'
)
tcr.clean <- tcr.clean[, ..cols.order]
setorderv(x=tcr.clean, cols='assay.method')
db.info[['iedb']] <- copy(tcr.clean)


cat('\n\n')
### ------------------------- VDJdb process ------------------------- ###
cat('### ------------------------- VDJdb process ------------------------- ###\n')

# ---> Files definition.
vdjdb.human.file <- paste0(vdjdb.data.path, '/vdjdb_human.tsv')

# ---> Loading data
tcr.info <- fread(file=vdjdb.human.file, na.strings=c('', 'NA', 'na'))

# ---> Data preprocessing

# @ Pairing information.
# For clonotypes that have chain pair info available, let's make them single entries. According to the DB developers: "The column complex.id contains pairing information. TRA and TRB that share an id in the list provided in this columns are known to be paired in one of the publications. If there is no information about paired record in database the complex.id field will be equal to 0 in the exported output data."

# Provide unique complex IDs to entries that have no pairing info available (i.e., whose complex ID is 0)
tcr.info[, complex.id:=as.character(complex.id)]
tcr.info[complex.id==0, complex.id:=paste(0, 1:.N, sep='.')]

# Label columns that are chain-specific.
chain.spc.cols <- c('CDR3', 'V', 'J', 'CDR3fix')
tmp.cols <- setdiff(x=colnames(tcr.info), y=c('Gene', chain.spc.cols))
tmp.data.1 <- unique(tcr.info[Gene=='TRA', ])
tmp.data.2 <- unique(tcr.info[Gene=='TRB', ])
tcr.info <- merge(x=tmp.data.1, y=tmp.data.2, by=tmp.cols, all=TRUE, suffixes=c('.tra', '.trb'), allow.cartesian=FALSE, sort=FALSE)

# @ Unique identifier for each DB entry.
tcr.info[, entry:=1:.N]

# @ Tidy columns.
# Function to spread the info from multiple variables provided into a single packed column.
tidy.data.cols.vdjdb <- function(input.data, packed.col){
    input.data[, tmp.col:=get(packed.col)]
    input.data[, tmp.col:=str_replace_all(string=tmp.col, pattern='\\{|\\}', replacement='')]
    tmp.data <- str_split(string=input.data[, tmp.col], pattern=', ')
    tmp.data <- lapply(X=tmp.data, FUN=function(x){
        tmp.data <- data.table(var=x)
        tmp.data <- separate(data=tmp.data, col='var', into=c('var', 'value'), sep=': ')
    })
    tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='entry')
    tmp.data[, `:=`(
        var=str_replace_all(string=var, pattern='"', replacement=''),
        value=str_replace_all(string=value, pattern='"', replacement='')
    )]
    tmp.data <- unique(tmp.data)
    tmp.data <- as.data.table(spread(data=tmp.data, key=var, value=value, fill=TRUE))
    tmp.data[tmp.data==""] <- NA
    colnames(tmp.data) <- paste(str_to_lower(string=packed.col), colnames(tmp.data), sep='.')
    input.data <- merge(x=input.data, y=tmp.data, by.x='entry', by.y=paste0(str_to_lower(string=packed.col), '.entry'), all.x=TRUE, sort=FALSE)
    cols.to.rm <- c('tmp.col', packed.col)
    tmp.cols <- setdiff(x=colnames(input.data), y=cols.to.rm)
    input.data <- input.data[, ..tmp.cols]
    return(input.data)
}
# Apply for each single packed column.
tcr.info <- tidy.data.cols.vdjdb(input.data=tcr.info, packed.col='Method')
tcr.info <- tidy.data.cols.vdjdb(input.data=tcr.info, packed.col='Meta')
# vdjdb.human <- tidy.data.cols.vdjdb(input.data=vdjdb.human, packed.col='CDR3fix')

# @ Remove columns w/ quite specific and hence unnecessary HLA information.
tmp.cols <- colnames(tcr.info)
tmp.cols <- tmp.cols[!str_detect(string=tmp.cols, pattern='^meta.HLA')]
tcr.info <- tcr.info[, ..tmp.cols]

# ---> Data exploration.

# @ Number of epitope/TCR pair entries for each type of chain pair
tmp.data <- tcr.info[,
    .(
        alpha.full=ifelse(
            test=!is.na(CDR3.tra),
            yes='TRA', no=NA
        ),
        beta.full=ifelse(
            test=!is.na(CDR3.trb),
            yes='TRB', no=NA
        )
    )
]
tmp.data <- tmp.data[, .N, by=.(alpha.full, beta.full)]
tmp.data[, tcr.type:=paste(alpha.full, beta.full, sep='/')]
tmp.data$tcr.type <- factor(x=tmp.data$tcr.type, levels=tmp.data$tcr.type)
setorderv(x=tmp.data, cols='tcr.type')
tmp.ggplot.1 <- ggplot(data=tmp.data, aes(x=tcr.type, y=N)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Number of DB entries for each type of chain pair', x='Chain pair type', y='Count') +
    theme_bw()

# @ Frequency of DB entries and unique epitopes recognized by unique alpha/beta clonotypes.
# A clonotype is defined as a set of alpha and beta CDR3 and TRAV, TRBV, TRAJ and TRBJ sequences.
tmp.data <- tcr.info[,
    .(
        alpha.chain=CDR3.tra,
        beta.chain=CDR3.trb,
        TRAV=V.tra,
        TRAJ=J.tra,
        TRBV=V.trb,
        TRBJ=J.trb,
        epitope=Epitope
    )
]
# Remove V and J gene allele information.
tmp.data[, TRAV:=str_replace(string=TRAV, pattern='\\*.+$', replacement='')]
tmp.data[, TRAJ:=str_replace(string=TRAJ, pattern='\\*.+$', replacement='')]
tmp.data[, TRBV:=str_replace(string=TRBV, pattern='\\*.+$', replacement='')]
tmp.data[, TRBJ:=str_replace(string=TRBV, pattern='\\*.+$', replacement='')]
# Identify total number of TCRs that react to each possible number of different epitopes.
tmp.data <- tmp.data[,
    .(entry.count=.N, epitope.count=uniqueN(epitope)),
    by=.(alpha.chain, beta.chain, TRAV, TRAJ, TRBV, TRBJ)
]
tmp.data[entry.count>5, entry.count:=5]
tmp.data[epitope.count>5, epitope.count:=5]
tmp.ggplot.2 <- ggplot(data=tmp.data, aes(x=entry.count)) +
    geom_histogram(bins=5, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Distribution of frequency of DB entries for unique clonotypes', x='Number of DB entries', y='Count', caption='TCRs that were not alpha/beta were discarded here\nNumber of entries was capped at 5.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))
tmp.ggplot.3 <- ggplot(data=tmp.data, aes(x=epitope.count)) +
    geom_histogram(bins=5, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Distribution of frequency of epitopes recognized by unique clonotypes', x='Number of recognized epitopes', y='Count', caption='TCRs that were not alpha/beta were discarded here\nNumber of epitopes was capped at 5.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))

# @ Distribution of frequency of epitopes per organism.
tmp.data <- tcr.info[,
    .(epitope.name=Epitope, epitope.gene=`Epitope gene`, organism=`Epitope species`)
]
tmp.data <- tmp.data[!is.na(epitope.name) & !is.na(organism), .(N=uniqueN(epitope.name)), by=.(organism)]
tmp.ggplot.4 <- ggplot(data=tmp.data, aes(x=N)) +
    geom_histogram(bins=10, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Distribution of frequency of epitopes per organism', x='Number of epitopes per organism', y='Count', caption='TCRs that were not alpha/beta were discarded here\nNumber of epitopes was capped at 10.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))

# @ Top 10 organisms w/ greatest number of epitopes assessed.
tmp.data <- tcr.info[,
    .(epitope.name=Epitope, epitope.gene=`Epitope gene`, organism=`Epitope species`)
]
tmp.data <- tmp.data[!is.na(epitope.name) & !is.na(organism), .(N=uniqueN(epitope.name)), by=.(organism)]
setorderv(x=tmp.data, cols='N', order=-1)
tmp.data <- tmp.data[1:10]
tmp.data$organism <- factor(x=tmp.data$organism, levels=tmp.data$organism)
tmp.ggplot.5 <- ggplot(data=tmp.data, aes(x=organism, y=N)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Top 10 organisms w/ greatest number of epitopes assessed', x='Organism', y='Count', caption='TCRs that were not alpha/beta were discarded here\nNumber of epitopes was capped at 10.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))

# @ Top 10 organisms w/ greatest number of DB entries.
tmp.data <- tcr.info[,
    .(epitope.name=Epitope, epitope.gene=`Epitope gene`, organism=`Epitope species`)
]
tmp.data <- tmp.data[!is.na(epitope.name) & !is.na(organism), .(N=.N), by=.(organism)]
setorderv(x=tmp.data, cols='N', order=-1)
tmp.data <- tmp.data[1:10]
tmp.data$organism <- factor(x=tmp.data$organism, levels=tmp.data$organism)
tmp.ggplot.6 <- ggplot(data=tmp.data, aes(x=organism, y=N)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Top 10 organisms w/ greatest number of DB entries', x='Organism', y='Count', caption='TCRs that were not alpha/beta were discarded here\nNumber of epitopes was capped at 10.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))

# @ Distribution of frequency of DB entries per epitope.
tmp.data <- tcr.info[,
    .(epitope.name=Epitope, epitope.gene=`Epitope gene`, organism=`Epitope species`)
]
tmp.data <- tmp.data[!is.na(epitope.name) & !is.na(organism), .(N=.N), by=.(epitope.name, organism)]
setorderv(x=tmp.data, cols='N', order=-1)
tmp.data[N>10, N:=10]
tmp.ggplot.7 <- ggplot(data=tmp.data, aes(x=N)) +
    geom_histogram(bins=10, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Distribution of frequency of epitopes per organism', x='Number of epitopes per organism', y='Count', caption='TCRs that were not alpha/beta were discarded here\nNumber of epitopes was capped at 10.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))

# @ Top 10 epitopes w/ greatest number of DB entries.
tmp.data <- tcr.info[,
    .(epitope.name=Epitope, epitope.gene=`Epitope gene`, organism=`Epitope species`)
]
tmp.data <- tmp.data[!is.na(epitope.name) & !is.na(organism), .(N=.N), by=.(epitope.name, organism)]
tmp.data[, full.epitope:=paste(organism, epitope.name, sep=' / ')]
setorderv(x=tmp.data, cols='N', order=-1)
tmp.data <- tmp.data[1:10]
tmp.ggplot.8 <- ggplot(data=tmp.data, aes(x=full.epitope, y=N)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Top 10 epitopes w/ greatest number of DB entries', x='Epitope', y='DB entries', caption='TCRs that were not alpha/beta were discarded here\nNumber of epitopes was capped at 10.') +
    theme_bw() + theme(axis.text.x=element_text(angle=45))

# ---> Retrieve info required for compilation.
# @ Pairs of TCRs and epitopes, along with their corresponding relevant information.
tcr.clean <- tcr.info[,
    .(
        alpha.chain=CDR3.tra,
        beta.chain=CDR3.trb,
        trav.full=V.tra,
        traj.full=J.tra,
        trbv.full=V.trb,
        trbj.full=J.trb,
        epitope.seq=Epitope,
        antigen=`Epitope gene`,
        epitope.organism=`Epitope species`,
        assay.method=method.identification,
        mhc=`MHC A`,
        donor.mhcs=meta.donor.MHC,
        pmid=str_replace(string=Reference, pattern='^PMID:', replacement='')
    )
]

# @ Visualization of the different number of DB entries we have per assay method.
tmp.data <- tcr.clean[, .(assay.method)]
tmp.data[is.na(assay.method), assay.method:='None']
tmp.data[str_detect(string=assay.method, pattern=';'), assay.method:='Multiple']
tmp.data.2 <- tmp.data[, .N, by=assay.method]
setorderv(x=tmp.data.2, cols='N', order=-1)
tmp.data$assay.method <- factor(x=tmp.data$assay.method, levels=tmp.data.2$assay.method)
tmp.ggplot.9 <- ggplot(data=tmp.data, aes(x=assay.method)) +
    geom_bar(width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Number of DB entries per assay method', x='Assay method', y='DB entries', caption='TCRs that were not alpha/beta were discarded here\nNumber of epitopes was capped at 10.') +
    theme_bw() + theme(axis.text.x=element_text(angle=45))

# @ Reference information.
# Visualization.
tmp.data <- tcr.clean[, .(pmid)]
tmp.data[is.na(pmid), pmid:='None']
tmp.data.2 <- tmp.data[, .N, by=pmid]
tmp.data.2[N>20, N:=20]
tmp.ggplot.10 <- ggplot(data=tmp.data.2, aes(x=N)) +
    geom_histogram(bins=10, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Distribution of the number of DB entries per PMID', x='Number of DB entries', y='Frequency', caption='TCRs that were not alpha/beta were discarded here\nNumber of DB entries was capped at 20.') +
    theme_bw() + theme(axis.text.x=element_text(angle=45))

# @ Clean MHC allele information.
# Basic gene information if provided.
tcr.clean[,
    mhc.gene:=str_replace(
        string=str_extract(string=mhc, pattern='HLA-[A-Z0-9]+\\**'),
        pattern='\\*', replacement=''
    )
]
tmp.check <- tcr.clean[!is.na(mhc.gene) & !mhc.gene %in% hla.allele.info[, hla.base.gene], .N]/tcr.clean[!is.na(mhc), .N]
if(tmp.check>0.05) stop('We failed to validate the HLA gene name for a fairly great amount of DB entries with HLA information.\n')
# Remove MHC base gene info if MHC gene not reported on IMGT.
tcr.clean[!is.na(mhc.gene) & !mhc.gene %in% hla.allele.info[, hla.base.gene], mhc.gene:=NA]
# Basic allele information if provided.
tcr.clean[
    !is.na(mhc.gene),
    mhc.allele:=str_replace(
        string=str_extract(string=mhc, pattern='HLA-[A-Z0-9]+\\*(\\d{2}:*){2}'),
        pattern=':$', replacement=''
    )
]
# Sanity check. Confirm validity of HLA alleles.
tmp.check <- tcr.clean[!is.na(mhc.allele) & !mhc.allele %in% hla.allele.info[, hla.base.allele], .N==0]
if(!tmp.check){
    warning('Some of the identified alleles in the DB cannot be confirmed based on IMGT.\n')
    # Discard MHC allele info because it doesn't seem fair.
    tcr.clean[
        !is.na(mhc.allele) & !mhc.allele %in% hla.allele.info[, hla.base.allele],
        mhc.allele:=NA
    ]
}

# @ Visualization of the different number of DB entries we have per HLA gene.
tmp.data <- tcr.clean[, .(mhc.gene)]
tmp.data[is.na(mhc.gene), mhc.gene:='None']
tmp.data.2 <- tmp.data[, .N, by=mhc.gene]
setorderv(x=tmp.data.2, cols='N', order=-1)
tmp.data$mhc.gene <- factor(x=tmp.data$mhc.gene, levels=tmp.data.2$mhc.gene)
tmp.ggplot.11 <- ggplot(data=tmp.data, aes(x=mhc.gene)) +
    geom_bar(width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Number of DB entries per HLA gene', x='HLA gene', y='DB entries', caption='TCRs that were not alpha/beta were discarded here\nNumber of epitopes was capped at 10.') +
    theme_bw() + theme(axis.text.x=element_text(angle=45))

# @ Clean VDJ gene usage info.
tcr.clean <- clean.vdj.info(tcr.data=tcr.clean, data.base='IEDB')

# ---> Output all plots.
tmp.file.name <- paste0(reports.path, '/BasicExploration_VDJdb.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot.1)
print(tmp.ggplot.2)
print(tmp.ggplot.3)
print(tmp.ggplot.4)
print(tmp.ggplot.5)
print(tmp.ggplot.6)
print(tmp.ggplot.7)
print(tmp.ggplot.8)
print(tmp.ggplot.9)
print(tmp.ggplot.10)
print(tmp.ggplot.11)
dev.off()

# ---> Final format of the table.
cols.order <- c(
    'alpha.chain', 'beta.chain',
    'trav.full', 'traj.full', 'trbv.full', 'trbj.full',
    'trav.base', 'traj.base', 'trbv.base', 'trbj.base',
    'epitope.organism', 'antigen', 'epitope.seq',
    'assay.method',
    'mhc', 'mhc.gene', 'mhc.allele',
    'pmid',
    'donor.mhcs'
)
tcr.clean <- tcr.clean[, ..cols.order]
setorderv(x=tcr.clean, cols='assay.method')
db.info[['vdjdb']] <- copy(tcr.clean)


cat('\n\n')
### ------------------------- McPAS process ------------------------- ###
cat('### ------------------------- McPAS process ------------------------- ###\n')

# ---> Files definition.
mcpas.info.file <- paste0(mcpas.data.path, '/McPAS-TCR.csv')

# ---> Loading data
tcr.info <- fread(file=mcpas.info.file, na.strings=c('', 'NA', 'na'))

# ---> Data preprocessing

# @ Unique entries for TCRs reported to have been found in Homo sapiens and which have associated chain and epitope sequences. Also, at least one piece of info related to the epitope has to exist in a valid entry.
tcr.info <- unique(tcr.info)
tcr.info <- tcr.info[!is.na(Species) & Species=='Human']

# ---> Data preprocessing for epitope info.
# @ DB entries w/ and w/out epitope information.
tcr.info[, epitope.info.provided:=FALSE]
tcr.info[
    # !is.na(Antigen.protein) | !is.na(Protein.ID), # In the end, I realized that you cannot unambiguously find a protein on UniProt based on the protein description/name, only based on  the protein ID.
    !is.na(Protein.ID),
    epitope.info.provided:=TRUE
]
tmp.ggplot.1 <- ggplot(data=tcr.info, aes(x=epitope.info.provided)) +
    geom_bar(width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Number of DB entries w/ and w/out epitope information', x='Is at least one piece of epitope info provided?', y='Count') +
    theme_bw()
# Filter out entries which lack epitope info.
tcr.info <- tcr.info[!is.na(Protein.ID)]

# @ Imputation of antigen organism origin info.
# Create a list of PIDs | protein names to retrieve further info from UniProt manually.
tmp.data <- tcr.info[, .(protein.access=Protein.ID)]
tmp.data <- unique(tmp.data)
tmp.file.name <- paste0(mcpas.data.path, '/AgProts_MCPas.csv')
fwrite(file=tmp.file.name, na=NA, x=tmp.data, col.names=FALSE)
# Load protein info retrieved manually from UniProt. Preprocess it as well.
this.file.name <- paste0(mcpas.data.path, '/AgProts_UniProt.tsv')
uniprot.info <- fread(file=this.file.name)
col.names <- c(
    Protein.ID='Entry',
    uniprot.prot.name='Protein names',
    uniprot.gene.name='Gene Names',
    uniprot.organism='Organism'
)
uniprot.info <- uniprot.info[, ..col.names]
colnames(uniprot.info) <- names(col.names)
uniprot.info[,
    uniprot.gene.name:=str_replace(
        string=str_replace_all(
            string=uniprot.gene.name,
            pattern='\\s+', replacement=';'
        ),
        pattern=';;(;)*', replacement=''
    )
]
# Merge w/ rest of database info.
tmp.check <- tcr.info[, .SD[Protein.ID %in% uniprot.info[, Protein.ID], .N]/.N]
if(tmp.check < 0.99) stop('Unexpected error. We have to capture protein info for majority of the DB entries to be able to continue with the process.\n')
tcr.info <- merge(x=tcr.info, y=uniprot.info, by='Protein.ID', all=FALSE, sort=FALSE)

# @ Assay method description.
# Keep only entries of interest.
tmp.vals <- c(
    `1`='pMHC tetramers',
    `2.1`='In vitro stimulation with a peptide',
    `2.2`='In vitro stimulation with a protein',
    `2.3`='In vitro stimulation with a pathogen',
    `2.4`='In vitro stimulation with tumor cells'
) # 3.0 option was excluded because we believe this type of isolation isn't reliable.
tcr.info[, Antigen.identification.method:=as.character(Antigen.identification.method)]
tcr.info <- tcr.info[Antigen.identification.method %in% names(tmp.vals), ]
tcr.info[, Antigen.identification.method:=tmp.vals[Antigen.identification.method]]

# ---> Data exploration.

# @ Number of epitope/TCR pair entries for each type of chain pair
tmp.data <- tcr.info[,
    .(
        alpha.full=ifelse(
            test=!is.na(CDR3.alpha.aa),
            yes='TRA', no=NA
        ),
        beta.full=ifelse(
            test=!is.na(CDR3.beta.aa),
            yes='TRB', no=NA
        )
    )
]
tmp.data <- tmp.data[, .N, by=.(alpha.full, beta.full)]
tmp.data[, tcr.type:=paste(alpha.full, beta.full, sep='/')]
tmp.data$tcr.type <- factor(x=tmp.data$tcr.type, levels=tmp.data$tcr.type)
setorderv(x=tmp.data, cols='tcr.type')
tmp.ggplot.2 <- ggplot(data=tmp.data, aes(x=tcr.type, y=N)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Number of DB entries for each type of chain pair', x='Chain pair type', y='Count', caption='DB entries w/out epitope info were discarded here.') +
    theme_bw()

# @ Frequency of DB entries and unique epitopes recognized by unique alpha/beta clonotypes.
# A clonotype is defined as a set of alpha and beta CDR3 and TRAV, TRBV, TRAJ and TRBJ sequences.
tmp.data <- tcr.info[
    !is.na(CDR3.alpha.aa) | !is.na(CDR3.beta.aa),
    .(
        alpha.chain=CDR3.alpha.aa,
        beta.chain=CDR3.beta.aa,
        TRAV=TRAV,
        TRAJ=TRAJ,
        TRBV=TRBV,
        TRBJ=TRBJ,
        epitope=Epitope.peptide
    )
]
tmp.data <- unique(tmp.data)
# Remove V and J gene allele information.
tmp.data[, TRAV:=str_replace(string=TRAV, pattern='\\*.+$', replacement='')]
tmp.data[, TRAJ:=str_replace(string=TRAJ, pattern='\\*.+$', replacement='')]
tmp.data[, TRBV:=str_replace(string=TRBV, pattern='\\*.+$', replacement='')]
tmp.data[, TRBJ:=str_replace(string=TRBV, pattern='\\*.+$', replacement='')]
# Identify total number of TCRs that react to each possible number of different epitopes.
tmp.data <- tmp.data[,
    .(entry.count=.N, epitope.count=uniqueN(epitope)),
    by=.(alpha.chain, beta.chain, TRAV, TRAJ, TRBV, TRBJ)
]
tmp.data[entry.count>10, entry.count:=10]
tmp.data[epitope.count>10, epitope.count:=10]
tmp.ggplot.3 <- ggplot(data=tmp.data, aes(x=entry.count)) +
    geom_histogram(bins=3, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Distribution of frequency of DB entries for unique clonotypes', x='Number of DB entries', y='Count', caption='DB entries w/out epitope or TCR chain info were discarded here.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))
tmp.ggplot.4 <- ggplot(data=tmp.data, aes(x=epitope.count)) +
    geom_histogram(bins=10, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Distribution of frequency of epitopes recognized by unique clonotypes', x='Number of recognized epitopes', y='Count', caption='DB entries w/out epitope or TCR chain info were discarded here.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))

# @ Distribution of frequency of epitopes per organism.
tmp.data <- tcr.info[
    !is.na(CDR3.alpha.aa) | !is.na(CDR3.beta.aa),
    .(epitope.name=Epitope.peptide, epitope.gene=`uniprot.gene.name`, organism=`uniprot.organism`)
]
tmp.data <- tmp.data[!is.na(epitope.name) & !is.na(organism), .(N=uniqueN(epitope.name)), by=.(organism)]
tmp.ggplot.5 <- ggplot(data=tmp.data, aes(x=N)) +
    geom_histogram(bins=10, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Distribution of frequency of epitopes per organism', x='Number of epitopes per organism', y='Count', caption='DB entries w/out epitope or TCR chain info were discarded here.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))

# @ Top 10 organisms w/ greatest number of epitopes assessed.
tmp.data <- tcr.info[
    !is.na(CDR3.alpha.aa) | !is.na(CDR3.beta.aa),
    .(epitope.name=Epitope.peptide, epitope.gene=`uniprot.gene.name`, organism=`uniprot.organism`)
]
tmp.data <- tmp.data[!is.na(epitope.name) & !is.na(organism), .(N=uniqueN(epitope.name)), by=.(organism)]
setorderv(x=tmp.data, cols='N', order=-1)
tmp.data <- tmp.data[1:10]
tmp.ggplot.6 <- ggplot(data=tmp.data, aes(x=organism, y=N)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Top 10 organisms w/ greatest number of epitopes assessed', x='Organism', y='Count', caption='DB entries w/out epitope or TCR chain info were discarded here.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))

# @ Top 10 organisms w/ greatest number of DB entries
tmp.data <- tcr.info[
    !is.na(CDR3.alpha.aa) | !is.na(CDR3.beta.aa),
    .(epitope.name=Epitope.peptide, epitope.gene=`uniprot.gene.name`, organism=`uniprot.organism`)
]
tmp.data <- tmp.data[!is.na(epitope.name) & !is.na(organism), .(N=.N), by=.(organism)]
setorderv(x=tmp.data, cols='N', order=-1)
tmp.data <- tmp.data[1:10]
tmp.ggplot.7 <- ggplot(data=tmp.data, aes(x=organism, y=N)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Top 10 organisms w/ greatest number of DB entries', x='Organism', y='Count', caption='DB entries w/out epitope or TCR chain info were discarded here.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))

# @ Distribution of frequency of DB entries per epitope.
tmp.data <- tcr.info[
    !is.na(CDR3.alpha.aa) | !is.na(CDR3.beta.aa),
    .(epitope.name=Epitope.peptide, epitope.gene=`uniprot.gene.name`, organism=`uniprot.organism`)
]
tmp.data <- tmp.data[!is.na(epitope.name) & !is.na(organism), .(N=.N), by=.(epitope.name, organism)]
setorderv(x=tmp.data, cols='N', order=-1)
tmp.ggplot.8 <- ggplot(data=tmp.data, aes(x=N)) +
    geom_histogram(bins=10, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Distribution of frequency of epitopes per organism', x='Number of epitopes per organism', y='Count', caption='DB entries w/out epitope or TCR chain info were discarded here.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))

# @ Top 10 epitopes w/ greatest number of DB entries.
tmp.data <- tcr.info[
    !is.na(CDR3.alpha.aa) | !is.na(CDR3.beta.aa),
    .(epitope.name=Epitope.peptide, epitope.gene=`uniprot.gene.name`, organism=`uniprot.organism`)
]
tmp.data <- tmp.data[!is.na(epitope.name) & !is.na(organism), .(N=.N), by=.(epitope.name, organism)]
tmp.data[, full.epitope:=paste(organism, epitope.name, sep=' / ')]
setorderv(x=tmp.data, cols='N', order=-1)
tmp.data <- tmp.data[1:10]
tmp.ggplot.9 <- ggplot(data=tmp.data, aes(x=full.epitope, y=N)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Top 10 epitopes w/ greatest number of DB entries', x='Epitope', y='DB entries', caption='TCRs that were not alpha/beta were discarded here\nNumber of epitopes was capped at 10.') +
    theme_bw() + theme(axis.text.x=element_text(angle=45))

# ---> Retrieve info required for compilation.
# @ Pairs of TCRs and epitopes, along with their corresponding relevant information.
tcr.clean <- tcr.info[
    !is.na(CDR3.alpha.aa) | !is.na(CDR3.beta.aa),
    .(
        alpha.chain=CDR3.alpha.aa,
        beta.chain=CDR3.beta.aa,
        trav.full=TRAV,
        traj.full=TRAJ,
        trbv.full=TRBV,
        trbj.full=TRBJ,
        epitope.seq=Epitope.peptide,
        antigen=uniprot.gene.name,
        epitope.organism=uniprot.organism,
        assay.method=Antigen.identification.method,
        mhc=MHC,
        pmid=PubMed.ID,
        pathology=Pathology,
        pathology.id=Pathology.Mesh.ID
    )
]

# @ Visualization of the different number of DB entries we have per assay method.
tmp.data <- tcr.clean[, .(assay.method)]
tmp.data[is.na(assay.method), assay.method:='None']
# tmp.data[str_detect(string=assay.method, pattern=';'), assay.method:='Multiple']
tmp.data.2 <- tmp.data[, .N, by=assay.method]
setorderv(x=tmp.data.2, cols='N', order=-1)
tmp.data$assay.method <- factor(x=tmp.data$assay.method, levels=tmp.data.2$assay.method)
tmp.ggplot.10 <- ggplot(data=tmp.data, aes(x=assay.method)) +
    geom_bar(width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Number of DB entries per assay method', x='Assay method', y='DB entries', caption='DB entries w/out epitope or TCR chain info were discarded here.') +
    theme_bw() + theme(axis.text.x=element_text(angle=45))

# @ Reference information.
# Visualization.
tmp.data <- tcr.clean[, .(pmid)]
tmp.data[is.na(pmid), pmid:='None']
tmp.data.2 <- tmp.data[, .N, by=pmid]
tmp.data.2[N>20, N:=20]
tmp.ggplot.11 <- ggplot(data=tmp.data.2, aes(x=N)) +
    geom_histogram(bins=10, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Distribution of the number of DB entries per PMID', x='Number of DB entries', y='Frequency', caption='TCRs that were not alpha/beta were discarded here\nNumber of DB entries was capped at 20.') +
    theme_bw() + theme(axis.text.x=element_text(angle=45))

# @ Clean MHC allele information.
# Specific changes.
tmp.vals <- c('DR3*02:02', 'DRB1*04-01', 'DRB1*15:03')
tcr.clean[mhc%in%tmp.vals, mhc:=paste0('HLA-', mhc)]
# Specific changes for alleles HLA-A1 and -A2, which will be taken as HLA-A*01 and -A*02, respectively.
tcr.clean[mhc=='HLA-A1', mhc:='HLA-A*01']
tcr.clean[mhc=='HLA-A2', mhc:='HLA-A*02']
# Basic gene information if provided.
tcr.clean[,
    mhc.gene:=str_replace(
        string=str_extract(string=mhc, pattern='HLA-[A-Z0-9]+\\**'),
        pattern='\\*', replacement=''
    )
]
tmp.check <- tcr.clean[!is.na(mhc.gene) & !mhc.gene %in% hla.allele.info[, hla.base.gene], .N]/tcr.clean[!is.na(mhc), .N]
if(tmp.check>0.05) warning('We failed to validate the HLA gene name for a fairly great amount of DB entries with HLA information.\n')
# Basic allele information if provided.
tcr.clean[
    !is.na(mhc.gene),
    mhc.allele:=str_replace(
        string=str_extract(string=mhc, pattern='HLA-[A-Z0-9]+\\*(\\d{2}:*){2}'),
        pattern=':$', replacement=''
    )
]
# Sanity check. Confirm validity of HLA alleles.
tmp.check <- tcr.clean[!is.na(mhc.allele) & !mhc.allele %in% hla.allele.info[, hla.base.allele], .N==0]
if(!tmp.check){
    warning('Some of the identified alleles in the DB cannot be confirmed based on IMGT.\n')
    # Discard MHC allele info because it doesn't seem fair.
    tcr.clean[
        !is.na(mhc.allele) & !mhc.allele %in% hla.allele.info[, hla.base.allele],
        mhc.allele:=NA
    ]
}

# @ Visualization of the different number of DB entries we have per HLA gene.
tmp.data <- tcr.clean[, .(mhc.gene)]
tmp.data[is.na(mhc.gene), mhc.gene:='None']
tmp.data.2 <- tmp.data[, .N, by=mhc.gene]
setorderv(x=tmp.data.2, cols='N', order=-1)
tmp.data$mhc.gene <- factor(x=tmp.data$mhc.gene, levels=tmp.data.2$mhc.gene)
tmp.ggplot.12 <- ggplot(data=tmp.data, aes(x=mhc.gene)) +
    geom_bar(width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Number of DB entries per HLA gene', x='HLA gene', y='DB entries', caption='TCRs that were not alpha/beta were discarded here\nNumber of epitopes was capped at 10.') +
    theme_bw() + theme(axis.text.x=element_text(angle=45))

# @ Clean VDJ gene usage info.
tcr.clean <- clean.vdj.info(tcr.data=tcr.clean, data.base='MCPas')

# ---> Output all plots.
tmp.file.name <- paste0(reports.path, '/BasicExploration_MCPas.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot.1)
print(tmp.ggplot.2)
print(tmp.ggplot.3)
print(tmp.ggplot.4)
print(tmp.ggplot.5)
print(tmp.ggplot.6)
print(tmp.ggplot.7)
print(tmp.ggplot.8)
print(tmp.ggplot.9)
print(tmp.ggplot.10)
print(tmp.ggplot.11)
print(tmp.ggplot.12)
dev.off()

# ---> Final format of the table.
cols.order <- c(
    'alpha.chain', 'beta.chain',
    'trav.full', 'traj.full', 'trbv.full', 'trbj.full',
    'trav.base', 'traj.base', 'trbv.base', 'trbj.base',
    'epitope.organism', 'antigen', 'epitope.seq',
    'assay.method',
    'mhc', 'mhc.gene', 'mhc.allele',
    'pmid',
    'pathology', 'pathology.id'
)
tcr.clean <- tcr.clean[, ..cols.order]
setorderv(x=tcr.clean, cols='assay.method')
db.info[['mcpas']] <- copy(tcr.clean)


cat('\n\n')
### ------------------- Reference harmonization -------------------- ###
cat('### ------------------- Reference harmonization -------------------- ###\n')

# Combine data from multiple DBs.
full.ref <- rbindlist(l=db.info, use.names=TRUE, fill=TRUE, idcol='data.base')
# Remove entries from IEDB that do not contain NA values.
full.ref[epitope.organism=='', epitope.organism:=NA]
full.ref <- full.ref[!(is.na(alpha.chain) & is.na(beta.chain)) & !(is.na(epitope.organism) | is.na(epitope.seq))]

# ---> Organism consensus values.
full.ref[, unique(epitope.organism)]
full.ref[, uniqueN(data.base), by=epitope.organism]
full.ref[, consensus.organism:='']

# @ Function to define consensus.
set.consens.org <- function(these.data, org.pttn, org.consens=NULL, do.apply=TRUE){
    tmp.vals <- these.data[
        consensus.organism=='' & str_detect(string=epitope.organism, pattern=org.pttn),
        .(N=.N),
        by=.(epitope.organism, data.base)
    ]
    if(!do.apply){
        tmp.vals <- apply(X=tmp.vals, MARGIN=1, FUN=function(x) paste0(x, collapse='\t'))
        cat(paste0(tmp.vals, collapse='\n'), '\n')
        return(NULL)
    }else{
        if(is.null(org.consens)) stop('Unexpected organism consensus value.\n')
        tmp.vals <- tmp.vals[, unique(epitope.organism)]
        these.data[
            epitope.organism %in% tmp.vals,
            consensus.organism:=org.consens
        ]
        return(these.data)
    }
}

# @ -- Eukaryotes
# Human
tmp.pttn <- 'Homo'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='Homo sapiens')
# Mouse
tmp.pttn <- 'musculus'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='Mus musculus')
# Rat
tmp.pttn <- 'Rattus norvegicus'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='Rattus norvegicus')
# Yeast
tmp.pttn <- 'Saccharomyces'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='Saccharomyces cerevisiae')
# Wheat
tmp.pttn <- '^Triticum|Wheat'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='Triticum aestivum')
# Selaginella moellendorffii
tmp.pttn <- '^Selaginella'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='Selaginella moellendorffii')

# @ -- Herpesviridae
# HSV
tmp.pttn <- 'Human herpesvirus 1|Human herpesvirus 2|HSV-2'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='HSV')
# VZV (varicella-zoster virus)
tmp.pttn <- 'Human herpesvirus 3'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='VZV')
# EBV
tmp.pttn <- 'herpesvirus 4|EBV'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='EBV')
# CMV
tmp.pttn <- 'Cytomegalovirus|CMV|herpesvirus 5'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='CMV')
# Non-human herpesviruses
tmp.pttn <- 'Retroperitoneal fibromatosis-associated herpesvirus|Murid betaherpesvirus 1 \\(Mouse cytomegalovirus 1\\)'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='Non-human herpesviruses')

# @ -- Respiratory viruses (including orthomyxoviridae, pneumoviridae, coronaviridae)
# IAV
tmp.pttn <- 'Influenza A|InfluenzaA|H1N1 subtype \\(H1N1\\)'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='IAV')
# IBV
tmp.pttn <- 'Influenza B virus'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='IBV')
# RSV
tmp.pttn <- 'Human respiratory syncytial virus'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='RSV')
# MERS-CoV
tmp.pttn <- 'MERS coronavirus'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='MERS-CoV')
# SARS-CoV
tmp.pttn <- 'SARS coronavirus BJ01|SARS-CoV1|SARS coronavirus Urbani'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='SARS-CoV')
# SARS-CoV-2
tmp.pttn <- 'CoV2|CoV-2|Severe acute respiratory syndrome coronavirus 2 Wuhan/Hu-1/2019'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='SARS-CoV-2')
# CCC. See: https://www.cdc.gov/coronavirus/general-information.html
tmp.pttn <- '229E|NL63|HKU1|OC43|CCC'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='CCC')
full.ref[epitope.organism=='NL63-related bat coronavirus', consensus.organism:=''] # Correction after inexact pattern for indicated entry.
# Non-human coronavirus
tmp.pttn <- 'NL63-related bat coronavirus|SARS coronavirus Tor2|^Chaerephon bat coronavirus'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='Non-human coronavirus')

# @ -- Flaviviridae
# DENV
tmp.pttn <- '^DENV|^Dengue|^dengue'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='DENV')
# YFV
tmp.pttn <- '^Yellow fever virus|YFV'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='YFV')

# @ -- Hepadnaviridae 
# HBV
tmp.pttn <- 'HBV|Hepatitis B'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='HBV')
# HCV
tmp.pttn <- 'HCV|Hepatitis C'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='HCV')
# HEV
tmp.pttn <- 'HEV|Hepatitis E'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='HEV')

# @ -- Retroviridae, papillomaviridae and polyomaviridae
# Human T-lymphotropic virus 1 (HTLV-I)
tmp.pttn <- 'Human T-cell leukemia virus 1|HTLV-1|Human T-cell leukemia virus type I'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='HTLV-1')
# HIV
tmp.pttn <- 'HIV|Human immunodeficiency virus 1'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='HIV-1')
# HPV
tmp.pttn <- 'HPV|papillomavirus'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='HPV')
# MCV
tmp.pttn <- 'Merkel cell polyomavirus|MCPyV'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='MCPyV')
# Human polyomavirus 2
full.ref[epitope.organism=='JC polyomavirus (Human polyomavirus (type JC))', consensus.organism:='Human polyomavirus 2']

# @ -- Bacteria
# Mtb
tmp.pttn <- 'tuberculosis'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='Mycobacterium tuberculosis')
# E. coli
tmp.pttn <- 'Escherichia coli|E.Coli'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='Escherichia coli')
# Salmonella - 6 different species, see https://en.wikipedia.org/wiki/Salmonella_enterica
tmp.pttn <- 'Salmonella enterica'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='Salmonella enterica')
# Pseudomonas aeruginosa
tmp.pttn <- 'Pseudomonas aeruginosa|PseudomonasAeruginosa'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='Pseudomonas aeruginosa')
# Pseudomonas fluorescens
tmp.pttn <- 'Pseudomonas fluorescens|PseudomonasFluorescens'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='Pseudomonas fluorescens')
# Akkermansia muciniphila
tmp.pttn <- 'Akkermansia muciniphila'
set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, do.apply=FALSE)
full.ref <- set.consens.org(these.data=full.ref, org.pttn=tmp.pttn, org.consens='Akkermansia muciniphila')
# Streptomyces kanamyceticus
full.ref[epitope.organism=='StreptomycesKanamyceticus', consensus.organism:='Streptomyces kanamyceticus']

# @ -- Other unclear organisms.
full.ref[epitope.organism=='AKR (endogenous) murine leukemia virus (AKR murine leukemia virus)', consensus.organism:='AKR murine leukemia virus']
full.ref[epitope.organism=='Chaerephon bat coronavirus/Kenya/KY22/2006 (Bat coronavirus BtKY22/Chaereph
on sp./Kenya/2006)', consensus.organism:='Bat coronavirus']

# @ -- General changes.
# Remove info within brackets for remaining entries w/out consensus info.
# First explore all possible values to confirm that it is ok to continue with this process (ok in the sense that the remaining string will make sense and is unambiguous compared with the consensus values we've already called)
# full.ref[consensus.organism=='' & str_detect(string=epitope.organism, pattern='\\(|\\)'), sort(unique(epitope.organism))]
full.ref[
    consensus.organism=='' & str_detect(string=epitope.organism, pattern='\\(|\\)'),
    consensus.organism:=str_replace(
        string=epitope.organism,
        pattern='\\s+\\([^\\(\\)]+\\)',
        replacement=''
    )
]
full.ref[
    str_detect(string=consensus.organism, pattern='\\(|\\)'),
    consensus.organism:=str_replace(
        string=consensus.organism,
        pattern='\\s+\\([^\\(\\)]+\\)',
        replacement=''
    )
]

# @ -- General changes.
# All remaining values look good enough.
# full.ref[consensus.organism=='', uniqueN(epitope.organism)]
# full.ref[consensus.organism=='', sort(unique(epitope.organism))]
full.ref[consensus.organism=='', consensus.organism:=epitope.organism]

# @ -- Manual inspection
# full.ref[, uniqueN(consensus.organism), by=epitope.organism][V1>1] # Sanity check.
tmp.data <- unique(full.ref[, .(epitope.organism, consensus.organism)])
setorderv(x=tmp.data, cols='consensus.organism')
tmp.file.name <- paste0(reports.path, '/ConsensusOrganismInspection.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA)

# ---> Unique identifiers.
# Unique identifiers for TCRs.
tmp.cols <- c('beta.chain', 'alpha.chain', 'trbv.base', 'trbj.base', 'trav.base', 'traj.base')
tmp.data <- full.ref[, ..tmp.cols]
tmp.data <- unique(tmp.data)
setorderv(x=tmp.data, cols=tmp.cols)
tmp.data[, tcr.id:=paste0(
    paste0('000000', collapse=''),
    1:.N
)]
tmp.data[, tcr.id:=str_extract(string=tcr.id, pattern='\\d{6}$')]
tmp.data[, tcr.id:=paste0('RTCR', tcr.id)]
full.ref <- merge(x=full.ref, y=tmp.data, by=tmp.cols, all=TRUE)

# Unique identifiers for epitopes.
tmp.cols <- c('epitope.organism', 'epitope.seq')
tmp.data <- full.ref[, ..tmp.cols]
tmp.data <- unique(tmp.data)
setorderv(x=tmp.data, cols=tmp.cols)
tmp.data[, epitope.id:=paste0(
    paste0('0000', collapse=''),
    1:.N
)]
tmp.data[, epitope.id:=str_extract(string=epitope.id, pattern='\\d{4}$')]
tmp.data[, epitope.id:=paste0('RE', epitope.id)]
full.ref <- merge(x=full.ref, y=tmp.data, by=tmp.cols, all=TRUE)

# ---> Final TCR specificity at the organism level
tmp.data <- full.ref[,
    .(
        specificity.ags.tag=paste0(
            unique(sort(consensus.organism)),
            collapse=';'
        )
    ),
    by=.(tcr.id)
]
tmp.data[,
    specificity.class.tag:=specificity.ags.tag
]
# Explore types of ag combinations for cross-reactive TCRs.
tmp.data.2 <- tmp.data[
    str_detect(string=specificity.class.tag, pattern=';'),
    .N,
    by=specificity.class.tag
]
setorderv(x=tmp.data.2, cols='N', order=-1)
# tmp.data.2[N>5]
tmp.file.name <- paste0(reports.path, '/IndAgsForCross-reactiveStuff.csv')
fwrite(file=tmp.file.name, x=tmp.data.2, na=NA, quote=TRUE)
# Define cross-reactive specificity class and merge w/ rest of info.
tmp.data[
    str_detect(string=specificity.class.tag, pattern=';'),
    specificity.class.tag:='Cross-reactive'
]
full.ref <- merge(x=full.ref, y=tmp.data, by='tcr.id')

# ---> Define experimental category based on assay type.
# Below are the definition of the possible categories.
#       1. Weak support only: everything that's NOT supported by tetramer-like sorting, cloning by limiting dilution or X-ray crystallography.
#       2. Strong support only: everything that IS supported by tetramer-like sorting, cloning by limiting dilution or X-ray crystallography.

# Manual inspection of all possible assays.
tmp.data <- full.ref[, .N, by=.(assay.method)]
setorderv(x=tmp.data, cols='N', order=-1)
# tmp.data.2[N>5]
tmp.file.name <- paste0(reports.path, '/MethodInspection.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)
# Assignment of a simplified assay description based on manual inspection.
tmp.data[
    str_detect(
        string=assay.method,
        pattern='tetramer|pentamer|streptamer|dextramer|multimer|pelimer'
    ),
    simp.assay.method:='Tetramer-like sorting'
]
tmp.data[
    is.na(simp.assay.method) & (
        assay.method %like% 'limiting-dilution-cloning' | assay.method %like% 'cloning-by-limiting-dilution' | assay.method %like% 'limited-dilution-cloning' | assay.method %like% 'limiting-diffusion-cloning'
    ),
    simp.assay.method:='Cloning by limiting dilution'
]
tmp.data[
    assay.method=='x-ray crystallography',
    simp.assay.method:='X-ray crystallography'
]
tmp.data[
    is.na(simp.assay.method),
    simp.assay.method:='Other or unspecified'
]
# Assign experimental category (see definition right above in this same code chunk).
tmp.data[, exp.category:=1]
tmp.vals <- c('Tetramer-like sorting', 'Cloning by limiting dilution', 'X-ray crystallography')
tmp.data[simp.assay.method %in% tmp.vals, exp.category:=2]
# Update table for further manual inspection.
tmp.file.name <- paste0(reports.path, '/MethodInspection.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)
# Merge w/ rest of info.
tmp.data <- tmp.data[, .(exp.category, simp.assay.method, assay.method)]
full.ref <- merge(x=full.ref, y=tmp.data, by='assay.method')


# ---> Exploration of all results.

# @ TCR overlap between databases.
tmp.vals <- full.ref[, unique(data.base)]
tmp.data <- lapply(X=tmp.vals, FUN=function(tmp.val){
    tmp.data <- unique(full.ref[data.base==tmp.val, unique(tcr.id)])
    return(tmp.data)
})
names(tmp.data) <- tmp.vals
tmp.cols <- c(`MCPas`='#594FE8', `VDJdb`='#E8884F', `IEDB`='#5A6958')
# As venn diagram
tmp.file.name <- paste0(reports.path, '/TCROvlp.tiff')
venn.diagram(filename=tmp.file.name, x=tmp.data, fill=tmp.cols)
# As UpSet plot.
tmp.ggplot.1 <- upset(data=fromList(tmp.data), main.bar.color='lightblue')

# @ TCR/epitope pair overlap between databases.
tmp.vals <- full.ref[, unique(data.base)]
tmp.data <- lapply(X=tmp.vals, FUN=function(tmp.val){
    tmp.data <- unique(full.ref[data.base==tmp.val, paste(tcr.id, epitope.id, sep='-')])
    return(tmp.data)
})
names(tmp.data) <- tmp.vals
tmp.cols <- c(`MCPas`='#594FE8', `VDJdb`='#E8884F', `IEDB`='#5A6958')
# As venn diagram
tmp.file.name <- paste0(reports.path, '/TCREpitopePairOvlp.tiff')
venn.diagram(filename=tmp.file.name, x=tmp.data, fill=tmp.cols)
# As UpSet plot.
tmp.ggplot.2 <- upset(data=fromList(tmp.data), main.bar.color='lightblue')

# @ Remove venn diagram creation process residuals.
tmp.files <- list.files(path=reports.path, pattern='.+tiff.+\\.log$', full.names=TRUE)
file.remove(tmp.files)

# @ Number of epitope/TCR pair entries for each type of chain pair
tmp.data <- full.ref[,
    .(
        alpha.full=ifelse(
            test=!is.na(alpha.chain),
            yes='TRA', no=NA
        ),
        beta.full=ifelse(
            test=!is.na(beta.chain),
            yes='TRB', no=NA
        )
    )
]
tmp.data <- tmp.data[, .N, by=.(alpha.full, beta.full)]
tmp.data[, tcr.type:=paste(alpha.full, beta.full, sep='/')]
setorderv(x=tmp.data, cols='tcr.type')
tmp.data$tcr.type <- factor(x=tmp.data$tcr.type, levels=tmp.data$tcr.type)
tmp.ggplot.3 <- ggplot(data=tmp.data, aes(x=tcr.type, y=N)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Number of DB entries for each type of chain pair', x='Chain pair type', y='Count') +
    theme_bw()

# @ Frequency of DB entries and unique epitopes recognized by unique alpha/beta clonotypes.
tmp.data <- full.ref[,
    .(entry.count=.N, epitope.count=uniqueN(epitope.id)),
    by=.(tcr.id)
]
tmp.data[entry.count>5, entry.count:=5]
tmp.data[epitope.count>5, epitope.count:=5]
tmp.ggplot.4 <- ggplot(data=tmp.data, aes(x=entry.count)) +
    geom_histogram(bins=5, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Frequency of DB entries for unique clonotypes', x='Number of DB entries', y='Count', caption='Number of entries capped at 5.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))
tmp.ggplot.5 <- ggplot(data=tmp.data, aes(x=epitope.count)) +
    geom_histogram(bins=5, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='frequency of epitopes recognized by unique clonotypes', x='Number of recognized epitopes', y='Count', caption='Number of epitopes capped at 5.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))

# @ Distribution of frequency of epitopes per organism.
tmp.data <- full.ref[, .(N=uniqueN(epitope.seq)), by=.(consensus.organism)]
tmp.data[N>10, N:=10]
tmp.ggplot.6 <- ggplot(data=tmp.data, aes(x=N)) +
    geom_histogram(bins=10, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Frequency of epitopes per organism', x='Number of epitopes per organism', y='Count', caption='Number of epitopes capped at 10.') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))

# @ Top 10 organisms w/ greatest number of epitopes assessed.
tmp.data <- full.ref[, .(N=uniqueN(epitope.seq)), by=.(consensus.organism)]
setorderv(x=tmp.data, cols='N', order=-1)
tmp.data <- tmp.data[1:10]
tmp.data$consensus.organism <- factor(x=tmp.data$consensus.organism, levels=tmp.data$consensus.organism)
tmp.ggplot.7 <- ggplot(data=tmp.data, aes(x=consensus.organism, y=N)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Top 10 organisms w/ greatest number of epitopes assessed', x='Organism', y='Count') +
    theme_bw() + theme(axis.text.x=element_text(angle=90))

# @ Top 10 organisms w/ greatest number of DB entries.
tmp.data <- full.ref[, .(N=.N), by=.(consensus.organism)]
# tmp.data[, full.epitope:=paste(organism, epitope.name, sep=' / ')]
setorderv(x=tmp.data, cols='N', order=-1)
tmp.data <- tmp.data[1:10]
tmp.data$consensus.organism <- factor(x=tmp.data$consensus.organism, levels=tmp.data$consensus.organism)
tmp.ggplot.8 <- ggplot(data=tmp.data, aes(x=consensus.organism, y=N)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Top 10 organisms w/ greatest number of DB entries', x='Organism', y='Number of DB entries') +
    theme_bw() + theme(axis.text.x=element_text(angle=45))

# @ Top 10 ag specificities w/ greatest number of DB entries.
tmp.data <- full.ref[, .(N=.N), by=.(specificity.class.tag)]
# tmp.data[, full.epitope:=paste(organism, epitope.name, sep=' / ')]
setorderv(x=tmp.data, cols='N', order=-1)
tmp.data <- tmp.data[1:10]
tmp.data$specificity.class.tag <- factor(x=tmp.data$specificity.class.tag, levels=tmp.data$specificity.class.tag)
tmp.ggplot.9 <- ggplot(data=tmp.data, aes(x=specificity.class.tag, y=N)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Top 10 ag specificities w/ greatest number of DB entries', x='Antigen specificity', y='Number of DB entries') +
    theme_bw() + theme(axis.text.x=element_text(angle=45))

# @ Top 10 epitopes w/ greatest number of DB entries.
tmp.data <- full.ref[, .(N=.N), by=.(consensus.organism, epitope.seq)]
tmp.data[, full.epitope:=paste(consensus.organism, epitope.seq, sep=' / ')]
setorderv(x=tmp.data, cols='N', order=-1)
tmp.data <- tmp.data[1:10]
tmp.data$full.epitope <- factor(x=tmp.data$full.epitope, levels=tmp.data$full.epitope)
tmp.ggplot.10 <- ggplot(data=tmp.data, aes(x=full.epitope, y=N)) +
    geom_bar(stat='identity', width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Top 10 epitopes w/ greatest number of DB entries', x='Epitope', y='Number of DB entries') +
    theme_bw() + theme(axis.text.x=element_text(angle=45))

# @ MHC restriction.
tmp.data <- full.ref[, .(mhc.gene)]
tmp.data[is.na(mhc.gene), mhc.gene:='None']
tmp.data.2 <- tmp.data[, .N, by=mhc.gene]
setorderv(x=tmp.data.2, cols='N', order=-1)
tmp.data$mhc.gene <- factor(x=tmp.data$mhc.gene, levels=tmp.data.2$mhc.gene)
tmp.ggplot.11 <- ggplot(data=tmp.data, aes(x=mhc.gene)) +
    geom_bar(width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Number of DB entries per HLA gene', x='HLA gene', y='Number of DB entries') +
    theme_bw() + theme(axis.text.x=element_text(angle=45))

# @ Methodology, simplified methods
tmp.data <- full.ref[, .(simp.assay.method)]
tmp.data.2 <- tmp.data[, .N, by=simp.assay.method]
setorderv(x=tmp.data.2, cols='N', order=-1)
tmp.data$simp.assay.method <- factor(x=tmp.data$simp.assay.method, levels=tmp.data.2$simp.assay.method)
tmp.ggplot.12 <- ggplot(data=tmp.data, aes(x=simp.assay.method)) +
    geom_bar(width=0.7, fill='lightblue', color='black') +
    scale_y_continuous(expand=c(0, NA)) +
    labs(title='Number of DB entries per assay', x='Assay', y='Number of DB entries') +
    theme_bw() + theme(axis.text.x=element_text(angle=45))

# @ Output all plots.
tmp.file.name <- paste0(reports.path, '/BasicExploration_FullRef.pdf')
pdf(file=tmp.file.name)
print(tmp.ggplot.1)
print(tmp.ggplot.2)
print(tmp.ggplot.3)
print(tmp.ggplot.4)
print(tmp.ggplot.5)
print(tmp.ggplot.6)
print(tmp.ggplot.7)
print(tmp.ggplot.8)
print(tmp.ggplot.9)
print(tmp.ggplot.10)
print(tmp.ggplot.11)
print(tmp.ggplot.12)
dev.off()

# ---> Report final reference.
# Column order.
tmp.cols <- c(
    'data.base', 'tcr.id',
    'beta.chain', 'alpha.chain', 'trbv.base', 'trbj.base', 'trav.base', 'traj.base', 'trbv.full', 'trbj.full', 'trav.full', 'traj.full',
    'specificity.class.tag', 'specificity.ags.tag', 'consensus.organism',
    'epitope.id', 'epitope.seq', 'epitope.organism', 'antigen',
    'mhc', 'mhc.gene', 'mhc.allele', 'assay.mhcs', 'donor.mhcs',
    'exp.category', 'simp.assay.method', 'assay.method',
    'pmid', 'pathology.id', 'pathology'
)
tmp.check <- length(setdiff(x=colnames(full.ref), y=tmp.cols)) == 0
if(!tmp.check) stop('Something went wrong while sorting columns in final reference file.\n')
full.ref <- full.ref[, ..tmp.cols]
# Full format.
tmp.file.name <- paste0(reports.path, '/AgSpcTCR_CompiledRef.csv')
fwrite(file=tmp.file.name, x=full.ref, na=NA, quote=TRUE)
# Clustering by similarity-ready format.
# Include only entries for well-represented organisms. Definition of "well-represented organisms" is still loose here.
tmp.orgs <- full.ref[, .N, by=consensus.organism][N>100, unique(consensus.organism)]
tmp.data <- full.ref[
    consensus.organism %in% tmp.orgs,
    .(
        clonotype.tag=tcr.id,
        donor.id.tag='REF-1',
        cdr3a.aa.seq=alpha.chain,
        cdr3b.aa.seq=beta.chain,
        cdr3a.nt.seq=NA, cdr3b.nt.seq=NA,
        tra.v=trav.base,
        tra.j=traj.base,
        trb.v=trbv.base,
        trb.j=trbj.base,
        specificity.class.tag,
        exp.category
    )
]
tmp.file.name <- paste0(reports.path, '/AgSpcTCR_CompiledRef_ClusteringReadyFormat.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA, quote=TRUE)


cat('\n\n')
### --------------------------------- END -------------------------------- ###
# ---> Print session info.
sessionInfo()

cat('PROGRAM FINISHED!\n\n')