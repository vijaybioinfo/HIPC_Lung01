###########    -----   Get data to start for the    ------    ###########
###########    ------------- BOOST project  --------------    ###########

# Version.
# Full version: 0.1

# Whole version: 0.1
# --> Subversion: 1
#   * First subversion.


### -------------------------- Description -------------------------- ###
# Script to get the basic data to get the final figures for the BOOST project.


### --------------------------- Libraries --------------------------- ###
library(Seurat)
library(stringr)
library(data.table)
library(tidyr)
library(gtools)


### --------------------------- Functions --------------------------- ###
coll.version <- 0.1
gen.data.path <- '/mnt/BioAdHoc/Groups/vd-sette/BOOST/paper_developments/BOOST_1/paper_items'
coll.file <- paste0(gen.data.path, '/jobs_scripts/functions_collection.', coll.version, '.R') # Collection borrowed from the DICE Tissue project.
file.exists(coll.file)
source(coll.file)


### ----------------------- General Arguments ----------------------- ###

# ---> General definitions.
# @ Seed.
set.seed(seed=1)
# @ Cluster labels for each dataset.
obj.extended.names <- c(
  cd4.all='CD4_Full',
  cd8.all='CD8_Full'
)
clust.labs <- c(
  cd4.all='RNA_snn_res.0.4',
  cd8.all='RNA_snn_res.0.1'
)
# ---> Path definitions.
# Module signatures.
# module.sigs.path <- paste0(gen.data.path, '/module_signatures')
# Donor info path
donor.data.path <- '/mnt/BioAdHoc/Groups/vd-sette/BOOST/paper_developments/BOOST_1/donors_metadata'
# ---> File definitions.
# Seurat objects.
obj.files <- c(
  cd4.all='/mnt/BioAdHoc/Groups/vd-sette/BOOST/seurat_analysis/BOOST_CD4/BOOST_Round-4_CD4_Subset-A/BOOST_Round-4_CD4_Subset-A_09-07-2023_qc-std_var-30_pc-30_hto-all_cohort-all_harmony-seq.batch.tag/seurat_objects/SeuratObjectForPrjBOOST_Round-4_CD4_Subset-A_WithArgs_NoPCs_30.RDS',
  cd8.all='/mnt/BioAdHoc/Groups/vd-sette/BOOST/seurat_analysis/BOOST_CD8/BOOST_Round-4_CD8/BOOST_Round-4_CD8_09-06-2023_qc-std_var-30_pc-30_hto-all_cohort-all_harmony-seq.batch.tag/seurat_objects/SeuratObjectForPrjBOOST_Round-4_CD8_WithArgs_NoPCs_30.RDS'
)
# ---> CD4RE status for donors.
cd4re.info.file <- paste0(donor.data.path, '/DonorCD4REClassification.0.1.csv')
# ---> Check directories and files.
if(!all(dir.exists(gen.data.path))) stop(paste0('Following paths must be already defined:\n', gen.data.path, '\n'))
reports.path <- gen.data.path
essential.files <- c(obj.files, cd4re.info.file)
essential.files <- essential.files[!unlist(lapply(X=essential.files, FUN=file.exists))]
if(length(essential.files) > 0) stop(paste0('Next essential files are not appropriately defined -make sure you\'ve got their names right-:\n', paste0(essential.files, collapse='\n'), '\n'))


### ------------------------- Data Loading -------------------------- ###

# ---> Seurat obejcts.
obj.list <- lapply(X=obj.files, FUN=readRDS)

# ---> Donors' CD4RE status
cd4re.info <- fread(file=cd4re.info.file)
cd4re.info[, donor.id:=str_replace(string=donor.id, pattern='\\s+', replacement='.')]


### ---------------------- Data preprocessing ----------------------- ###


### ------------------------ Backup samples ------------------------- ###

# ---> Fetch backup samples from unprocessed object.
tmp.data <- lapply(X=names(obj.list), FUN=function(cell.type){
    tmp.data <- as.data.table(obj.list[[cell.type]]@meta.data)
    tmp.data <- tmp.data[
        !is.na(hashtag.tag) & hashtag.tag!='Negative' & is.na(dose.2.time.point.tag) & is.na(dose.3.time.point.tag) & is.na(dose.4.time.point.tag),
        .(cell.count.backup=.N),
        by=.(ext.donor.id.tag, donor.id.tag, exp.sample.id.tag)
    ]
    return(tmp.data)
})
names(tmp.data) <- names(obj.list)
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, idcol='data.set')
tmp.file.name <- paste0(reports.path, '/BackupSampleInfo.csv')
fwrite(file=tmp.file.name, x=tmp.data, na=NA)

# ---> Determine number of cells for santandard samples that are backupped.
# This is the same as the one above except that it was modified by adding the standard sample info manually.
tmp.file.name <- paste0(donor.data.path, '/BackupSampleInfo_1.csv')
bkup.samples <- fread(file=tmp.file.name)
for(cell.type in names(obj.list)){
    tmp.data.1 <- bkup.samples[data.set==cell.type, standard.sample.id]
    tmp.data.2 <- as.data.table(obj.list[[cell.type]]@meta.data)
    tmp.data.2 <- tmp.data.2[
        !is.na(ext.donor.id.tag) & ext.donor.id.tag %in% tmp.data.1,
        .(cell.count.standard=.N),
        by=ext.donor.id.tag
    ]
    tmp.data.2[, data.set:=cell.type]
    bkup.samples <- merge(x=bkup.samples, y=tmp.data.2, by.x=c('data.set', 'standard.sample.id'), by.y=c('data.set', 'ext.donor.id.tag'), all=TRUE)
}
bkup.samples[!is.na(cell.count.standard.x), cell.count.standard:=cell.count.standard.x]
bkup.samples[!is.na(cell.count.standard.y), cell.count.standard:=cell.count.standard.y]
bkup.samples[, `:=`(cell.count.standard.x=NULL, cell.count.standard.y=NULL)]
# Cell count ratio Ri
bkup.samples[, cell.count.ratio:=cell.count.backup/cell.count.standard]
# Save to make final decision.
col.order <- c('data.set', 'donor.id.tag', 'exp.sample.id.tag', 'standard.sample.id', 'ext.donor.id.tag', 'cell.count.standard', 'cell.count.backup', 'cell.count.ratio', 'dose.cohorts.standard', 'days.from.last.dose.standard', 'days.from.last.dose.backup', 'months.from.last.dose.standard', 'months.from.last.dose.backup')
bkup.samples <- bkup.samples[, ..col.order]
tmp.file.name <- paste0(donor.data.path, '/BackupSampleInfo_2.csv')
fwrite(file=tmp.file.name, x=bkup.samples, na=NA)

# ---> Apply specific criteria to determine how to proceed for each pair of standard and backup samples.
# We defined these criteria in an xlsx file that can be found in the same folder as the csv file (same basename).
tmp.file.name <- paste0(donor.data.path, '/BackupSampleInfo_3.csv')
bkup.samples <- fread(file=tmp.file.name)

# ---> Action per case of standard/backup sample pair.

# @ Combine sample pairs.
for(cell.type in names(obj.list)){
    cat(cell.type, '\n')
    # Pairs to which the rule has to be applied.
    tmp.data.1 <- bkup.samples[
        data.set==cell.type & rule=='combine',
        .(sample.id.backup=ext.donor.id.tag, sample.id.standard=standard.sample.id)
    ]
    # Current metadata
    tmp.data.2 <- as.data.table(obj.list[[cell.type]]@meta.data)
    tmp.data.2[, barcode:=Cells(obj.list[[cell.type]])]
    # Sanity check. Confirm that the constant tags are the same between standard and backup sample pairs for each pair.
    constant.tags <- c(
        'donor.id.tag', 'gender.tag', 'ethnicity.tag', # Donor demographics and ID.
        'peptide.pool.tag', 'geo.location.site.tag',
        'dose.1.type.tag', 'dose.2.type.tag', 'dose.3.type.tag', 'dose.4.type.tag', 'dose.5.type.tag'
    )
    tmp.cols <- c('ext.donor.id.tag', constant.tags)
    tmp.data <- unique(tmp.data.2[, ..tmp.cols])
    to.check <- merge(x=tmp.data.1, y=tmp.data, by.x='sample.id.backup', by.y='ext.donor.id.tag', all.x=TRUE)
    to.check <- merge(x=to.check, y=tmp.data, by.x='sample.id.standard', by.y='ext.donor.id.tag', all.x=TRUE, suffixes=c('|bkup', '|std'))
    to.check <- as.data.table(gather(data=to.check, key='constant.tag', value='value', -`sample.id.standard`, -`sample.id.backup`))
    to.check <- as.data.table(separate(data=to.check, col=constant.tag, into=c('constant.tag', 'sample'), sep='\\|'))
    to.check <- as.data.table(spread(data=to.check, key=sample, value=value))
    to.check <- to.check[, .SD[!is.na(bkup), all(bkup==std)] & .SD[is.na(bkup), all(is.na(std))]]
    if(!to.check) stop('Unexpected error. Values are different for tags that are meant to be constant for a given donor between some backup and standard sample pairs.\n')
    # Proceed with the transfer of information.
    transfer.tags <- c(
        'dose.2.time.point.tag', 'dose.3.time.point.tag', 'dose.4.time.point.tag',
        'cohort.dose.2.flag.tag', 'cohort.dose.3.flag.tag', 'cohort.dose.4.flag.tag'
    )
    for(idx in 1:tmp.data.1[, .N]){
        cat('\t', idx, '\n')
        sample.id.std <- tmp.data.1[idx, sample.id.standard]
        sample.id.bkup <- tmp.data.1[idx, sample.id.backup]
        for(transfer.tag in transfer.tags){
            cat('\t\t', transfer.tag, '\n')
            rep.val <- unique(
                obj.list[[cell.type]]@meta.data[
                    !is.na(obj.list[[cell.type]]@meta.data[, 'ext.donor.id.tag']) & obj.list[[cell.type]]@meta.data[, 'ext.donor.id.tag']==sample.id.std,
                    transfer.tag
                ]
            )
            if(length(rep.val)>1) stop('Unexpected error\n')
            obj.list[[cell.type]]@meta.data[
                !is.na(obj.list[[cell.type]]@meta.data[, 'ext.donor.id.tag']) & obj.list[[cell.type]]@meta.data[, 'ext.donor.id.tag']==sample.id.bkup,
                transfer.tag
            ] <-  rep.val
        }
    }
}

# Perform sample replacement for samples whose backup is a better fit than the current standard.
# First for samples whose standard version is not even part of the seurat object yet.
# Information for standard sample "6673-LC-2" is replaced by the information of backup sample "6673-LC-1" (in reality, this happened due to a mislabeling issue, yet the code fits here so that the issue can be fixed).
# Information for standard sample "5610-VA-1" is replaced by the information of backup sample "5610-VA-2".
new.vals <- list(
    `5610-VA-1`=list(
        'sample.id.bkup'='5610-VA-2',
        'dose.2.time.point.tag'=0,
        'dose.3.time.point.tag'=NA,
        'dose.4.time.point.tag'=NA,
        'cohort.dose.2.flag.tag'='In cohort',
        'cohort.dose.3.flag.tag'='Other',
        'cohort.dose.4.flag.tag'='Other'
    ),
    `1203-VA-8`=list(
        'sample.id.bkup'='1203-VA-7',
        'dose.2.time.point.tag'=1,
        'dose.3.time.point.tag'=NA,
        'dose.4.time.point.tag'=NA,
        'cohort.dose.2.flag.tag'='In cohort',
        'cohort.dose.3.flag.tag'='Other',
        'cohort.dose.4.flag.tag'='Other'
    ),
    `6673-LC-1`=list(
        'sample.id.bkup'='6673-LC-2',
        'exp.sample.id.tag'='6673-LC-1', 'age.tag'=79, # because of the mislabelling issue.
        'blood.draw.date.tag'='2022-12-14',
        "gender.tag"='M', "ethnicity.tag"='White', "geo.location.site.tag"='LJI',
        "dose.1.type.tag"='Pfizer', "dose.2.type.tag"='Pfizer', "dose.3.type.tag"='Pfizer', "dose.4.type.tag"='Pfizer', "dose.5.type.tag"=NA, 
        "cv.vacc.date.1.tag"='2021-01-30', "cv.vacc.date.2.tag"='2021-02-17', "cv.vacc.date.3.tag"='2021-10-30', "cv.vacc.date.4.tag"='2022-06-22', "cv.vacc.date.5.tag"=NA,
        "max.shot.no.tag"=4,
        'dose.2.time.point.tag'=NA,
        'dose.3.time.point.tag'=NA,
        'dose.4.time.point.tag'=6,
        'cohort.dose.2.flag.tag'='Other',
        'cohort.dose.3.flag.tag'='Other',
        'cohort.dose.4.flag.tag'='In cohort'
    )
)
for(sample.id.std in names(new.vals)){
    sample.id.bkup <- new.vals[[sample.id.std]][['sample.id.bkup']]
    transfer.tags <- setdiff(x=names(new.vals[[sample.id.std]]), y='sample.id.bkup')
    for(cell.type in names(obj.list)){
        for(transfer.tag in transfer.tags){
            cat('\t\t', transfer.tag, '\n')
            rep.val <- new.vals[[sample.id.std]][[transfer.tag]]
            if(length(rep.val)!=1) stop('Unexpected error\n')
            # Transfer values from standard to backup sample.
            obj.list[[cell.type]]@meta.data[
                !is.na(obj.list[[cell.type]]@meta.data[, 'ext.donor.id.tag']) & obj.list[[cell.type]]@meta.data[, 'ext.donor.id.tag']==sample.id.bkup,
                transfer.tag
            ] <-  rep.val
            # Remove timepoint status for standard sample (i.e., to confiscate its "standard" status).
            obj.list[[cell.type]]@meta.data[
                !is.na(obj.list[[cell.type]]@meta.data[, 'ext.donor.id.tag']) & obj.list[[cell.type]]@meta.data[, 'ext.donor.id.tag']==sample.id.std,
                transfer.tag
            ] <-  NA
        }
    }
}
# Very last specific step needed to complete the correction of the mislabelling issue
for(cell.type in names(obj.list)){
    obj.list[[cell.type]]@meta.data[
        !is.na(obj.list[[cell.type]]@meta.data[, 'ext.donor.id.tag']) & obj.list[[cell.type]]@meta.data[, 'ext.donor.id.tag']=='6673-LC-2',
        'ext.donor.id.tag'
    ] <-  '6673-LC-1'
}

# Confirmation that only samples that remain without cohort assignment are the "backup" samples that were ignored or the standard samples that were replaced by the backup sample. Cross-compare with the "BackUpSample3" sheet.
# tmp.data <- as.data.table(obj.list[['cd4.all']]@meta.data)
# tmp.data[!is.na(hashtag.tag) & hashtag.tag!='Negative' & is.na(dose.2.time.point.tag) & is.na(dose.3.time.point.tag) & is.na(dose.4.time.point.tag), sort(unique(ext.donor.id.tag))] 
# tmp.data[!is.na(hashtag.tag) & hashtag.tag!='Negative', .N]
# tmp.data <- as.data.table(obj.list[['cd8.all']]@meta.data)
# tmp.data[!is.na(hashtag.tag) & hashtag.tag!='Negative' & is.na(dose.2.time.point.tag) & is.na(dose.3.time.point.tag) & is.na(dose.4.time.point.tag), sort(unique(ext.donor.id.tag))] 
# tmp.data[!is.na(hashtag.tag) & hashtag.tag!='Negative', .N]


### ---------------------- CD4RE and BTI status --------------------- ###

# ---> CD4RE status per sample.
for(cell.type in names(obj.list)){
    # Fetch basic data.
    meta.data <- as.data.table(obj.list[[cell.type]]@meta.data)
    meta.data[, barcode:=Cells(obj.list[[cell.type]])]
    # Create alternative donor ID based on cohort and time point (one per cohort)
    cohort.info <- list(
        `D2`=c('cohort.dose.2.flag.tag', 'dose.2.time.point.tag'),
        `D3`=c('cohort.dose.3.flag.tag', 'dose.3.time.point.tag'),
        `D4`=c('cohort.dose.4.flag.tag', 'dose.4.time.point.tag')
    )
    for(cohort in names(cohort.info)){
        meta.data[, dose.flag.tag:=get(cohort.info[[cohort]][1])]
        meta.data[, dose.time.point.tag:=get(cohort.info[[cohort]][2])]
        meta.data[
            !is.na(dose.flag.tag) & dose.flag.tag=='In cohort',
            tmp.donor.id:=paste(
                donor.id.tag, cohort,
                ifelse(
                    test=dose.time.point.tag>9,
                    yes=dose.time.point.tag,
                    no=paste0('0', dose.time.point.tag)
                ),
                sep='.'
            )
        ]
        tmp.data <- cd4re.info[donor.id %like% cohort, .(donor.id, cd4re.si.tag, cd4re.status.tag)]
        colnames(tmp.data)[colnames(tmp.data)!='donor.id'] <- paste0(cohort, '|', colnames(tmp.data)[colnames(tmp.data)!='donor.id'])
        meta.data <- merge(x=meta.data, y=tmp.data, by.x='tmp.donor.id', by.y='donor.id', all.x=TRUE, all.y=FALSE, sort=FALSE)
        meta.data[, tmp.donor.id:=NULL]
    }
    # Consensus CD4RE status across cohorts per sample.
    tmp.cols <- c('barcode', paste0(rep(x=names(cohort.info), each=2), c('|cd4re.si.tag', '|cd4re.status.tag')))
    tmp.data <- meta.data[, ..tmp.cols]
    tmp.data <- as.data.table(gather(data=tmp.data, key='variable', value='value', -`barcode`))
    tmp.data <- tmp.data[!is.na(value)]
    tmp.data <- as.data.table(separate(data=tmp.data, col=variable, into=c('dose', 'variable'), sep='\\|'))
    tmp.data[, dose:=NULL]
    tmp.data <- unique(tmp.data)
    tmp.check <- tmp.data[, .N, by=barcode][N>2, .N==0] # Sanity check, after rmeoving the dose information, there should be single entries (for each of two CD4RE status-related variables) for each barcode.
    if(!tmp.check) stop('Unexpected error. Multiple different CD4RE variable values defined for a single sample.\n')
    tmp.data <- as.data.table(spread(data=tmp.data, key=variable, value=value))
    meta.data <- merge(x=meta.data, y=tmp.data, by='barcode', all.x=TRUE, all.y=FALSE, sort=FALSE)
    # Find samples that were left without any CD4RE info.
    tmp.vals <- meta.data[is.na(cd4re.status.tag), unique(ext.donor.id.tag)]
    tmp.warn <- paste0('Following samples were left without CD4RE info:\n', paste0(tmp.vals, collapse=', '), '\n')
    if(length(tmp.vals)) warning(tmp.warn)
    # Save tidy CD4RE info to metadata.
    tmp.check <- all(Cells(obj.list[[cell.type]]) == meta.data[, barcode])
    if(!tmp.check) stop('Unexpected error, something went wrong.\n')
    obj.list[[cell.type]]@meta.data[, 'cd4re.si.tag'] <- meta.data[, cd4re.si.tag]
    obj.list[[cell.type]]@meta.data[, 'cd4re.status.tag'] <- meta.data[, cd4re.status.tag]
}

# ---> Exploration.
# Status does not depend on cell type. We take both metadata to perform this part of the analysis.
meta.data <- lapply(X=names(obj.list), FUN=function(x) as.data.table(obj.list[[x]]@meta.data))
names(meta.data) <- names(obj.list)
meta.data <- rbindlist(l=meta.data, use.names=TRUE, idcol='cell.type')
tmp.cols <- c('donor.id.tag', 'ext.donor.id.tag', 'dose.2.time.point.tag', 'dose.3.time.point.tag', 'dose.4.time.point.tag', 'cd4re.status.tag')
tmp.data <- unique(meta.data[!is.na(hashtag.tag) & hashtag.tag!='Negative' & !is.na(cd4re.status.tag), ..tmp.cols])
tmp.data <- as.data.table(gather(data=tmp.data, key='dose', value='time.point', -`donor.id.tag`, -`ext.donor.id.tag`, -`cd4re.status.tag`))
tmp.data <- tmp.data[!is.na(time.point)]
tmp.data[, dose:=str_replace(string=dose, pattern='.time.point.tag', replacement='')]
tmp.data[, cd4re.status.tag:=ifelse(test=cd4re.status.tag=='positive', yes=1, no=0)]
tmp.data[, ext.donor.id.tag:=NULL]
tmp.data <- unique(tmp.data)
tmp.vals <- paste0('dose.', 2:4)
tmp.data <- lapply(X=tmp.vals, FUN=function(x){
    tmp.data <- tmp.data[dose==x]
    tmp.data <- spread(data=tmp.data, key=time.point, value=cd4re.status.tag)
    return(tmp.data)
})
tmp.data <- rbindlist(l=tmp.data, use.names=TRUE, fill=TRUE)
# Clustering.
tmp.cols <- setdiff(x=colnames(tmp.data), y=c('donor.id.tag', 'dose'))
clust.data <- tmp.data[, ..tmp.cols]
clust.data[is.na(clust.data)] <- -1
dist.mat <- dist(as.matrix(clust.data), method='euclidean')
hclust.obj <- hclust(dist.mat, method='average')
# Data to plot.
data.to.plot <- as.data.frame(tmp.data[, ..tmp.cols])
row.names(data.to.plot) <- tmp.data[, paste(donor.id.tag, dose, sep='/')]
# Set color scale and breaks for heatmap.
col.breaks <- seq(from=min(range(data.to.plot, na.rm=TRUE)), to=max(range(data.to.plot, na.rm=TRUE)), length.out=100)
mid.point <- which.min(abs(col.breaks - 0))
hmap.col.scale.1 <- colorRampPalette(c('blue', 'black'))(mid.point)
hmap.col.scale.2 <- colorRampPalette(c('black', 'red'))(100-(mid.point+1))
hmap.col.scale <- c(hmap.col.scale.1, hmap.col.scale.2)
# Set rows metadata.
row.meta <- as.data.frame(tmp.data[, .(dose=str_replace(string=dose, pattern='dose\\.', replacement='Dose '))])
row.names(row.meta) <- tmp.data[, paste(donor.id.tag, dose, sep='/')]
# Plot
tmp.file.name <- paste0(reports.path, '/CD4REAlongTimeInspection.pdf')
pheatmap(
    mat=data.to.plot,
    color=hmap.col.scale, breaks=col.breaks,
    scale='none', 
    cluster_rows=hclust.obj, cluster_cols=FALSE,
    annotation_row=row.meta, # annotation_colors=ann.colors,
    show_rownames=FALSE, show_colnames=TRUE,
    filename=tmp.file.name
)


### --------------------------- Other tags -------------------------- ###

# ---> Early and late samples.
# See definition elsewhere.
es.time.points <- c(1, 3, 6, 9) # 2-dose cohort.
ls.time.points <- c(3, 6, 9, 12) # 4-dose cohort.
# Define tag.
for(cell.type in names(obj.list)){
    # Fetch basic data.
    meta.data <- as.data.table(obj.list[[cell.type]]@meta.data)
    meta.data[, barcode:=Cells(obj.list[[cell.type]])]
    # Define a tag independently for early samples and late samples.
    meta.data[, `:=`(es.tag=FALSE, ls.tag=FALSE)]
    meta.data[
        !is.na(dose.2.time.point.tag) & dose.2.time.point.tag %in% es.time.points,
        es.tag:=TRUE
    ]
    meta.data[
        !is.na(dose.4.time.point.tag) & dose.4.time.point.tag %in% ls.time.points,
        ls.tag:=TRUE
    ]
    # Confirm that early sample and late sample statuses are mutually exclusive.
    tmp.check <- meta.data[ls.tag & es.tag, .N==0]
    if(!tmp.check) stop('Unexpected error.\n')
    # Unique early vs late sample tag.
    meta.data[es.tag==TRUE, evl.sample.tag:='Early']
    meta.data[ls.tag==TRUE, evl.sample.tag:='Late']
    # Save tidy info to metadata.
    tmp.check <- all(Cells(obj.list[[cell.type]]) == meta.data[, barcode])
    if(!tmp.check) stop('Unexpected error, something went wrong.\n')
    obj.list[[cell.type]]@meta.data[, 'evl.sample.tag'] <- meta.data[, evl.sample.tag]
}

# ---> Unique cohort tag.
# Define a unique cohort for each experimental sample. This is the same regardless of cell type.
meta.data <- lapply(X=names(obj.list), FUN=function(x) as.data.table(obj.list[[x]]@meta.data))
names(meta.data) <- names(obj.list)
meta.data <- rbindlist(l=meta.data, use.names=TRUE, idcol='cell.type')
tmp.cols <- c('donor.id.tag', 'ext.donor.id.tag', 'dose.2.time.point.tag', 'dose.3.time.point.tag', 'dose.4.time.point.tag')
tmp.data <- unique(meta.data[!is.na(hashtag.tag) & hashtag.tag!='Negative', ..tmp.cols])
tmp.data <- as.data.table(gather(data=tmp.data, key='dose', value='time.point', -`donor.id.tag`, -`ext.donor.id.tag`))
tmp.data <- tmp.data[!is.na(time.point)]
tmp.data[, dose:=str_extract(string=dose, pattern='\\d')]
# Find and inspect experimental samples that belong to multiple dose cohorts.
to.check <- tmp.data[, .(cohort.count=uniqueN(dose)), by=ext.donor.id.tag][cohort.count>1, ext.donor.id.tag]
to.check <- tmp.data[ext.donor.id.tag %in% to.check]
to.check <- spread(data=to.check, key=dose, value=time.point)
# This immediately proves that, when disregarding the 0 time point from any cohort, experimental samples can be unambiguously assigned to a single dose cohort. Let's complete the assignment by disregarding the time point 0 from any kind of cohort.
tmp.data <- tmp.data[
    time.point!=0,
    .(
        ext.donor.id.tag,
        dose.cohort.tag=paste('dose', dose, sep='.'),
        dose.cohort.time.point.tag=paste('dose', dose, 'tp', time.point, sep='.')
    )
]
to.check <- tmp.data[, .N==uniqueN(ext.donor.id.tag)]
if(!to.check) stop('Encountered ambiguous cohort assignments to experimental samples.\n')
# Merge info into metadata.
for(cell.type in names(obj.list)){
    # Fetch basic data.
    meta.data <- as.data.table(obj.list[[cell.type]]@meta.data)
    meta.data[, barcode:=Cells(obj.list[[cell.type]])]
    # Define a tag independently for early samples and late samples.
    meta.data <- merge(x=meta.data, y=tmp.data, by='ext.donor.id.tag', all.x=TRUE, all.y=FALSE, sort=FALSE)
    # Save tidy info to metadata.
    tmp.check <- all(Cells(obj.list[[cell.type]]) == meta.data[, barcode])
    if(!tmp.check) stop('Unexpected error, something went wrong.\n')
    obj.list[[cell.type]]@meta.data[, 'dose.cohort.tag'] <- meta.data[, dose.cohort.tag]
    obj.list[[cell.type]]@meta.data[, 'dose.cohort.time.point.tag'] <- meta.data[, dose.cohort.time.point.tag]
}


### ------------------------- Final details ------------------------- ###

# ---> Specific removal of cell clusters.
# ---> HLTY CD4 (EXAMPLE FROM PREV. PROJECT)
# tmp.lab <- 'hlty.cd4'
# # Clusters to remove: 8
# tmp.remove <- '8'
# tmp.data <- obj.list[[tmp.lab]]@meta.data[, clust.labs[tmp.lab]]
# tmp.cells <- Cells(obj.list[[tmp.lab]])[!tmp.data %in% tmp.remove]
# obj.list[[tmp.lab]] <- subset(x=obj.list[[tmp.lab]], cells=tmp.cells)

# ---> Save seurat objects.
lapply(X=names(obj.list), FUN=function(tmp.obj){
    tmp.file.name <- paste0(reports.path, '/SeuratObj_', obj.extended.names[tmp.obj], '.RDS')
    if(!file.exists(tmp.file.name)) saveRDS(file=tmp.file.name, object=obj.list[[tmp.obj]]) else cat(paste0('File already exists for object: ', tmp.obj, '\n'))
})
