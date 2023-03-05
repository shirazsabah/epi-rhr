###########################         Import, remove duplicates and compress NJR        ###########################

pacman::p_load(pacman, data.table, rio, tidyverse, logger)

setwd(data_dir)
setDTthreads(10)

## Import
# rHR
rh <- import("NJR/HipAllRevisions.txt", quote="")
rh <- as.data.table(rh)

# pHR
ph <- import("NJR/HipPrimaryOutcomes.txt", quote="")
ph <- as.data.table(ph)

log_info('Initial NJR import: {ph %>% count()} pHA records & {rh %>% count()} rHA records')

## Change all names to lower case
rh <- rh %>% janitor::clean_names()
ph <- ph %>% janitor::clean_names()

## Sort out dates
rh[, op_date := as.Date(fasttime::fastPOSIXct(op_date))]
ph[, primary_op_date := as.Date(fasttime::fastPOSIXct(primary_op_date))]
ph[, revision_date := as.Date(fasttime::fastPOSIXct(revision_date))]

## Drop 'cat_no', 'details' and 'text' columns (by setting them to NULL)
rh[, grep("details", names(rh)) := NULL]
rh[, grep("text", names(rh)) := NULL]
ph[, grep("details", names(ph)) := NULL]
ph[, grep("cat_no", names(ph)) := NULL]
ph[, grep("text", names(ph)) := NULL]

## Replace "NULL", "" or "na" with NAs 
# Need to exclude date columns from this

# Revision dataset
p_load(car)
.cols <- setdiff(colnames(rh), "op_date")
rh[,(.cols):=lapply(.SD, recode, '"NULL"=NA'), .SDcols = .cols]
rh[,(.cols):=lapply(.SD, recode, '""=NA'), .SDcols = .cols]
rh[,(.cols):=lapply(.SD, recode, '"na"=NA'), .SDcols = .cols]

# Primary dataset
.cols <- setdiff(colnames(ph), c("primary_op_date", "revision_date"))
ph[,(.cols):=lapply(.SD, recode, '"NULL"=NA'), .SDcols = .cols]
ph[,(.cols):=lapply(.SD, recode, '""=NA'), .SDcols = .cols]
ph[,(.cols):=lapply(.SD, recode, '"na"=NA'), .SDcols = .cols]

## Convert all character fields to factor

# Revision dataset
changeCols <- colnames(rh)[which(as.vector(rh[,lapply(.SD, class)]) == "character")]
rh[,(changeCols):= lapply(.SD, as.factor), .SDcols = changeCols]

# Primary dataset
changeCols <- colnames(ph)[which(as.vector(ph[,lapply(.SD, class)]) == "character")]
ph[,(changeCols):= lapply(.SD, as.factor), .SDcols = changeCols]

# Convert numeric fields
# Except for those which need to remain numeric
# Revision dataset
changeCols <- colnames(rh)[which(as.vector(rh[,lapply(.SD, class)]) == "numeric")]
changeCols <- setdiff(changeCols, "bmi")
rh[,(changeCols):= lapply(.SD, as.factor), .SDcols = changeCols]

# Primary dataset
changeCols <- colnames(ph)[which(as.vector(ph[,lapply(.SD, class)]) == "numeric")]
changeCols <- setdiff(changeCols, c("bmi", "primary_to_outcome_years"))
ph[,(changeCols):= lapply(.SD, as.factor), .SDcols = changeCols]

# Remove duplicates based on ALL variables
phslim <- unique(ph)
rhslim <- unique(rh)

# Attrition data for flowchart
no_phr_supp <- ph[,.N]
no_rhr_supp <- rh[,.N]
no_uniq_phr_supp <- phslim[,.N]
no_uniq_rhr_supp <- rhslim[,.N]

# Save the environment for quickly re-loading the compressed NJR with all fields
rm(list=setdiff(ls(), c("phslim","rhslim", "no_phr_supp", "no_rhr_supp", "no_uniq_phr_supp", "no_uniq_rhr_supp")))

save.image("R_IMAGES/epi-rh-njr.RData")
save.image("R_IMAGES/epi-rh-njr-backup.RData")