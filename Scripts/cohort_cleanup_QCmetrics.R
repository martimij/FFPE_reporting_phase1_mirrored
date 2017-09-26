# Martina Mijuskovic
# FFPE project - small variant reporting
# Cohort cleanup and QC metrics
# Sept 2017

library(dplyr)
library(ggplot2)
library(data.table)

##### Clean cohort and QC metrics data ##### 

### Load upload report

# Look into unfiltered upload report
upload_full <- read.table(paste0("./Data/upload_report.", today, ".txt"), sep = "\t")
colnames(upload_full) <- as.character(fread(paste0("./Data/upload_report.", today, ".txt"), skip = 14, nrows = 1, header = F))
upload_full$Path <- as.character(upload_full$Path)
table(upload_full$Status, upload_full$Type)  # some have unknown status, checking those
table((upload_full %>% filter(Status == "unknown", Type != "rare disease") %>% pull(`Delivery Version`)))  # "Unknown" are non-V4 versions
# V1 V1.5   V2   V4 
# 112  101  350    0

# Restrict upload report to cancer, V4 (but NOT qc_passed; some samples that shouldn't fail still have "failed" status)
upload <- read.table(paste0("./Data/upload_report.", today, ".txt"), sep = "\t")
colnames(upload) <- as.character(fread(paste0("./Data/upload_report.", today, ".txt"), skip = 14, nrows = 1, header = F))
#upload <- upload %>% filter(`Delivery Version` == "V4", Status == "qc_passed", Type %in% c("cancer germline", "cancer tumour"))
upload <- upload %>% filter(`Delivery Version` == "V4", Type %in% c("cancer germline", "cancer tumour"))
upload$Path <- as.character(upload$Path)
table(upload$Status)  # 40 have qc_failed status


### Load cohort and QC data

# Load full list of FFPE plates
FFPE_plates <- read.table("./Data/FFPE_main_program_plates.txt")
FFPE_plates <- as.character(FFPE_plates$V1)

# Load Alona's QC metrics table
QC <- read.csv("./Data/pre-seq_metrics.FFPE.seq_metrics.csv")
dim(QC)  # 245
sum(duplicated(QC$WELL_ID))
summary(QC)  # 11 samples seem to have NAs for most fields
QC %>% filter(is.na(TUMOUR_PURITY))  # these 11 fail contamination in ContEst and ConPair
QC_failed <- as.character(QC %>% filter(is.na(TUMOUR_PURITY)) %>% pull(WELL_ID))
# Remove QC failed samples from the list
QC <- QC %>% filter(!WELL_ID %in% QC_failed)
dim(QC) # 234 samples left that pass contamination QC (NOTE that 3 samples still have 1-3% contmaination but don't fail)
summary(QC)
sum(is.na(QC))
# Examine missing values
QC %>% filter(TUMOUR_TYPE == "N/A") %>% select(WELL_ID, TUMOUR_TYPE, COLLECTION_DATE)  # 50 samples missing tumour type and collection date
# Write table with missing data
write.table((QC %>% filter(TUMOUR_TYPE == "N/A") %>% select(WELL_ID, TUMOUR_TYPE, COLLECTION_DATE)), file = "./Data/FFPE_missing_tumourType_collDate.tsv", sep = "\t", quote = F, col.names = T, row.names = F)


# Find FFPE plate samples in the (cancer, clean) upload report
upload$Platekey <- as.character(upload$Platekey)
upload$Plate <- sapply(upload$Platekey, function(x){strsplit(x, split = "-")[[1]][1] })
FFPE_list <- upload %>% filter(Plate %in% FFPE_plates)
dim(FFPE_list)  # 246
sum(duplicated(FFPE_list$Platekey))  # 1 duplicated (LP3000140-DNA_A03)
FFPE_list %>% filter(Platekey == FFPE_list[duplicated(FFPE_list$Platekey),]$Platekey)  # One version has a delivery problem, removing (ID BE00000001)
FFPE_list <- FFPE_list %>% filter(DeliveryID != "BE00000001")
table(FFPE_list$Status, exclude = NULL)  # 5 samples are QC failed, 240 pass
FFPE_list %>% filter(Status == "qc_failed") %>% select(Platekey, DeliveryID, `Delivery Date`, `Delivery Version`, Status)


# Add matching germlines from the upload report, check versions and QC status
FFPE_list$Path <- as.character(FFPE_list$Path )
FFPE_list$Platekey_normal <- sapply(FFPE_list$Path, function(x){strsplit(strsplit(x, split = "/")[[1]][6], split = "Normal")[[1]][2]})
dim(upload %>% filter(Platekey %in% FFPE_list$Platekey_normal))  # 246
FFPE_list$Platekey_normal[!FFPE_list$Platekey_normal %in% upload$Platekey]  # 0 (no normals missing from upload)

# Look at matching normals in the full upload report (2 normals are qc_failed)
table((upload_full %>% filter(Platekey %in% FFPE_list$Platekey_normal) %>% pull(`Delivery Version`)), (upload_full %>% filter(Platekey %in% FFPE_list$Platekey_normal) %>% pull(Status)))
upload_full %>% filter(Platekey %in% FFPE_list$Platekey_normal, `Delivery Version` == "V4", Status == "qc_failed")
normals_failed <- upload_full %>% filter(Platekey %in% FFPE_list$Platekey_normal, `Delivery Version` == "V4", Status == "qc_failed") %>% pull(Platekey)
FFPE_list[FFPE_list$Platekey_normal %in% normals_failed,]
upload_full %>% filter(Platekey %in% normals_failed) 
# LP3000069-DNA_E03 QC-passed when delivered in V2, but failed when it was lifted to V4 (LP3000069-DNA_E03); Bertha error: "Minimum covered positions requirement not met: Positions covered at >= 15x is 94.5%. Must be >= 95%."
# LP3000115-DNA_G11 QC failed due to invalid BAM (DeliveryID CANCT40004) but then re-delivered and passed (HX01751975) ??? the odd thing is that HX01751975 is delivered earlier

# # Exclude FFPE sample with matching normal LP3000069-DNA_E03  ---> update: NOT excluding this, needs Bertha QC intake re-run
# FFPE_list <- FFPE_list %>% filter(Platekey_normal != "LP3000069-DNA_E03")
# dim(FFPE_list) 

# Compare clean FFPE list with Alona's QC metrics table containing samples that pass contamination QC
FFPE_list %>% filter(Platekey %in% QC_failed)  # 11/11 samples that fail contamination are still here
# Remove contaminated samples from the FFPE list
dim(FFPE_list)  # 245
FFPE_list <- FFPE_list %>%  filter(!Platekey %in% QC_failed)  
dim(FFPE_list)  # 234 (check if they are all in Alona's QC metrics list)
QC %>% filter(!WELL_ID %in% FFPE_list$Platekey) # all there


##### Final FFPE cohort list ##### 

### Total of 234 samples that have no QC or other issues (confirmation pending Bertha QC intake re-run)

# Add QC data to FFPE list
FFPE_list <- full_join(FFPE_list, QC, by = c("Platekey" = "WELL_ID"))
dim(FFPE_list) # 234
write.csv(FFPE_list, file = "./Data/Clean_FFPE_samplelist.csv", quote = F, row.names = F)


##### Collect FFPE + normal data from catalog ##### 

# Get germline and tumour QC data and make sure germline is not contaminated








##### Summary analysis of QC metrics ##### 

summary(QC)

# Examine contaminated samples (i.e. non "Pass")
QC %>% filter(TUMOUR_CONTAMINATION_CONTEST != "Pass") %>% select(WELL_ID, TUMOUR_CONTAMINATION_CONTEST, TUMOUR_CONTAMINATION_CONPAIR, TUMOUR_PURITY, TUMOUR_TYPE, COSMIC_COV_LT30X, MEDIAN_COV)

# Establish median +/- 2 SD values for FFPE cohort
QC %>% summarise(AT_DROP_MED = median(AT_DROP), AT_DROP_SD = sd(AT_DROP), GC_DROP_MED = median(GC_DROP), GC_DROP_SD = sd(GC_DROP))
