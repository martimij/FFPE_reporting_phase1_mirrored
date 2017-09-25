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

# Restrict upload report to cancer, V4 and qc_passed
upload <- read.table(paste0("./Data/upload_report.", today, ".txt"), sep = "\t")
colnames(upload) <- as.character(fread(paste0("./Data/upload_report.", today, ".txt"), skip = 14, nrows = 1, header = F))
upload <- upload %>% filter(`Delivery Version` == "V4", Status == "qc_passed", Type %in% c("cancer germline", "cancer tumour"))
upload$Path <- as.character(upload$Path)


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
FFPE_list <- upload %>% filter(Plate %in% FFPE_plates)  # 240 FFPE samples
sum(duplicated(FFPE_list$Platekey))  # 0

# Add matching germlines from the upload report, check versions and QC status
FFPE_list$Path <- as.character(FFPE_list$Path )
FFPE_list$Platekey_normal <- sapply(FFPE_list$Path, function(x){strsplit(strsplit(x, split = "/")[[1]][6], split = "Normal")[[1]][2]})
dim(upload %>% filter(Platekey %in% FFPE_list$Platekey_normal))  # 239 (one matching normal missing)
FFPE_list$Platekey_normal[!FFPE_list$Platekey_normal %in% upload$Platekey]  # "LP3000069-DNA_E03" missing

# Look at matching normals in the full upload report (2 qc failed, check if one is the missing one)
table((upload_full %>% filter(Platekey %in% FFPE_list$Platekey_normal) %>% pull(`Delivery Version`)), (upload_full %>% filter(Platekey %in% FFPE_list$Platekey_normal) %>% pull(Status)))
upload_full %>% filter(Platekey %in% FFPE_list$Platekey_normal, `Delivery Version` == "V4", Status == "qc_failed")
normals_failed <- upload_full %>% filter(Platekey %in% FFPE_list$Platekey_normal, `Delivery Version` == "V4", Status == "qc_failed") %>% pull(Platekey)
FFPE_list[FFPE_list$Platekey_normal %in% normals_failed,]
upload_full %>% filter(Platekey == "LP3000069-DNA_E03") # QC-passed when delivered in V2, but failed when it was lifted to V4 (LP3000069-DNA_E03); Bertha error: "Minimum covered positions requirement not met: Positions covered at >= 15x is 94.5%. Must be >= 95%."

# Exclude FFPE sample with matching normal LP3000069-DNA_E03
FFPE_list <- FFPE_list %>% filter(Platekey_normal != "LP3000069-DNA_E03")
dim(FFPE_list)  # 239 samples total

# Compare clean FFPE list with Alona's QC metrics table containing samples that pass contamination QC
FFPE_list %>% filter(Platekey %in% QC_failed)  # 10/11 samples that fail contamination are still here
# Remove contaminated samples from the FFPE list
dim(FFPE_list)  # 239
FFPE_list <- FFPE_list %>%  filter(!Platekey %in% QC_failed)  
dim(FFPE_list)  # 229  # Alona's list has 234; where are the missing 5 samples?
missing_FFPE <- QC %>% filter(!WELL_ID %in% FFPE_list$Platekey) %>% pull(WELL_ID)
upload_full %>% filter(Platekey %in% missing_FFPE) # These 5 were lifted from V2 to V4
upload_full$Platekey <- as.character(upload_full$Platekey)
table((upload_full %>% filter(Platekey %in% missing_FFPE) %>% pull(`Delivery Version`)), (upload_full %>% filter(Platekey %in% missing_FFPE) %>% pull(Platekey)))
table((upload_full %>% filter(Platekey %in% missing_FFPE) %>% pull(`Delivery Version`)), (upload_full %>% filter(Platekey %in% missing_FFPE) %>% pull(Platekey)), (upload_full %>% filter(Platekey %in% missing_FFPE) %>% pull(Status)))
# 4 fail QC in V4 but pass in V2 and one (LP3000074-DNA_D01) is the one where matching normal fails


##### Final FFPE cohort list ##### 

### Total of 229 samples that have no QC or other issues

# Add QC data to FFPE list
FFPE_list <- inner_join(FFPE_list, QC, by = c("Platekey" = "WELL_ID"))
dim(FFPE_list) # 229
write.csv(FFPE_list, file = "./Data/Clean_FFPE_samplelist.csv", quote = F, row.names = F)




##### Summary analysis of QC metrics ##### 

summary(QC)

# Examine contaminated samples (i.e. non "Pass")
QC %>% filter(TUMOUR_CONTAMINATION_CONTEST != "Pass") %>% select(WELL_ID, TUMOUR_CONTAMINATION_CONTEST, TUMOUR_CONTAMINATION_CONPAIR, TUMOUR_PURITY, TUMOUR_TYPE, COSMIC_COV_LT30X, MEDIAN_COV)

# Establish median +/- 2 SD values for FFPE cohort
QC %>% summarise(AT_DROP_MED = median(AT_DROP), AT_DROP_SD = sd(AT_DROP), GC_DROP_MED = median(GC_DROP), GC_DROP_SD = sd(GC_DROP))
