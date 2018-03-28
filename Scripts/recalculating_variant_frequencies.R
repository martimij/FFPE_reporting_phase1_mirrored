# Martina Mijuskovic
# Small variant flagging
# Re-calculating variant frequencies in 3 cancer cohorts (main program samples)
# March 2018

library(dplyr)
library(VariantAnnotation)
library(ggplot2)
library(VennDiagram)
library(scales)
library(ensembldb)
library(data.table)
library(reshape)
#library(jsonlite)

today <- Sys.Date()

########## Get upload report ########## 

system(paste0("wget ", "https://upload-reports.gel.zone/upload_report.", today, ".txt"))
upload_23Mar2018 <- read.table(paste0("./upload_report.", today, ".txt"), sep = "\t")
colnames(upload_23Mar2018) <- as.character(fread(paste0("./upload_report.", today, ".txt"), skip = 14, nrows = 1, header = F))
upload_23Mar2018$Path <- as.character(upload_23Mar2018$Path)
dim(upload_23Mar2018)  #51063
table(upload_23Mar2018$Status, upload_23Mar2018$Type)  # some have unknown status, checking those
#                   cancer germline cancer tumour rare disease unknown
#                               268           295         2193       0
# awaiting_delivery               1            11            9       0
# delete_permitted                0             0            1       0
# delivered                       0             0            0       1
# md5_passed                     19            17            2       0
# qc_failed                      29            26           24       1
# qc_passed                    5643          5456        36880     180
# upload_failed                   1             0            5       1
table((upload_23Mar2018 %>% dplyr::filter(Status == "", Type != "rare disease") %>% pull(`Delivery Version`)))  # "Unknown" are non-V4 versions
# unknown      V1    V1.5      V2      V4 
# 0     112     101     350       0 

# Restrict upload report to cancer, V4 and qc_passed
upload <- upload_23Mar2018 %>% dplyr::filter(`Delivery Version` == "V4", Type %in% c("cancer germline", "cancer tumour"), Status == "qc_passed")
table(upload$Type)

# Subset to new samples only (excluding any potential samples from the same patient that were already processed)
FF_list$SAMPLE_WELL_ID_normal <- sapply(1:length(FF_list$SAMPLE_WELL_ID), function(x){
  strsplit(FF_list$Path[x], "_Normal")[[1]][2]
})
dim(upload %>% dplyr::filter(Platekey %in% FFPE_list$Platekey))  # 244
dim(upload %>% dplyr::filter(Platekey %in% FF_list$SAMPLE_WELL_ID))  # 1061
dim(upload %>% dplyr::filter(Platekey %in% FFPE_list$Platekey_normal)) # 245
dim(upload %>% dplyr::filter(Platekey %in% FF_list$SAMPLE_WELL_ID_normal)) # 1057

upload_subset <- upload %>% dplyr::filter(!(Platekey %in% FFPE_list$Platekey), !(Platekey %in% FF_list$SAMPLE_WELL_ID), !(Platekey %in% FFPE_list$Platekey_normal), !(Platekey %in% FF_list$SAMPLE_WELL_ID_normal))
dim(upload_subset)  # 6443
table(upload_subset$Type)  # 3103 potential new samples in the total cancer cohort (in addition to previous 1306)
# cancer germline   cancer tumour    rare disease         unknown 
#           3340            3103               0               0 


########## Get library prep info ########## 

# Read library prep info from stats.json delivered by Illumina (only for tumour, not GL)

cancer_tumour <- upload_subset %>% dplyr::filter(Type == "cancer tumour")
dim(cancer_tumour)  # 3103
cancer_tumour$STATS_JSON_PATH <- sapply(1:dim(cancer_tumour)[1], function(x){
  paste0(cancer_tumour$Path[x], "/", cancer_tumour$Platekey[x], ".stats.json")
})

# Write the file
write.csv(cancer_tumour, file = "./Data/new_tumour_samples_23Mar2019.csv", quote = F, col.names = T, row.names = F)

### HPC - read sample list and get library prep info 
library(dplyr)
library(jsonlite)
cancer_tumour <- read.csv("/home/mmijuskovic/small_variant_freq_March2018/new_tumour_samples_23Mar2019.csv")
cancer_tumour$STATS_JSON_PATH <- as.character(cancer_tumour$STATS_JSON_PATH)
cancer_tumour$LIBRARY_TYPE <- sapply(1:dim(cancer_tumour)[1], function(x){
  if(file.exists(cancer_tumour$STATS_JSON_PATH[x])){  # check if file exists
    fromJSON(cancer_tumour$STATS_JSON_PATH[x])$sample_library_type
    #fromJSON("/genomes/by_date/2018-03-21/CX03283736/CancerLP3000771-DNA_B03_NormalLP3000770-DNA_B03/LP3000771-DNA_B03.stats.json")$sample_library_type
  }
  else { "" }
})

# Write table including library prep info
write.csv(cancer_tumour, file = "/home/mmijuskovic/small_variant_freq_March2018/new_tumour_samples_wLibraryPrep_23Mar2019.csv", quote = F, row.names = F, col.names = T)



### Local: read library prep info table

cancer_tumour <- read.csv("./Data/new_tumour_samples_wLibraryPrep_23Mar2019.csv")

# Overview of library types
as.data.frame(table(cancer_tumour$LIBRARY_TYPE)) # some unknown, some just "pcr", some missing JSONs

### Library type cleanup
cancer_tumour$LIBRARY_TYPE <- as.character(cancer_tumour$LIBRARY_TYPE)
cancer_tumour[cancer_tumour$LIBRARY_TYPE == "pcr_ffpe",]$LIBRARY_TYPE <- "FFPE"
cancer_tumour[cancer_tumour$LIBRARY_TYPE == "TruSeq FFPE High Throughput",]$LIBRARY_TYPE <- "FFPE"
cancer_tumour[cancer_tumour$LIBRARY_TYPE == "pcrfree",]$LIBRARY_TYPE <- "TruSeq PCR-Free"
cancer_tumour[cancer_tumour$LIBRARY_TYPE == "TruSeq PCR-Free High Throughput",]$LIBRARY_TYPE <- "TruSeq PCR-Free"
cancer_tumour[cancer_tumour$LIBRARY_TYPE == "TruSeq Nano High Throughput",]$LIBRARY_TYPE <- "TruSeq Nano"

# Check 1 "Unknown" sample
cancer_tumour %>% dplyr::filter(LIBRARY_TYPE == "Unknown") # Manual check of JSON shows "sample_library_type": "Unknown,TruSeq PCR-Free High Throughput"
cancer_tumour[cancer_tumour$Platekey == "LP3000539-DNA_E05",]$LIBRARY_TYPE <- "TruSeq PCR-Free"

# Check 25 "pcr" samples
cancer_tumour %>% dplyr::filter(LIBRARY_TYPE == "pcr") # Manual checks show no more information but "pcr"; exclude that they are FFPE
FFPE_list_new %>% dplyr::filter(Platekey %in% cancer_tumour[cancer_tumour$LIBRARY_TYPE == "pcr",]$Platekey)  # None on the FFPE plates
# Check if any samples are on the new FFPE plate LP3000606 and if correct library type is extracted
cancer_tumour[grepl("LP3000606", cancer_tumour$Path),]  # Yes, all correct library type
# Save sample names with "pcr" sample_library_type to check later against Illumina's manifests; 
cancer_tumour_pcr_libType <- cancer_tumour %>% dplyr::filter(LIBRARY_TYPE == "pcr")
# Set to nano
cancer_tumour[cancer_tumour$LIBRARY_TYPE == "pcr",]$LIBRARY_TYPE <- "TruSeq Nano"


# Look at samples with missing data for library type
cancer_tumour %>% dplyr::filter(LIBRARY_TYPE == "") # 455
# How many are lifted in-house (CANCP DeliveryID) --> all of them apparently, but deliveries range up to Feb 2018
dim(cancer_tumour[cancer_tumour$LIBRARY_TYPE == "",][grepl("CANCP", cancer_tumour[cancer_tumour$LIBRARY_TYPE == "",]$DeliveryID),]) # 454
cancer_tumour[cancer_tumour$LIBRARY_TYPE == "",][!grepl("CANCP", cancer_tumour[cancer_tumour$LIBRARY_TYPE == "",]$DeliveryID),]
sum(grepl("CANCP", cancer_tumour[cancer_tumour$LIBRARY_TYPE == "",]$DeliveryID)) # 454
# Checking upload report samples for some of them
upload_23Mar2018 %>% dplyr::filter(Platekey == "LP2000283-DNA_D02")  # 3 deliveries, V1.5, V2, V4 (Feb 2018), likely a pilot sample
upload_23Mar2018 %>% dplyr::filter(Platekey == "LP2000287-DNA_E02")  # 3 deliveries, V1.5, V2, V4 (Feb 2018), likely a pilot sample
upload_23Mar2018 %>% dplyr::filter(Platekey == "LP2000297-DNA_F01")  # 2 deliveries, V2, V4 (2017)
upload_23Mar2018 %>% dplyr::filter(Platekey == "LP2000939-DNA_H02")  # 2 deliveries, V2, V4 (2016)
# Checking delivery dates for unknown samples
cancer_tumour$Delivery.Date <- as.character(cancer_tumour$Delivery.Date)
table(cancer_tumour[cancer_tumour$LIBRARY_TYPE == "",]$Delivery.Date)  # Most on 2016-11-23, 2017-01-26, 2017-02-20, 2017-11-03




########## FFPE cleanup ########## 


########## FF cleanup ########## 

# Remove previously merged samples

# Remove known low quality samples

