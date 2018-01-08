# Martina Mijuskovic
# FFPE project - small variant reporting
# Calculating and analysing recurrent FFPE variants (main program samples)
# Sept 2017

library(dplyr)
library(VariantAnnotation)
library(ggplot2)
library(VennDiagram)
library(scales)
library(ensembldb)
library(data.table)

today <- Sys.Date()


##### Get VCF paths (FFPE) ##### 

# Get the list of FFPE paths (main program, 232 non-excluded samples)
table(FFPE_list$Excluded, FFPE_list$Exclusion_reason) # exclude sample that still hasn't been reprocessed (LP3000074-DNA_D06)
samples_for_recurr <- FFPE_list %>% filter(Exclusion_reason == "")
# Write the sample list (to be transferred to HPC)
write.csv(samples_for_recurr, file = "./Data/samples_for_recurrent_calc.csv", quote = F, row.names = F)

##### Get VCF paths for normalization (run on LSF)

# Read the table of FFPE samples for reporting/recurrent variant calculation
samples_for_recurr <- read.csv("/home/mmijuskovic/small_variant_freq/FFPE/samples_for_recurrent_calc.csv", header = T)

# Get VCF paths (note that there are 2 somatic SNV VCF paths for 26 FFPE trios, they were normalized before)
paths <- unlist(sapply(samples_for_recurr$Path, function(x){
  command <- paste("find", paste0(x, "/SomaticVariations"), "-iname *.somatic.vcf.gz", sep = " ")
  system(command, intern = T)
}))
paths <- paths[!grepl("PASS.duprem.atomic.left.split.somatic.vcf.gz", paths)]
samples_for_recurr$SNV_VCF_path <- as.character(paths)

# Write the table with VCF paths
write.csv(samples_for_recurr, file = "/home/mmijuskovic/small_variant_freq/FFPE/samples_for_recurrent_calc_withVCFpaths.csv", quote = F, row.names = F)

# Write sample names and VCF paths into the same table
write.table((samples_for_recurr %>% select(Platekey, SNV_VCF_path)), file = "/home/mmijuskovic/small_variant_freq/FFPE/samples_for_recurrent_calc_normalization_input.csv", sep = ",", col.names = F, quote = F, row.names = F)

# Split input table into 8 files (for 8 parallel jobs, 29 samples each)
sapply(0:7, function(x){
  # Get sample index
  i <- (x*29)+1
  j <- i+28
  z <- x+1
  # Write input table containing (next) 29 samples
  tabl <- (samples_for_recurr %>% select(Platekey, SNV_VCF_path))[i:j,]
  write.table(tabl, file = paste0("/home/mmijuskovic/small_variant_freq/FFPE/norm_VCFs_input/samples_for_recurrent_calc_normalization_input_", z, ".csv"), sep = ",", col.names = F, quote = F, row.names = F)
  
})


##### Get VCF paths for adding GT on normalized VCFs (run on LSF)

# Read the table of FFPE samples for reporting/recurrent variant calculation (includes original VCF paths)
samples_for_recurr <- read.csv("/home/mmijuskovic/small_variant_freq/FFPE/samples_for_recurrent_calc_withVCFpaths.csv", header = T)

# Add normalized VCF paths
samples_for_recurr$norm_SNV_VCF_path <- paste0("/home/mmijuskovic/small_variant_freq/FFPE/norm_VCFs/", samples_for_recurr$Platekey, ".somatic.duprem.left.split.vcf.gz")
  
# Write out the table of normalized VCF paths
write.csv(samples_for_recurr, file = "/home/mmijuskovic/small_variant_freq/FFPE/samples_for_recurrent_calc_withNormVCFpaths.csv", quote = F, row.names = F)

# Split input table into 8 files (for 8 parallel jobs, 29 samples each)
sapply(0:7, function(x){
  # Get sample index
  i <- (x*29)+1
  j <- i+28
  z <- x+1
  # Write input table containing (next) 29 samples
  tabl <- (samples_for_recurr %>% select(Platekey, norm_SNV_VCF_path))[i:j,]
  write.table(tabl, file = paste0("/home/mmijuskovic/small_variant_freq/FFPE/addGT_input/samples_for_recurrent_calc_addGT_input_", z, ".csv"), sep = ",", col.names = F, quote = F, row.names = F)
  
})


##### Get VCF paths for merging somatic VCFs

vcf_list <- paste0("/home/mmijuskovic/small_variant_freq/FFPE/GT_VCFs/", samples_for_recurr$Platekey, ".GT.duprem.left.split.vcf.gz")
write.table(vcf_list, file = "./Data/FFPE_mainProgram_2017.txt", col.names = F, quote = F, row.names = F)







##### Get VCF paths (FF) ##### 

# Read the list of samples, subset to FF pass contamination QC
FF_list <- read.csv("./Data/ready_to_upload.09-24-17.csv")

# Subset to FF, pass contamination
FF_list <- FF_list %>% filter(SAMPLE_TYPE == "FF", TUMOUR_CONTAMINATION == "Pass")
  
# Sanity check
table(FF_list$PASS_TO_SEQ, exclude = NULL)
table(FF_list$LIBRARY_TYPE, exclude = NULL)
table(FF_list$LIBRARY_TYPE, FF_list$CENTER_CODE, exclude = NULL)
table(FF_list$CENTER_CODE, exclude = NULL)
sum(duplicated(FF_list$SAMPLE_WELL_ID))  # 0

# Check QC metrics by GMC
ggplot(FF_list, aes(x=CENTER_CODE, y=AT_DROP, fill = CENTER_CODE)) + geom_boxplot()
ggplot(FF_list, aes(x=CENTER_CODE, y=GC_DROP, fill = CENTER_CODE)) + geom_boxplot()
ggplot(FF_list, aes(x=CENTER_CODE, y=COVERAGE_HOMOGENEITY, fill = CENTER_CODE)) + geom_boxplot()
ggplot(FF_list, aes(x=CENTER_CODE, y=AV_FRAGMENT_SIZE_BP, fill = CENTER_CODE)) + geom_boxplot()
ggplot(FF_list, aes(x=CENTER_CODE, y=MAPPING_RATE_PER, fill = CENTER_CODE)) + geom_boxplot()
ggplot(FF_list, aes(x=CENTER_CODE, y=CHIMERIC_PER, fill = CENTER_CODE)) + geom_boxplot()


# By mutational burden
ggplot(FF_list, aes(x=CENTER_CODE, y=SNV, fill = CENTER_CODE)) + geom_boxplot()
ggplot(FF_list, aes(x=CENTER_CODE, y=INDEL, fill = CENTER_CODE)) + geom_boxplot()



### Add BAM paths

# Check for samples in upload reprot ("upload" object contains only cancer, only V4)
dim(FF_list)  # 1070
dim(upload %>% filter(Platekey %in% FF_list$SAMPLE_WELL_ID)) # 1071
sum(duplicated((upload %>% filter(Platekey %in% FF_list$SAMPLE_WELL_ID) %>% pull(Platekey))))  # 1 duplicated (LP3000354-DNA_C01)
dup <- upload[upload$Platekey %in% FF_list$SAMPLE_WELL_ID,][duplicated(upload[upload$Platekey %in% FF_list$SAMPLE_WELL_ID,]$Platekey),]$Platekey
upload %>% filter(Platekey %in% dup)  
  
# Excluding sample LP3000354-DNA_C01 (two V4 entries in upload report, one says it's germline, other tumour, strange delivery IDs: CF02206571, CF02206572)  
FF_list <- FF_list %>% filter(SAMPLE_WELL_ID != "LP3000354-DNA_C01")
# Check
dim(FF_list)  # 1069
dim(upload %>% filter(Platekey %in% FF_list$SAMPLE_WELL_ID)) # 1069
table((upload %>% filter(Platekey %in% FF_list$SAMPLE_WELL_ID) %>% pull(Status)), (upload %>% filter(Platekey %in% FF_list$SAMPLE_WELL_ID) %>% pull(`Delivery Version`)), exclude = NULL)

# Exclude 8 QC failed samples (upload report status)
ff_failed <- upload %>% filter(Platekey %in% FF_list$SAMPLE_WELL_ID, Status != "qc_passed") %>% pull(Platekey)
FF_list <- FF_list %>% filter(!SAMPLE_WELL_ID %in% ff_failed)
# Check
dim(FF_list)  # 1061
dim(upload %>% filter(Platekey %in% FF_list$SAMPLE_WELL_ID)) # 1061

# Add upload report info to the FF_list
FF_list <- left_join(FF_list, upload, by = c("SAMPLE_WELL_ID" = "Platekey"))

# Sanity check
sum(duplicated(FF_list$SAMPLE_WELL_ID))  # 0
summary(FF_list)

# Write the sample list (to be transferred to HPC)
write.csv(FF_list, file = "./Data/samples_for_recurrent_calc_FF_all.csv", quote = F, row.names = F)



##### Get VCF paths for normalization (run on LSF)

# Read the table of FF samples for reporting/recurrent variant calculation
samples_for_recurr <- read.csv("/home/mmijuskovic/small_variant_freq/FF/samples_for_recurrent_calc_FF_all.csv", header = T)

# Get VCF paths (note that some were normalized before)
paths <- unlist(sapply(samples_for_recurr$Path, function(x){
  command <- paste("find", paste0(x, "/SomaticVariations"), "-iname *.somatic.vcf.gz", sep = " ")
  system(command, intern = T)
}))
paths <- paths[!grepl("PASS.duprem.atomic.left.split.somatic.vcf.gz", paths)]
samples_for_recurr$SNV_VCF_path <- as.character(paths)

# Write the table with VCF paths
write.csv(samples_for_recurr, file = "/home/mmijuskovic/small_variant_freq/FF/samples_for_recurrent_calc_FF_all_withVCFpaths.csv", quote = F, row.names = F)

# Write sample names and VCF paths into the same table
write.table((samples_for_recurr %>% select(SAMPLE_WELL_ID, SNV_VCF_path)), file = "/home/mmijuskovic/small_variant_freq/FF/samples_for_recurrent_calc_FF_all_normalization_input.csv", sep = ",", col.names = F, quote = F, row.names = F)

# Split input table into 35 files (for 35 parallel jobs, 30 samples each)
sapply(0:34, function(x){
  # Get sample index
  i <- (x*30)+1
  j <- i+29
  z <- x+1
  # Write input table containing (next) 30 samples
  tabl <- (samples_for_recurr %>% select(SAMPLE_WELL_ID, SNV_VCF_path))[i:j,]
  write.table(tabl, file = paste0("/home/mmijuskovic/small_variant_freq/FF/norm_VCFs_input/samples_for_recurrent_calc_normalization_input_", z, ".csv"), sep = ",", col.names = F, quote = F, row.names = F)
  
})

# Write the final table of 11 samples as input table No 36 for the job #36
i <- 1051
j <- 1061
z <- 36
tabl_fin <- (samples_for_recurr %>% select(SAMPLE_WELL_ID, SNV_VCF_path))[i:j,]
write.table(tabl_fin, file = paste0("/home/mmijuskovic/small_variant_freq/FF/norm_VCFs_input/samples_for_recurrent_calc_normalization_input_", z, ".csv"), sep = ",", col.names = F, quote = F, row.names = F)






##### Get VCF paths for adding GT on normalized VCFs (run on LSF)

# Read the table of FFPE samples for reporting/recurrent variant calculation (includes original VCF paths)
samples_for_recurr <- read.csv("/home/mmijuskovic/small_variant_freq/FF/samples_for_recurrent_calc_FF_all_withVCFpaths.csv", header = T)

# Add normalized VCF paths
samples_for_recurr$norm_SNV_VCF_path <- paste0("/home/mmijuskovic/small_variant_freq/FF/norm_VCFs/", samples_for_recurr$SAMPLE_WELL_ID, ".somatic.duprem.left.split.vcf.gz")

# Write out the table of normalized VCF paths
write.csv(samples_for_recurr, file = "/home/mmijuskovic/small_variant_freq/FF/samples_for_recurrent_calc_FF_all_withNormVCFpaths.csv", quote = F, row.names = F)

# Split input table into 35 files (for 8 parallel jobs, 29 samples each)
sapply(0:34, function(x){
  # Get sample index
  i <- (x*30)+1
  j <- i+29
  z <- x+1
  # Write input table containing (next) 30 samples
  tabl <- (samples_for_recurr %>% select(SAMPLE_WELL_ID, norm_SNV_VCF_path))[i:j,]
  write.table(tabl, file = paste0("/home/mmijuskovic/small_variant_freq/FF/addGT_input/samples_for_recurrent_calc_addGT_input_", z, ".csv"), sep = ",", col.names = F, quote = F, row.names = F)
  
})

# Write the final table of 11 samples as input table No 36 for the job #36
i <- 1051
j <- 1061
z <- 36
tabl_fin <- (samples_for_recurr %>% select(SAMPLE_WELL_ID, norm_SNV_VCF_path))[i:j,]
write.table(tabl_fin, file = paste0("/home/mmijuskovic/small_variant_freq/FF/addGT_input/samples_for_recurrent_calc_addGT_input_", z, ".csv"), sep = ",", col.names = F, quote = F, row.names = F)




##### Get VCF paths for merging somatic VCFs (all)

vcf_list <- paste0("/home/mmijuskovic/small_variant_freq/FF/GT_VCFs/", samples_for_recurr$SAMPLE_WELL_ID, ".GT.duprem.left.split.vcf.gz")
#write.table(vcf_list, file = "./Data/FF_all_mainProgram_2017.txt", col.names = F, quote = F, row.names = F)
write.table(vcf_list, file = "/home/mmijuskovic/small_variant_freq/FF/FF_all_mainProgram_2017.txt", col.names = F, quote = F, row.names = F)


##### Get VCF paths for merging somatic VCFs (by library prep)

table(FF_list$LIBRARY_TYPE, exclude = NULL)

vcf_list_nano <- paste0("/home/mmijuskovic/small_variant_freq/FF/GT_VCFs/", FF_list[FF_list$LIBRARY_TYPE == "TruSeq Nano",]$SAMPLE_WELL_ID, ".GT.duprem.left.split.vcf.gz")
write.table(vcf_list_nano, file = "./Data/FF_nano_mainProgram_2017.txt", col.names = F, quote = F, row.names = F) 
  
vcf_list_PCRfree <- paste0("/home/mmijuskovic/small_variant_freq/FF/GT_VCFs/", FF_list[FF_list$LIBRARY_TYPE == "TruSeq PCR-Free",]$SAMPLE_WELL_ID, ".GT.duprem.left.split.vcf.gz")
write.table(vcf_list_PCRfree, file = "./Data/FF_PCRfree_mainProgram_2017.txt", col.names = F, quote = F, row.names = F)


### Write FF PCR-free paths for SV tiering

temp <- FF_list %>% filter(LIBRARY_TYPE == "TruSeq PCR-Free", !TumourType %in% "UNKNOWN", !is.na(TumourType)) %>% select(SAMPLE_WELL_ID, Path, TumourType)
write.table(temp, file = "/Users/MartinaMijuskovic/cancer_SV_pipeline_dev/Data/FF_PCRfree_wClinical_forSVtiering.tsv", sep = "\t", row.names = F, col.names = T, quote = F)



##### Check for low quality samples ##### 

# Jan 5 2018

### FF LOW QUALITY

# List of low quality FF samples (detected by PCA, reviewed or to-be-reviewed manually)
low_qual <- read.csv("./Data/Low_qual_cancer_samples.csv")
dim(FF_list[FF_list$SAMPLE_WELL_ID %in% low_qual$WELL.ID,])  # 51 of FF (13 of these are nano)
# Add low qual flag
FF_list$LOW_QUAL <- 0
FF_list[FF_list$SAMPLE_WELL_ID %in% low_qual$WELL.ID,]$LOW_QUAL <- 1
# Review (nothing major discovered, see notebook)
FF_list %>% filter(LOW_QUAL ==1) %>% dplyr::select(SAMPLE_WELL_ID, CENTER_CODE, SAMPLE_TYPE, LIBRARY_TYPE, TUMOUR_PURITY, `BAM Date`, TumourType, AT_DROP, COVERAGE_HOMOGENEITY, CHIMERIC_PER, AV_FRAGMENT_SIZE_BP, MAPPING_RATE_PER, SNV, INDEL)
# Check outcome (2 nano and 5 PCR-free samples need to be removed and VF re-calculated for production)
low_qual %>% filter(WELL.ID %in% FF_list[FF_list$LOW_QUAL == 1,]$SAMPLE_WELL_ID) %>% dplyr::select(WELL.ID, CENTER, LIBRARY.PREP, PROBLEM.REPORTED.TO.GMC, COMMENTS, STATE)
to_remove_nano <- low_qual %>% filter(WELL.ID %in% FF_list[FF_list$LOW_QUAL == 1,]$SAMPLE_WELL_ID, STATE == "FAILED", LIBRARY.PREP == "Nano") %>% pull(WELL.ID)  # 8 (6%)
to_remove_PCRfree <- low_qual %>% filter(WELL.ID %in% FF_list[FF_list$LOW_QUAL == 1,]$SAMPLE_WELL_ID, STATE == "FAILED", LIBRARY.PREP == "PCR-Free") %>% pull(WELL.ID)  # 14 (1.5%)
  
# Flag FF samples to be removed
FF_list$TO_REMOVE <- 0
FF_list[FF_list$SAMPLE_WELL_ID %in% c(as.character(to_remove_nano), as.character(to_remove_PCRfree)),]$TO_REMOVE <- 1

# Review samples to be removed
low_qual %>% filter(WELL.ID %in% c(as.character(to_remove_nano), as.character(to_remove_PCRfree))) %>% dplyr::select(WELL.ID, CENTER, LIBRARY.PREP, AT.DROP, COVERAGE.HOMOGENEITY, CHIMERIC.PER, AV.FRAGMENT.SIZE.BP, DESEASE.TYPE, PROBLEM.REPORTED.TO.GMC,COMMENTS,STATE)

# Write the clean list of FF samples
write.table(FF_list, file = "./Data/Clean_FF_sample_list_Jan2018.csv", quote = F, row.names = F, col.names = T, sep = ",")

# Write the list of samples that need to be removed from FF nano and FF PCR-free
write.table((FF_list %>% filter(LIBRARY_TYPE == "TruSeq Nano", TO_REMOVE == 1) %>% dplyr::select(SAMPLE_WELL_ID)), file = "./Data/FFnano_to_remove.txt", sep = "/t", row.names = F, col.names = F, quote = F)
write.table((FF_list %>% filter(LIBRARY_TYPE == "TruSeq PCR-Free", TO_REMOVE == 1) %>% dplyr::select(SAMPLE_WELL_ID)), file = "./Data/FFPCRfree_to_remove.txt", sep = "/t", row.names = F, col.names = F, quote = F)



### CONTAMINTATION

# List of contaminated samples (ContEst, ConPair, VerifyBamID, Array concordance)
contamin <- read.csv("./Data/Contamin_cancer_samples.csv")

# Check for samples in contaminated list (NO ISSUE FOUND)
dim(FF_list)  # 1061
table(FF_list$LIBRARY_TYPE)
FF_list[FF_list$SAMPLE_WELL_ID %in% contamin$WELL.ID,]  # none of FF
dim(FFPE_list)  # 245
FFPE_list_clean <- FFPE_list %>% filter(Exclusion_reason == "")
dim(FFPE_list_clean)
FFPE_list_clean[FFPE_list_clean$Platekey %in% contamin$WELL.ID,]
FFPE_contamin <- FFPE_list_clean[FFPE_list_clean$Platekey %in% contamin$WELL.ID,]$Platekey
contamin %>% filter(WELL.ID %in% FFPE_contamin)  # Low contamination of 3 FFPE samples, ok to keep


########## Read merged VCFs (PASS, coding regions) ########## 

###### FFPE

# Read merged VCF
ffpe_vcf <- readVcf("/Users/MartinaMijuskovic/FFPE_reporting_phase1/Data/mergedVCFs/FFPE_mainProgram_2017.PASSonly.merged.AF.info_only.coding_only.sorted.vcf.gz")
ffpe_info <- as.data.frame(info(ffpe_vcf))

# Add variant info
ffpe_info$ID <- names(ranges(ffpe_vcf))
ffpe_info$CHR <- as.character(seqnames(ffpe_vcf))
ffpe_info$POS <- ffpe_vcf@rowRanges@ranges@start
ffpe_info$REF <- as.data.frame(ref(ffpe_vcf))$x
ffpe_info$ALT <- as.data.frame(alt(ffpe_vcf))$value

# Make key
ffpe_info$KEY <- sapply(1:dim(ffpe_info)[x], function(x){
  paste(ffpe_info$CHR[x], ffpe_info$POS[x], ffpe_info$REF[x], ffpe_info$ALT[x], sep = "_")
})

# Cleanup
rm(ffpe_vcf)



###### FF

# Read merged VCF (PCRfree)
ff_PCRfree_vcf <- readVcf("/Users/MartinaMijuskovic/FFPE_reporting_phase1/Data/mergedVCFs/FF_PCRfree_mainProgram_2017.PASSonly.merged.AF.info_only.coding_only.sorted.vcf.gz")
ff_nano_vcf  <- readVcf("/Users/MartinaMijuskovic/FFPE_reporting_phase1/Data/mergedVCFs/FF_nano_mainProgram_2017.PASSonly.merged.AF.info_only.coding_only.sorted.vcf.gz")
ff_PCRfree_info <-  as.data.frame(info(ff_PCRfree_vcf))
ff_nano_info <-  as.data.frame(info(ff_nano_vcf))


# Add variant info
ff_PCRfree_info$ID <- names(ranges(ff_PCRfree_vcf))
ff_PCRfree_info$CHR <- as.character(seqnames(ff_PCRfree_vcf))
ff_PCRfree_info$POS <- ff_PCRfree_vcf@rowRanges@ranges@start
ff_PCRfree_info$REF <- as.data.frame(ref(ff_PCRfree_vcf))$x
ff_PCRfree_info$ALT <- as.data.frame(alt(ff_PCRfree_vcf))$value

ff_nano_info$ID <- names(ranges(ff_nano_vcf))
ff_nano_info$CHR <- as.character(seqnames(ff_nano_vcf))
ff_nano_info$POS <- ff_nano_vcf@rowRanges@ranges@start
ff_nano_info$REF <- as.data.frame(ref(ff_nano_vcf))$x
ff_nano_info$ALT <- as.data.frame(alt(ff_nano_vcf))$value

# Make key
ff_PCRfree_info$KEY <- sapply(1:dim(ff_PCRfree_info)[x], function(x){
  paste(ff_PCRfree_info$CHR[x], ff_PCRfree_info$POS[x], ff_PCRfree_info$REF[x], ff_PCRfree_info$ALT[x], sep = "_")
})
ff_nano_info$KEY <- sapply(1:dim(ff_nano_info)[x], function(x){
  paste(ff_nano_info$CHR[x], ff_nano_info$POS[x], ff_nano_info$REF[x], ff_nano_info$ALT[x], sep = "_")
})

# Cleanup
rm(ff_PCRfree_vcf,ff_nano_vcf)



####### Sanity check

# Number of variants corresponding to expected (compared to vt peek)? -ok
dim(ffpe_info) # 586660
dim(ff_nano_info) # 279270
dim(ff_PCRfree_info) # 1659842

# Duplicate keys? (none expected)
sum(duplicated(ffpe_info$KEY))  # 216460
sum(duplicated(ffpe_info)) # 216460 (all duplicates also have full INFO duplicated)
length(unique(ffpe_info$KEY)) # 370200
sum(duplicated(ff_nano_info$KEY))  # 101689
sum(duplicated(ff_PCRfree_info$KEY))  # 636164


# Look at examples of duplicated variants (artefacts created during coding region subset; not present in the previous step)
ffpe_dups <- ffpe_info$KEY[duplicated(ffpe_info$KEY)]
head(ffpe_info %>% filter(KEY %in% ffpe_dups))

# Check if all duplicates are exact copies of the same row
sum(duplicated(ffpe_info)) # should be 216460

# Removing duplicates
ffpe_info <- ffpe_info[!duplicated(ffpe_info$KEY),]
ff_nano_info <- ff_nano_info[!duplicated(ff_nano_info$KEY),]
ff_PCRfree_info <- ff_PCRfree_info[!duplicated(ff_PCRfree_info$KEY),]

# Check
dim(ffpe_info) # 370200
dim(ff_nano_info) # 177581
dim(ff_PCRfree_info) # 1023678




########## Variant frequency plots and summary ########## 

# Calculate VF (variant frequency) as AC/total samples (FFPE: 232, FF-nano: 136, FF-PCRfree: 925)
ffpe_info$VF <- as.numeric(ffpe_info$AC) / 232
ff_nano_info$VF <- as.numeric(ff_nano_info$AC) / 136
ff_PCRfree_info$VF <- as.numeric(ff_PCRfree_info$AC) / 925

# Plot VF (FFPE only)
ggplot(ffpe_info, aes(x=VF)) +
  geom_histogram(bins = 50) +
  #scale_x_log10() +
  scale_y_log10() +
  blank

# Put all data together to compare VFs
ffpe_info$group <- "FFPE"
ff_nano_info$group <- "FF TruSeq nano"
ff_PCRfree_info$group <- "FF TruSeq PCRfree"
all_coding <- rbind(ffpe_info, ff_nano_info, ff_PCRfree_info)

# Plot VF (all)
# ggplot(all_coding, aes(x=VF, col = group)) +
#   geom_freqpoly(bins = 100) +
#   scale_x_log10(limits = c(0.0001,1)) +
#   scale_y_log10() +
#   #scale_x_continuous(limits = c(0,1)) +
#   blank

pdf(file = "./Plots/recurrent_variants/somatic_var_freq.pdf", width = 10, height = 5)
print(
ggplot(all_coding, aes(x=VF, col = group)) +
  geom_freqpoly(bins = 50) +
  #scale_x_log10(limits = c(0.0001,1)) +
  scale_y_log10() +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,by=0.1)) +
  labs(x = "Variant Frequency", y = "Variant Count (log)") +
  blank +
  bigger
)
dev.off()


### Distribution by variant type

# Add type (SNV, indel)
# Below is too slow...
# all_coding$VAR_TYPE <- sapply(1:dim(all_coding)[1], function(x){
#   if(nchar(all_coding$REF)[x]==1 & nchar(all_coding$ALT)[x]==1){ "SNV" }
#   else { "INDEL" }
# })

all_coding_recurr$VAR_TYPE <- sapply(1:dim(all_coding_recurr)[1], function(x){
  if(nchar(all_coding_recurr$REF)[x]==1 & nchar(all_coding_recurr$ALT)[x]==1){ "SNV" }
  else { "INDEL" }
})

all_coding$REF_length <- nchar(all_coding$REF)
all_coding$ALT_length <- nchar(all_coding$ALT)
all_coding$VAR_TYPE <- all_coding$REF_length + all_coding$ALT_length
all_coding$VAR_TYPE <- as.character(all_coding$VAR_TYPE)
all_coding[all_coding$VAR_TYPE != "2",]$VAR_TYPE <- "INDEL"
all_coding[all_coding$VAR_TYPE != "INDEL",]$VAR_TYPE <- "SNV"



# Summary of recurrent by variant type and group
table(all_coding_recurr$group, all_coding_recurr$VAR_TYPE)
table(all_coding$group, all_coding$VAR_TYPE)

# Plot frequency distribution by variant type
pdf(file = "./Plots/recurrent_variants/somatic_var_freq_varType_FFPE.pdf", width = 10, height = 5)
print(
  ggplot((all_coding %>% filter(group == "FFPE")), aes(x=VF, col = VAR_TYPE)) +
    geom_freqpoly(bins = 50) +
    #scale_x_log10(limits = c(0.0001,1)) +
    scale_y_log10() +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,by=0.1)) +
    labs(x = "Variant Frequency", y = "Variant Count (log)") +
    blank +
    bigger
)
dev.off()

pdf(file = "./Plots/recurrent_variants/somatic_var_freq_varType_FFpcrfree.pdf", width = 10, height = 5)
print(
  ggplot((all_coding %>% filter(group == "FF TruSeq PCRfree")), aes(x=VF, col = VAR_TYPE)) +
    geom_freqpoly(bins = 50) +
    #scale_x_log10(limits = c(0.0001,1)) +
    scale_y_log10() +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,by=0.1)) +
    labs(x = "Variant Frequency", y = "Variant Count (log)") +
    blank +
    bigger
)
dev.off()

pdf(file = "./Plots/recurrent_variants/somatic_var_freq_varType_FFnano.pdf", width = 10, height = 5)
print(
  ggplot((all_coding %>% filter(group == "FF TruSeq nano")), aes(x=VF, col = VAR_TYPE)) +
    geom_freqpoly(bins = 50) +
    #scale_x_log10(limits = c(0.0001,1)) +
    scale_y_log10() +
    scale_x_continuous(limits = c(0,1), breaks = seq(0,1,by=0.1)) +
    labs(x = "Variant Frequency", y = "Variant Count (log)") +
    blank +
    bigger
)
dev.off()


### Summary of coding variants by frequency bin

# Total variants per group
all_coding %>% group_by(group) %>% summarise(TOTAL_VARIANTS = n())

# Variants per group at VF>10%
table((all_coding$VF >= 0.1), all_coding$group)
# VF < 5%
table((all_coding$VF < 0.05), all_coding$group)


### Barplot (by group) of VF bins

# Add bins
# Add VAF bin to variants table (note that nano sample size is 136 so lowest VF is 1/136=0.7%)
all_coding$VF_BIN <- sapply(1:dim(all_coding)[1], function(x){
  if (all_coding$VF[x] >= 0.1) {">10%"}
  else if (all_coding$VF[x] >= 0.05) {"5-10%"}
  else if (all_coding$VF[x] >= 0.01) {"1-5%"}
  else {"<1%"}
})
table(all_coding$VF_BIN, all_coding$group)
table(all_coding$VF_BIN, all_coding$group, all_coding$simpleRepeat_overlap)

# Flag private variants
all_coding <- all_coding %>% mutate(private = case_when(group == "FFPE" & VF < 0.008 ~ 1,
                                                        group == "FF TruSeq nano" & VF < 0.01 ~ 1,
                                                        group == "FF TruSeq PCRfree" & VF < 0.002 ~ 1,
                                                        TRUE ~ 0))
table(all_coding$group, all_coding$private)
table(all_coding$private, all_coding$VF_BIN, all_coding$group)

# Make new variable for plotting
all_coding <- all_coding %>% mutate(`Variant Frequency` = case_when(private == 0 ~ VF_BIN, TRUE ~ "private"))
table(all_coding$group, all_coding$`Variant Frequency`)
table(all_coding$VAR_TYPE, all_coding$group)
table(all_coding$group, all_coding$`Variant Frequency`, all_coding$VAR_TYPE)


# Barplot by VF bins
all_coding_summary <- as.data.frame(table(all_coding$group, all_coding$`Variant Frequency`))
# Reorder the levels for plotting
# levels(all_coding_summary$Var2)
# levels(all_coding_summary$Var2) <- c("<1%" ,    ">10%" ,   "1-5%" ,   "5-10%"  , "private") 
all_coding_summary$Var2 <- factor(all_coding_summary$Var2, levels = rev(c("private", "<1%", "1-5%", "5-10%", ">10%")))

pdf(file = paste0("./Plots/recurrent_variants/withPrivate_byGroup.pdf"))
print(
ggplot(all_coding_summary, aes(x=Var1, y=Freq, fill = Var2))  +
  geom_bar(position = "fill",stat = "identity") + 
  scale_y_continuous(labels = percent_format()) + 
  bigger +
  tiltedX +
  labs(x="",y="") +
  blank
)
dev.off()

# Same barplot, but omitting private variants
all_coding_summary_nonPrivate <- as.data.frame(table(all_coding[all_coding$private == 0,]$group, all_coding[all_coding$private == 0,]$`Variant Frequency`))
# Reorder the levels for plotting
all_coding_summary_nonPrivate$Var2 <- factor(all_coding_summary_nonPrivate$Var2, levels = rev(c("<1%", "1-5%", "5-10%", ">10%")))

pdf(file = paste0("./Plots/recurrent_variants/nonPrivate_byGroup.pdf"))
print(
ggplot(all_coding_summary_nonPrivate, aes(x=Var1, y=Freq, fill = Var2))  +
  geom_bar(position = "fill",stat = "identity") + 
  scale_y_continuous(labels = percent_format()) + 
  bigger +
  tiltedX +
  labs(x="",y="") +
  blank
)
dev.off()


# Putting private and <1% together

all_coding_summary2 <- as.data.frame(table(all_coding$group, all_coding$VF_BIN))
all_coding_summary2$Var2 <- factor(all_coding_summary2$Var2, levels = rev(c("<1%", "1-5%", "5-10%", ">10%")))

pdf(file = paste0("./Plots/recurrent_variants/VF_bins_byGroup.pdf"))
print(
  ggplot(all_coding_summary2, aes(x=Var1, y=Freq, fill = Var2))  +
    geom_bar(position = "fill",stat = "identity") + 
    scale_y_continuous(labels = percent_format()) + 
    bigger +
    tiltedX +
    labs(x="",y="") +
    blank
)
dev.off()





########## Analysis by tumour type ########## 

# FFPE tumour type distribution
as.data.frame(table(FFPE_list[FFPE_list$Exclusion_reason == "",]$TUMOUR_TYPE, exclude = NULL))

# FF PCR-free tumour type distribution  (no tumour type available)
#as.data.frame(table(FF_list[FF_list$LIBRARY_TYPE == "TruSeq PCR-Free",]$TUMOUR_TYPE, exclude = NULL))

# FF nano tumour type distribution (no tumour type available)
#as.data.frame(table(FF_list[FF_list$LIBRARY_TYPE == "TruSeq Nano",]$TUMOUR_TYPE, exclude = NULL))




########## Recurrent variants summary and overlap ########## 

# Subset variants to VF >= 10%
dim(all_coding)
all_coding_recurr <- all_coding %>% filter(VF >= 0.1)
dim(all_coding_recurr)

# Recurrent variants by group
table(all_coding_recurr$group)

##### Overlap of recurrent variants
pcrfree_keys <- all_coding_recurr %>% filter(group == "FF TruSeq PCRfree") %>% pull(KEY)  # 3
nano_keys <- all_coding_recurr %>% filter(group == "FF TruSeq nano") %>% pull(KEY)  # 1405
ffpe_keys <- all_coding_recurr %>% filter(group == "FFPE") %>% pull(KEY)  # 1456


### FF PCR-free recurrent variants

sum(pcrfree_keys %in% nano_keys)  # 2
sum(pcrfree_keys %in% ffpe_keys)  # 0

# Recurrent variants in FF PCR-free
all_coding_recurr %>% filter(group == "FF TruSeq PCRfree") %>% dplyr::select(KEY, ID, VF, AF1000G, CSQT)
all_coding_recurr %>% filter(group == "FF TruSeq PCRfree", KEY %in% nano_keys) %>% dplyr::select(KEY, ID, VF, AF1000G, CSQT)


### Overlap of FFPE and FF nano recurrent variants

length(ffpe_keys)  # 1456
length(nano_keys) # 1405
sum(ffpe_keys %in% nano_keys) # 846 OVERLAP





### Recurrent variants observed in the germline

# Allele frequency in 1000G > 1%
table(all_coding_recurr[all_coding_recurr$AF1000G > 0.01,]$group, exclude = NULL)

nano_keys_leaks <- all_coding_recurr %>% filter(group == "FF TruSeq nano", AF1000G > 0.01) %>% pull(KEY) 
ffpe_keys_leaks <- all_coding_recurr %>% filter(group == "FFPE", AF1000G > 0.01) %>% pull(KEY)

sum(nano_keys_leaks %in% ffpe_keys_leaks )  # 115 overlap


### Overlap between recurrent variants in FFPE and FF nano
grid.newpage()
pdf(file = "./Plots/recurrent_variants/NanoVsFFPE_venn.pdf")
print(
  draw.pairwise.venn(length(nano_keys), length(ffpe_keys), sum(ffpe_keys %in% nano_keys), c("FF nano", "FFPE"), fill = c("hotpink", "turquoise3"), lwd = 1, cat.fontfamily = rep("ArialMT", 2), fontfamily = rep("ArialMT", 3), alpha = rep(0.4, 2), cex = rep(1.5, 3), cat.cex = rep(1.2, 2))
)
dev.off()
grid.newpage()



### Recurrent FFPE variants NOT observed in FF cohorts

# KEYs of FFPE recurrent variants (>5%)
recurr_FFPE <- all_coding %>% filter(group == "FFPE", VF >= 0.05) %>% pull(KEY)  # 5447
# NUmber of recurrent FFPE variants NOT in other cohorts  
sum(!recurr_FFPE %in% all_coding[all_coding$group != "FFPE",]$KEY)  # 103
recurr_FFPE_specific <- recurr_FFPE[!recurr_FFPE %in% all_coding[all_coding$group != "FFPE",]$KEY]
# List FFPE-specific recurrent variants
all_coding %>% filter(group == "FFPE", KEY %in% recurr_FFPE_specific) %>% dplyr::select(KEY, VF_BIN, rs_ID, simpleRepeat_overlap, VAR_TYPE, AF1000G)
# Various summaries of FFPE-specific recurrent variants
table(all_coding %>% filter(group == "FFPE", KEY %in% recurr_FFPE_specific) %>% pull(rs_ID))
table(all_coding %>% filter(group == "FFPE", KEY %in% recurr_FFPE_specific) %>% pull(simpleRepeat_overlap))
table(all_coding %>% filter(group == "FFPE", KEY %in% recurr_FFPE_specific) %>% pull(VAR_TYPE))



### How many recurrent FF nano and FFPE variants are not observed in FF PCR-free at all

# KEYs of FF nano recurrent variants (>10%)
recurr_FF_nano <- all_coding %>% filter(group == "FF TruSeq nano", VF >= 0.1) %>% pull(KEY)  # 1405

# How many of recurrent FF nano are observed in FF PCR free
sum(recurr_FF_nano %in% all_coding[all_coding$group == "FF TruSeq PCRfree",]$KEY)  # 1244 observed in PCR-free

# Frequency at which these are present in different cohorts
pdf(file = "./Plots/recurrent_variants/FF_nano_recurr_otherCohorts.pdf")
print(
ggplot((all_coding %>% filter(KEY %in% recurr_FF_nano)), aes(x=VF, fill = group)) +
  geom_histogram(alpha = 0.5, position = "identity") +
  blank +
  bigger
)
dev.off()

# KEYs of FFPE recurrent variants (>10%)
recurr_FFPE_2 <- all_coding %>% filter(group == "FFPE", VF >= 0.1) %>% pull(KEY)  # 1456

# How many of recurrent FF nano are observed in FF PCR free
sum(recurr_FFPE_2 %in% all_coding[all_coding$group == "FF TruSeq PCRfree",]$KEY)  # 1316 observed in PCR-free

# Frequency at which these are present in different cohorts
pdf(file = "./Plots/recurrent_variants/FFPE_recurr_otherCohorts.pdf")
print(
ggplot((all_coding %>% filter(KEY %in% recurr_FFPE_2)), aes(x=VF, fill = group)) +
  geom_histogram(alpha = 0.5, position = "identity") +
  blank +
  bigger
)
dev.off()


###### Comparison to FFPE trio analysis ###### 

### Is the variant observed as recurrent (>10%) in the FFPE trio analysis?

# Load FFPE trio recurrent variant table (non-unique, all entries there)
trio_recurr <- read.table("/Users/MartinaMijuskovic/FFPE_trio_analysis/Data/SNV/var_recurrent_FF_FFPE_BOTH.tsv", header = T, sep = "\t")
dim(trio_recurr)
head(trio_recurr)
table(trio_recurr$RECURR_TYPE)

# Subset to variants recurrent in FFPE 
trio_recurr_ffpe <- trio_recurr %>% filter(SAMPLE_TYPE == "FFPE")
trio_recurr_ffpe_summary <- lapply(unique(trio_recurr_ffpe$KEY), function(x){
  data.frame(KEY = x, NUM_OBS_FFPE = dim(trio_recurr_ffpe[trio_recurr_ffpe$KEY == x,])[1])
}
)
trio_recurr_ffpe_summary <- bind_rows(trio_recurr_ffpe_summary)

dim(trio_recurr_ffpe_summary) # 433
length(unique(trio_recurr_ffpe_summary$KEY))  # 433
dim(trio_recurr_ffpe_summary %>% filter(NUM_OBS_FFPE >= 6))  # 197  (in agreement to FFPE trio analysis)


# Sanity check (compare merged VCF variants to FFPE trio analysis)
recurr_trio_keys <- paste0("chr", (trio_recurr_ffpe_summary %>% filter(NUM_OBS_FFPE >= 6) %>% pull(KEY)))  # 197
full_ffpe_cohort_recurr_keys <- all_coding_recurr %>% filter(group == "FFPE", VF >= 0.1) %>% pull(KEY) # 1456
sum(recurr_trio_keys %in% full_ffpe_cohort_recurr_keys)  # 109/197 overlap (recurrent at 10%+ in full cohort and trios)
non_overl_ffpe_recurr_keys <- recurr_trio_keys[!recurr_trio_keys %in% full_ffpe_cohort_recurr_keys]  # 88 don't overlap
# Are FFPE trio recurrent keys that are not >10% in the full cohort observed in full cohort? If yes, what frequencies (5-10%)?
sum(non_overl_ffpe_recurr_keys %in% all_coding$KEY)  # 88/88 observed

# Summary of VF for non-overlapping variants in the full FFPE cohort
summary(all_coding %>% filter(group == "FFPE", KEY %in% non_overl_ffpe_recurr_keys) %>% pull(VF))

# Plot frequency distribution of FFPE trio recurrent variants that are not recurrent in full FFPE cohort
pdf(file = "./Plots/recurrent_variants/ffpe_trio_recurr_not_observedInFullCohort.pdf", width = 10, height = 5)
print(
ggplot(all_coding[all_coding$KEY %in% non_overl_ffpe_recurr_keys,], aes(x=VF, col = group)) +
  geom_freqpoly(bins = 30) +
  #scale_x_log10(limits = c(0.0001,1)) +
  #scale_y_log10() +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,by=0.1)) +
  #labs(x = "Variant Frequency", y = "Variant Count (log)") +
  labs(x = "Variant Frequency", y = "Variant Count") +
  blank +
  bigger
)
dev.off()









###### Overlap with simple repeats ###### 

### Create bed files (0-based pos) with all variants to calculate TRF simple repeat overlap

# Correct (direct) TRF overlap (for deletions, use bed file with deletion width)
# Change to bed 0-based format
all_coding$length <- sapply(1:dim(all_coding)[1], function(x) { max(all_coding$REF_length[x], all_coding$ALT_length[x]) - 1 })
all_coding$POS <- as.numeric(all_coding$POS)

all_coding <- all_coding %>% mutate(start_bed = case_when(length == 0 ~ (POS - 1), 
                                                    length != 0 ~ POS),
                              end_bed = case_when(length == 0 ~ POS,
                                                  length != 0 ~ (POS + length)))
all_coding_bed <- all_coding %>% dplyr::select(CHR, start_bed, end_bed, KEY)

write.table(all_coding_bed, file = "./Data/all_coding.bed", quote = F, row.names = F, col.names = F, sep = "\t")


# Call bedtools to find  overlaps with TRF simple repeats
system(paste("bedtools coverage -a /Users/MartinaMijuskovic/FFPE_reporting_phase1/Data/all_coding.bed -b /Users/MartinaMijuskovic/FFPE/simpleRepeat.hg38.bed > /Users/MartinaMijuskovic/FFPE_reporting_phase1/Data/all_coding_simpleRepeat_DIRECToverlap.bed"), intern = T)
simpleRepeat_overlap <- read.table("/Users/MartinaMijuskovic/FFPE_reporting_phase1/Data/all_coding_simpleRepeat_DIRECToverlap.bed", sep = "\t")
names(simpleRepeat_overlap) <- c("CHR", "START", "END", "KEY", "NumOverlap", "BPoverlap", "BPTotal", "PCT")
# Add new flag
simpleRepeat_keys <- unique(as.character(simpleRepeat_overlap %>% filter(NumOverlap != 0) %>% .$KEY))  # 76714 keys
all_coding$simpleRepeat_overlap <- 0
all_coding[(all_coding$KEY %in% simpleRepeat_keys),]$simpleRepeat_overlap <- 1
table(all_coding$simpleRepeat_overlap, exclude = NULL)


### Proportion of recurrent variants overlapping simple repeats

table(all_coding$VF > 0.1, all_coding$group, all_coding$simpleRepeat_overlap, exclude = NULL)  # too few >10% variants in FF PCR-free for comparison



# Correct factor order
table(all_coding$VF_BIN)
all_coding$VF_BIN <- factor(all_coding$VF_BIN, levels = c(">10%", "5-10%", "1-5%", "<1%"))
table(all_coding$VF_BIN)

table(all_coding$group)
all_coding$group <- factor(all_coding$group, levels = c("FF TruSeq PCRfree", "FF TruSeq nano", "FFPE"))
table(all_coding$group)

# Sort data
all_coding <- with(all_coding, all_coding[order(VF_BIN, simpleRepeat_overlap, group),])

# Simple repeat overlap summary
table(all_coding$VF_BIN, all_coding$simpleRepeat_overlap, all_coding$group)
#all_coding %>% group_by(group, VF_BIN, simpleRepeat_overlap) %>% summarise(n())
table(all_coding$VF > 0.01, all_coding$simpleRepeat_overlap, all_coding$group)


### Barplot colored by simple repeat overlap by VF bins, by group

pdf(file = "./Plots/recurrent_variants/simpleRep_overlap_bins.pdf", width = 10, height = 5)
print(
  ggplot(all_coding, aes(x=group, fill = factor(simpleRepeat_overlap))) +
    #geom_bar(stat = "identity") +
    geom_bar() +
    facet_grid(~VF_BIN) +
    #blank +
    tiltedX +
    labs(fill="") +
    bigger
  )  
dev.off()

# Barplot (excluding <1% variants)
pdf(file = "./Plots/recurrent_variants/simpleRep_overlap_bins_recurrOnly.pdf", width = 10, height = 5)
print(
  ggplot((all_coding %>% filter(VF_BIN != "<1%")), aes(x=group, fill = factor(simpleRepeat_overlap))) +
    #geom_bar(stat = "identity") +
    geom_bar() +
    facet_grid(~VF_BIN) +
    #blank +
    tiltedX +
    labs(fill="", x="") +
    bigger
)  
dev.off()


###### Tiering into Domains 1-3  ###### 

# Get the list of ENSEMBL transcripts with tiers (GO v1.5)
all_GO_tr <- read.table("./Data/transcripts_v1p5/GENOMONCOLOGY_SOLID_TUMOUR.tsv", header = T, sep = "\t") 

# Domain 1 transcripts
domain1_tr <- unique(all_GO_tr$transcript_ID)
length(domain1_tr) # 113

# Cancer census genes
cancer_census_tr <- read.table("./Data/transcripts_v1p5/CANCER_CENSUS_GENES.tsv", header = T, sep = "\t")

# Domain 2 transcripts
domain2_tr <- unique(cancer_census_tr$transcript_ID)
length(domain2_tr)  # 566



### Function to annotate small variants POS with ENSEMBL transcripts (assumes chr names have "chr")
annGenes_smallVar <- function(vcf_info, transcripts){
  
  require(ensembldb)
  require(dplyr)
  
  # Create number-only chromosome names
  vcf_info$CHR_ <- vcf_info$CHR
  vcf_info$CHR <- sub(vcf_info$CHR, pattern = "chr", replacement = "")
  
  txdb_pth <- '/Users/MartinaMijuskovic/cancer_SV_pipeline_dev/Homo_sapiens.GRCh38.84.sqlite'
  if (!file.exists(txdb_pth)) {
    #txdb_pth <- ensDbFromGtf(gtf="../Ensembl_db/Homo_sapiens.GRCh38.84.gtf.gz")
    print("No Homo_sapiens.GRCh38.84.sqlite found. See https://blog.liang2.tw/posts/2016/05/biocondutor-ensembl-reference/ for installation")
  }
  txdb <- EnsDb(txdb_pth)
  
  # Get transcript ranges from Ensembl transcript database
  transcripts_gr <- transcripts(txdb, filter=TxidFilter(transcripts))
  
  # Flag if start overlaps with any transcripts
  #vcf_info$all_transcript_related_start <- as.numeric(overlapsAny(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 1)), transcripts_gr))
  vcf_info$POS_in_transcr <- as.numeric(overlapsAny(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$POS, width = 1)), transcripts_gr))
  
  # # Flag if end overlaps with any transcripts
  # vcf_info$all_transcript_related_end <- sapply(1:dim(vcf_info)[1], function(x){
  #   if(is.na(vcf_info$END[x])){
  #     return(NA)
  #   } 
  #   else {
  #     if(overlapsAny(GRanges(seqnames=vcf_info$CHR[x], ranges=IRanges(start = vcf_info$END[x], width = 1)), transcripts_gr)) {
  #       return(1)
  #     }
  #     else {
  #       return(0)
  #     }
  #   }
  # })
  

  
  # Initiate gene annotations
  vcf_info$ann_start <- ""
  #vcf_info$ann_end  <- ""
  
  # Find exact overlaps of start and end position with transcripts
  #overlaps_start <- findOverlaps(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$START, width = 100000)), transcripts_gr)  # testing
  overlaps_start <- findOverlaps(GRanges(seqnames=vcf_info$CHR, ranges=IRanges(start = vcf_info$POS, width = 1)), transcripts_gr)
  #overlaps_end <- findOverlaps(GRanges(seqnames=vcf_info[!is.na(vcf_info$all_transcript_related_end),]$CHR, ranges=IRanges(start = vcf_info[!is.na(vcf_info$all_transcript_related_end),]$END, width = 1)), transcripts_gr)
  
  # Add start position overlaps to table
  if (dim(vcf_info[queryHits(overlaps_start),])[1] > 0) {
    vcf_info[queryHits(overlaps_start),]$ann_start <- sapply(seq(length(overlaps_start)), function(x){
      i <- queryHits(overlaps_start)[x]
      j <- subjectHits(overlaps_start)[x]
      vcf_info[i,]$ann_start <- transcripts_gr@elementMetadata$tx_name[j]
    })
  }
  
  # # Add end position overlaps to table --- FIXED
  # if (dim(vcf_info[queryHits(overlaps_end),])[1] > 0) {
  #   vcf_info[!is.na(vcf_info$all_transcript_related_end),][queryHits(overlaps_end),]$ann_end <- sapply(seq(length(overlaps_end)), function(x){
  #     i <- queryHits(overlaps_end)[x]
  #     j <- subjectHits(overlaps_end)[x]
  #     vcf_info[!is.na(vcf_info$all_transcript_related_end),][i,]$ann_end <- transcripts_gr@elementMetadata$tx_name[j]
  #   })
  # }
  
  return(vcf_info)  
}


# Annotate all with Domain 1
recurr_coding_dom1 <- annGenes_smallVar((all_coding %>% filter(VF_BIN != "<1%")), domain1_tr)
names(recurr_coding_dom1)[63] <- "Domain 1"
names(recurr_coding_dom1)[64] <- "Domain1_tr"
  
# Annotate all with Domain 2
recurr_coding_dom2 <- annGenes_smallVar(recurr_coding_dom1, domain2_tr)
recurr_coding_dom12 <- recurr_coding_dom2
rm(recurr_coding_dom1, recurr_coding_dom2)
names(recurr_coding_dom12)[65] <- "Domain 2"
names(recurr_coding_dom12)[66] <- "Domain6_tr"

# Create domain variable
recurr_coding_dom12$Domain <- ""
recurr_coding_dom12[recurr_coding_dom12$`Domain 1` == 1,]$Domain <- "Domain 1"
recurr_coding_dom12[recurr_coding_dom12$`Domain 2` == 1,]$Domain <- "Domain 2"

###### Analysis by domain (>1% VF only) ###### 

# Summary by group and domain
table(recurr_coding_dom12$Domain, recurr_coding_dom12$group)

# Summary by group, domain and VF
table(recurr_coding_dom12$Domain, recurr_coding_dom12$VF_BIN, recurr_coding_dom12$group)

# Summary by type
table(recurr_coding_dom12$Domain, recurr_coding_dom12$VAR_TYPE, recurr_coding_dom12$group)



### Flag if variant has an assigned rs number in the VCF ID field (given by Strelka)
# This means variant has been observed before in the germline (likely germline leaks or low quality variants)

all_coding$rs_ID <- grepl("rs", all_coding$ID)
recurr_coding_dom12$rs_ID <- grepl("rs", recurr_coding_dom12$ID)

# Recurrent variants by "rs" ID
table(recurr_coding_dom12$Domain, recurr_coding_dom12$group, recurr_coding_dom12$rs_ID)
table(recurr_coding_dom12$group, recurr_coding_dom12$rs_ID)

table(all_coding$group, all_coding$VF_BIN, all_coding$rs_ID)

# Manually examine PCR free 5-10% variants that don't overlap simple repeats not have the rs ID
recurr_coding_dom12 %>% filter(group == "FF TruSeq PCRfree", VF >= 0.05, simpleRepeat_overlap == 0, rs_ID == 0)  # 90




###### Summary by tumour type ###### 

# Add tumour types to FF list
tumour_type_data <- read.csv("../cancer_SV_pipeline_dev/Data/cancer_samples.csv", header = T)

# Check if all FF samples are there
sum(FF_list$SAMPLE_WELL_ID %in% tumour_type_data$sampleId)  # 696 there
sum(!FF_list$SAMPLE_WELL_ID %in% tumour_type_data$sampleId)  # 365 missing

missing_tumour_type <- FF_list$SAMPLE_WELL_ID[!FF_list$SAMPLE_WELL_ID %in% tumour_type_data$sampleId]

# Summary of FF tumour types (I have all FFPEs)
table((tumour_type_data %>% filter(sampleId %in% FF_list$SAMPLE_WELL_ID) %>% pull(tumorType)), exclude = NULL)

# Clean data (use defined dictionary for tumour types)
# Dictionary (cancer data model 4.0: ADULT_GLIOMA, BLADDER, BREAST, CARCINOMA_OF_UNKNOWN_PRIMARY, CHILDHOOD, COLORECTAL, ENDOMETRIAL_CARCINOMA, HAEMONC, HEPATOPANCREATOBILIARY, LUNG, MALIGNANT_MELANOMA, NASOPHARYNGEAL, ORAL_OROPHARYNGEAL, OVARIAN, PROSTATE, RENAL, SARCOMA, SINONASAL, TESTICULAR_GERM_CELL_TUMOURS, UPPER_GASTROINTESTINAL, NON_HODGKINS_B_CELL_LYMPHOMA_LOW_MOD_GRADE, CLASSICAL_HODGKINS, NODULAR_LYMPHOCYTE_PREDOMINANT_HODGKINS, T_CELL_LYMPHOMA
table(tumour_type_data$tumorType)

tumour_type_data <- tumour_type_data %>% mutate(TumourType = case_when(
  tumorType == "Adult Glioma" ~ "ADULT_GLIOMA",
  tumorType == "adult_glioma" ~ "ADULT_GLIOMA",
  tumorType == "bladder" ~ "BLADDER",
  tumorType == "Bladder" ~ "BLADDER",
  tumorType == "Breast" ~ "BREAST",
  tumorType == "breast" ~ "BREAST",
  tumorType == "colorectal" ~ "COLORECTAL",
  tumorType == "Colorectal" ~ "COLORECTAL",
  tumorType == "Endometrial Carcinoma" ~ "ENDOMETRIAL_CARCINOMA",
  tumorType == "Hepatopancreatobiliary" ~ "HEPATOPANCREATOBILIARY",
  tumorType == "lung" ~ "LUNG",
  tumorType == "Lung" ~ "LUNG",
  tumorType == "Malignant Melanoma" ~ "MALIGNANT_MELANOMA",
  tumorType == "malignant_melanoma" ~ "MALIGNANT_MELANOMA",
  tumorType == "ovarian" ~ "OVARIAN",
  tumorType == "Ovarian" ~ "OVARIAN",
  tumorType == "Prostate" ~ "PROSTATE",
  tumorType == "renal" ~ "RENAL",
  tumorType == "Renal" ~ "RENAL",
  tumorType == "Sarcoma" ~ "SARCOMA",
  tumorType == "Testicular Germ Cell Tumours" ~ "TESTICULAR_GERM_CELL_TUMOURS",
  tumorType == "Upper Gastrointestinal" ~ "UPPER_GASTROINTESTINAL",
  tumorType == "Unknown" ~ "UNKNOWN",
  tumorType == "" ~ "UNKNOWN"
  
))

# Check
table(tumour_type_data$tumorType, tumour_type_data$TumourType)


### Get tumour types from catalog

getTumourType <- function(SAMPLE_WELL_ID, sessionID, studyID="1000000038"){
  require(jsonlite)
  require(dplyr)
  #                  curl -X GET --header "Accept: application/json" --header "Authorization: Bearer " "https://opencgainternal.gel.zone/opencga/webservices/rest/v1/samples/LP2000907-DNA_A01/annotationsets?sid=eyJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJtbWlqdXNrb3ZpYyIsImF1ZCI6Ik9wZW5DR0EgdXNlcnMiLCJpYXQiOjE1MTMwOTg0NTYsImV4cCI6MTUxMzEwMDI1Nn0.JxLrfOwG7VzJQFfGkUGZ8d6DCrig94eGq5oHNl-RYms&study=1000000038&asMap=false"
  command <- paste0('curl -X GET --header "Accept: application/json" --header "Authorization: Bearer " "https://opencgainternal.gel.zone/opencga/webservices/rest/v1/samples/', SAMPLE_WELL_ID, '/annotationsets?sid=', sessionID, '&study=', studyID, '&asMap=false"')
  
  info_json <- fromJSON(system(command, intern = T), flatten = T)
  
  if (info_json$response$numResults == 0) {
    result <- data.frame(
      WELL_ID = SAMPLE_WELL_ID,
      TumourType = NA
    )
    return(result)  
  }
  
  else if (length(info_json$response$result[[1]]) == 0 ){
    result <- data.frame(
      WELL_ID = SAMPLE_WELL_ID,
      TumourType = NA
    )
    return(result)  
  }
  
  else {
    result <- data.frame(
      WELL_ID = SAMPLE_WELL_ID,
      TumourType = if(dim(info_json$response$result[[1]]$annotations[[1]][info_json$response$result[[1]]$annotations[[1]]$name == "diseaseType",])[1] == 0){
        NA
      }
                  else{
                    info_json$response$result[[1]]$annotations[[1]][info_json$response$result[[1]]$annotations[[1]]$name == "diseaseType",]$value
                  }
    )
    return(result)  
  }
  
}

sessionID <- "eyJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJtbWlqdXNrb3ZpYyIsImF1ZCI6Ik9wZW5DR0EgdXNlcnMiLCJpYXQiOjE1MTMxMDA0MDcsImV4cCI6MTUxMzEwMjIwN30.DZiaG_ye9BWkdvUgHGNsHTBrJ-5XzIjvYxvEkB87Vl0"

# Test
#getTumourType("LP2000907-DNA_A01", sessionID)  # Works!

FF_tumour_types_catalog <- lapply(FF_list$SAMPLE_WELL_ID, getTumourType, sessionID)
FF_tumour_types_catalog <- as.data.frame(rbind_list(FF_tumour_types_catalog))

# Check results
table(FF_tumour_types_catalog$TumourType, exclude = NULL)  # 634 NA

# How many of missing ones are in catalog
sum(missing_tumour_type %in% FF_tumour_types_catalog[!is.na(FF_tumour_types_catalog$TumourType),]$WELL_ID)  # 161 additional ones found in catalog


# Test study *00036
FF_tumour_types_catalog_s36 <- lapply(FF_list$SAMPLE_WELL_ID, getTumourType, sessionID, studyID="1000000036")
FF_tumour_types_catalog_s36 <- as.data.frame(rbind_list(FF_tumour_types_catalog_s36))
table(FF_tumour_types_catalog_s36$TumourType, exclude = NULL)  # 1061 NA


### Add tumour types to FF list

# From Alona's table
FF_list$TumourType <- tumour_type_data[match(FF_list$SAMPLE_WELL_ID, tumour_type_data$sampleId),]$TumourType
table(FF_list$TumourType, exclude = NULL)  # 365 NAs
# From catalog (study 38)
sum(FF_list[is.na(FF_list$TumourType),]$SAMPLE_WELL_ID %in% FF_tumour_types_catalog[!is.na(FF_tumour_types_catalog$TumourType),]$WELL_ID) # 161 NAs can be added
FF_list$TumourType_cat <- FF_tumour_types_catalog[match(FF_list$SAMPLE_WELL_ID, FF_tumour_types_catalog$WELL_ID),]$TumourType
# Check
table(FF_list$TumourType, FF_list$TumourType_cat, exclude = NULL)
# Consolidate
FF_list[is.na(FF_list$TumourType),]$TumourType <- FF_list[is.na(FF_list$TumourType),]$TumourType_cat
table(FF_list$TumourType, FF_list$TumourType_cat, exclude = NULL)

### PCT of tumour types per group

# Clean FFPE tumour types
table(FFPE_list$TUMOUR_TYPE, exclude = NULL)
FFPE_list <- FFPE_list %>% mutate(TumourType = case_when(
  TUMOUR_TYPE == "Adult Glioma" ~ "ADULT_GLIOMA",
  TUMOUR_TYPE == "adult_glioma" ~ "ADULT_GLIOMA",
  TUMOUR_TYPE == "bladder" ~ "BLADDER",
  TUMOUR_TYPE == "Bladder" ~ "BLADDER",
  TUMOUR_TYPE == "Breast" ~ "BREAST",
  TUMOUR_TYPE == "breast" ~ "BREAST",
  TUMOUR_TYPE == "colorectal" ~ "COLORECTAL",
  TUMOUR_TYPE == "Colorectal" ~ "COLORECTAL",
  TUMOUR_TYPE == "Endometrial Carcinoma" ~ "ENDOMETRIAL_CARCINOMA",
  TUMOUR_TYPE == "Hepatopancreatobiliary" ~ "HEPATOPANCREATOBILIARY",
  TUMOUR_TYPE == "lung" ~ "LUNG",
  TUMOUR_TYPE == "Lung" ~ "LUNG",
  TUMOUR_TYPE == "Malignant Melanoma" ~ "MALIGNANT_MELANOMA",
  TUMOUR_TYPE == "malignant_melanoma" ~ "MALIGNANT_MELANOMA",
  TUMOUR_TYPE == "ovarian" ~ "OVARIAN",
  TUMOUR_TYPE == "Ovarian" ~ "OVARIAN",
  TUMOUR_TYPE == "Prostate" ~ "PROSTATE",
  TUMOUR_TYPE == "renal" ~ "RENAL",
  TUMOUR_TYPE == "Renal" ~ "RENAL",
  TUMOUR_TYPE == "Sarcoma" ~ "SARCOMA",
  TUMOUR_TYPE == "Testicular Germ Cell Tumours" ~ "TESTICULAR_GERM_CELL_TUMOURS",
  TUMOUR_TYPE == "Upper Gastrointestinal" ~ "UPPER_GASTROINTESTINAL",
  TUMOUR_TYPE == "Unknown" ~ "UNKNOWN",
  TUMOUR_TYPE == "" ~ "UNKNOWN"
  
))
# Check
table(FFPE_list$TumourType, FFPE_list$TUMOUR_TYPE, exclude = NULL)

# Put data together
all_cancer_cohorts_byTumourType <- rbind((FF_list %>% dplyr::select(SAMPLE_WELL_ID, LIBRARY_TYPE, TumourType)), (FFPE_list %>% mutate(SAMPLE_WELL_ID = Platekey, LIBRARY_TYPE = "FFPE") %>% dplyr::select(SAMPLE_WELL_ID, LIBRARY_TYPE, TumourType)))
table(all_cancer_cohorts_byTumourType$LIBRARY_TYPE, all_cancer_cohorts_byTumourType$TumourType, exclude = NULL)

# 11 FFPE tumour types are unknown, checking in catalog
sessionID <- "eyJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJtbWlqdXNrb3ZpYyIsImF1ZCI6Ik9wZW5DR0EgdXNlcnMiLCJpYXQiOjE1MTMxNzI2MDcsImV4cCI6MTUxMzE3NDQwN30.wYt6SErQcendjzPzun2ufWQ0zf-csXH7fQKFIEk3rqA"
FFPE_list_TumourType_catalog <- lapply(FFPE_list[FFPE_list$TumourType == "UNKNOWN",]$Platekey, getTumourType, sessionID)
FFPE_list_TumourType_catalog <- as.data.frame(rbind_all(FFPE_list_TumourType_catalog))
table(FFPE_list_TumourType_catalog$TumourType, exclude = NULL)  # 4 found in catalog

# Adding FFPE tumour types found in catalog
FFPE_list[FFPE_list$TumourType == "UNKNOWN",]$TumourType <- FFPE_list_TumourType_catalog[match(FFPE_list[FFPE_list$TumourType == "UNKNOWN",]$Platekey, FFPE_list_TumourType_catalog$WELL_ID),]$TumourType
table(FFPE_list$TumourType, exclude = NULL)

# Redo full table
all_cancer_cohorts_byTumourType <- rbind((FF_list %>% dplyr::select(SAMPLE_WELL_ID, LIBRARY_TYPE, TumourType)), (FFPE_list %>% mutate(SAMPLE_WELL_ID = Platekey, LIBRARY_TYPE = "FFPE") %>% dplyr::select(SAMPLE_WELL_ID, LIBRARY_TYPE, TumourType)))

# Summarize
as.data.frame(table(all_cancer_cohorts_byTumourType$LIBRARY_TYPE, all_cancer_cohorts_byTumourType$TumourType, exclude = NULL))
all_cancer_cohorts_byTumourType[is.na(all_cancer_cohorts_byTumourType$TumourType),]$TumourType <- "UNKNOWN"
all_cancer_cohorts_byTumourType$LIBRARY_TYPE <- as.character(all_cancer_cohorts_byTumourType$LIBRARY_TYPE)
TumourType_summary <- as.data.frame(table(all_cancer_cohorts_byTumourType$LIBRARY_TYPE, all_cancer_cohorts_byTumourType$TumourType))
TumourType_summary <- TumourType_summary %>% mutate(PCT_TOTAL = case_when(
  Var1 == "FFPE" ~ Freq/232*100,
  Var1 == "TruSeq Nano" ~ Freq/136*100,
  Var1 == "TruSeq PCR-Free" ~ Freq/925*100
))


# Plot proportion of each tumour type in respective cohorts
# Reorder groups
TumourType_summary$Var1 <- factor(TumourType_summary$Var1, levels = c("TruSeq Nano", "TruSeq PCR-Free", "FFPE"))

pdf(file = "./Plots/recurrent_variants/byTumourType.pdf", width = 10, height = 5)
print(
ggplot(TumourType_summary, aes(x=Var2, y=PCT_TOTAL, fill = Var1)) +
  geom_col(position = 'dodge') +
  tiltedX +
  labs(x="", y="% total") +
  blank
)
dev.off()







