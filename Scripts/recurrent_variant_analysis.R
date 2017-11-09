# Martina Mijuskovic
# FFPE project - small variant reporting
# Calculating and analysing recurrent FFPE variants (main program samples)
# Sept 2017

library(dplyr)

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
table(FF_list$CENTER_CODE, exclude = NULL)
sum(duplicated(FF_list$SAMPLE_WELL_ID))  # 0

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




##### Get VCF paths for merging somatic VCFs

vcf_list <- paste0("/home/mmijuskovic/small_variant_freq/FF/GT_VCFs/", samples_for_recurr$SAMPLE_WELL_ID, ".GT.duprem.left.split.vcf.gz")
#write.table(vcf_list, file = "./Data/FF_all_mainProgram_2017.txt", col.names = F, quote = F, row.names = F)
write.table(vcf_list, file = "/home/mmijuskovic/small_variant_freq/FF/FF_all_mainProgram_2017.txt", col.names = F, quote = F, row.names = F)





