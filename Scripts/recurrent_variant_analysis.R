# Martina Mijuskovic
# FFPE project - small variant reporting
# Calculating and analysing recurrent FFPE variants (main program samples)
# Sept 2017

library(dplyr)
library(VariantAnnotation)
library(ggplot2)

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



########## Variant overlap and frequency plots ########## 

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








