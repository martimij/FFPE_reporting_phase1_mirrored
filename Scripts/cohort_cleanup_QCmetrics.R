# Martina Mijuskovic
# FFPE project - small variant reporting
# Cohort cleanup and QC metrics
# Sept 2017

library(dplyr)
library(ggplot2)
library(data.table)
library(jsonlite)

today <- Sys.Date()

### Helper objects for plotting
# Blank theme
blank <-  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.title=element_blank())
# Black regression line (linear model)
regr_line <- geom_smooth(method = "lm", se = F, aes(group = 1), linetype = 2, col = "black", size = 0.5)
# Bigger plot text
bigger <- theme(legend.text=element_text(size=15), legend.title = element_text(size=15), axis.title = element_text(size=15), axis.text = element_text(size=15))
# Tilted x-axis labels
tiltedX <- theme(axis.text.x=element_text(angle=45,hjust=1))



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
# Write table with tumour type missing data
write.table((QC %>% filter(TUMOUR_TYPE == "N/A") %>% select(PATIENT_ID, WELL_ID, TUMOUR_TYPE, COLLECTION_DATE)), file = "./Data/FFPE_missing_tumourType_collDate.tsv", sep = "\t", quote = F, col.names = T, row.names = F)


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


#####  FFPE cohort list ##### 

### Total of 234 samples that have no QC or other issues (confirmation pending Bertha QC intake re-run and VerifyBamID conntamination check for germline)

# Add QC data to FFPE list
FFPE_list <- full_join(FFPE_list, QC, by = c("Platekey" = "WELL_ID"))
dim(FFPE_list) # 234
write.csv(FFPE_list, file = "./Data/Clean_FFPE_samplelist.csv", quote = F, row.names = F)


##### Collect FFPE + normal data from catalog ##### 

# Get germline and tumour QC data and make sure germline is not contaminated, add germline contamination data to FFPE_list, rewrite

# Function from DataFromBerthaCat01.R (fixed for new URL)
getBerthaQCmetrics <- function(SAMPLE_WELL_ID, sessionID, studyID="1000000036"){
  require(jsonlite)
  require(dplyr)
  # Create Catalog 1.0 search command line
  #command <- paste0('curl -X GET --header "Accept: application/json" "http://bio-prod-opencgainternal-tomcat-01.gel.zone:8080/opencga/webservices/rest/v1/files/search?sid=', sessionID, '&study=', studyID, '&name=', SAMPLE_WELL_ID, '.bam&skipCount=false&limit=1"')
  # FIXED (Sep 26 2017)
  command <- paste0('curl -X GET --header "Accept: application/json" "https://opencgainternal.gel.zone/opencga/webservices/rest/v1/files/search?sid=', sessionID, '&study=', studyID, '&name=', SAMPLE_WELL_ID, '.bam&skipCount=false&lazy=true"')
  #curl -X GET --header "Accept: application/json" "https://opencgainternal.gel.zone/opencga/webservices/rest/v1/files/search?sid=dILUsdxVTzLd0mOB3EWS&study=1000000036&name=LP2000907-DNA_A01.bam&skipCount=false&lazy=true"
  
  # Read in the json file containing the QC metrics for the specific sample
  bam_json <- fromJSON(system(command, intern = T), flatten = T)
  # Check if sample data is there (bam file registered in the catalog or not?)
  if (bam_json$response$numResults == 0) {
    result <- data.frame(
      WELL_ID = SAMPLE_WELL_ID,
      SAMPLE_TYPE = NA,
      TUMOUR_PURITY = NA,
      GC_DROP = NA,
      AT_DROP = NA,
      COVERAGE_HOMOGENEITY = NA,
      CHIMERIC_PER = NA, # corrected
      SOFT_CLIPPED_BASES_PER = NA, # added
      DEAMINATION_MISMATCHES_PER = NA, # missing
      AV_FRAGMENT_SIZE_BP = NA, 
      MED_FRAGMENT_LENGTH_ILLUMINA = NA,
      MAPPING_RATE_PER =  NA, 
      GbQ30NoDupsNoClip = NA, 
      perc_bases_ge_15x_mapQ_ge11 = NA, 
      SNVs =  NA, 
      INDELs =  NA, 
      SVs =  NA, 
      CNVs = NA, 
      SNV_LOG = NA, 
      INDEL_LOG = NA, 
      DUPL_RATE_ILLUMINA = NA,
      NON_SPATIAL_DUPL_RATE_ILLUMINA = NA,
      COSMIC_COV_LT30X = NA, # calculated with different method in Bertha
      MEDIAN_COV = NA,
      DIVERSITY_ILLUMINA = NA, 
      CONTAMINATION_PCT_VERIFYBAMID = NA,
      CONTAMINATION_PCT_ILLUMINA = NA, 
      DISCORDANCE_PCT = NA,
      BERTHA_NUM_TIMES_PROCESSED = NA,
      BERTHA_VERSION = NA
    )
    return(result)  
  }
  
  # Figure out if the sample is tumor or germline
  tumor <- sum(grepl("CANCER", names(bam_json$response$result[[1]]))) != 0
  
  # Number of times sample has been processed by Bertha
  n <- dim(bam_json$response$result[[1]]$attributes.processes[[1]])[1]
  
  # Extract QC metrics into a table (based on sample type)
  if (tumor) {
    result <- data.frame(
      WELL_ID = SAMPLE_WELL_ID,
      SAMPLE_TYPE = "TUMOR",
      #
      if (is.null(bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.ESTIMATED_PURITY)){
        TUMOUR_PURITY = NA
      }
      else {
        TUMOUR_PURITY = bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.ESTIMATED_PURITY
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.AT_GC_DROP.gc_drop)){
        GC_DROP = NA
      }
      else {
        GC_DROP = bam_json$response$result[[1]]$stats.AT_GC_DROP.gc_drop
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.AT_GC_DROP.at_drop)){
        AT_DROP = NA
      }
      else {
        AT_DROP = bam_json$response$result[[1]]$stats.AT_GC_DROP.at_drop
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.WHOLE_GENOME_COVERAGE.coverageSummary))
      {
        COVERAGE_HOMOGENEITY = NA
      }
      else {
        COVERAGE_HOMOGENEITY = round(bam_json$response$result[[1]]$stats.WHOLE_GENOME_COVERAGE.coverageSummary[[1]][bam_json$response$result[[1]]$stats.WHOLE_GENOME_COVERAGE.coverageSummary[[1]]$scope == "autosomes",]$localRMSD, 2)
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.SAMTOOLS_STATS_ALL.SAMTOOLS_READS_MAPPED_AND_PAIRED))
      {
        CHIMERIC_PER = NA
      }
      else {
        CHIMERIC_PER = round(bam_json$response$result[[1]]$stats.SAMTOOLS_STATS_ALL.SAMTOOLS_PAIRS_ON_DIFFERENT_CHROMOSOMES*200 / bam_json$response$result[[1]]$stats.SAMTOOLS_STATS_ALL.SAMTOOLS_READS_MAPPED_AND_PAIRED, 2)
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.TUMOR_PERCENT_SOFT_CLIPPED_BASES))
      {
        SOFT_CLIPPED_BASES_PER = NA
      }
      else {
        SOFT_CLIPPED_BASES_PER = bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.TUMOR_PERCENT_SOFT_CLIPPED_BASES
      },  
      #
      DEAMINATION_MISMATCHES_PER = NA, # missing
      #
      if (is.null(bam_json$response$result[[1]]$stats.SAMTOOLS_STATS_ALL.SAMTOOLS_INSERT_SIZE_AVERAGE))
      {
        AV_FRAGMENT_SIZE_BP = NA
      }
      else {
        AV_FRAGMENT_SIZE_BP = bam_json$response$result[[1]]$stats.SAMTOOLS_STATS_ALL.SAMTOOLS_INSERT_SIZE_AVERAGE
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.tumour_stats.FRAGMENT_LENGTH_MEDIAN)){
        MED_FRAGMENT_LENGTH_ILLUMINA = NA
      }
      else {
        MED_FRAGMENT_LENGTH_ILLUMINA = bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.tumour_stats.FRAGMENT_LENGTH_MEDIAN 
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.SAMTOOLS_STATS_ALL.SAMTOOLS_RAW_TOTAL_SEQUENCES))
      {
        MAPPING_RATE_PER = NA
      }
      else {
        MAPPING_RATE_PER =  round(bam_json$response$result[[1]]$stats.SAMTOOLS_STATS_ALL.SAMTOOLS_READS_MAPPED*100 / bam_json$response$result[[1]]$stats.SAMTOOLS_STATS_ALL.SAMTOOLS_RAW_TOTAL_SEQUENCES, 2)
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.GEL_METRICS.GbQ30NoDupsNoClip)){
        GbQ30NoDupsNoClip = NA
      }
      else {
        GbQ30NoDupsNoClip = round(bam_json$response$result[[1]]$stats.GEL_METRICS.GbQ30NoDupsNoClip, 2)
      }, 
      #
      if (is.null(bam_json$response$result[[1]]$stats.GEL_METRICS.perc_bases_ge_15x_mapQ_ge11)){
        perc_bases_ge_15x_mapQ_ge11 = NA
      }
      else {
        perc_bases_ge_15x_mapQ_ge11 = bam_json$response$result[[1]]$stats.GEL_METRICS.perc_bases_ge_15x_mapQ_ge11
      }, 
      #
      if (is.null(bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_SNVS)){
        SNVs = NA
      }
      else {
        SNVs =   bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_SNVS
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_INDELS)){
        INDELs = NA
      }
      else {
        INDELs =  bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_INDELS
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_SV_BREAKENDS)){
        SVs = NA
      }
      else {
        SVs = (bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_SV_BREAKENDS / 2) +
          bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_SV_DELETIONS +
          bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_SV_INSERTIONS +
          bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_SV_INVERSIONS +
          bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_SV_TANDEM_DUPLICATIONS
      }, 
      #
      if (is.null(bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_CNVS)){
        CNVs = NA
      }
      else {
        CNVs = bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_CNVS
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_SNVS)){
        SNV_LOG = NA
      }
      else {
        SNV_LOG = log(bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_SNVS, 10)
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_INDELS)){
        INDEL_LOG = NA
      }
      else {
        INDEL_LOG = log(bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.pair_stats.SOMATIC_INDELS, 10)
      },
      # 
      if (is.null(bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.tumour_stats.PERCENT_DUPLICATE_ALIGNED_READS)){
        DUPL_RATE_ILLUMINA = NA
      }
      else {
        DUPL_RATE_ILLUMINA = bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.tumour_stats.PERCENT_DUPLICATE_ALIGNED_READS
      },
      
      if (is.null(bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.tumour_stats.PERCENT_NON_SPATIAL_DUPLICATE_READ_PAIRS)){
        NON_SPATIAL_DUPL_RATE_ILLUMINA = NA
      }
      else {
        NON_SPATIAL_DUPL_RATE_ILLUMINA = bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.tumour_stats.PERCENT_NON_SPATIAL_DUPLICATE_READ_PAIRS
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.VARIANTS_COVERAGE.coverageSummary[[1]])){
        COSMIC_COV_LT30X = NA
      }
      else {
        #COSMIC_COV_LT30X = round((1 - bam_json$response$result[[1]]$stats.VARIANTS_COVERAGE.coverageSummary[[1]][bam_json$response$result[[1]]$stats.VARIANTS_COVERAGE.coverageSummary[[1]]$scope == "allchrs",]$gte30x)*100, 2)
        COSMIC_COV_LT30X = round((1 - bam_json$response$result[[1]]$stats.VARIANTS_COVERAGE.coverageSummary[[1]][bam_json$response$result[[1]]$stats.VARIANTS_COVERAGE.coverageSummary[[1]]$scope == "autosomes",]$gte30x)*100, 2)
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.WHOLE_GENOME_COVERAGE.coverageSummary))
      {
        MEDIAN_COV = NA
      }
      else {
        MEDIAN_COV = round(bam_json$response$result[[1]]$stats.WHOLE_GENOME_COVERAGE.coverageSummary[[1]][bam_json$response$result[[1]]$stats.WHOLE_GENOME_COVERAGE.coverageSummary[[1]]$scope == "autosomes",]$med, 2)
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.tumour_stats.DIVERSITY)){
        DIVERSITY_ILLUMINA = NA
      }
      else {
        DIVERSITY_ILLUMINA = bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT_CANCER.tumour_stats.DIVERSITY
      }, 
      #
      CONTAMINATION_PCT_VERIFYBAMID = NA,
      CONTAMINATION_PCT_ILLUMINA = NA, 
      DISCORDANCE_PCT = NA,
      #
      if (is.null(n)){
        BERTHA_NUM_TIMES_PROCESSED = NA
      }
      else {
        BERTHA_NUM_TIMES_PROCESSED = n
      },
      #
      if (is.null(bam_json$response$result[[1]]$attributes.processes[[1]]$version[n])){
        BERTHA_VERSION = NA
      }
      else {
        BERTHA_VERSION = bam_json$response$result[[1]]$attributes.processes[[1]]$version[n]
      }
    )
  }
  
  else { # germline
    result <- data.frame(
      WELL_ID = SAMPLE_WELL_ID,
      SAMPLE_TYPE = "GL",
      TUMOUR_PURITY = NA,
      #
      if (is.null(bam_json$response$result[[1]]$stats.AT_GC_DROP.gc_drop)){
        GC_DROP = NA
      }
      else {
        GC_DROP = bam_json$response$result[[1]]$stats.AT_GC_DROP.gc_drop
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.AT_GC_DROP.at_drop)){
        AT_DROP = NA
      }
      else {
        AT_DROP = bam_json$response$result[[1]]$stats.AT_GC_DROP.at_drop
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.WHOLE_GENOME_COVERAGE.coverageSummary[[1]]))  # changing this to test the 1st element of the list
      {
        COVERAGE_HOMOGENEITY = NA
      }
      else { 
        COVERAGE_HOMOGENEITY = round(bam_json$response$result[[1]]$stats.WHOLE_GENOME_COVERAGE.coverageSummary[[1]][bam_json$response$result[[1]]$stats.WHOLE_GENOME_COVERAGE.coverageSummary[[1]]$scope == "autosomes",]$localRMSD, 2)
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.SAMTOOLS_STATS_ALL.SAMTOOLS_READS_MAPPED_AND_PAIRED))
      {
        CHIMERIC_PER = NA
      }
      else {
        CHIMERIC_PER = round(bam_json$response$result[[1]]$stats.SAMTOOLS_STATS_ALL.SAMTOOLS_PAIRS_ON_DIFFERENT_CHROMOSOMES*200 / bam_json$response$result[[1]]$stats.SAMTOOLS_STATS_ALL.SAMTOOLS_READS_MAPPED_AND_PAIRED, 2)
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT.PERCENT_SOFT_CLIPPED_BASES))
      {
        SOFT_CLIPPED_BASES_PER = NA
      }
      else {
        SOFT_CLIPPED_BASES_PER = bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT.PERCENT_SOFT_CLIPPED_BASES
      },  
      #
      DEAMINATION_MISMATCHES_PER = NA, # missing
      #
      if (is.null(bam_json$response$result[[1]]$stats.SAMTOOLS_STATS_ALL.SAMTOOLS_INSERT_SIZE_AVERAGE))
      {
        AV_FRAGMENT_SIZE_BP = NA
      }
      else {
        AV_FRAGMENT_SIZE_BP = bam_json$response$result[[1]]$stats.SAMTOOLS_STATS_ALL.SAMTOOLS_INSERT_SIZE_AVERAGE
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT.FRAGMENT_LENGTH_MEDIAN)){
        MED_FRAGMENT_LENGTH_ILLUMINA = NA
      }
      else {
        MED_FRAGMENT_LENGTH_ILLUMINA = bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT.FRAGMENT_LENGTH_MEDIAN
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.SAMTOOLS_STATS_ALL.SAMTOOLS_RAW_TOTAL_SEQUENCES))
      {
        MAPPING_RATE_PER = NA
      }
      else {
        MAPPING_RATE_PER =  round(bam_json$response$result[[1]]$stats.SAMTOOLS_STATS_ALL.SAMTOOLS_READS_MAPPED*100 / bam_json$response$result[[1]]$stats.SAMTOOLS_STATS_ALL.SAMTOOLS_RAW_TOTAL_SEQUENCES, 2)
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.GEL_METRICS.GbQ30NoDupsNoClip)){
        GbQ30NoDupsNoClip = NA
      }
      else {
        GbQ30NoDupsNoClip = round(bam_json$response$result[[1]]$stats.GEL_METRICS.GbQ30NoDupsNoClip, 2)
      }, 
      #
      if (is.null(bam_json$response$result[[1]]$stats.GEL_METRICS.perc_bases_ge_15x_mapQ_ge11)){
        perc_bases_ge_15x_mapQ_ge11 = NA
      }
      else {
        perc_bases_ge_15x_mapQ_ge11 = bam_json$response$result[[1]]$stats.GEL_METRICS.perc_bases_ge_15x_mapQ_ge11
      }, 
      
      SNVs = NA, 
      INDELs =  NA, 
      SVs = NA, 
      CNVs = NA, 
      SNV_LOG = NA, 
      INDEL_LOG = NA,
      #
      if (is.null(bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT.PERCENT_DUPLICATE_ALIGNED_READS)){
        DUPL_RATE_ILLUMINA = NA
      }
      else {
        DUPL_RATE_ILLUMINA = bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT.PERCENT_DUPLICATE_ALIGNED_READS
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT.PERCENT_NON_SPATIAL_DUPLICATE_READ_PAIRS)){
        NON_SPATIAL_DUPL_RATE_ILLUMINA = NA
      }
      else {
        NON_SPATIAL_DUPL_RATE_ILLUMINA = bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT.PERCENT_NON_SPATIAL_DUPLICATE_READ_PAIRS
      },
      #
      COSMIC_COV_LT30X = NA,
      #
      if (is.null(bam_json$response$result[[1]]$stats.WHOLE_GENOME_COVERAGE.coverageSummary[[1]])) # changing this to test the 1st element of the list
      {
        MEDIAN_COV = NA
      }
      else { 
        MEDIAN_COV = round(bam_json$response$result[[1]]$stats.WHOLE_GENOME_COVERAGE.coverageSummary[[1]][bam_json$response$result[[1]]$stats.WHOLE_GENOME_COVERAGE.coverageSummary[[1]]$scope == "autosomes",]$med, 2)
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT.DIVERSITY)){
        DIVERSITY_ILLUMINA = NA
      }
      else {
        DIVERSITY_ILLUMINA = bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT.DIVERSITY
      }, 
      #
      if (is.null(bam_json$response$result[[1]]$stats.VERIFY_BAM_ID.FREEMIX)){
        CONTAMINATION_PCT_VERIFYBAMID = NA
      }
      else {
        CONTAMINATION_PCT_VERIFYBAMID = round( (bam_json$response$result[[1]]$stats.VERIFY_BAM_ID.FREEMIX * 100), 2)
      },
      #
      if (is.null(bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT.CONTAMINATION)){
        CONTAMINATION_PCT_ILLUMINA = NA
      }
      else {
        CONTAMINATION_PCT_ILLUMINA = bam_json$response$result[[1]]$stats.ILLUMINA_SUMMARY_REPORT.CONTAMINATION
      }, 
      # Added discordance percentage (germline only)
      if (is.null(bam_json$response$result[[1]]$stats.ARRAY_CONCORDANCE.numberOfSites)){
        DISCORDANCE_PCT = NA
      }
      else {
        DISCORDANCE_PCT = round( (bam_json$response$result[[1]]$stats.ARRAY_CONCORDANCE.numberOfDiscordantSites * 100 / bam_json$response$result[[1]]$stats.ARRAY_CONCORDANCE.numberOfSites), 2)
      },
      #
      if (is.null(n)){
        BERTHA_NUM_TIMES_PROCESSED = NA
      }
      else {
        BERTHA_NUM_TIMES_PROCESSED = n
      },
      #
      if (is.null(bam_json$response$result[[1]]$attributes.processes[[1]]$version[n])){
        BERTHA_VERSION = NA
      }
      else {
        BERTHA_VERSION = bam_json$response$result[[1]]$attributes.processes[[1]]$version[n]
      }
    ) 
  }
  names(result) <- c("WELL_ID", "SAMPLE_TYPE", "TUMOUR_PURITY", "GC_DROP", "AT_DROP", "COVERAGE_HOMOGENEITY", "CHIMERIC_PER", "SOFT_CLIPPED_BASES_PER", "DEAMINATION_MISMATCHES_PER", "AV_FRAGMENT_SIZE_BP",  "MED_FRAGMENT_LENGTH_ILLUMINA", "MAPPING_RATE_PER", "GbQ30NoDupsNoClip", "perc_bases_ge_15x_mapQ_ge11", "SNVs", "INDELs", "SVs", "CNVs ", "SNV_LOG", "INDEL_LOG", "DUPL_RATE_ILLUMINA", "NON_SPATIAL_DUPL_RATE_ILLUMINA", "COSMIC_COV_LT30X", "MEDIAN_COV", "DIVERSITY_ILLUMINA", "CONTAMINATION_PCT_VERIFYBAMID", "CONTAMINATION_PCT_ILLUMINA", "DISCORDANCE_PCT", "BERTHA_NUM_TIMES_PROCESSED", "BERTHA_VERSION")
  return(result)
}

### Collect all data from catalog
all_ids <- c(FFPE_list$Platekey, FFPE_list$Platekey_normal)
QC_catalog <- bind_rows(lapply(all_ids, getBerthaQCmetrics, sessionID = "dILUsdxVTzLd0mOB3EWS"))  # Re-running Sept 27 2017, Oct 5 2017

# Check that all tumour and germline data is there
dim(QC_catalog)  # 469  (one extra entry, possibly known sample LP3000115-DNA_G11 that was used for testing and has 2 catalog entries)
table(QC_catalog$SAMPLE_TYPE, exclude = NULL)  # 1 NA
QC_catalog[duplicated(QC_catalog$WELL_ID),]$WELL_ID  # LP3000115-DNA_G11 duplicated
QC_catalog %>% filter(is.na(SAMPLE_TYPE))  # LP3000067-DNA_H01 
# Identifying sample LP3000067-DNA_H01
FFPE_list %>% filter(Platekey == "LP3000067-DNA_H01" | Platekey_normal == "LP3000067-DNA_H01")  # germline of LP3000074-DNA_F05
upload_full %>% filter(Platekey == "LP3000067-DNA_H01")  # Processed with Bertha 1.1.8 but all fields missing from catalog
# Removing extra (incomplete) entry of LP3000115-DNA_G11
QC_catalog <- QC_catalog %>% filter(!(WELL_ID == "LP3000115-DNA_G11" & is.na(GC_DROP)))

# Sanity check of fields & missing data assessment (specifically, checking if tumours have GbQ30NoDupsNoClip and germlines GbQ30NoDupsNoClip, perc_bases_ge_15x_mapQ_ge11 and CONTAMINATION_PCT_VERIFYBAMID)
summary(QC_catalog)

# Make a list of germline samples with missing VerifyBamID data, check if it's already reported in JIRA (URGENTLY needed)
gl_missing_verifybamid <- as.character(QC_catalog %>% filter(SAMPLE_TYPE != "TUMOR", is.na(CONTAMINATION_PCT_VERIFYBAMID)) %>% pull(WELL_ID)) 
length(gl_missing_verifybamid)  # 154
table(QC_catalog[QC_catalog$WELL_ID %in% gl_missing_verifybamid,]$BERTHA_VERSION, exclude = NULL)  # all Bertha v 1.0
# Check whether these were reported already as mising in JIRA ticket BERTHA-356
reported_missing <- read.table("./Data/missing.contamination.V4.2017-08-10.txt", sep = "\t", header = T)
table(reported_missing$Type, exclude = NULL)  # only 9 cancer germlines
sum(gl_missing_verifybamid %in% reported_missing$Platekey)  # none reported here

# Write a list with missing GL verifybamid with delivery IDs
missing_verifybamid_list <- upload_full %>% filter(Platekey %in% gl_missing_verifybamid, `Delivery Version` == "V4") %>% select(Platekey, DeliveryID, `Delivery Version`)
sum(duplicated(missing_verifybamid_list$Platekey))  # 0
sum(duplicated(missing_verifybamid_list$DeliveryID))  # 0
write.csv(missing_verifybamid_list, file = "./Data/missing_verifybamid_list.csv", quote = F, row.names = F)

# # Collect FFPE tumour data separately from catalog
# QC_catalog_ffpe <- bind_rows(lapply(FFPE_list$Platekey, getBerthaQCmetrics, sessionID = "dILUsdxVTzLd0mOB3EWS"))
# 
# # Sanity check of FFPE samples
# summary(QC_catalog_ffpe)
# table(QC_catalog_ffpe$BERTHA_VERSION, exclude = NULL) #  6 different versions: 1.0 ,1.4.2, 1.5.0, 1.7.0, 1.7.1,  <NA>
# # Check if the 1.7+ versions have complete data
# summary(QC_catalog_ffpe[(QC_catalog_ffpe$BERTHA_VERSION %in% c("1.7.0", "1.7.1")),])  # 51 samples, CONTAMINATION_PCT_ILLUMINA and DISCORDANCE_PCT missing (ok)
# summary(QC_catalog_ffpe[(QC_catalog_ffpe$BERTHA_VERSION == "1.4.2"),]) # 16 samples, lots of missing metrics
# summary(QC_catalog_ffpe[(QC_catalog_ffpe$BERTHA_VERSION == "1.5.0"),])  # 5 samples, lots of missing metrics
# # Look up samples with missing Bertha version
# QC_catalog_ffpe %>% filter(is.na(BERTHA_VERSION))  # 4 missing; Bertha 1.3.2, 1.3.1 
# # Bottom line: intake QC should be re-run for all FFPE version 1.5.0 and under due to missing stats and metrics (not urgent)
# 
# 
# # Collect matching NORMAL data separately from catalog
# QC_catalog_gl <- bind_rows(lapply(FFPE_list$Platekey_normal, getBerthaQCmetrics, sessionID = "dILUsdxVTzLd0mOB3EWS"))  # error in sample # 115
# # Debugging
# getBerthaQCmetrics(FFPE_list$Platekey_normal[115], sessionID = "dILUsdxVTzLd0mOB3EWS") 
# # Bertha v 1.1.8, error bc 1st element of WHOLE_GENOME_COVERAGE list is NULL
# # MOre issues with this sample - it NEVER completed intake QC successfully
# 
# # Checking germline data again after changing the function so WHOLE_GENOME_COVERAGE won't fail
# QC_catalog_gl <- bind_rows(lapply(FFPE_list$Platekey_normal, getBerthaQCmetrics, sessionID = "dILUsdxVTzLd0mOB3EWS"))  # OK, but extra sample
# QC_catalog_gl[QC_catalog_gl$WELL_ID == QC_catalog_gl[duplicated(QC_catalog_gl$WELL_ID),]$WELL_ID,]  # problematic sample LP3000115-DNA_G11 listed twice (two deliveries in V4)
# # Remove second copy of the problematic sample
# QC_catalog_gl <- QC_catalog_gl %>% filter(!(WELL_ID == "LP3000115-DNA_G11" & is.na(AV_FRAGMENT_SIZE_BP)))
# 
# # Sanity check
# summary(QC_catalog_gl) # CONTAMINATION_PCT_VERIFYBAMID missing for 156 samples
# table(QC_catalog_gl$BERTHA_VERSION, exclude = NULL) #  8 different versions: 1.0 ,1.4.0, 1.4.2, 1.5.0, 1.6.0, 1.7.0, 1.7.1,  <NA>
# # Samples that miss VerifyBamID
# summary(QC_catalog_gl %>% filter(is.na(CONTAMINATION_PCT_VERIFYBAMID)))  # These must re-run
# table(QC_catalog_gl[is.na(QC_catalog_gl$CONTAMINATION_PCT_VERIFYBAMID),]$BERTHA_VERSION, exclude = NULL) # 155 on Bertha 1.0, one NA 
# # Samples that have VerifyBamID
# summary(QC_catalog_gl %>% filter(!is.na(CONTAMINATION_PCT_VERIFYBAMID)))  # 10 missing coverage homog, median cov, 2 missing bertha version
# table((QC_catalog_gl %>% filter(!is.na(CONTAMINATION_PCT_VERIFYBAMID)) %>% pull(BERTHA_VERSION)), exclude = NULL)
# 
# # Check if samples with missing VerifyBamID are already reported (with correct DeliveryID) in JIRA
# missing_verifybamid <- QC_catalog_gl %>% filter(is.na(CONTAMINATION_PCT_VERIFYBAMID)) %>% pull(WELL_ID)
# missing_cont_info <- c("LP2000773-DNA_H03","LP2000773-DNA_H01","LP3000364-DNA_D01","LP3000364-DNA_H09","LP2000762-DNA_E05","LP2000762-DNA_D09","LP2000773-DNA_C04")
# sum(missing_verifybamid %in% missing_cont_info) # none here



### Collect all data from catalog after Bertha re-run (Oct 5 2017)
all_ids <- c(FFPE_list$Platekey, FFPE_list$Platekey_normal)
QC_catalog <- bind_rows(lapply(all_ids, getBerthaQCmetrics, sessionID = "eyJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJtbWlqdXNrb3ZpYyIsImF1ZCI6Ik9wZW5DR0EgdXNlcnMiLCJpYXQiOjE1MDcyODg2NjAsImV4cCI6MTUwNzI5MDQ2MH0.sMcfyhV5SeTbiB608bG-_ESHcI3I2nXa5PTev8FUdsg"))  # Oct 6 2017

dim(QC_catalog)  # 469

table(QC_catalog$SAMPLE_TYPE, exclude = NULL)  # 1 NA
QC_catalog[duplicated(QC_catalog$WELL_ID),]$WELL_ID  # LP3000115-DNA_G11 duplicated
QC_catalog %>% filter(is.na(SAMPLE_TYPE))  # LP3000067-DNA_H01 
# Identifying sample LP3000067-DNA_H01
FFPE_list %>% filter(Platekey == "LP3000067-DNA_H01" | Platekey_normal == "LP3000067-DNA_H01")  # germline of LP3000074-DNA_F05
upload_full %>% filter(Platekey == "LP3000067-DNA_H01")  # Processed with Bertha 1.1.8 but all fields missing from catalog
# Removing extra (incomplete) entry of LP3000115-DNA_G11
QC_catalog <- QC_catalog %>% filter(!(WELL_ID == "LP3000115-DNA_G11" & is.na(GC_DROP)))

summary(QC_catalog) 
QC_catalog[is.na(QC_catalog$GbQ30NoDupsNoClip),]


# Look for any still missing VerifyBamID contamination for GL
sum(is.na(QC_catalog[QC_catalog$SAMPLE_TYPE == "GL",]$CONTAMINATION_PCT_VERIFYBAMID))  # 1 missing
summary((QC_catalog %>% filter(SAMPLE_TYPE == "GL") %>% pull(CONTAMINATION_PCT_VERIFYBAMID)))
QC_catalog %>% filter(SAMPLE_TYPE == "GL", CONTAMINATION_PCT_VERIFYBAMID > 2) %>% select(WELL_ID, CONTAMINATION_PCT_VERIFYBAMID, CONTAMINATION_PCT_ILLUMINA) # 1 contaminated (LP3000279-DNA_B05)

# VerifyBamID from catalog to the FFPE list
FFPE_list$CONTAMINATION_PCT_VERIFYBAMID <- QC_catalog[match(FFPE_list$Platekey_normal, QC_catalog$WELL_ID),]$CONTAMINATION_PCT_VERIFYBAMID
FFPE_list[is.na(FFPE_list$CONTAMINATION_PCT_VERIFYBAMID),] # sample LP3000067-DNA_H01 (GL) still missing contamination info, fails integrity check in Bertha; excluding it
FFPE_list$Excluded <- 0
FFPE_list[FFPE_list$Platekey_normal %in% c("LP3000067-DNA_H01", "LP3000279-DNA_B05"),]$Excluded <- 1
FFPE_list$Exclusion_reason <- ""

# Check completeness of catalog data for tumour samples and QC metrics (after 7 with missing data were re-run)
summary(QC_catalog[QC_catalog$SAMPLE_TYPE == "TUMOR",])  
missing_ffpe <- read.table("./Data/tumour_ffpe_missingQC.txt")
QC_catalog %>% filter(WELL_ID %in% missing_ffpe$V1)  # data still missing, abandoning catalog reading as a source of info

# Add exclusion reasons to the final FFPE manifest
FFPE_list$Exclusion_reason <- ""

# Add missing tumour types and collection date to the FFPE list
missing_clinical <- read.csv("./Data/missing_tumour_types_completed.csv")
FFPE_list[FFPE_list$TUMOUR_TYPE == "N/A",]$TUMOUR_TYPE <- missing_clinical[match(FFPE_list[FFPE_list$TUMOUR_TYPE == "N/A",]$Platekey, missing_clinical$WELL_ID),]$TUMOUR_TYPE
FFPE_list$COLLECTION_DATE <- as.character(FFPE_list$COLLECTION_DATE)
FFPE_list[is.na(FFPE_list$COLLECTION_DATE ),]$COLLECTION_DATE <- ""
missing_clinical$COLLECTION_DATE <- as.Date(missing_clinical$COLLECTION_DATE, format = "%d/%m/%y")
#FFPE_list[FFPE_list$COLLECTION_DATE == "",]$COLLECTION_DATE <- missing_clinical[match(FFPE_list[FFPE_list$COLLECTION_DATE == "",]$Platekey, missing_clinical$WELL_ID),]$COLLECTION_DATE



############ FINAL FFPE cohort cleanup ############ 

# New QC list from Alona (Oct 5 2017) that includes GL samples
QC_all <- read.csv("./Data/ready_to_upload.09-24-17.csv")
dim(QC_all)  # 3505
summary(QC_all)
table(QC_all$SAMPLE_TYPE, QC_all$LIBRARY_TYPE, exclude = NULL)  # some samples have type "FFPE" but no library prep info (pilot FFPE?)
sum(duplicated(QC_all$SAMPLE_WELL_ID))  # 0

# Make new, clean FFPE list
FFPE_list <- upload %>% filter(Plate %in% FFPE_plates)
dim(FFPE_list)  # 246
sum(duplicated(FFPE_list$Platekey))  # 1 duplicated (LP3000140-DNA_A03)
FFPE_list %>% filter(Platekey == FFPE_list[duplicated(FFPE_list$Platekey),]$Platekey)  # One version has a delivery problem, removing (ID BE00000001)
FFPE_list <- FFPE_list %>% filter(DeliveryID != "BE00000001")
table(FFPE_list$Status, exclude = NULL)  # 5 samples are QC failed previously, one of them still fails (invalid BAM: LP3000074-DNA_D06)
FFPE_list %>% filter(Status == "qc_failed") %>% select(Platekey, DeliveryID, `Delivery Date`, `Delivery Version`, Status)


# Add matching germlines from the upload report, check versions and QC status
FFPE_list$Path <- as.character(FFPE_list$Path )
FFPE_list$Platekey_normal <- sapply(FFPE_list$Path, function(x){strsplit(strsplit(x, split = "/")[[1]][6], split = "Normal")[[1]][2]})
dim(upload %>% filter(Platekey %in% FFPE_list$Platekey_normal))  # 246
FFPE_list$Platekey_normal[!FFPE_list$Platekey_normal %in% upload$Platekey]  # 0 (no normals missing from upload)

# Add exclusion variable to the final FFPE list
FFPE_list$Excluded <- 0
FFPE_list$Exclusion_reason <- ""

# Exclude the sample with invalid tumour BAM
FFPE_list[FFPE_list$Platekey == "LP3000074-DNA_D06",]$Excluded <- 1
FFPE_list[FFPE_list$Platekey == "LP3000074-DNA_D06",]$Exclusion_reason <- "Invalid tumour BAM"

# Look at matching normals in the full upload report (2 normals are qc_failed)
table((upload_full %>% filter(Platekey %in% FFPE_list$Platekey_normal) %>% pull(`Delivery Version`)), (upload_full %>% filter(Platekey %in% FFPE_list$Platekey_normal) %>% pull(Status)))
upload_full %>% filter(Platekey %in% FFPE_list$Platekey_normal, `Delivery Version` == "V4", Status == "qc_failed")
normals_failed <- upload_full %>% filter(Platekey %in% FFPE_list$Platekey_normal, `Delivery Version` == "V4", Status == "qc_failed") %>% pull(Platekey)
FFPE_list[FFPE_list$Platekey_normal %in% normals_failed,]
upload_full %>% filter(Platekey %in% normals_failed) 
# LP3000069-DNA_E03 QC-passed when delivered in V2, but failed when it was lifted to V4 (LP3000069-DNA_E03); Bertha error: "Minimum covered positions requirement not met: Positions covered at >= 15x is 94.5%. Must be >= 95%."
# This is now processed with "nometrics_1-8-0" workflow so it doesn't calculate intake metrics (since the delivery was already accepted); keeping this sample
# LP3000115-DNA_G11 this sample was used for testing (CANCT40004 ID, small BAM) and should be ignored

# Add tumour QC data to the FFPE list from upload report
dim(FFPE_list) # 245  14
FFPE_list <- left_join(FFPE_list, QC_all, by = c("Platekey" = "SAMPLE_WELL_ID"))
dim(FFPE_list) # 245  43

# Add germline QC data to the FFPE list
QC_gl <- QC_all %>% filter(SAMPLE_WELL_ID %in% FFPE_list$Platekey_normal)
QC_gl <- QC_gl %>% select(-(SNV), -(INDEL), -(TUMOUR_PURITY), -(PASS_TO_ANALYSIS), -(TUMOUR_CONTAMINATION), -(DUPLICATE_CHECK), -(SEX_CHECK), -(CENTER_CODE), -(SAMPLE_TYPE), -(PATIENT_ID))
names(QC_gl) <- paste0(names(QC_gl), "_GL")
FFPE_list <- left_join(FFPE_list, QC_gl, by = c("Platekey_normal" = "SAMPLE_WELL_ID_GL"))

# Check if all germline passes QC + add exclusions
# perc_bases_ge_15x_mapQ_ge11 > 95%
# GbQ30NoDupsNoClip > 85x10^9 (bases with Q>=30)

summary(FFPE_list$perc_bases_ge_15x_mapQ_ge11_GL) # 1 NA
summary(FFPE_list$GbQ30NoDupsNoClip_GL)  # 1 NA
FFPE_list[is.na(FFPE_list$perc_bases_ge_15x_mapQ_ge11_GL),] # Sample with NA is LP3000279-DNA_B05 (contaminated)

# Add exclusion
FFPE_list[FFPE_list$Platekey_normal == "LP3000279-DNA_B05",]$Excluded <- 1
FFPE_list[FFPE_list$Platekey_normal == "LP3000279-DNA_B05",]$Exclusion_reason <- "GL contaminated"


# Check if all tumour passes QC + add exclusions
# 212.5x10^9 bases with Q>=30
summary(FFPE_list$GbQ30NoDupsNoClip)  # 11 NA
FFPE_list[is.na(FFPE_list$GbQ30NoDupsNoClip),] # 11 NAs fail tumour contamination

# Add exclusion
tumour_contamin <- FFPE_list %>% filter(TUMOUR_CONTAMINATION == "Fail") %>% pull(Platekey)
FFPE_list[FFPE_list$Platekey %in% tumour_contamin,]$Excluded <- 1
FFPE_list[FFPE_list$Platekey %in% tumour_contamin,]$Exclusion_reason <- "TUMOUR contaminated"

# Add tumour type and collection date (from original QC, adding missing data)
missing_clinical <- read.csv("./Data/missing_tumour_types_completed.csv")
missing_clinical$COLLECTION_DATE <- as.Date(missing_clinical$COLLECTION_DATE, format = "%d/%m/%y")
QC$COLLECTION_DATE <- as.Date(QC$COLLECTION_DATE)
QC[QC$TUMOUR_TYPE == "N/A",]$TUMOUR_TYPE <- missing_clinical[match(QC[QC$TUMOUR_TYPE == "N/A",]$WELL_ID, missing_clinical$WELL_ID),]$TUMOUR_TYPE
table(QC$TUMOUR_TYPE)
QC[QC$TUMOUR_TYPE == "",]$COLLECTION_DATE <- "1000-01-01"  # Make a placeholder for samples with tumour contaminated so I can merge the collection dates accurately
sum(is.na(QC$COLLECTION_DATE))  # 50
QC[is.na(QC$COLLECTION_DATE),]$COLLECTION_DATE <- missing_clinical[match(QC[is.na(QC$COLLECTION_DATE),]$WELL_ID, missing_clinical$WELL_ID),]$COLLECTION_DATE
QC[QC$TUMOUR_TYPE == "",]$COLLECTION_DATE <- NA
# Merge with FFPE list
FFPE_list$TUMOUR_TYPE <- QC[match(FFPE_list$Platekey, QC$WELL_ID),]$TUMOUR_TYPE
FFPE_list$COLLECTION_DATE <- QC[match(FFPE_list$Platekey, QC$WELL_ID),]$COLLECTION_DATE

# Sanity check
summary(FFPE_list[FFPE_list$Excluded == 0,])

# LP3000074-DNA_D06 (tumour sample that was remapped - BAM needs fixing mate info, but apparently ok)

### Add last sample to the cohort
FFPE_list[FFPE_list$Platekey == "LP3000074-DNA_D06",]$Exclusion_reason <- "Invalid tumour BAM fixed"
FFPE_list[FFPE_list$Platekey == "LP3000074-DNA_D06",]$Excluded <- 0

# Change vector type
FFPE_list$PATIENT_ID <- as.character(FFPE_list$PATIENT_ID)

##### Final FFPE cohort list 
write.csv(FFPE_list, file = "./Data/Clean_FFPE_samplelist.csv", quote = F, row.names = F)  # Note that "Status" field is not updated after cleanup (all pass QC)


##### Sample level QC metrics summary ##### 

# Distribution of tumour types
FFPE_list$TUMOUR_TYPE <- as.character(FFPE_list$TUMOUR_TYPE)
as.data.frame(table(FFPE_list[FFPE_list$Excluded == 0,]$TUMOUR_TYPE))

### Calculate mean +/- 2 SD for QC metrics

QC_summary <- FFPE_list %>% filter(Excluded == 0) %>% summarise(AT_drop_mean = mean(AT_DROP), 
                                                                AT_drop_sd = sd(AT_DROP), 
                                                                GC_drop_mean = mean(GC_DROP), 
                                                                GC_drop_sd = sd(GC_DROP),
                                                                COVERAGE_HOMOGENEITY_mean = mean(COVERAGE_HOMOGENEITY),
                                                                COVERAGE_HOMOGENEITY_sd = sd(COVERAGE_HOMOGENEITY),
                                                                MAPPING_RATE_PER_mean = mean(MAPPING_RATE_PER),
                                                                MAPPING_RATE_PER_sd = sd(MAPPING_RATE_PER),
                                                                AV_FRAGMENT_SIZE_BP_mean = mean(AV_FRAGMENT_SIZE_BP),
                                                                AV_FRAGMENT_SIZE_BP_sd = sd(AV_FRAGMENT_SIZE_BP)
                                                                )

### Plot metrics, color outliers (mean +/- 2 SD) in red

# Collect samples outside of mean +/- 2 SD range in any QC_summary metrics, color red in plots
QC_outliers <- unique(c((FFPE_list %>% filter(AT_DROP > (QC_summary$AT_drop_mean+2*QC_summary$AT_drop_sd)) %>% pull(Platekey)),
                 (FFPE_list %>% filter(GC_DROP < (QC_summary$GC_drop_mean-2*QC_summary$GC_drop_sd)) %>% pull(Platekey)),
                 (FFPE_list %>% filter(COVERAGE_HOMOGENEITY > (QC_summary$COVERAGE_HOMOGENEITY_mean+2*QC_summary$COVERAGE_HOMOGENEITY_sd)) %>% pull(Platekey)),
                 (FFPE_list %>% filter(MAPPING_RATE_PER < (QC_summary$MAPPING_RATE_PER_mean-2*QC_summary$MAPPING_RATE_PER_sd)) %>% pull(Platekey)),
                 (FFPE_list %>% filter(AV_FRAGMENT_SIZE_BP > (QC_summary$AV_FRAGMENT_SIZE_BP_mean+2*QC_summary$AV_FRAGMENT_SIZE_BP_sd)) %>% pull(Platekey))))

FFPE_list$QC_outlier_2sd <- 0
FFPE_list[FFPE_list$Platekey %in% QC_outliers,]$QC_outlier_2sd <- 1


# AT drop
pdf(file = "./Plots/cohort_QC/AT_DROP.pdf")
print(ggplot(FFPE_list[FFPE_list$Excluded == 0,], aes(x="", y=AT_DROP, color = (QC_outlier_2sd + 1))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter() +
  labs(x = "", y = "AT dropout")  + 
  bigger +
  tiltedX +
  blank +
  theme(legend.position="none")
  )
dev.off()

# GC drop
pdf(file = "./Plots/cohort_QC/GC_DROP.pdf")
print(ggplot(FFPE_list[FFPE_list$Excluded == 0,], aes(x="", y=GC_DROP, color = (QC_outlier_2sd+1))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter() +
  labs(x = "", y = "GC dropout")  + 
  bigger + 
  tiltedX +
  blank +
  theme(legend.position="none"))
dev.off()

# Coverage unevenness
pdf(file = "./Plots/cohort_QC/COVERAGE_HOMOGENEITY.pdf")
print(ggplot(FFPE_list[FFPE_list$Excluded == 0,], aes(x="", y=COVERAGE_HOMOGENEITY, color = (QC_outlier_2sd+1))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter() +
  labs(x = "", y = "Unevenness of coverage")  + 
  bigger + 
  tiltedX +
  blank +
  theme(legend.position="none"))
dev.off()


# Mapping rate
pdf(file = "./Plots/cohort_QC/MAPPING_RATE_PER.pdf")
print(ggplot(FFPE_list[FFPE_list$Excluded == 0,], aes(x="", y=MAPPING_RATE_PER, color = (QC_outlier_2sd+1))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter() +
  labs(x = "", y = "Mapping rate (%)")  + 
  bigger + 
  tiltedX +
  blank +
  theme(legend.position="none"))
dev.off()

# Fragment size
pdf(file = "./Plots/cohort_QC/AV_FRAGMENT_SIZE_BP.pdf")
print(ggplot(FFPE_list[FFPE_list$Excluded == 0,], aes(x="", y=AV_FRAGMENT_SIZE_BP, color = (QC_outlier_2sd+1))) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter() +
  labs(x = "", y = "Average fragment size")  + 
  bigger + 
  tiltedX +
  blank +
  theme(legend.position="none"))
dev.off()


# # Plot mapping rate vs fragment size
# ggplot(FFPE_list[FFPE_list$Excluded == 0,], aes(x=AV_FRAGMENT_SIZE_BP, y=MAPPING_RATE_PER, color =(QC_outlier_3sd+1))) + 
#   geom_point() +
#   regr_line +
#   blank




  