# Martina Mijuskovic
# Download data from Bertha catalog 1.0 & compare QC metrics calculations with Alona's
# Apr 27 2017

library(dplyr)
library(jsonlite)
library(data.table)


############ Final function (corrected May 31 2017) #############

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
  
  # NEW: Get only the result with gelStatus == READY ----> UPDATE: not sure this is necessary
  
  
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
      if (is.null(bam_json$response$result[[1]]$stats.WHOLE_GENOME_COVERAGE.coverageSummary[[1]]))  # changed
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
      if (is.null(bam_json$response$result[[1]]$stats.WHOLE_GENOME_COVERAGE.coverageSummary[[1]]))  # changed
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
      if (is.null(bam_json$response$result[[1]]$stats.WHOLE_GENOME_COVERAGE.coverageSummary)[[1]])  # changed
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
      if (is.null(bam_json$response$result[[1]]$stats.WHOLE_GENOME_COVERAGE.coverageSummary[[1]]))  # changed
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


# Test
bam_json <- getBerthaQCmetrics("LP2000907-DNA_A01", sessionID = "dILUsdxVTzLd0mOB3EWS")

# Collect QC metrics from Catalog 1.0
#bertha <- bind_rows(lapply(completeIDs, getBerthaQCmetrics, sessionID = "TjE7iIBWF2ss3Dz2oZjC"))

