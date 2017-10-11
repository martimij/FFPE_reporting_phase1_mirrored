# Martina Mijuskovic
# FFPE project - small variant reporting
# Check cancer participant information in the catalog
# Sept 2017

library(dplyr)
library(ggplot2)
library(data.table)
library(jsonlite)

# Download all cancer participats currently in the catalog (Study: 1000000041, used for interpretation pipeline development by Antonio, everything there at the moment)

getCatalogParticipants <- function(PATIENT_ID, sessionID, studyID="1000000041"){
  require(jsonlite)
  require(dplyr)
  # Create Catalog 1.0 search command line
  #curl -X GET --header "Accept: application/json" --header "Authorization: Bearer " "https://opencgainternal.gel.zone/opencga/webservices/rest/v1/individuals/search?sid=eyJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJiZXJ0aGEiLCJhdWQiOiJPcGVuQ0dBIHVzZXJzIiwiaWF0IjoxNTA3NjQwNjU5LCJleHAiOjE1MDc2NDI0NTl9.kGXQJE-xRLzL4DfB1dwHwyiUF9gkfUJBmjtkm6SzTYw&study=1000000041&skipCount=false&limit=1"
  #curl -X GET --header "Accept: application/json" --header "Authorization: Bearer " "https://opencgainternal.gel.zone/opencga/webservices/rest/v1/individuals/search?sid=eyJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJiZXJ0aGEiLCJhdWQiOiJPcGVuQ0dBIHVzZXJzIiwiaWF0IjoxNTA3NjQwNjU5LCJleHAiOjE1MDc2NDI0NTl9.kGXQJE-xRLzL4DfB1dwHwyiUF9gkfUJBmjtkm6SzTYw&study=1000000041&name=217000038&skipCount=false&limit=1"

  command <- paste0('curl -X GET --header "Accept: application/json" --header "Authorization: Bearer " "https://opencgainternal.gel.zone/opencga/webservices/rest/v1/individuals/search?sid=', sessionID, '&study=', studyID, '&name=', PATIENT_ID, '&skipCount=false&limit=1"')
  # Read in the json file containing the QC metrics for the specific sample
  bam_json <- fromJSON(system(command, intern = T), flatten = T)
  
  if (bam_json$response$numResults == 0) {
    result <- data.frame(
      PATIENT_ID = PATIENT_ID,
      Platekey_normal = NA,
      Platekey = NA)
    return(result)
  }
  result <- data.frame(PATIENT_ID = bam_json$response$result[[1]]$name, 
                       Platekey_normal = bam_json$response$result[[1]]$annotationSets[[1]]$annotations[[2]]$value[1],
                       Platekey = bam_json$response$result[[1]]$annotationSets[[1]]$annotations[[2]]$value[2])  # Needs fixing in case FF and FFPE exist for this patient
  return(result)
  }

###  Check if all FFPE participants are registered 

### Collect all data from catalog ### 
participants_catalog <- bind_rows(lapply(FFPE_list$PATIENT_ID, getCatalogParticipants, sessionID = "eyJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJiZXJ0aGEiLCJhdWQiOiJPcGVuQ0dBIHVzZXJzIiwiaWF0IjoxNTA3NjQyODE0LCJleHAiOjE1MDc2NDQ2MTR9.ez-COsAYmyIvPtWTwr3HSzehD3K_NyuhwrRPn2GFBkU")) 

# Number of missing samples
sum(is.na(participants_catalog$Platekey))  # 50

# Check whether the missing samples are those that were missing tumour type info
missing_clinical$PATIENT_ID <- FFPE_list[match(missing_clinical$WELL_ID, FFPE_list$Platekey),]$PATIENT_ID
sum((participants_catalog %>% filter(is.na(participants_catalog$Platekey)) %>% pull(PATIENT_ID)) %in% missing_clinical$PATIENT_ID)  # 50

# Sanity check of catalog vs FFPE list info
FFPE_list %>% filter(PATIENT_ID %in% head(participants_catalog$PATIENT_ID)) %>% select(PATIENT_ID, Platekey_normal, Platekey)  # DOES NOT MATCH - make sure to read FFPE, not FF

### Look for missing data in LabKey ###

temp <- read.csv("./Data/all_FFPE_labkey.csv") # table created by filtering for "FFPE tumour"

sum(FFPE_list$Platekey %in% temp$PlateWell)  # only 164/245 samples there
sum(FFPE_list$PATIENT_ID %in% temp$Participant.ID)  # 170
sum(missing_clinical$WELL_ID %in% temp$PlateWell)  # missing 50 samples all there

temp2 <- read.csv("./")

### Create clinical JSONs for 50 missing samples from Catalog ### 

# Create "table_matched"
table_matched <- missing_clinical %>% select(PATIENT_ID, WELL_ID)
table_matched$Platekey_normal <- FFPE_list[match(table_matched$WELL_ID, FFPE_list$Platekey),]$Platekey_normal
table_matched <- table_matched %>% select(PATIENT_ID, Platekey_normal, WELL_ID)
names(table_matched) <- c("gelId","germlineSampleId","tumorSampleId")

# Create "table_participant"
# gelId,center,centerPatientId,labkeyParticipantId,sex,programmeConsent,primaryFindingConsent,secondaryFindingConsent,carrierStatusConsent,assignedICD10,sampleDiagnosis


# Create "table_samples"
# gelId,sampleId,labId,gelPhase,sampleType,sampleDiagnosis,tumorType,tumorSubType,preservationMethod,phase,method,tumorContent,grade,tnm_stage_version,tmn_stage_grouping,cellularity,sampleSource,clinic_sample_date_time




