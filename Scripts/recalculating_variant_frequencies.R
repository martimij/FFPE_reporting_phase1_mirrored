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

today <- Sys.Date()


########## List of samples to add ########## 

### FFPE

# FFPE_list is the original full list (232 non-exluded)
# FFPE_list_new contains the updated list

# FF (nano and PCR-free)



