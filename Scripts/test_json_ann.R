# Read subset of the self-annotation JSON created by CellBase annotation
# 
# VCF used for self-annotation:
# /genomes/bertha-test/resources/bertha/data/GRCh38Decoy/annotations/index_somatic_variants/cancer_mainProgram_2017.merged.AF.sorted.vcf.gz
#
# command to create subset JSON:
# head -75168 /home/mmijuskovic/small_variant_freq/allCohortsMergedVCF/test/cancer_mainProgram_2017.merged.AF.sorted.self_annotation.json > \
# /home/mmijuskovic/small_variant_freq/allCohortsMergedVCF/test/cancer_mainProgram_2017.merged.AF.sorted.self_annotation.subset.json
#

library(dplyr)
library(jsonlite)

#j <- fromJSON("./Tests/cancer_mainProgram_2017.merged.AF.sorted.self_annotation.subset.json")  # doesn't work, since apparently this file is 1 JSON per line, rather than 1 JSON

js <- lapply(readLines("./Tests/cancer_mainProgram_2017.merged.AF.sorted.self_annotation.subset.json"), fromJSON, flatten = T)

# Annotation
js[[1]]$annotation$additionalAttributes$somatic_agg_vcf
names(js[[1]]$annotation$additionalAttributes$somatic_agg_vcf[[1]])

# Original INFO
as.data.frame(js[[1]]$studies$files)
names(js[[1]]$studies$files[[1]]) 

# Find lines where the INFO doesn't match the annotation (either one missing)

js_df <- data.frame(
  CHR = sapply(1:length(js), function(x){
  js[[x]]$chromosome
}),
POS = sapply(1:length(js), function(x){
  js[[x]]$start
}),
REF = sapply(1:length(js), function(x){
  js[[x]]$reference
}),
ALT = sapply(1:length(js), function(x){
  js[[x]]$alternate
}),
PCRfree_INFO = sapply(1:length(js), function(x){
  "attributes.AF_FFpcrfree" %in% names(js[[x]]$studies$files[[1]])
}),
nano_INFO = sapply(1:length(js), function(x){
  "attributes.AF_FFnano" %in% names(js[[x]]$studies$files[[1]])
}),  
FFPE_INFO = sapply(1:length(js), function(x){
  "attributes.AF_FFPE" %in% names(js[[x]]$studies$files[[1]])
}),
PCRfree_ann = sapply(1:length(js), function(x){
  "AF_FFpcrfree" %in% names(js[[x]]$annotation$additionalAttributes$somatic_agg_vcf[[1]])
}),
nano_ann = sapply(1:length(js), function(x){
  "AF_FFnano" %in% names(js[[x]]$annotation$additionalAttributes$somatic_agg_vcf[[1]])
}),
FFPE_ann = sapply(1:length(js), function(x){
  "AF_FFPE" %in% names(js[[x]]$annotation$additionalAttributes$somatic_agg_vcf[[1]])
})
)

# Find exceptions

js_df %>% filter((!(PCRfree_INFO + PCRfree_ann) %in% c(0,2)))
dim(js_df)  # 75168
dim(js_df %>% filter((!(PCRfree_INFO + PCRfree_ann) %in% c(0,2))))  # 10224

head(js_df %>% filter((!(nano_INFO + nano_ann) %in% c(0,2))))
dim(js_df %>% filter((!(nano_INFO + nano_ann) %in% c(0,2))))  # 3048
dim(js_df %>% filter(((nano_INFO + nano_ann) %in% c(0,2))))  # 72120

head(js_df %>% filter((!(FFPE_INFO + FFPE_ann) %in% c(0,2))))
dim(js_df %>% filter((!(FFPE_INFO + FFPE_ann) %in% c(0,2))))  # 5120
dim(js_df %>% filter(((FFPE_INFO + FFPE_ann) %in% c(0,2))))  # 70048


# Proper annotations
head(js_df %>% filter(((nano_INFO + nano_ann) %in% c(0,2))),20)  # looks like all SNVs are proper and indels not (does it have to do with "skip_normalize" option in CellBase? -YES)


#### Re-test after annotation WITHOUT "skip_normalize" option
# NOTE: "normalize" in this case meant just changing from VCF format for indels to format that doesn't include the base previous to indel

js_norm <- lapply(readLines("./Tests/cancer_mainProgram_2017.merged.AF.sorted.self_annotation_normSubset.json"), fromJSON, flatten = T)


# Find lines where the INFO doesn't match the annotation (either one missing)

js_norm_df <- data.frame(
  CHR = sapply(1:length(js_norm), function(x){
    js_norm[[x]]$chromosome
  }),
  POS = sapply(1:length(js_norm), function(x){
    js_norm[[x]]$start
  }),
  REF = sapply(1:length(js_norm), function(x){
    js_norm[[x]]$reference
  }),
  ALT = sapply(1:length(js_norm), function(x){
    js_norm[[x]]$alternate
  }),
  PCRfree_INFO = sapply(1:length(js_norm), function(x){
    "attributes.AF_FFpcrfree" %in% names(js_norm[[x]]$studies$files[[1]])
  }),
  nano_INFO = sapply(1:length(js_norm), function(x){
    "attributes.AF_FFnano" %in% names(js_norm[[x]]$studies$files[[1]])
  }),  
  FFPE_INFO = sapply(1:length(js_norm), function(x){
    "attributes.AF_FFPE" %in% names(js_norm[[x]]$studies$files[[1]])
  }),
  PCRfree_ann = sapply(1:length(js_norm), function(x){
    "AF_FFpcrfree" %in% names(js_norm[[x]]$annotation$additionalAttributes$somatic_agg_vcf[[1]])
  }),
  nano_ann = sapply(1:length(js_norm), function(x){
    "AF_FFnano" %in% names(js_norm[[x]]$annotation$additionalAttributes$somatic_agg_vcf[[1]])
  }),
  FFPE_ann = sapply(1:length(js_norm), function(x){
    "AF_FFPE" %in% names(js_norm[[x]]$annotation$additionalAttributes$somatic_agg_vcf[[1]])
  })
)

# Find exceptions
dim(js_norm_df %>% filter((!(nano_INFO + nano_ann) %in% c(0,2))))  # 20  (wrong)
dim(js_norm_df %>% filter(((nano_INFO + nano_ann) %in% c(0,2))))  # 75148 (right)
js_norm_df %>% filter((!(nano_INFO + nano_ann) %in% c(0,2)))


dim(js_norm_df %>% filter((!(FFPE_INFO + FFPE_ann) %in% c(0,2))))  # 14 (wrong)
dim(js_norm_df %>% filter(((FFPE_INFO + FFPE_ann) %in% c(0,2))))  # 75154 (right)
js_norm_df %>% filter((!(FFPE_INFO + FFPE_ann) %in% c(0,2)))

dim(js_norm_df %>% filter((!(PCRfree_INFO + PCRfree_ann) %in% c(0,2))))  # 21  (wrong)
dim(js_norm_df %>% filter(((PCRfree_INFO + PCRfree_ann) %in% c(0,2))))  # 75147 (right)
js_norm_df %>% filter((!(PCRfree_INFO + PCRfree_ann) %in% c(0,2)))

# All variants annotated wrong in this subset
js_norm_df %>% filter((!(nano_INFO + nano_ann) %in% c(0,2)) | (!(FFPE_INFO + FFPE_ann) %in% c(0,2)) | (!(PCRfree_INFO + PCRfree_ann) %in% c(0,2)))






