#!/bin/bash
#BSUB -J mergedVCF_test
#BSUB -n 1
#BSUB -q bio
#BSUB -o /home/mmijuskovic/bsub_outputs/output.%J
#BSUB -e /home/mmijuskovic/bsub_outputs/error.%J

module load bcftools 
module load vt/0.577


#### VCF aggregation test log

touch /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt


#### FFPE 

echo "********** FFPE before INFO modification **********" &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt
vt peek /home/mmijuskovic/small_variant_freq/FFPE/mergedVCF/final_for_production/FFPE_mainProgram_2017.PASSonly.merged.AF.info_only.vcf.gz &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt

echo "********** FFPE after INFO modification **********" &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt
vt peek /home/mmijuskovic/small_variant_freq/FFPE/mergedVCF/final_for_production/FFPE_mainProgram_2017.PASSonly.merged.AF.info_only.vcf.gz.renamed.vcf.gz &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt

#### FF nano

echo "********** FF nano before INFO modification **********" &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt
vt peek /home/mmijuskovic/small_variant_freq/FF/mergedVCF/final_for_production/fixed_for_production/FF_nano_mainProgram_2017.clean.PASSonly.merged.AF.info_only.vcf.gz &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt

echo "********** FF nano after INFO modification **********" &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt
vt peek /home/mmijuskovic/small_variant_freq/FF/mergedVCF/final_for_production/fixed_for_production/FF_nano_mainProgram_2017.clean.PASSonly.merged.AF.info_only.vcf.gz.renamed.vcf.gz &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt

#### FF PCR free

echo "********** FF PCRfree before INFO modification **********" &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt
vt peek /home/mmijuskovic/small_variant_freq/FF/mergedVCF/final_for_production/fixed_for_production/FF_PCRfree_mainProgram_2017.clean.PASSonly.merged.AF.info_only.vcf.gz &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt

echo "********** FF PCRfree after INFO modification **********" &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt
vt peek /home/mmijuskovic/small_variant_freq/FF/mergedVCF/final_for_production/fixed_for_production/FF_PCRfree_mainProgram_2017.clean.PASSonly.merged.AF.info_only.vcf.gz.renamed.vcf.gz &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt



#### Final merged VCF

echo "********** Merged VCF header modifications **********" &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt
bcftools view /home/mmijuskovic/small_variant_freq/allCohortsMergedVCF/cancer_mainProgram_2017.merged.AF.sorted.vcf.gz -h | grep cohort &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt

echo "********** AF annotations in merged VCF **********" &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt

echo "Number of FFPE annotations"  &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt
bcftools view /home/mmijuskovic/small_variant_freq/allCohortsMergedVCF/cancer_mainProgram_2017.merged.AF.sorted.vcf.gz -H | grep AF_FFPE | wc -l &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt

echo "Number of FF nano annotations"  &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt
bcftools view /home/mmijuskovic/small_variant_freq/allCohortsMergedVCF/cancer_mainProgram_2017.merged.AF.sorted.vcf.gz -H | grep AF_FFnano | wc -l &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt

echo "Number of FF PCR free annotations"  &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt
bcftools view /home/mmijuskovic/small_variant_freq/allCohortsMergedVCF/cancer_mainProgram_2017.merged.AF.sorted.vcf.gz -H | grep AF_FFpcrfree | wc -l &>> /home/mmijuskovic/small_variant_freq/all_cohorts_merge_test.txt







