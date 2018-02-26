#!/bin/bash
#BSUB -J merge3cohortVCFs
#BSUB -n 1
#BSUB -q bio
#BSUB -o /home/mmijuskovic/bsub_outputs/output.%J
#BSUB -e /home/mmijuskovic/bsub_outputs/error.%J
 
# Merges 3 cohort info_only VCFs
# Sorts merged VCF
# bgzips, tabixes modified VCF


# Set the environment
source /etc/profile.d/modules.sh
source /lsf/prod/conf/profile.lsf

module load bcftools 
module load vt/0.577

# Enter working directory
echo "Entering working directory /home/mmijuskovic/small_variant_freq/allCohortsMergedVCF"
cd /home/mmijuskovic/small_variant_freq/allCohortsMergedVCF

# Merge 3 cancer cohort VCFs
echo "Merging 3 cohort VCFs .... "
bcftools merge -m none \
/home/mmijuskovic/small_variant_freq/FF/mergedVCF/final_for_production/fixed_for_production/FF_PCRfree_mainProgram_2017.clean.PASSonly.merged.AF.info_only.vcf.gz.renamed.vcf.gz \
/home/mmijuskovic/small_variant_freq/FF/mergedVCF/final_for_production/fixed_for_production/FF_nano_mainProgram_2017.clean.PASSonly.merged.AF.info_only.vcf.gz.renamed.vcf.gz \
/home/mmijuskovic/small_variant_freq/FFPE/mergedVCF/final_for_production/FFPE_mainProgram_2017.PASSonly.merged.AF.info_only.vcf.gz.renamed.vcf.gz \
-o cancer_mainProgram_2017.merged.AF.vcf

# Sort VCF
echo "Sorting merged VCF .... "
vt sort cancer_mainProgram_2017.merged.AF.vcf -o cancer_mainProgram_2017.merged.AF.sorted.vcf

# bgzip and index
echo "Indexing .... "
bgzip cancer_mainProgram_2017.merged.AF.sorted.vcf
tabix cancer_mainProgram_2017.merged.AF.sorted.vcf.gz

