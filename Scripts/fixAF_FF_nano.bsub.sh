#!/bin/bash
#BSUB -J fixAF_FF_nano_VCF
#BSUB -n 1
#BSUB -q bio
#BSUB -o /home/mmijuskovic/bsub_outputs/output.%J
#BSUB -e /home/mmijuskovic/bsub_outputs/error.%J
 
# Set the environment
source /etc/profile.d/modules.sh
source /lsf/prod/conf/profile.lsf

module load bcftools 
module load vt/0.577

# Working directory where fixed VCF will be written
workdir="/home/mmijuskovic/small_variant_freq/FF/mergedVCF/final_for_production/fixed_for_production/"

# Merged VCF to be fixed (contains genotypes)
vcf_input="/home/mmijuskovic/small_variant_freq/FF/mergedVCF/final_for_production/FF_nano_mainProgram_2017.PASSonly.merged.sorted.vcf.gz"

# List of samples to remove
exclude_samples="/home/mmijuskovic/small_variant_freq/FF/mergedVCF/final_for_production/fixed_for_production/FFnano_to_remove.txt"

# Call merging script
/home/mmijuskovic/small_variant_freq/fixAF_inMergedVCF.sh $workdir $vcf_input $exclude_samples
