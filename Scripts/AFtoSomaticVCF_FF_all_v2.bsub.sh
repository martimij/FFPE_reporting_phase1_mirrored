#!/bin/bash
#BSUB -J addAF_merged_FF_VCF
#BSUB -n 1
#BSUB -q bio
#BSUB -o /home/mmijuskovic/bsub_outputs/output.%J
#BSUB -e /home/mmijuskovic/bsub_outputs/error.%J
 
# Set the environment
source /etc/profile.d/modules.sh
source /lsf/prod/conf/profile.lsf

module load bcftools 
module load vt/0.577

# Working directory where merged VCF will be written
workdir="/home/mmijuskovic/small_variant_freq/FF/mergedVCF/"

# List of VCFs to be merged (need to have GT, PL values)
vcf_list="/home/mmijuskovic/small_variant_freq/FF/mergedVCF/FF_all_mainProgram_2017.txt"

# Call merging script
/home/mmijuskovic/small_variant_freq/AFtoSomaticVCF_v2.sh $workdir $vcf_list
