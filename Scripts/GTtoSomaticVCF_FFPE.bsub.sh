#!/bin/bash
#BSUB -J addGT_FFPE_VCFs[1-8]
#BSUB -n 1
#BSUB -q bio
#BSUB -o /home/mmijuskovic/bsub_outputs/output.%I.%J
#BSUB -e /home/mmijuskovic/bsub_outputs/error.%I.%J
 
# Set the environment
source /etc/profile.d/modules.sh
source /lsf/prod/conf/profile.lsf

module load bcftools 
module load vt/0.577

# Working directory where modified VCFs will be written
workdir="/home/mmijuskovic/small_variant_freq/FFPE/GT_VCFs/"

# Read the CSV input file in the format "sample_name,VCF_path" (29 samples/job)
IFS=","
while read f1 f2
do
        # Read the input line
        echo "Sample name is    : $f1"
        echo "VCF path is		: $f2"
        
        # Call script for adding GT values
        /home/mmijuskovic/small_variant_freq/GTtoSomaticVCF.sh $workdir $f2 
        
done < /home/mmijuskovic/small_variant_freq/FFPE/addGT_input/samples_for_recurrent_calc_addGT_input_${LSB_JOBINDEX}.csv


