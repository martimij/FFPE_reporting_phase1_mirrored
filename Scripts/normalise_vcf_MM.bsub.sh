#!/bin/bash
#BSUB -J norm_FFPE_VCFs[1-8]
#BSUB -n 1
#BSUB -q bio
#BSUB -o /home/mmijuskovic/bsub_outputs/output.%I.%J
#BSUB -e /home/mmijuskovic/bsub_outputs/error.%I.%J

# Set the environment
source /etc/profile.d/modules.sh
source /lsf/prod/conf/profile.lsf

# MM CHANGE: Instead of loading the test environment, load python module that should
# contain newest vcf_dedupper tool
#source /genomes/software/src/test-venvs/cancer-test/bin/activate  

module load python/2.7.12
module load bcftools/1.3
module load htslib/1.3
module load vt/ee9a751

# Working directory where normalized VCFs will be written
workdir="/home/mmijuskovic/small_variant_freq/FFPE/norm_VCFs/"

# Read the CSV input file in the format "sample_name,VCF_path" (29 samples/job)
IFS=","
while read f1 f2
do
        # Read the input line
        echo "Sample name is    : $f1"
        echo "VCF path is		: $f2"
        
        # Call normalization script
        /home/mmijuskovic/small_variant_freq/normalise_vcf_MM.sh $workdir $f1 $f2 
        
done < /home/mmijuskovic/small_variant_freq/FFPE/norm_VCFs_input/samples_for_recurrent_calc_normalization_input_${LSB_JOBINDEX}.csv


