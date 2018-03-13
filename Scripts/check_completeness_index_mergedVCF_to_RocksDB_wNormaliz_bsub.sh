#!/bin/bash
#BSUB -J mergedVCF_index_test
#BSUB -n 1
#BSUB -q bio
#BSUB -o /home/mmijuskovic/bsub_outputs/output.%J
#BSUB -e /home/mmijuskovic/bsub_outputs/error.%J

module load cellbase/v4.5.1 

# Create log file
logfile=/home/mmijuskovic/small_variant_freq/allCohortsMergedVCF/test/mergedVCF_index_test_log.txt

# Annotate the original merged VCF with the same file to check that index was created correctly
# CHANGE (Mar 12 2018: remove "skip_normalize")
cellbase.sh variant-annotation --species hsapiens --assembly GRCh38 --local \
--custom-file /genomes/bertha-test/resources/bertha/data/GRCh38Decoy/annotations/index_somatic_variants/cancer_mainProgram_2017.merged.AF.sorted.vcf.gz \
--custom-file-fields AF_FFPE,AF_FFpcrfree,AF_FFnano \
--custom-file-id somatic_agg_vcf \
--include cytoband \
-i /genomes/bertha-test/resources/bertha/data/GRCh38Decoy/annotations/index_somatic_variants/cancer_mainProgram_2017.merged.AF.sorted.vcf.gz \
-o /home/mmijuskovic/small_variant_freq/allCohortsMergedVCF/test/cancer_mainProgram_2017.merged.AF.sorted.self_annotation.json


# Examine the output json to check if all variants are annotated
echo "Count of AF_FFPE" >> $logfile
grep AF_FFPE /home/mmijuskovic/small_variant_freq/allCohortsMergedVCF/test/cancer_mainProgram_2017.merged.AF.sorted.self_annotation.json | wc -l &>> $logfile

echo "Count of AF_FFpcrfree" >> $logfile
grep AF_FFpcrfree /home/mmijuskovic/small_variant_freq/allCohortsMergedVCF/test/cancer_mainProgram_2017.merged.AF.sorted.self_annotation.json | wc -l &>> $logfile

echo "Count of AF_FFnano" >> $logfile
grep AF_FFnano /home/mmijuskovic/small_variant_freq/allCohortsMergedVCF/test/cancer_mainProgram_2017.merged.AF.sorted.self_annotation.json | wc -l &>> $logfile

