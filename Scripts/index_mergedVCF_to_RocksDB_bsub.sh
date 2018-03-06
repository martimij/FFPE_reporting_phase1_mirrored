#!/bin/bash
#BSUB -J mergedVCF_index
#BSUB -n 1
#BSUB -q bio
#BSUB -o /home/mmijuskovic/bsub_outputs/output.%J
#BSUB -e /home/mmijuskovic/bsub_outputs/error.%J

module load cellbase/v4.5.1 

# Annotate a small VCF in order to trigger indexing
cellbase.sh variant-annotation --species hsapiens --assembly GRCh38 --local \
--skip-normalize \
--custom-file /genomes/bertha-test/resources/bertha/data/GRCh38Decoy/annotations/index_somatic_variants/cancer_mainProgram_2017.merged.AF.sorted.vcf.gz \
--custom-file-fields AF_FFPE,AF_FFpcrfree,AF_FFnano \
--custom-file-id somatic_agg_vcf \
-i /home/mmijuskovic/small_variant_freq/FFPE/norm_VCFs/LP3000156-DNA_H03.somatic.duprem.left.split.vcf.gz \
-o /home/mmijuskovic/small_variant_freq/FFPE/test/test_VF_LP3000156-DNA_H03_annotation.json



