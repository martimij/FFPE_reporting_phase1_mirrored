#!/bin/bash
#BSUB -J mergedVCF_new_annot_test
#BSUB -n 1
#BSUB -q bio
#BSUB -o /home/mmijuskovic/bsub_outputs/output.%J
#BSUB -e /home/mmijuskovic/bsub_outputs/error.%J


# Test the CellBase fix to skip decomposing variants
# Annotate full merged VCF with itself to check completeness

/home/jlopez/cellbase-SNAPSHOT/bin/cellbase.sh variant-annotation --species hsapiens --assembly GRCh38 --local \
--skip-decompose \
--custom-file /home/mmijuskovic/small_variant_freq/allCohortsMergedVCF/cancer_mainProgram_2017.merged.AF.sorted.vcf.gz \
--custom-file-fields AF_FFPE,AF_FFpcrfree,AF_FFnano \
--custom-file-id somatic_agg_vcf \
--include cytoband \
-i /home/mmijuskovic/small_variant_freq/allCohortsMergedVCF/cancer_mainProgram_2017.merged.AF.sorted.vcf.gz \
-o /home/mmijuskovic/small_variant_freq/cellbase_decomp_fix_test/cancer_mainProgram_2017.merged.AF.sorted.vcf_newIndex_annotation.json