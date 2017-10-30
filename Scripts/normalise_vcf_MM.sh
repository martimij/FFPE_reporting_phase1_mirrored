#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Rename arguments (workdir is where the output goes; sample_name is LP sample ID; 
# input is VCF path
workdir=$1
sample_name=$2
input=$3

# Set reference
reference="/genomes/bertha-test/resources/bertha/data/GRCh38Decoy/reference/GRCh38Decoy_no_alt_no_ambiguous.fa"

echo "$(date)"
echo "Sample ${sample_name}: Normalising variants"


cd $workdir

echo "Filtering...."

bcftools view -e 'ALT="."' -Oz -o ${sample_name}.somatic.filter.vcf.gz ${input}
bcftools index ${sample_name}.somatic.filter.vcf.gz

# split records with multiple alternate alleles into multiple biallelic records
# an e.g. 1/2 genotype will be split to 1/. and ./1
# You must specify the -s ("smart") option in order for INFO and FORMAT fields of type A and R to be retained and decomposed appropriately.

echo "Decomposing...."

vt decompose -s ${sample_name}.somatic.filter.vcf.gz -o ${sample_name}.somatic.split.vcf.gz
bcftools index ${sample_name}.somatic.split.vcf.gz

echo "Normalizing...."
# left-align indels and trim redundant bases
vt normalize -w 10000 -r $reference ${sample_name}.somatic.split.vcf.gz -o ${sample_name}.somatic.left.split.vcf.gz
bcftools index ${sample_name}.somatic.left.split.vcf.gz

# decompose complex variants to allelic primitives
# this will only operate on variants where REF and ALT alleles have the same length unless the -a option is specified
# not specifying -a for now
#vt decompose_blocksub -p ${sample_name}.somatic.left.split.vcf.gz -o ${sample_name}.somatic.atomic.left.split.vcf.gz
#bcftools index ${sample_name}.somatic.atomic.left.split.vcf.gz

echo "Deduplicating...."
# remove duplicates
# this program requires all fields of records to be identical before they are considered duplicates
# cf vt uniq, which defines duplicates based only on CHROM, POS, REF, ALT
vcf_dedupper --input-vcf ${sample_name}.somatic.left.split.vcf.gz --sort --variant-caller strelka --selection-method af --equality-mode 1 --sample-name ${sample_name} | bgzip --stdout > ${sample_name}.somatic.duprem.left.split.vcf.gz
bcftools index ${sample_name}.somatic.duprem.left.split.vcf.gz
tabix -p vcf ${sample_name}.somatic.duprem.left.split.vcf.gz

echo "Cleanup...."
# This script extracts somatic calls only, adds genotype, replaces old sample name in vcf header, reindexes
# Old file is backed up as "somatic.duprem.left.split.vcf.orig.gz"
/genomes/software/apps/InterpretationWorkflowCancer/bash/cleanup_somaticVCF.sh -v ${sample_name}.somatic.duprem.left.split.vcf.gz -c ${sample_name}
