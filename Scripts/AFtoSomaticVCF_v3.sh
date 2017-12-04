#!/bin/bash -e
set -e

### Pre-processing somatic Strelka VCF for allele frequency calculations ###
#
# Merges the VCF files from the input list (must contain GT and PL values)
# Sorts merged VCF
# Adds AC, AN, NS, AF to INFO field
# Removes samples
# bgzips, tabixes modified VCF
# Cleanup of intermediate files
#
# NEEDS: module load bcftools vt/0.577
#
# MMijuskovic, Nov 3 2017
# 
# CHANGES:
# Nov 30 2017 (v2): 
# 1) modified merging to take only PASS filter variants into account (if filter is not PASS
# in none of the samples, there is no output; if it's in at least one of the samples, the 
# non-PASS genotypes will be set as missing) 
# 2) creating a subset VCF that contains coding regions only (bed file with regions hard coded)
#
# Dec 4 2017 (v3):
# 1) not removing original merged sorted VCF in the cleanup stage
# 2) indexing original merged sorted VCF
# 3) changing bcftools merge to set ./. genotypes to 0/0 so AF can be calculated correctly downstream

# FIX:
# Logfile has basename in addition to date and is placed in the working directory


#######################################
# SETUP
#######################################


# Check input variables
if [ -n "$1" ];then
	echo Working directory is "$1" 
	else
	echo "Usage: $0 <working directory> <VCF list>"
	exit 1
fi

if [ -n "$2" ];then
	echo Input VCF list is "$2"
	else
	echo "Usage: $0 <working directory> <VCF list>"
	exit 1
fi


today=$(date '+%y%m%d')
echo $today


# Rename input variables
working_dir=$1
vcf_list=$2

# Bed file with coding regions (GENCODE v24 from UCSC table browser with exons + 8 bp splice regions)
coding_bed=/home/mmijuskovic/small_variant_freq/GENCODEv24.hg38.bed
echo "Using $coding_bed file for subsetting"

# Enter working directory
cd $working_dir



# Get basename from VCF list filename
base_name="$(echo $vcf_list | awk 'BEGIN{FS="/"} END{print $NF}' | awk 'BEGIN{FS="."} END{print $1}')"

# Start a log file
logfile="$today"-"$base_name"_AFtoSomaticVCF_v2_log.txt
touch $logfile
echo "Logfile is $logfile"
echo "$0" >> "$logfile"
echo $today >> "$logfile"
echo "Working directory is: $working_dir" >> "$logfile"


# Check that VCF list file exists
if [ ! -f $vcf_list ]; then
	echo "VCF list does not exist"
	exit 1
else
	echo "Using list of VCF files: $vcf_list" >> "$logfile"
fi

# Check that coding region bed file exists
if [ ! -f $coding_bed ]; then
	echo "Bed file does not exist"
	exit 1
else
	echo "Using $coding_bed file for subsetting" >> "$logfile"
fi



#######################################
# MERGE VCFs
#######################################

# Read VCF paths from the input list

echo "Files in VCF list:" >> "$logfile"
for i in $(cat $vcf_list); do echo $i >> "$logfile"; done;

# Make a variable with all VCF paths from the list	
vcf_paths=$(for i in $(cat $vcf_list); do echo $i; done)

# Merge VCFs using only PASS filter variants and set missing genotypes to reference ("-0")
echo "Merging VCFs ...PASS filter only.... " >> "$logfile"
#bcftools merge -m none -f PASS $vcf_paths -o "$base_name".PASSonly.merged.vcf
bcftools merge -m none -f PASS -0 $vcf_paths -o "$base_name".PASSonly.merged.vcf

# Sort VCF
echo "Sorting merged VCF .... " >> "$logfile"
vt sort "$base_name".PASSonly.merged.vcf -o "$base_name".PASSonly.merged.sorted.vcf


#######################################
# Add AF to INFO field of merged VCF
#######################################

# Calculate AF
echo "Calculating AF in merged VCF .... " >> "$logfile"
vt estimate -e AF "$base_name".PASSonly.merged.sorted.vcf > "$base_name".PASSonly.merged.AF.vcf

# Remove samples (individual genotypes)
echo "Removing samples .... " >> "$logfile"
bcftools view -G "$base_name".PASSonly.merged.AF.vcf > "$base_name".PASSonly.merged.AF.info_only.vcf

# bgzip and index
echo "Indexing .... " >> "$logfile"
bgzip "$base_name".PASSonly.merged.AF.info_only.vcf
tabix "$base_name".PASSonly.merged.AF.info_only.vcf.gz


#######################################
# Subset to coding regions only
#######################################

# Subset the merged file (no samples) to coding regions only
echo "Subsetting for coding regions .... " >> "$logfile"
echo "Bed file used for subsetting is $coding_bed" >> "$logfile"
bcftools view -R $coding_bed "$base_name".PASSonly.merged.AF.info_only.vcf.gz > "$base_name".PASSonly.merged.AF.info_only.coding_only.vcf

# Sort, bgzip and index the subsetted file
echo "Sorting subsetted VCF.... " >> "$logfile"
vt sort "$base_name".PASSonly.merged.AF.info_only.coding_only.vcf -o "$base_name".PASSonly.merged.AF.info_only.coding_only.sorted.vcf

echo "Indexing .... " >> "$logfile"
bgzip "$base_name".PASSonly.merged.AF.info_only.coding_only.sorted.vcf
tabix "$base_name".PASSonly.merged.AF.info_only.coding_only.sorted.vcf.gz

#######################################
# CLEANUP
#######################################

echo "Cleaning up .... " >> "$logfile"
rm "$base_name".PASSonly.merged.vcf
#rm "$base_name".PASSonly.merged.sorted.vcf
rm "$base_name".PASSonly.merged.AF.info_only.coding_only.vcf

# bgzip and index original sorted merge
echo "Indexing original merged sorted VCF.... " >> "$logfile"
bgzip "$base_name".PASSonly.merged.sorted.vcf
tabix "$base_name".PASSonly.merged.sorted.vcf.gz








