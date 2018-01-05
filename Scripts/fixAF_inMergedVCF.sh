#!/bin/bash -e
set -e

### Remove failed samples from merged VCF and calculate AF ###
#
# Removes samples from a list
# Adds AC, AN, NS, AF to INFO field
# Removes genotypes
# bgzips, tabixes modified VCF
# Cleanup of intermediate files
#
# NEEDS: module load bcftools vt/0.577
#
# MMijuskovic, Jan 5 2018



#######################################
# SETUP
#######################################


# Check input variables
if [ -n "$1" ];then
	echo Working directory is "$1" 
	else
	echo "Usage: $0 <working directory> <merged VCF with genotypes>"
	exit 1
fi

if [ -n "$2" ];then
	echo Input VCF is "$2"
	else
	echo "Usage: $0 <working directory> <merged VCF with genotypes>"
	exit 1
fi


today=$(date '+%y%m%d')
echo $today


# Rename input variables
working_dir=$1
vcf_input=$2

# Bed file with coding regions (GENCODE v24 from UCSC table browser with exons + 8 bp splice regions)
coding_bed=/home/mmijuskovic/small_variant_freq/GENCODEv24.hg38.bed
echo "Using $coding_bed file for subsetting"

# Enter working directory
cd $working_dir



# Get basename from VCF list filename
base_name="$(echo $vcf_input | awk 'BEGIN{FS="/"} END{print $NF}' | awk 'BEGIN{FS="."} END{print $1}')"

# Start a log file
logfile="$today"-"$base_name"_AFtoSomaticVCF_v2_log.txt
touch $logfile
echo "Logfile is $logfile"
echo "$0" >> "$logfile"
echo $today >> "$logfile"
echo "Working directory is: $working_dir" >> "$logfile"


# Check that VCF list file exists
if [ ! -f $vcf_input ]; then
	echo "VCF input file does not exist"
	exit 1
else
	echo "Using VCF input file: $vcf_input" >> "$logfile"
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
# 
# echo "Files in VCF list:" >> "$logfile"
# for i in $(cat $vcf_list); do echo $i >> "$logfile"; done;

# Make a variable with all VCF paths from the list	
# vcf_paths=$(for i in $(cat $vcf_list); do echo $i; done)
# 
# Merge VCFs using only PASS filter variants and set missing genotypes to reference ("-0")
# echo "Merging VCFs ...PASS filter only.... " >> "$logfile"
# bcftools merge -m none -f PASS -0 $vcf_paths -o "$base_name".PASSonly.merged.vcf
# 
# Sort VCF
# echo "Sorting merged VCF .... " >> "$logfile"
# vt sort "$base_name".PASSonly.merged.vcf -o "$base_name".PASSonly.merged.sorted.vcf
# 

###########################################
# Remove failed samples from the merged VCF
###########################################

# Remove samples


# Sort VCF




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








