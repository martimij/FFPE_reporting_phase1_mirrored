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

# Start a log file
logfile="$today"_AFtoSomaticVCF_log.txt
touch $logfile
echo "Logfile is $logfile"
echo "$0" >> "$logfile"
echo $today >> "$logfile"

# Enter working directory
cd $working_dir
echo "Working directory is: $working_dir" >> "$logfile"

# Check that VCF list file exists
if [ ! -f $vcf_list ]; then
	echo "VCF list does not exist"
	exit 1
else
	echo "Using list of VCF files: $vcf_list" >> "$logfile"
fi


# Get basename from VCF list filename
base_name="$(echo $vcf_list | awk 'BEGIN{FS="/"} END{print $NF}' | awk 'BEGIN{FS="."} END{print $1}')"



#######################################
# MERGE VCFs
#######################################

# Read VCF paths from the input list

echo "Files in VCF list:" >> "$logfile"
for i in $(cat $vcf_list); do echo $i >> "$logfile"; done;

# Make a variable with all VCF paths from the list	
vcf_paths=$(for i in $(cat $vcf_list); do echo $i; done)

# Merge VCFs
echo "Merging VCFs .... " >> "$logfile"
bcftools merge -m none $vcf_paths -o "$base_name".merged.vcf

# Sort VCF
echo "Sorting merged VCF .... " >> "$logfile"
vt sort "$base_name".merged.vcf -o "$base_name".merged.sorted.vcf



#######################################
# Add AF to INFO field of merged VCF
#######################################

# Calculate AF
echo "Calculating AF in merged VCF .... " >> "$logfile"
vt estimate -e AF "$base_name".merged.sorted.vcf > "$base_name".merged.AF.vcf

# Remove samples (individual genotypes)
echo "Removing samples .... " >> "$logfile"
bcftools view -G "$base_name".merged.AF.vcf > "$base_name".merged.AF.info_only.vcf

# bgzip and index
echo "Indexing .... " >> "$logfile"
bgzip "$base_name".merged.AF.info_only.vcf
tabix "$base_name".merged.AF.info_only.vcf.gz



#######################################
# CLEANUP
#######################################

echo "Cleaning up .... " >> "$logfile"
rm "$base_name".merged.vcf
rm "$base_name".merged.sorted.vcf










