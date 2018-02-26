#!/bin/bash -e
set -e

### Renames AF, AC, AN, NS according to the cohort name ###
#
# Modifies header of input VCF with INFO field descriptions (eg. AF_FFPE instead of AF)
# Changes the INFO field names (eg. instead of AF: AF_FFPE or AF_FFpcrfree or AF_FFnano)
#
# NEEDS: module load bcftools 
#
# MMijuskovic, Feb 26 2018



#######################################
# SETUP
#######################################

today=$(date '+%y%m%d')
echo $today

# Check input variables (expects working direcory, VCF path and a string used to modify the fields)

if [ -n "$1" ];then
	echo Working directory is "$1" 
	else
	echo "Usage: $0 <working directory> <VCF path> <cohort_name>"
	exit 1
fi

if [ -n "$2" ];then
	echo Input VCF is "$2"
	else
	echo "Usage: $0 <working directory> <VCF path> <cohort_name>"
	exit 1
fi

if [ -n "$3" ];then
	echo INFO fields AN AC AF NS will be modified using the cohort name "$3"
	else
	echo "Usage: $0 <working directory> <VCF path> <cohort_name>"
	exit 1
fi



# Rename input variables
working_dir=$1
vcf_path=$2
cohort_name=$3

# Get file name from VCF path
file_name="$(echo $vcf_path | awk 'BEGIN{FS="/"} END{print $NF}' | awk 'BEGIN{FS=".vcf"} END{print $1}')"


# Start a log file
logfile="$today"_"$file_name"_rename_AF_AC_AN_NS_log.txt
touch $logfile
echo "Logfile is $logfile"
echo "$0" >> "$logfile"
echo $today >> "$logfile"

# Enter working directory
cd $working_dir
echo "Working directory is: $working_dir" >> "$logfile"

# Check that VCF file exists
if [ ! -f $vcf_path ]; then
	echo "VCF does not exist"
	exit 1
else
	echo "Modifying VCF: $vcf_path" >> "$logfile"
fi





#######################################
# MODIFY VCF
#######################################

# Change AC, AF, AN, NS INFO field names and write the modified VCF without header
echo "Modifying INFO fields... " >> "$logfile"
bcftools view $vcf_path -H | sed -r "s/;AF=/;AF_${cohort_name}=/g" | sed -r "s/;AC=/;AC_${cohort_name}=/g" | sed -r "s/;AN=/;AN_${cohort_name}=/g" | sed -r "s/;NS=/;NS_${cohort_name}=/g" > "$vcf_path".renamed.noHeader.vcf

# Extract and modify header from the original VCF (some VCFs have different descr of AN and AC)
echo "Modifying header .... " >> "$logfile"
bcftools view $vcf_path -h | sed -r "s/INFO=<ID=AF,/INFO=<ID=AF_${cohort_name},/g" | sed -r "s/Alternate Allele Frequency/Alternate Allele Frequency in the ${cohort_name} cohort/g" \
| sed -r "s/INFO=<ID=AC,/INFO=<ID=AC_${cohort_name},/g" | sed -r "s/Alternate Allele Counts/Alternate Allele Counts in the ${cohort_name} cohort/g" \
| sed -r "s/Allele count in genotypes/Allele count in genotypes in the ${cohort_name} cohort/g" \
| sed -r "s/INFO=<ID=AN,/INFO=<ID=AN_${cohort_name},/g" | sed -r "s/Total Number Allele Counts/Total Number Allele Counts in the ${cohort_name} cohort/g" \
| sed -r "s/Total number of alleles in called genotypes/Total number of alleles in called genotypes in the ${cohort_name} cohort/g" \
| sed -r "s/INFO=<ID=NS,/INFO=<ID=NS_${cohort_name},/g" | sed -r "s/Number of Samples With Data/Number of Samples With Data in the ${cohort_name} cohort/g" \
> "$vcf_path".header

# Attach modified header to modified VCF (bcftools reheader doesn't work with headerless VCF)
cat "$vcf_path".header "$vcf_path".renamed.noHeader.vcf > "$vcf_path".renamed.vcf 

# bgzip and index
echo "Indexing .... " >> "$logfile"
bgzip "$vcf_path".renamed.vcf
tabix "$vcf_path".renamed.vcf.gz



#######################################
# CLEANUP
#######################################

echo "Cleaning up .... " >> "$logfile"
rm "$vcf_path".renamed.noHeader.vcf
rm "$vcf_path".header











