#!/bin/bash -e
set -e

### Pre-processing somatic Strelka VCF for allele frequency calculations ###
#
# Adds GT (0/1) and PL (200,0,200) FORMAT fields to a one-sample VCF
# Modifies VCF header to include GT, PL declarations
# Sorts, bgzips, tabixes modified VCF
# Cleanup of intermediate files
#
# NEEDS: module load bcftools vt/0.577
#
# MMijuskovic, Nov 3 2017



#######################################
# SETUP
#######################################

today=$(date '+%y%m%d')
echo $today

# Check input variables
if [ -n "$1" ];then
	echo Working directory is "$1" 
	else
	echo "Usage: $0 <working directory> <VCF path>"
	exit 1
fi

if [ -n "$2" ];then
	echo Input VCF is "$2"
	else
	echo "Usage: $0 <working directory> <VCF path>"
	exit 1
fi


# Rename input variables
working_dir=$1
vcf_path=$2

# Start a log file
logfile="$today"_GTtoSomaticVCF_log.txt
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


# Get sample name from VCF path
sample_name="$(echo $vcf_path | awk 'BEGIN{FS="/"} END{print $NF}' | awk 'BEGIN{FS="."} END{print $1}')"




#######################################
# MODIFY VCF
#######################################

# Add GT,PL values and write modified VCF without header
echo "Adding GT and PL values .... " >> "$logfile"
bcftools view $vcf_path -H | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8, "GT:PL:"$9, "0/1:200,0,200:"$10}' > "$sample_name".noHead.GT.duprem.left.split.vcf

# Extract and modify header from the original VCF
echo "Modifying header .... " >> "$logfile"
bcftools view $vcf_path -h > "$vcf_path".header
# Add GT to header
sed '/FORMAT=<ID=DP,/i ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (0/1 means allele is observed in the tumour)">' "$vcf_path".header > "$vcf_path".header.mod1
# Add PL to header
sed '/FORMAT=<ID=DP,/i ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification. Manually set to HET.">' "$vcf_path".header.mod1 > "$vcf_path".header.mod2

# Attach modified header to modified VCF (bcftools reheader doesn't work with headerless VCF)
cat "$vcf_path".header.mod2 "$sample_name".noHead.GT.duprem.left.split.vcf > "$sample_name".withHead.GT.duprem.left.split.vcf 

# Sort VCF
echo "Sorting VCF .... " >> "$logfile"
vt sort "$sample_name".withHead.GT.duprem.left.split.vcf -o "$sample_name".GT.duprem.left.split.vcf

# bgzip and index
echo "Indexing .... " >> "$logfile"
bgzip "$sample_name".GT.duprem.left.split.vcf
tabix "$sample_name".GT.duprem.left.split.vcf.gz



#######################################
# CLEANUP
#######################################

echo "Cleaning up .... " >> "$logfile"
rm "$sample_name".noHead.GT.duprem.left.split.vcf
rm "$vcf_path".header
rm "$vcf_path".header.mod1
rm "$vcf_path".header.mod2
rm "$sample_name".withHead.GT.duprem.left.split.vcf










