#!/bin/bash -e
# Filename: cleanup_somaticVCF.sh
# Function: for given vcf.gz file,
#	    1)extract somatic calls only
#           2)add genotype for untyped sample
#           3)remove trailing prefix from sample name
#           4)calculate md5
# Usage	: $0 -v <CancerPoo-NormalBah.vcf.gz> -c <cancer_sample_name>
# note	: update md5sum.txt at the calling script
# 20150922 Mikyung Jang	1st working ver
#-----------------------------------------------------

declare -r SCRIPT_NAME=$(basename "$BASH_SOURCE" .sh)

## exit the shell(default status code: 1) after printing the message to stderr
bail() {
    echo -ne "$1" >&2
    exit ${2-1}
} 

## help message
declare -r HELP_MSG="Usage: $SCRIPT_NAME [OPTION]...
Description: this script extract somatic calls only, add genotype, replaces old sample name in vcf header, reindex
Options:
  -v	vcf file
  -c	cancer sample
  -b    bcftools path
  -t    tabix path
  -h	display this help and exit
"

## print the usage and exit the shell(default status code: 2)
usage() {
    declare status=2
    if [[ "$1" =~ ^[0-9]+$ ]]; then
        status=$1
        shift
    fi
    bail "${1}$HELP_MSG" $status
}

while getopts ":v:c:b:t:" opt; do
    case $opt in
        v) VCF_FILE=$OPTARG ;;
        c) CANCER_SAMPLE=$OPTARG ;;
        b) BCFTOOLS=$OPTARG ;;
        t) TABIX=$OPTARG ;;
        h) usage 0	;;
        \?) usage "Invalid option: -$OPTARG \n" ;;
    esac
done
[[ "$#" -lt 4 ]] && usage

if [ -z $BCFTOOLS ] || [ -z $TABIX ]; then
 source /etc/profile.d/modules.sh
 module load htslib
 module load bcftools
 BCFTOOLS=`which bcftools`
 TABIX=`which tabix`
fi
FBASE=${VCF_FILE%.vcf.gz}
OLD_FILE=$FBASE.vcf.orig.gz
OLD_INDEX=$FBASE.vcf.orig.gz.tbi

# 1. check existing orig file
# 1.1 if orig, IN<-orig
if [ -r $OLD_FILE ] && [ -r $OLD_INDEX ]; then
    cp $OLD_FILE $FBASE.ori.vcf.gz
    cp $OLD_INDEX $FBASE.ori.vcf.gz.tbi
    INFILE=$FBASE.ori.vcf.gz
    isOrig=1
elif [ ! -r $VCF_FILE ]; then
    exit "cannot read input file $VCF_FILE"
else
    INFILE=$VCF_FILE
fi
OUTFILE=$FBASE.single.vcf.gz

# 2. make new header
# 2.1 extract sample name from file and input param
#     - if sample <2, exit
#     - else get err_sample_name
$BCFTOOLS query -l $INFILE > $INFILE.spl
Spl_count=`wc -l $INFILE.spl|awk '{print $1}' `
if [ "$Spl_count" -lt 2 ]; then
    echo "sample count $Spl_count, assume no SAMPLE"
    cat $INFILE.spl; rm $INFILE.spl
    exit 0
fi
C_Sample_Old=`grep $CANCER_SAMPLE $INFILE.spl`

# 2.2 extract header - remove last 3 line, eg.
##$BCFTOOLS_viewVersion=1.2-dirty+htslib-1.2.1
##$BCFTOOLS_viewCommand=view -h CancerLP2000836-DNA_F01_NormalLP2000834-DNA_F01.somatic.SV.vcf.gz
#CHROM	POS ID	REF ALT	QUAL	FILTER	INFO	FORMAT	NormalID_LP2000834-DNA_F01  TumorID_LP2000836-DNA_F01
$BCFTOOLS view -h $INFILE |head -n -3 > $FBASE.header

# 2.3 add FORMAT/GT
# 2.4 replace sample name
FORMAT="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%b\n"
printf "$FORMAT" $CANCER_SAMPLE >> $FBASE.header

# 3. generate new vcf file
# 3.1)extract somatic calls only
# 3.2)reheader
# 3.3)clear faulty AC and AN
$BCFTOOLS view -s $C_Sample_Old -Ov $INFILE |
    $BCFTOOLS reheader -h $FBASE.header /dev/stdin |
    sed -e 's/;AC=0//' -e 's/;AN=0//' |
    $BCFTOOLS convert -Oz -o $OUTFILE /dev/stdin

if [ -n "$isOrig" ]; then
    rm $INFILE $INFILE.tbi
else
    echo "old file is backed up as $OLD_FILE" #--debug
    mv $VCF_FILE $OLD_FILE
    mv $VCF_FILE.tbi $OLD_INDEX
fi
mv $OUTFILE $VCF_FILE
rm $FBASE.header $INFILE.spl

# 4. index new vcf file
if [ -e $VCF_FILE.tbi ]; then rm $VCF_FILE.tbi; fi
$TABIX -p vcf $VCF_FILE
#chown pipeline.pipeline $VCF_FILE $VCF_FILE.tbi  #--activate when +permission

# 5. generate md5 for new vcf file
md5sum $VCF_FILE > $VCF_FILE.md5
md5sum $VCF_FILE.tbi >> $VCF_FILE.md5

