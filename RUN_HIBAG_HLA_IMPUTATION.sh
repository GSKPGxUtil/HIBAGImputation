#!/bin/csh -f

######################################################################################################
#
# RUN_HIBAG_HLA_IMPUTATION.sh
# Author: Judong Shen (judong.x.shen@gsk.com)
#
# DESCRIPTION: This script runs HIBAG imputation of classical HLA alleles using GWAS data.
# 
# INPUTS: 
# 1. Plink GWAS dataset
# 2. Ethnicity data including SUBJID and Ethnicity columns
# 
# DEPENDENCIES: 
# 1. PLINK (1.07)
# 2. CheckSNPOverlap.R (R scripts to check SNP overlapping betweem GWAS data in MHC region and the classifiers)
# 3. RaceSUBJID.R (R scripts to extract SUBJID based on unique race information from the Ethnicity file)
# 4. HLAImputation.R (R scripts to impute HLA alleles)
# 5. ResultSummary.R (R scripts to merge and summarize the imputed HLA alleles)
# 6. ResultConvert.R (R scripts to convert the merged imputed HLA alleles to the .info and .dose data)
#
# USAGE: ./RUN_HIBAG_HLA_IMPUTATION.sh PLINKDATAFILE (.bed/.bim/.fam) SUBJECTEthnitiyFILE (.txt)
# For example, under /GWD/bioinfo/projects/statgen/HIBAG_Classifiers/HIBAG_Pipeline/Data
#               export R_LIBS_USER=/home/jxs62889/R/library
#               nohup /GWD/bioinfo/projects/statgen/HIBAG_Classifiers/HIBAG_Pipeline/Scripts/RUN_HIBAG_HLA_IMPUTATION.sh pgx1111 Subject_Ethnitiy.txt > pgx1111_HIBAG.log &
#
######################################################################################################

echo "Start imputation work..."
date

set PLINK="/GWD/bioinfo/apps/bin/plink"
set RDIR="/GWD/bioinfo/tools/bin"
set SCRIPTDIR="/GWD/appbase/projects/statgen/GXapp/HIBAGImputation"
set R3="/GWD/appbase/projects/statgen/R3.0.0/R-3.0.0/bin/R"
setenv R_LIBS_USER /GWD/appbase/projects/statgen/GXapp/R-packages

# CHECK FOR DEPENDENCIES
if (! -e $PLINK) then
    echo "Please install PLINK (http://pngu.mgh.harvard.edu/~purcell/plink/download.shtml) first."; exit 1
else if (! -e $SCRIPTDIR/RUN_HIBAG_HLA_IMPUTATION.sh) then
    echo "Please copy RUN_HIBAG_HLA_IMPUTATION.sh into this the script directory."; exit 1
else if (! -e $SCRIPTDIR/CheckSNPOverlap.R) then
    echo "Please copy CheckSNPOverlap.R into this the script directory."; exit 1
else if (! -e $SCRIPTDIR/RaceSUBJID.R) then
    echo "Please copy RaceSUBJID.R into this the script directory."; exit 1
else if (! -e $SCRIPTDIR/HLAImputation.R) then
    echo "Please copy HLAImputation.R into this the script directory."; exit 1
else if (! -e $SCRIPTDIR/ResultSummary.R) then
    echo "Please copy ResultSummary.R into this the script directory."; exit 1
else if (! -e $SCRIPTDIR/ResultConvert.R) then
    echo "Please copy ResultConvert.R into this the script directory."; exit 1
else if (! -e $SCRIPTDIR/ExtractMatchingAlleles.R) then
    echo "Please copy ExtractMatchingAlleles.R into this the script directory."; exit 1
endif

if ($#argv != 2) then
    echo "USAGE: /GWD/bioinfo/projects/statgen/HIBAG_Classifiers/HIBAG_Pipeline/Scripts/RUN_HIBAG_HLA_IMPUTATION.sh DATA (.bed/.bim/.fam) ETHNICITY (.txt)"; exit 1
endif

# INPUTS
set INPUT=$1
set ETHNICITY=$2

# Functions to run
set SubsetToMHC=1
set CheckSNPOverlap=1
set RaceSUBJID=1
set EXTRACT_MHC=0
set EXTRACT_MATCHING_ALLELES=1
set IMPUTE=1
set SUMMARY=1
set CONVERT=1

if ($SubsetToMHC) then
    printf "Subsetting data to MHC ...\n"
    $PLINK --bfile $INPUT --chr 6 --from-bp 25651242 --to-bp 33544122 --make-bed --out $INPUT.MHC
    set INPUT=$INPUT.MHC
    echo "Subsetting data to MHC done!"
    echo "  "
endif

date

if ($CheckSNPOverlap) then
    printf "Check SNP overlapping between GWAS data and classifiers from the MHC ...\n"
    mkdir -p Results_CheckSNPOverlap
    $RDIR/R64-2.14.0 --vanilla --slave --args $INPUT $ETHNICITY < $SCRIPTDIR/CheckSNPOverlap.R
    echo "CheckSNPOverlap Done!"
    echo "  "
endif

date

if ($RaceSUBJID) then
    printf "Extract SUBJID based on unique race information from the Ethnicity file ...\n"
    mkdir -p ProcessedData
    $RDIR/R64-2.14.0 --vanilla --slave --args $ETHNICITY < $SCRIPTDIR/RaceSUBJID.R
    echo "RaceSUBJID Done!"
    echo "  "
endif

#Why doing this if advise user to do ahead? Risk to have these hard coded
#Currently hg18 coord lifted to hg19 but HIBAG authors used wider region for hg19
if ($EXTRACT_MHC) then
    printf "Extracting SNPs from the MHC for each race group ...\n"
    #    for file in ./ProcessedData/*.txt ; do
    #    $PLINK --bfile $INPUT --chr 6 --from-bp 25651263 --to-bp 33426849 --keep ./ProcessedData/$file --make-bed --out ./ProcessedData/$INPUT.$file.MHC
    #    done
    if (-e ./ProcessedData/European.txt) then
         $PLINK --bfile $INPUT --chr 6 --from-bp 25651263 --to-bp 33426849 --keep ./ProcessedData/European.txt --make-bed --out ./ProcessedData/$INPUT.European.MHC
    endif
    if (-e ./ProcessedData/Asian.txt) then
         $PLINK --bfile $INPUT --chr 6 --from-bp 25651263 --to-bp 33426849 --keep ./ProcessedData/Asian.txt --make-bed --out ./ProcessedData/$INPUT.Asian.MHC
    endif
    if (-e ./ProcessedData/Hispanic.txt) then
         $PLINK --bfile $INPUT --chr 6 --from-bp 25651263 --to-bp 33426849 --keep ./ProcessedData/Hispanic.txt --make-bed --out ./ProcessedData/$INPUT.Hispanic.MHC
    endif
    if (-e ./ProcessedData/African.txt) then
         $PLINK --bfile $INPUT --chr 6 --from-bp 25651263 --to-bp 33426849 --keep ./ProcessedData/African.txt --make-bed --out ./ProcessedData/$INPUT.African.MHC
    endif
    if (-e ./ProcessedData/Other.txt) then
         $PLINK --bfile $INPUT --chr 6 --from-bp 25651263 --to-bp 33426849 --keep ./ProcessedData/Other.txt --make-bed --out ./ProcessedData/$INPUT.Broad.MHC
    endif
    echo "EXTRACT_MHC Done!"
    echo "  "
endif

if ($EXTRACT_MATCHING_ALLELES) then
    printf "Extract matching alleles based on selected classifiers ...\n"
    $RDIR/R64-2.14.0 --vanilla --slave --args $INPUT < $SCRIPTDIR/ExtractMatchingAlleles.R
    echo "ExtractMatchingAlleles Done!"
    echo "  "
endif

date

if ($IMPUTE) then
    printf "Start to impute HLA alleles for each race ..."
    mkdir -p Results_ImputedHLAAlleles
    $RDIR/R64-2.14.0 --vanilla --slave --args $INPUT $ETHNICITY < $SCRIPTDIR/HLAImputation.R
    echo "IMPUTE Done!"
    echo " "
endif

date

if ($SUMMARY) then
    printf "Start to merge and summarize the imputed HLA alleles ..."
    mkdir -p Results_ImputedHLAAlleles_Summary
    $R3 --vanilla --slave --args $INPUT $ETHNICITY < $SCRIPTDIR/ResultSummary.R
    echo "SUMMARY Done!"
    echo "  "
endif

if ($CONVERT) then
    printf "Start to convert the imputed HLA data ..."
    mkdir -p Results_ImputedHLAAlleles_Converted
    $RDIR/R64-2.14.0 --vanilla --slave --args $INPUT $ETHNICITY < $SCRIPTDIR/ResultConvert.R
    gzip ./Results_ImputedHLAAlleles_Converted/*.*
    echo "CONVERT Done!"
    echo "  "
endif

echo "All imputation work done!"
date

