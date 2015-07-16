**HIBAG HLA Imputation Pipeline**

Judong Shen

The HIBAG HLA imputation pipeline is a tool that can automate the task of HLA genotype imputation using [HIBAG R package](http://www.bioconductor.org/packages/release/bioc/html/HIBAG.html). HIBAG is a state of the art software package for imputing HLA types using SNP data. It combines the concepts of attribute bagging, an ensemble classifier method, with haplotype inference for SNPs and HLA types. Attribute bagging is a technique which improves the accuracy and stability of classifier ensembles using bootstrap aggregating and random variable selection. HIBAG can be used by researchers with [published parameter estimates](http://www.biostat.washington.edu/~bsweir/HIBAG/#estimates) instead of requiring access to large training sample datasets. The published parameter estimates are stored in the published classifiers which were built from various genotyping arrays. All the genotype arrays have ancestry-specific models for European, Asian, Hispanic or African ancestry and most also have a mixed-ancestry model (“Broad”).

HIBAG HLA imputation pipeline is implemented by using shell scripts and R scripts. The shell scripts call all the R scripts and PLINK commands to automate the entire HLA imputation processes. The details of the HIBAG HLA imputation pipeline are as follows:
* A commented script called RUN_HIBAG_HLA_IMPUTATION.sh that can be executed to completely run the HLA imputation workflow. 
* Input of RUN_HIBAG_HLA_IMPUTATION.sh:
  * A set of GWAS data or xMHC GWAS data including the SNP genotype data (with rsids) in xMHC region (hg19, chr6: 025759242-033534827) in PLINK format (.bed/bim/fam).
  * A tab-delimited text file including exactly two columns “SUBJID” and “Ethnicity”. The subject ids included in this file should be exactly the same as those in the PLINK .fam file. The “Ethnicity” column can contain up to five ancestries: European, Asian, Hispanic, African and Other. 
* The shell scripts call all the R scripts and PLINK commands to automate the entire HLA imputation processes. 
  * CheckSNPOverlap: compare the SNPs included in the GWAS data and those included in all the available classifiers
  * RaceSUBJID: generate the subset of subject ids for each ethnicity and then extract its corresponding subset of GWAS or SNP data in xMHC region respectively
  * HLAImputation: impute HLA alleles for each of the seven HLA loci (HLA-A, -B, -C, -DRB1, -DQA1, -DQB1 and-DPB1) and each ethnicity respectively by using the corresponding optimal classifiers identified in CheckSNPOverlap step and the corresponding subset GWAS or SNP data generated in the RaceSUBJID step.
  * ResultSummary: merge and summarize the imputed HLA alleles. 
  * ResultConvert: recode the best guessed HLA genotypes into binary format for each unique HLA alleles and recode the carrier genotypes as 0/1/2 and 0/1 for additive and dominant genetics model respectively, which does not consider the imputation uncertainty. In addition, also calculate the dosage and other information such as allele frequency and Rsq values etc. for each unique HLA allele additive and dominant genetics model respectively by incorporating the imputation uncertainty.
