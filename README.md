# HIBAG HLA Imputation Workflow

Judong Shen
Andrew Slater

This workflow implements HIBAG<sup>1</sup> to impute 4-digit classical HLA alleles from SNV genotypes in the xMHC region and convert the probabilities to a binary-expanded set of doses in [minimac](http://genome.sph.umich.edu/wiki/Minimac) format.  HIBAG was developed in a collaboration between GSK and the University of Washington which maintains the R package and hosts a series of pre-fit classification models on their [website](http://www.biostat.washington.edu/~bsweir/HIBAG).

Currently, the pre-fit models are all trained from a single reference dataset of individuals with both classical HLA genotypes (determined by direct assaying) and SNV genotypes from arrays of the Illumina 1M class (unclear which specific version(s)). To train each model, the SNV genotypes in the reference dataset were subset to the variants on the array of interest that are polymorphic in the ancestry group of interest. For example, the Asian model for the Affymetrix Genome-Wide Human SNP Array 5.0 was trained using the Asian individuals from the reference dataset, removing monomorphic SNVs in this subset of individuals and removing SNVs not assayed by the Affymetrix Genome-Wide Human SNP Array 5.0.

1) Zheng X, Shen J, Cox C, Wakefield J, Ehm M, Nelson M, Weir BS. HIBAG – HLA Genotype Imputation with Attribute Bagging. Pharmacogenomics Journal (2013). [doi: 10.1038/tpj.2013.18](http://dx.doi.org/10.1038/tpj.2013.18).

### High-level workflow overview

This workflow consists of a csh driver script which calls R scripts to perform the following steps:

1. Check the SNV overlap between the dataset to be imputed and each of the pre-fit models.
2. For each ancestry group and HLA locus, select the best pre-fit model from the results of (1).
3. For each ancestry group and HLA locus, predict the HLA genotypes using the model selected in (2).
4. Plot the probability distributions from (3).
5. Transform the genotype-level probabilities to binary-expanded doses and format as a pair of minimac dose / info files.


### Details of workflow - how to run

#### Assumptions / pre-requisites / how to prepare data to be imputed
* Data is plink binary format on GRCh37
  * The variants in the the xMHC region are represented in the bim file with '6' in column 1 and the GRCh37 coordinate in column 3.
  * The reference dataset used for training the pre-fit models did not contain indels so no need to worry about their representation (e.g. VCF conventions) as they won't be used.
  * Strand also does not matter as the reference dataset used for training was not resolved to any particular strand - a "hard alignment" will be done by matching alleles and dropping ambiguous 'A/T' and 'C/G' SNVs.
  * SNP name also does not matter as the "hard alignment" will be done by coordinate.
* An ancestry map file is included
  * Two space-delimited columns with headers "SUBJID" and "Ethnicity"
  * SUBJID should match the FID and IID in the fam file (i.e. FID must match IID)
  * Valid Ethnicity values are European, Asian, Hispanic, African, and Broad where Broad refers to the pre-fit models trained with all ancestry groups and recommended for best performance in mixed / unknown ancestry groups.
  * Current recommendation is to use self-reported ancestry as follows:
    * If ethnicity is Hispanic, ancestry is Hispanic. Otherwise, use this table to map race to ancestry:
Race | Ancestry
------------ | ------------
African American/African Heritage | African
American Indian or Alaskan Native | Broad
Asian - Central/South Asian Heritage | Broad
Asian - East Asian Heritage | Asian
Asian - Japanese Heritage | Asian
Asian - South East Asian Heritage | Asian
Native Hawaiian or Other Pacific Islander | Broad
White - Arabic/North African Heritage | European
White - White/Caucasian/European Heritage | European

#### Running workflow
* Call the driver script. Here is example command if your plink dataset is named PGxNNN.bed and your ancestry map file is named ancestry.txt and both are in the current directory. nohup is recommended as it will take several hours to run. In this example, stdout and stderr are re-directed to files in the current directory (re-running will overwrite these files).
```
  nohup /GWD/appbase/projects/statgen/GXapp/HIBAGImputation/RUN_HIBAG_HLA_IMPUTATION.sh PGxNNN ancestry.txt >myrun.out 2>myrun.err
```
Alternatively, to submit to SGE with e-mail notification
```
  qsub -N PGxNNN -q dl580 -b y -l mt=5G -m e -cwd \
  -e myrun.err \
  -o myrun.out \
  /GWD/appbase/projects/statgen/GXapp/HIBAGImputation/RUN_HIBAG_HLA_IMPUTATION.sh \
  PGxNNN ancestry.txt
```

#### What the workflow is doing - how to read outputs

