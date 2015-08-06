# HIBAG HLA Imputation Workflow

Judong Shen & Andrew Slater

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
  * The driver script will subset the data to the xMHC before starting so no need to do this yourself.
  * CAUTION: Do not use any datasets downstream of the whole-genome imputation "alignment" step (imputation-Makefile-v2) as this will remove all variants not in the reference haplotypes an many of the xMHC variants are missing from these haplotypes due to difficulty in sequencing this region.
* An ancestry map file is included
  * Two space-delimited columns with headers "SUBJID" and "Ethnicity"
  * SUBJID should match the FID and IID in the fam file (i.e. FID must match IID)
  * Valid Ethnicity values are European, Asian, Hispanic, African, and Broad where Broad refers to the pre-fit models trained with all ancestry groups and recommended for best performance in mixed / unknown ancestry groups.
  * Current recommendation is to use self-reported ancestry as follows:
    * If ethnicity is Hispanic, ancestry is Hispanic. Otherwise, use this table to map race to ancestry:

          Race | Ancestry
          ---- | --------
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
Call the driver script. Below is an example command if your plink dataset is named PGxNNN.bed, your ancestry map file is named ancestry.txt and both are in the current directory. nohup is recommended as it will take several hours to run. In this example, stdout and stderr are re-directed to files in the current directory (re-running will overwrite these files).
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
The workflow proceeds sequentially (cost/benefit of adding parallel computing support is unclear) as follows:

1. Subset the data to the xMHC region creating a plink dataset in the same directory using the same name with '.MHC' appended. The boundaries of this region are defined within the HIBAG package and for GRCh37 are 6:25651242-33544122. Note, this does not correspond to a lift-over of the NCBI36 region referenced in the original HIBAG paper. Future upgrades to GRCh38 should look for the region definition within the HIBAG package and not rely on a lifting of these coordinates.
2. For each ancestry group present in the data per the ancestry map file, iterate over the relevant pre-fit models in /GWD/appbase/projects/RD-MDD-GX_PUBLIC/HIBAG_Classifiers (downloaded July 2, 2015) where the name of the ancestry group is present in the file name. For each model and locus:
  1. Determine the coordinate and alleles of the SNVs and their contribution to the model.
  2. Remove any ambiguous 'A/T' and 'C/G' SNVs.
  3. Remove any SNVs whose alleles do not match those observed in the data for the same coordinate.
  4. Using SNVs remaining, summarize overlap with data and append to Results_CheckSNPOverlap/comparison.txt.
    * TO DO - describe metrics reported in this file
  5. Append list of names of overlapping SNVs in data to Results_CheckSNPOverlap/[Model File Name].extract.IDs.txt.
3. Using the overlap summary in Results_CheckSNPOverlap/comparison.txt, for each ancestry group and locus, select the optimal model by ranking on both mean.accuracy (higher = better rank) and sum.miss.pctl (lower = better rank) and summing these two ranks to get an overall rank. Ties in this overall rank are broken by selecting max pct.model.in.data. Selected models are written to Results_CheckSNPOverlap/SelectedClassifiers.txt.
4. Parse the ancestry map file into a series of plink --keep files, one for each ancestry group in ProcessedData/.
5. For each ancestry group and locus, create a subset of the data in ProcessedData/ be restricting to: 
  * the individuals within the group using the --keep file in ProcessedData/
  * the SNVs in the selected model using the --extract file in Results_CheckSNPOverlap/ 
6. For each ancestry group and locus, use the relevant data subset in ProcessedData/ and the selected model to impute HLA genotypes to Results_ImputedHLAAlleles/.
7. Plot the distribution of the probabilities of the "best guess" genotypes (highest probability) by locus across all ancestry groups to Results_ImputedHLAAlleles_Summary/PosteriorProbabilityPlot_allSubjects.pdf and by locus and ancestry group to Results_ImputedHLAAlleles_Summary/PosteriorProbabilityPlot_ByRace.pdf.
8. Convert the probabilities to binary-expanded additive doses in minimac format:
  * Initialize info file by creating a binary expanded marker for each allele:
    * SNP = [Coord]_[Locus]*[allele]     where Coord is the midpoint of the locus
    * Al1 = [allele]
    * Al2 = X
  * Sum the probabilities for each allele over all possible genotypes counting the homozygote genotype twice to obtain an additive dose (i.e. scaled from 0-2 like minimac).
  * Calculate Freq1 (half the mean of the doses), MAF and Rsq to finish the info file (formula for calculating AvgCall is unknown).
  * gzip files in Results_ImputedHLAAlleles_Converted/.
