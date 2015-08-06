# HIBAG HLA Imputation Workflow

Judong Shen
Andrew Slater

This workflow implements HIBAG<sup>1</sup> to impute 4-digit classical HLA alleles from SNV genotypes in the xMHC region and convert the probabilities to a binary-expanded set of doses in [minimac](http://genome.sph.umich.edu/wiki/Minimac) format.  HIBAG was developed in a collaboration between GSK and the University of Washington which maintains the R package and hosts a series of pre-fit classification models on their [website](www.biostat.washington.edu/~bsweir/HIBAG).

Currently, the pre-fit models are all trained from a single reference dataset of individuals with both classical HLA genotypes (determined by direct assaying) and SNV genotypes from arrays of the Illumina 1M class (unclear which specific version(s)). To train each model, the SNV genotypes in the reference dataset were subset to the variants on the array of interest that are polymorphic in the ancestry group of interest. For example, the Asian model for the Affymetrix Genome-Wide Human SNP Array 5.0 was trained using the Asian individuals from the reference dataset, removing monomorphic SNVs in this subset of individuals and removing SNVs not assayed by the Affymetrix Genome-Wide Human SNP Array 5.0.

### High-level workflow overview

This workflow consists or a csh driver script which calls R scripts to perform the following steps:

1. Check the SNV overlap between the dataset to be imputed and each of the pre-fit models.
2. For each ancestry group and HLA locus, select the best pre-fit model from the results of (1).
3. For each ancestry group and HLA locus, predict the HLA genotypes using the model selected in (2).
4. Plot the probability distributions from (3).
5. Transform the genotype-level probabilities to binary-expanded doses and format as a pair of minimac dose / info files.

1) Zheng X, Shen J, Cox C, Wakefield J, Ehm M, Nelson M, Weir BS. HIBAG – HLA Genotype Imputation with Attribute Bagging. Pharmacogenomics Journal (2013). doi: 10.1038/tpj.2013.18.
