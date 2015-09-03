######################################################################################################
#
# ExtractMatchingAlleles.R
# Author: Andrew Slater (andrew.j.slater@gsk.com)
#
######################################################################################################

extract.matching <- function(plink.root){
  
  #Read selected classifiers
  classifiers <- read.table("./Results_CheckSNPOverlap/SelectedClassifiers.txt", as.is=TRUE, head=TRUE, sep="\t")
  
  #For each ancestry group and locus, plink extract IDs of markers with overlapping alleles and subset to subjects in group
  for (i in 1:nrow(classifiers)) {
    plink.success <- system(paste("/GWD/bioinfo/apps/bin/plink --bfile ", plink.root, " --extract ./Results_CheckSNPOverlap/", classifiers[i,]$Classifier, 
      ".extract.IDs.txt --keep ./ProcessedData/", classifiers[i,]$Ancestry,".txt --make-bed --out ./ProcessedData/", 
      basename(plink.root), ".", classifiers[i,]$Ancestry, ".", classifiers[i,]$HLA.locus, sep=""))
    if (plink.success != 0) {
      stop("plink extract failed")
    }
  }
}


myargs = commandArgs(TRUE)
#check all arguments specified
if (length(myargs) < 1) {
  cat("Invalid arguments, should be: \"--args [root of plink dataset] \"\n")
  q()
}
in.data <- myargs[1]
extract.matching(plink.root=in.data)
