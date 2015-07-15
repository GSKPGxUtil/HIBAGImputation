######################################################################################################
#
# HLAImputation.R
# Author: Judong Shen (judong.x.shen@gsk.com)
#
######################################################################################################

paste0 <- function(...) paste(..., sep="")  

hla.imp <- function(pgx,race,model,hla.id){

  cat(paste(pgx, race, model, sep=", "),"\n")
  # Load the published parameter estimates from classifier
  model.list <- get(load(model))
  
  #########################################################################
  # Import your PLINK BED file
  #
  bed.f <- paste("./ProcessedData/", pgx,".bed",sep="")
  fam.f <- paste("./ProcessedData/", pgx,".fam",sep="")
  bim.f <- paste("./ProcessedData/", pgx,".bim",sep="")
  
  geno <- hlaBED2Geno(bed.fn=bed.f, fam.fn=fam.f, bim.fn=bim.f, assembly=c("hg19"))
  summary(geno)
  cat("\n")
  
    # HLA imputation 
    model <- hlaModelFromObj(model.list[[hla.id]])
    perf = summary(model)
    cat("\n")

    # SNPs in the model
    cat("SNPs in the model:", "\n")
    # Not all SNPs in model are used
    length(perf$snp.position[perf$snp.hist > 0])

    # Intersecting SNPs in the model and the input geno data
    cat("SNP ID #: ",length(intersect(model$snp.id, geno$snp.id)), "\n")
    cat("SNP ID: ",(intersect(model$snp.id, geno$snp.id)), "\n")
    cat("SNP Position #: ",length(intersect(model$snp.position, geno$snp.position)), "\n")
    cat("SNP Position: ",(intersect(model$snp.position, geno$snp.position)), "\n")
    cat("\n")

    log.file <- paste("./Results_ImputedHLAAlleles/",pgx,"_","HLA-",hla.id,"_imputation.log",sep="")
    sink(log.file) # redirect console output to a file
    # best-guess genotypes
    pred.guess <- predict(model, geno, type="response", match.type="Position", allele.check=TRUE)
    cat("\n")
    # posterior probabilities
    pred.prob <- predict(model, geno, type="prob", match.type="Position", allele.check=TRUE)
    cat("\n")

    name1 <- paste("./Results_ImputedHLAAlleles/",pgx,"_","HLA-",hla.id,"_ImputationResults.RData",sep="")
    name2 <- paste("./Results_ImputedHLAAlleles/",pgx,"_","HLA-",hla.id,"_Imputation.prob.txt",sep="")
    name3 <- paste("./Results_ImputedHLAAlleles/",pgx,"_","HLA-",hla.id,"_Imputation.guess.txt",sep="")
    save(pred.guess, pred.prob, file=name1)
    write.table(pred.prob, file = name2,quote = F,col.names = T,na="", row.names = T, sep = "\t" )
    write.table(pred.guess$value, file = name3,quote = F,col.names = T,na="", row.names = F, sep = "\t" )
    sink()
  
  #Close model to reclaim memory
  hlaClose(model)

}

library(HIBAG)

myargs = commandArgs(TRUE)
#check all arguments specified
if (length(myargs) < 1) {
  cat("Invalid arguments, should be: \"--args [root of plink dataset] \"\n")
  q()
}
pgx0 <- myargs[1]

classifier.loc <- "/GWD/appbase/projects/RD-MDD-GX_PUBLIC/HIBAG_Classifiers"

classifiers <- read.table("./Results_CheckSNPOverlap/SelectedClassifiers.txt", as.is=TRUE, head=TRUE, sep="\t")


for (i in 1:nrow(classifiers)){
  pgx <- paste(basename(pgx0), ".", classifiers[i,]$Ancestry, ".", classifiers[i,]$HLA.locus, sep="")
  model.run <- paste(classifier.loc, classifiers[i,]$Classifier, sep="/")
  hla.imp(pgx=pgx, race=classifiers[i,]$Ancestry, model=model.run, hla.id=classifiers[i,]$HLA.locus)
}
