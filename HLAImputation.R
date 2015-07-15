######################################################################################################
#
# HLAImputation.R
# Author: Judong Shen (judong.x.shen@gsk.com)
#
######################################################################################################

paste0 <- function(...) paste(..., sep="")  

hla.imp <- function(pgx,race,model){
  library(HIBAG)
  cat(paste(pgx, race, model, sep=", "),"\n")
  # Load the published parameter estimates from European ancestry
  model.list <- get(load(model))
  
  #########################################################################
  # Import your PLINK BED file
  #
  bed.f <- paste("./ProcessedData/", pgx,".bed",sep="")
  fam.f <- paste("./ProcessedData/", pgx,".fam",sep="")
  bim.f <- paste("./ProcessedData/", pgx,".bim",sep="")
  geno <- hlaBED2Geno(bed.fn=bed.f, fam.fn=fam.f, bim.fn=bim.f)
  summary(geno)
  cat("\n")
  
  # HLA imputation for each HLA allele
  # hla.id <- "A"   # "B", "C", "DRB1", "DQA1", "DQB1", "DPB1"
  for (hla.id in c("A", "B", "C", "DRB1", "DQA1", "DQB1", "DPB1")){
    model <- hlaModelFromObj(model.list[[hla.id]])
    summary(model)
    cat("\n")
    
    # SNPs in the model
    cat("SNPs in the model:", "\n")
    head(model$snp.id)
    head(model$snp.position)
    
    # Intersecting SNPs in the model and the input geno data
    cat("SNP ID #: ",length(intersect(model$snp.id, geno$snp.id)), "\n")
    cat("SNP ID: ",(intersect(model$snp.id, geno$snp.id)), "\n")
    cat("SNP Position #: ",length(intersect(model$snp.position, geno$snp.position)), "\n")
    cat("SNP Position: ",(intersect(model$snp.position, geno$snp.position)), "\n")
    cat("\n")
    
    log.file <- paste("./Results_ImputedHLAAlleles/",pgx,"_","HLA-",hla.id,"_imputation.log",sep="")
    sink(log.file) # redirect console output to a file
    # best-guess genotypes
    pred.guess <- predict(model, geno, type="response")
    cat("\n")
    # posterior probabilities
    pred.prob <- predict(model, geno, type="prob")
    cat("\n")
    
    name1 <- paste("./Results_ImputedHLAAlleles/",pgx,"_","HLA-",hla.id,"_ImputationResults.RData",sep="")
    name2 <- paste("./Results_ImputedHLAAlleles/",pgx,"_","HLA-",hla.id,"_Imputation.prob.txt",sep="")
    name3 <- paste("./Results_ImputedHLAAlleles/",pgx,"_","HLA-",hla.id,"_Imputation.guess.txt",sep="")
    save(pred.guess, pred.prob, file=name1)
    write.table(pred.prob, file = name2,quote = F,col.names = T,na="", row.names = T, sep = "\t" )
    write.table(pred.guess$value, file = name3,quote = F,col.names = T,na="", row.names = F, sep = "\t" )
    sink()
  }
}

library(HIBAG)

myargs = commandArgs(TRUE)
#check all arguments specified
if (length(myargs) < 1) {
  cat("Invalid arguments, should be: \"--args [root of plink dataset] \"\n")
  q()
}
pgx0 <- myargs[1]
race.file <- myargs[2]
race.d <- read.table(race.file,na.strings = c("",".",NA),check.name=FALSE,as.is=T, head=T, sep="")
races <- unique(race.d$Ethnicity)
if ("Other" %in% races){
  races[races=="Other"] <- "Broad"
}
classifier.loc <- "/GWD/appbase/projects/RD-MDD-GX_PUBLIC/HIBAG_Classifiers"

for (race in races){
  pgx <- paste(pgx0,race,"MHC",sep=".")
  model.file <- "./Results_CheckSNPOverlap/SelectedClassifiers.txt"
  model <- read.table(model.file,na.strings = c("",".",NA),check.name=FALSE,as.is=T, head=T, sep="")
  model.run <- model$classifier[match(race, model$race)]
  model.run <- paste(classifier.loc, model.run, sep="/")
  hla.imp(pgx=pgx, race=race, model=model.run)
}





