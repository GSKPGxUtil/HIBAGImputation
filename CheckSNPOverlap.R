######################################################################################################
#
# CheckSNPOverlap.R
# Author: Judong Shen (judong.x.shen@gsk.com)
#
######################################################################################################

# overlap.check: checking the SNPs includes in the GWAS data and those included in all classifiers
# Make sure the GWAS data includes the SNPs with rsids, but not other forms of ids. 
overlap.check <- function(in.data, classifier.loc){
  bim.file <- paste(in.data, "bim",sep=".")
  bim.d <- read.table(bim.file,na.strings = c("",".",NA),check.name=FALSE,as.is=T, head=F, sep="")
  gwas.snps <- bim.d$V2
  n.data <- length(gwas.snps)
  
  classifiers <- list.files(classifier.loc)
  classifiers <- classifiers[grep("-HLA4-hg19.RData", classifiers)]
  i <- 0
  out <- list()
  for (classifier in classifiers){
    cat(classifier, "\n")
    i <- i + 1
    load(paste(classifier.loc, classifier,sep="/"))
    classifier1 <- gsub("-HLA4-hg19.RData","",classifier)
    classifier1 <- unlist(strsplit(classifier1,"-"))

    j <- 0
    out1 <- list()
    for (hla.id in c("A", "B", "C", "DRB1", "DQA1", "DQB1", "DPB1")){
      j <- j + 1
      if ( classifier %in%  classifiers[grep("^Illumina1M", classifiers)]){
        model.snps <- HLA4[[hla.id]]$snp.id
      } else {
        model.snps <- HLAModelList[[hla.id]]$snp.id
      }      

      n.model <- length(model.snps)
      n.overlap <- length(intersect(gwas.snps, model.snps))
      p.overlap <- n.overlap/n.model
      vec <- c(classifier,classifier1, hla.id, n.data, n.model, n.overlap, p.overlap)
      out1[[j]] <- vec
    }
    out1 <- do.call(rbind, out1)
    out[[i]] <- out1
  }
  out <- do.call(rbind, out)
  out <- data.frame(out)
  names(out) <- c("classifier","Platform", "Ethnicity", "Gene", "NumData","NumModel", "NumOverlap","Percent")
  out[,-c(1:4)] <- apply(out[,-c(1:4)],2,as.numeric)
  write.table(out, "./Results_CheckSNPOverlap/SNPOverlap.txt", quote = F, col.names = T, row.names = F, sep = "\t")
  out
}

# model.iden: identify classifiers
model.iden <- function(overlap.p, race.file){
  overlap.stat <- ddply(overlap.p, c("Platform"), summarize,
                 MeanPercent = max(Percent),
                 SDPercent = sd(Percent))
  good.Platform <- as.character(overlap.stat$Platform[which.max(overlap.stat$MeanPercent)])
  overlap.p1 <- subset(overlap.p, Platform==good.Platform)
  good.classifiers <- as.character(unique(overlap.p1$classifier))
  platforms <- unlist(lapply( strsplit(good.classifiers,"-"),function(x)x[[2]]))
  good.classifiers <- data.frame(classifier=good.classifiers, race=platforms)
  
  race.d <- read.table(race.file,na.strings = c("",".",NA),check.name=FALSE,as.is=T, head=T, sep="")
  races <- unique(race.d$Ethnicity)
  if ("Other" %in% races){
    races[races=="Other"] <- "Broad"
  }
  
  good.classifiers <- subset(good.classifiers, race %in% races)
  write.table(overlap.stat, "./Results_CheckSNPOverlap/SNPOverlap_byPlatform.txt", quote = F, col.names = T, row.names = F, sep = "\t")
  write.table(good.classifiers, "./Results_CheckSNPOverlap/SelectedClassifiers.txt", quote = F, col.names = T, row.names = F, sep = "\t")
  invisible(good.classifiers)
}


library(HIBAG)
library(plyr)

myargs = commandArgs(TRUE)
#check all arguments specified
if (length(myargs) < 1) {
  cat("Invalid arguments, should be: \"--args [PLINK Data] \"\n")
  q()
}
in.data <- myargs[1]
race.file <- myargs[2]
classifier.loc <- "/GWD/bioinfo/projects/statgen/HIBAG_Classifiers/HIBAG_allClassifiers"

overlap.res <- overlap.check(in.data=in.data, classifier.loc=classifier.loc)
model.iden(overlap.p=overlap.res, race.file=race.file)

