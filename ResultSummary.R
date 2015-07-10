######################################################################################################
#
# CheckSNPOverlap.R
# Author: Judong Shen (judong.x.shen@gsk.com)
# Merge and summarize the imputed HLA alleles
#
######################################################################################################

postprob.conv <- function(res){
  #in.data <- out
  no <- grep("prob", names(res))
  res1 <- res[,no]
  res2 <- sapply(c(seq(0.5,0.9,0.1),0.95), function(x){
    apply(res1, 2, function(y)sum(y>x))
  })
  res2 <- as.data.frame(res2)
  names(res2) <- as.character(c(seq(0.5,0.9,0.1),0.95))
  res2$locus <- rownames(res2)
  res2$locus <- gsub("_prob","", res2$locus)
  res2 <- res2[,c("locus", as.character(c(seq(0.5,0.9,0.1),0.95)))]
  res2$locus <- as.factor(res2$locus)
  rownames(res2) <- 1:nrow(res2)
  
  library(reshape2,lib.loc="/home/jxs62889/R/library")
  res2.lng <- melt(res2, id=c("locus"))
  res2.lng$perc <- res2.lng$value/nrow(res)
  res2.lng
}

postprob.plot <- function(res){
  #res <- out
  
  library(scales,lib.loc="/home/jxs62889/R/library")
  library(ggplot2,lib.loc="/home/jxs62889/R/library")
  library(reshape2,lib.loc="/home/jxs62889/R/library")

  all.data <- postprob.conv(res)
  p <- ggplot(aes(x=variable, y=perc, group=locus, colour=locus), data=all.data)
  p <- p + geom_line() +
    geom_point(size=3)+
    scale_color_discrete("HLA Locus") +
    scale_y_continuous(limits=c(0,1),
                       labels=percent) +
    labs(x="Posterior probability cutoff", y="Percent of subjects with Post. Prob. > Post. Prob. cutoff", title="All subjects")
  ggsave(file = "./Results_ImputedHLAAlleles_Summary/PosteriorProbabilityPlot_allSubjects.pdf",plot = p, width = 5, height = 5)
  
  res$Ancestry[res$Ancestry=="Broad"] <- "Other"
  races <- unique(res$Ancestry)
  if (length(races)>1){
    res.race.out <- list()
    k <- 0
    for (race in races){
      k <- k + 1 
      res.race <- subset(res, Ancestry==race)
      race.data <- postprob.conv(res.race)
      race.data$Ancestry <- race
      res.race.out[[k]] <- race.data
    }
    res.race.out <- do.call("rbind",res.race.out)
    
    p1 <- ggplot(aes(x=variable, y=perc, group=locus, colour=locus), data=res.race.out)
    p1 <- p1 + geom_line() +
      facet_grid(. ~ Ancestry) + 
      geom_point(size=3)+
      scale_color_discrete("HLA Locus") +
      scale_y_continuous(limits=c(0,1),
                         labels=percent) +
      labs(x="Posterior probability cutoff", y="Percent of subjects with Post. Prob. > Post. Prob. cutoff", title="")
    ggsave(file = "./Results_ImputedHLAAlleles_Summary/PosteriorProbabilityPlot_ByRace.pdf", plot=p1, w=10, h=5)
  } 
  
}

# res.merge: merge the imputed HLA alleles
res.merge <- function(pgx, race){
  race[race=="Other"] <- "Broad"
  out <- list()
  i <- 0
  for (r in race){
    i <- i + 1
    race.hla <- list()
      n <- 0
      for (hla.id in c("A", "B", "C", "DRB1", "DQA1", "DQB1", "DPB1")){
        file <- paste("./Results_ImputedHLAAlleles/",pgx,".",r,".MHC_","HLA-",hla.id,"_Imputation.guess.txt",sep="")
        hla.a <- read.delim(file, stringsAsFactors=F) 
        names(hla.a)[-1] <- paste("HLA-", hla.id, "_", c("a1","a2","prob"),sep="")
        
        hla.a$pgx <- pgx
        if (hla.id!="A"){
          n <- n + 1
          hla.a <- hla.a[,grep("^HLA",names(hla.a))]
          race.hla[[n]] <- hla.a
        } else {
          hla.a <- hla.a[,c("sample.id","pgx","HLA-A_a1","HLA-A_a2","HLA-A_prob")]
          hla.a1 <- hla.a
        }
      }
    race.hla <- do.call("cbind",race.hla)
    race.hla <- cbind(hla.a1, race.hla)    
    
    race.hla$Ancestry <- r
    #write.table(race.hla, file = paste("./Results_ImputedHLAAlleles_Summary/Imputed_HLAalleles_",r,sep=""),quote = F,col.names = T,na="", row.names = F, sep = "\t" )
    
    out[[i]] <- race.hla
  }
  out <- do.call("rbind",out)
  hla.c <- names(out)[grep("^HLA",names(out))]
  out <- out[,c("Ancestry","sample.id","pgx",hla.c)]
  out$Ancestry[out$Ancestry=="Broad"] <- "Other"
  write.table(out, file = "./Results_ImputedHLAAlleles_Summary/Imputed_HLAalleles_AllSubjects.txt",quote = F,col.names = T,na="", row.names = F, sep = "\t" )
  out
}

myargs = commandArgs(TRUE)
#check all arguments specified
if (length(myargs) < 2) {
  cat("Invalid arguments, should be: \"--args [root of plink dataset] [ancestry file] \"\n")
  q()
}
pgx0 <- myargs[1]
race.file <- myargs[2]
race.d <- read.table(race.file,na.strings = c("",".",NA),check.name=FALSE,as.is=T, head=T, sep="")
races <- unique(race.d$Ethnicity)
if ("Other" %in% races){
  races[races=="Other"] <- "Broad"
}

hla.res <- res.merge(pgx=pgx0, race=races)
postprob.plot(hla.res)


