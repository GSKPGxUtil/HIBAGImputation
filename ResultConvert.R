######################################################################################################
#
# ResultConvert.R
# Author: Judong Shen (judong.x.shen@gsk.com)
# Convert the best guessed HLA genotypes and its posterior prob. into the .info.gz and .dose.gz files
#
######################################################################################################

# info: generate the first 5 columns of .info data
info <- function(dose.data){
  alleles <- names(dose.data)
  alleles <- alleles[grep("HLA", alleles)]
  Al1 <- unlist(lapply( strsplit(alleles,"\\*"),function(x)x[[2]]))
  Al2 <- rep("X", length(alleles))
  Freq1 <- as.numeric(0.5*apply(dose.data[,grep("HLA", names(dose.data))],2,mean,na.rm=T))  
  MAF <- pmin(Freq1, 1-Freq1)
  info.out <- data.frame(SNP=alleles, Al1, Al2, Freq1, MAF)
  info.out
}

# dosage(): consider the uncertainty due to posterior prob. for both additive and dominant model
dosage <- function(pgx, race, model){
  file <- "./Results_ImputedHLAAlleles_Summary/Imputed_HLAalleles_AllSubjects.txt"
  res <- read.delim(file, stringsAsFactors=F) 
  res1 <- res[,grep("HLA",names(res))]
  res1 <- res1[,-grep("_prob",names(res1))]
  
  out.dose <- list()
  out.call <- list()
  i <- 0
  for (r in race){
    i <- i + 1
    dose.locus <- list()
    call.locus <- list()
    n <- 0
    for (hla.id in c("A", "B", "C", "DRB1", "DQA1", "DQB1", "DPB1")){
      n <- n + 1
      alleles <- unique(c(res1[,paste("HLA.",hla.id,"_a1",sep="")],res1[,paste("HLA.",hla.id,"_a2",sep="")]))
      prob.file <- paste("./Results_ImputedHLAAlleles/",pgx,".",r,".MHC_","HLA-",hla.id,"_Imputation.prob.txt",sep="")
      prob.res <- read.delim(prob.file, stringsAsFactors=F)
      allele.pairs <- rownames(prob.res)
      
      dosage <- list()
      call <- list()
      a <- 0
      for (allele in alleles){
        a <- a + 1
        allele.pairs.dose <- rep(0, length(allele.pairs))
        allele.pairs.dose[grep(allele,allele.pairs)] <- 1
        if (model=="Additive"){
          allele.pairs.dose[allele.pairs==paste(allele,allele,sep=".")] <- 2  
        }
        
        pr0 <- apply(prob.res,2,function(x){sum(x[allele.pairs.dose==0])})
        pr1 <- apply(prob.res,2,function(x){sum(x[allele.pairs.dose==1])})
        if (model=="Additive"){
          pr2 <- apply(prob.res,2,function(x){sum(x[allele.pairs.dose==2])})
          dosage[[a]] <- pr2*2 + pr1*1 + pr0*0
          call[[a]] <- max(pr2,pr1,pr0)
        } else {
          dosage[[a]] <- pr1*1 + pr0*0
          call[[a]] <- max(pr1,pr0)
        }
      }
      dosage <- do.call(rbind, dosage)
      dosage <- as.data.frame(t(dosage))
      colnames(dosage) <- paste("HLA-", hla.id,"*",alleles,sep="")
      dose.locus[[n]] <- dosage
      
      call <- do.call(rbind, call)
      call <- as.data.frame(t(call))
      colnames(call) <- paste("HLA-", hla.id,"*",alleles,sep="")
      call.locus[[n]] <- call      
      
    }
    dose.locus <- do.call(cbind, dose.locus)
    call.locus <- do.call(cbind, call.locus)
    
    sid <- rownames(dose.locus)
    sid <- paste(sid,sid,sep="->")
    sid <- data.frame(SAMPLEID=sid, DOSE="DOSE")
    dose.locus <- cbind(sid, dose.locus)
    
    info.locus <- info(dose.locus)
    info.locus$AvgCall <- as.numeric(call.locus)
    Freq1 <- info.locus$Freq1
    if (model=="Additive"){
      info.locus$Rsq <- as.numeric(apply(dose.locus[,grep("HLA",names(dose.locus))],2,var))/(2*Freq1*(1-Freq1))  
    } else {
      info.locus$Rsq <- as.numeric(apply(dose.locus[,grep("HLA",names(dose.locus))],2,var))/(2*Freq1*(1-2*Freq1))  
    }
    
    info.locus$Genotyped <- info.locus$LooRsq  <- info.locus$EmpR <- info.locus$EmpRsq <- info.locus$Dose1 <- info.locus$Dose2 <- "-"
    
    if (r=="Broad"){
      r1 <- "Other"
      dose.file <- paste("./Results_ImputedHLAAlleles_Converted/", "Imputed_HLAalleles_",r1,"_",model,".dose",sep="")
      info.file <- paste("./Results_ImputedHLAAlleles_Converted/", "Imputed_HLAalleles_",r1,"_",model,".info",sep="")
    } else {
      dose.file <- paste("./Results_ImputedHLAAlleles_Converted/", "Imputed_HLAalleles_",r,"_",model,".dose",sep="") 
      info.file <- paste("./Results_ImputedHLAAlleles_Converted/", "Imputed_HLAalleles_",r,"_",model,".info",sep="")
    }
    
    #write.table(dose.locus, file = dose.file,quote = F,col.names = T,na="", row.names = F, sep = "\t" ) 
    #write.table(info.locus, file = info.file,quote = F,col.names = T,na="", row.names = F, sep = "\t" ) 
    
    out.dose[[i]] <- dose.locus
    out.call[[i]] <- call.locus
  }
  out.dose <- do.call(rbind, out.dose)
  out.call <- do.call(rbind, out.call)
  out.call1 <- apply(out.call,2,max, na.rm=T)
  
  info.all <- info(out.dose)
  info.all$AvgCall <- as.numeric(out.call1)
  Freq1 <- info.all$Freq1
  if (model=="Additive"){
    info.all$Rsq <- as.numeric(apply(out.dose[,grep("HLA",names(out.dose))],2,var))/(2*Freq1*(1-Freq1))  
  } else {
    info.all$Rsq <- as.numeric(apply(out.dose[,grep("HLA",names(out.dose))],2,var))/(2*Freq1*(1-2*Freq1))
  }
  
  info.all$Genotyped <- info.all$LooRsq  <- info.all$EmpR <- info.all$EmpRsq <- info.all$Dose1 <- info.all$Dose2 <- "-"
  
  dose.file1 <- paste("./Results_ImputedHLAAlleles_Converted/Imputed_HLAalleles_AllSubjects_",model,".dose",sep="")
  write.table(out.dose, file = dose.file1,quote = F,col.names = T,na="", row.names = F, sep = "\t" )  
  info.file1 <- paste("./Results_ImputedHLAAlleles_Converted/Imputed_HLAalleles_AllSubjects_",model,".info",sep="")
  write.table(info.all, file = info.file1, quote = F,col.names = T,na="", row.names = F, sep = "\t" )  
  
  invisible(list(dose=out.dose, info=info.all))
}

# dosage2(): Consider two situations by considering or not considering the uncertainty due to posterior prob. for both additive and dominant model
dosage2 <- function(pgx, race, model, uncertainty){
  file <- "./Results_ImputedHLAAlleles_Summary/Imputed_HLAalleles_AllSubjects.txt"
  res <- read.delim(file, stringsAsFactors=F) 
  res1 <- res[,grep("HLA",names(res))]
  res1 <- res1[,-grep("_prob",names(res1))]
  
  out.dose <- list()
  out.call <- list()
  i <- 0
  for (r in race){
    #cat(r, "\n")
    i <- i + 1
    dose.locus <- list()
    call.locus <- list()
    n <- 0
    for (hla.id in c("A", "B", "C", "DRB1", "DQA1", "DQB1", "DPB1")){
      #cat(hla.id, "\n")
      n <- n + 1
      
      if (r=="Broad"){
        res1.locus <- res[res$Ancestry=="Other",paste("HLA.",hla.id,c("_a1","_a2"),sep="")]
      } else {
        res1.locus <- res[res$Ancestry==r,paste("HLA.",hla.id,c("_a1","_a2"),sep="")]      
      }

      alleles <- unique(c(res1[,paste("HLA.",hla.id,"_a1",sep="")],res1[,paste("HLA.",hla.id,"_a2",sep="")]))
      prob.file <- paste("./Results_ImputedHLAAlleles/",pgx,".",r,".MHC_","HLA-",hla.id,"_Imputation.prob.txt",sep="")
      prob.res <- read.delim(prob.file, stringsAsFactors=F)
      allele.pairs <- rownames(prob.res)
      
      dosage <- list()
      call <- list()
      a <- 0
      for (allele in alleles){
        a <- a + 1
        allele.pairs.dose <- rep(0, length(allele.pairs))
        allele.pairs.dose[grep(allele,allele.pairs)] <- 1
        if (model=="Additive"){
          allele.pairs.dose[allele.pairs==paste(allele,allele,sep=".")] <- 2  
        }
        
        pr0 <- apply(prob.res,2,function(x){sum(x[allele.pairs.dose==0])})
        pr1 <- apply(prob.res,2,function(x){sum(x[allele.pairs.dose==1])})
        if (model=="Additive"){
          pr2 <- apply(prob.res,2,function(x){sum(x[allele.pairs.dose==2])})
          if (uncertainty){
            dosage[[a]] <- pr2*2 + pr1*1 + pr0*0  
          } else {
            dose.tmp <- apply(res1.locus,1,function(x){sum(x %in% allele)})
            names(dose.tmp) <- names(pr0)
            dosage[[a]] <- dose.tmp  
          }
          call[[a]] <- max(pr2,pr1,pr0)
        } else {
          if (uncertainty){
            dosage[[a]] <- pr1*1 + pr0*0 
          } else {
            dose.tmp <- apply(res1.locus,1,function(x){sum(any(x %in% allele))})
            names(dose.tmp) <- names(pr0)
            dosage[[a]] <- dose.tmp 
          }          
          call[[a]] <- max(pr1,pr0)
        }
      }
      dosage <- do.call(rbind, dosage)
      dosage <- as.data.frame(t(dosage))
      colnames(dosage) <- paste("HLA-", hla.id,"*",alleles,sep="")
      dose.locus[[n]] <- dosage
      
      call <- do.call(rbind, call)
      call <- as.data.frame(t(call))
      colnames(call) <- paste("HLA-", hla.id,"*",alleles,sep="")
      call.locus[[n]] <- call      
      
    }
    dose.locus <- do.call(cbind, dose.locus)
    call.locus <- do.call(cbind, call.locus)
    
    sid <- rownames(dose.locus)
    sid <- paste(sid,sid,sep="->")
    sid <- data.frame(SAMPLEID=sid, DOSE="DOSE")
    dose.locus <- cbind(sid, dose.locus)
    
    info.locus <- info(dose.locus)
    info.locus$AvgCall <- as.numeric(call.locus)
    Freq1 <- info.locus$Freq1
    if (model=="Additive"){
      info.locus$Rsq <- as.numeric(apply(dose.locus[,grep("HLA",names(dose.locus))],2,var))/(2*Freq1*(1-Freq1))  
    } else {
      info.locus$Rsq <- as.numeric(apply(dose.locus[,grep("HLA",names(dose.locus))],2,var))/(2*Freq1*(1-2*Freq1))  
    }
    
    info.locus$Genotyped <- info.locus$LooRsq  <- info.locus$EmpR <- info.locus$EmpRsq <- info.locus$Dose1 <- info.locus$Dose2 <- "-"
    
    if (r=="Broad"){
      r1 <- "Other"
      if (uncertainty){
        dose.file <- paste("./Results_ImputedHLAAlleles_Converted/", "Imputed_HLAalleles_",r1,"_",model,".dose",sep="")
        info.file <- paste("./Results_ImputedHLAAlleles_Converted/", "Imputed_HLAalleles_",r1,"_",model,".info",sep="")
      } else {
        dose.file <- paste("./Results_ImputedHLAAlleles_Converted/", "Imputed_HLAalleles_",r1,"_",model,"_BestGuessed.dose",sep="")
        info.file <- paste("./Results_ImputedHLAAlleles_Converted/", "Imputed_HLAalleles_",r1,"_",model,"_BestGuessed.info",sep="")
      }
    } else {
      if (uncertainty){
        dose.file <- paste("./Results_ImputedHLAAlleles_Converted/", "Imputed_HLAalleles_",r,"_",model,".dose",sep="") 
        info.file <- paste("./Results_ImputedHLAAlleles_Converted/", "Imputed_HLAalleles_",r,"_",model,".info",sep="")
      } else {
        dose.file <- paste("./Results_ImputedHLAAlleles_Converted/", "Imputed_HLAalleles_",r,"_",model,"_BestGuessed.dose",sep="") 
        info.file <- paste("./Results_ImputedHLAAlleles_Converted/", "Imputed_HLAalleles_",r,"_",model,"_BestGuessed.info",sep="")
      }
    }
    
    #write.table(dose.locus, file = dose.file,quote = F,col.names = T,na="", row.names = F, sep = "\t" ) 
    #write.table(info.locus, file = info.file,quote = F,col.names = T,na="", row.names = F, sep = "\t" ) 
    
    out.dose[[i]] <- dose.locus
    out.call[[i]] <- call.locus
  }
  out.dose <- do.call(rbind, out.dose)
  out.call <- do.call(rbind, out.call)
  out.call1 <- apply(out.call,2,max, na.rm=T)
  
  info.all <- info(out.dose)
  info.all$AvgCall <- as.numeric(out.call1)
  Freq1 <- info.all$Freq1
  if (uncertainty){
    if (model=="Additive"){
      info.all$Rsq <- as.numeric(apply(out.dose[,grep("HLA",names(out.dose))],2,var))/(2*Freq1*(1-Freq1))  
    } else {
      info.all$Rsq <- as.numeric(apply(out.dose[,grep("HLA",names(out.dose))],2,var))/(2*Freq1*(1-2*Freq1))
    }
    info.all$Rsq <- pmin(info.all$Rsq, 1)
  } else {
    info.all$Rsq <- NA
  }
  
  info.all$Genotyped <- info.all$LooRsq  <- info.all$EmpR <- info.all$EmpRsq <- info.all$Dose1 <- info.all$Dose2 <- "-"
  
  if (uncertainty){
    info.file1 <- paste("./Results_ImputedHLAAlleles_Converted/Imputed_HLAalleles_AllSubjects_",model,".info",sep="")
    dose.file1 <- paste("./Results_ImputedHLAAlleles_Converted/Imputed_HLAalleles_AllSubjects_",model,".dose",sep="")
  } else {
    info.file1 <- paste("./Results_ImputedHLAAlleles_Converted/Imputed_HLAalleles_AllSubjects_",model,"_BestGuessed.info",sep="")
    dose.file1 <- paste("./Results_ImputedHLAAlleles_Converted/Imputed_HLAalleles_AllSubjects_",model,"_BestGuessed.dose",sep="")
  }
  write.table(info.all, file = info.file1, quote = F,col.names = T,na="", row.names = F, sep = "\t" )  
  write.table(out.dose, file = dose.file1, quote = F,col.names = F,na="", row.names = F, sep = "\t" )
  invisible(list(dose=out.dose, info=info.all))
}


myargs = commandArgs(TRUE)
#check all arguments specified
if (length(myargs) < 1) {
  cat("Invalid arguments, should be: \"--args [PLINK Data] \"\n")
  q()
}
pgx0 <- myargs[1]
race.file <- myargs[2]
race.d <- read.table(race.file,na.strings = c("",".",NA),check.name=FALSE,as.is=T, head=T, sep="")
races <- unique(race.d$Ethnicity)
if ("Other" %in% races){
  races[races=="Other"] <- "Broad"
}

dosage2(pgx=pgx0, race=races, model="Additive", uncertainty=F)
dosage2(pgx=pgx0, race=races, model="Dominant", uncertainty=F)
dosage2(pgx=pgx0, race=races, model="Additive", uncertainty=T)
dosage2(pgx=pgx0, race=races, model="Dominant", uncertainty=T)
