######################################################################################################
#
# ResultConvert.R
# Author: Judong Shen (judong.x.shen@gsk.com)
# Convert the best guessed HLA genotypes and its posterior prob. into the .info.gz and .dose.gz files
#
######################################################################################################

# info: generate the first 5 columns of .info data as expected in minimac output
# dose.data is dataframe of minimac-like dose file (first column is sample ID, second column is "DOSE")
info <- function(dose.data){
  #Determine alleles from column names - excluding first 2 columns
  alleles <- names(dose.data)
  alleles <- alleles[grep("HLA", alleles)]
  #Split locus and allele and keep allele only (locus will be in marker name - "SNP" column)
  Al1 <- unlist(lapply( strsplit(alleles,"\\*"),function(x)x[[2]]))
  #Dummy allele 2
  Al2 <- rep("X", length(alleles))
  #Calculate frequency of allele as average of doses (divide by 2 since diploid) rounding to 5 to match minimac output
  Freq1 <- round(as.numeric(0.5*apply(dose.data[,grep("HLA", names(dose.data))],2,mean,na.rm=T)),5)  
  #Calculate MAF based on Freq1, rounding to 5 to match minimac output
  MAF <- round(pmin(Freq1, 1-Freq1),5)
  #Return dataframe with first 5 columns of minimac info file
  info.out <- data.frame(SNP=alleles, Al1, Al2, Freq1, MAF,stringsAsFactors=FALSE)
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

# dosage3(): trimmed down version of dosage2 to calculate just quantitative doses under additive model
dosage3 <- function() {
	#List of probability files for all ancestry groups & loci
	prob.files = list.files(path="./Results_ImputedHLAAlleles",pattern=".prob.txt$",full.names=TRUE)
	
	#For each locus, join probabilities across ancestry groups
	for (hla.id in c("A", "B", "C", "DRB1", "DQA1", "DQB1", "DPB1")){
		#cat('processing locus: ',hla.id,'\n')
		#Subset of probability files for this locus
		prob.files.locus = prob.files[grep(paste('HLA-',hla.id,sep=""),prob.files)]
		#Read first file
		if (length(prob.files.locus) > 0) {
			prob.df = read.table(prob.files.locus[1],row.names=NULL)
		}
		#Remove first file from list
		prob.files.locus = prob.files.locus[-1]
		#For each subsequent file, read and merge by genotype (first column)
		while (length(prob.files.locus) > 0) {
			prob.df2 = read.table(prob.files.locus[1],row.names=NULL)
			prob.dfm = merge(prob.df,prob.df2,by = 1, all=TRUE)
			prob.df = prob.dfm
			prob.files.locus = prob.files.locus[-1]
		}

		#Unique alleles observed across all ancestry groups
		alleles = unique(unlist(strsplit(prob.df[,1],'/')))
		
		#Initialize dataframe of samples (rows) x binary-expanded allele doses (columns) to match minimac output
		sid = names(prob.df)[-1]
		sid = paste(sid,sid,sep="->")
		dose.locus = data.frame(sid = sid, dose = "DOSE", stringsAsFactors=FALSE)
		
		for (allele in alleles) {
			prob.allele = apply(prob.df[,-1],2,function(x) {
				#Sum probabilities of all genotypes with allele and add homozygote genotype so it counts twice
				#Rounding to 3 to match minimac output
				round(sum(sum(x[grep(allele,prob.df[,1])],na.rm = TRUE),x[grep(paste(allele,allele,sep='/'),prob.df[,1])],na.rm = TRUE),3)
			})
		#Append doses for this allele to dose.locus dataframe
		dose.locus = cbind(dose.locus,data.frame(prob.allele))
		}
		
		#Update column names of dose.locus dataframe to include alleles
		locus.coord = switch(hla.id, A="6:29911349:HLA-A",
			B="6:31323307:HLA-B",
			C="6:31238217:HLA-C",
			DPB1="6:33049341:HLA-DPB1",
			DQA1="6:32605398:HLA-DQA1",
			DQB1="6:32631702:HLA-DQB1",
			DRB1="6:32552086:HLA-DRB1", paste("HLA-",hla.id,sep=""))
		names(dose.locus) = c("sid","dose",paste(locus.coord,"*",alleles,sep=""))
		
		if (exists("dose")) {
			#Append this locus' doses to dose dataframe
			dose = merge(dose, dose.locus, by=c(1,2), all=TRUE)
		} else {
			#Use this first dose.locus to start dose dataframe
			dose = dose.locus
		}
	#cat('dim dose: ',dim(dose),'\n')
	}
	
	#Dataframe of binary-expanded alleles (rows) x metrics (columns) to match expected minimac output
	info.df = info(dose)
	
	#Unclear from minimac documentation how to calculate AvgCall
	info.df$AvgCall = "-"
	
	#Calculate rsq - need to confirm this formula
	xbar <- colMeans(dose[,-c(1,2)]/2)
	x2bar <- colMeans((dose[,-c(1,2)]/2)^2)
	#Rounding to 5 to match minimac output
	info.df$Rsq = round((x2bar - xbar^2)/(xbar*(1-xbar)),5)
	#Replacing NaNs due to division by 0 when xbar is 0 with 0
	info.df$Rsq[is.na(info.df$Rsq)] = 0
	
	#Set all other info columns to blank as irrelevant (comparisons of assayed genotypes to imputed)
	info.df = cbind(info.df, data.frame(Genotyped = "-", LooRsq = "-", EmpR = "-", EmpRsq = "-", Dose1 = "-", Dose2 = "-",stringsAsFactors=FALSE))
	
	#Write results
	info.file <- paste("./Results_ImputedHLAAlleles_Converted/Imputed_HLAalleles_AllSubjects_Additive.info",sep="")
	dose.file <- paste("./Results_ImputedHLAAlleles_Converted/Imputed_HLAalleles_AllSubjects_Additive.dose",sep="")

	#Prevent scientific notation:
	options(scipen=999)
	write.table(info.df, file = info.file, quote = F,col.names = T,na="", row.names = F, sep = "\t" )  
	write.table(dose, file = dose.file, quote = F,col.names = F,na="", row.names = F, sep = "\t" )
}
# dosage2(): Consider two situations by considering or not considering the uncertainty due to posterior prob. for both additive and dominant model
dosage2 <- function(pgx, race, model, uncertainty){
  #collated best-guess genotypes across all ancestry groups & loci
  file <- "./Results_ImputedHLAAlleles_Summary/Imputed_HLAalleles_AllSubjects.txt"
  res <- read.delim(file, stringsAsFactors=F) 
  #Only keep all genotypes
  res1 <- res[,grep("HLA",names(res))]
  res1 <- res1[,-grep("_prob",names(res1))]
  
  #Initialize lists to store coded doses for each race
  out.dose <- list()
  out.call <- list()
  i <- 0
  for (r in race){
    #cat(r, "\n")
    i <- i + 1
    #Initialize lists to store coded doses for each locus
    dose.locus <- list()
    call.locus <- list()
    n <- 0
    for (hla.id in c("A", "B", "C", "DRB1", "DQA1", "DQB1", "DPB1")){
      #cat(hla.id, "\n")
      n <- n + 1
      
      #Best-guess genotypes for this ancestry group and locus (sample order ?)
      if (r=="Broad"){
        res1.locus <- res[res$Ancestry=="Other",paste("HLA.",hla.id,c("_a1","_a2"),sep="")]
      } else {
        
        res1.locus <- res[res$Ancestry==r,paste("HLA.",hla.id,c("_a1","_a2"),sep="")]      
      }

      #Determine unique alleles observed for binary expansion
      #Missing alleles not represented in best-guess genotypes - when tabulating binary-expanded doses based on probabilities, won't sum to 2 as expected
      alleles <- unique(c(res1[,paste("HLA.",hla.id,"_a1",sep="")],res1[,paste("HLA.",hla.id,"_a2",sep="")]))
      #Table of probabilities of each genotype (rows) for each sample (columns)
      #Note, columns won't sum to 0 as expected because of precision of small numbers
      prob.file <- paste("./Results_ImputedHLAAlleles/",pgx,".",r,".",hla.id,"_","HLA-",hla.id,"_Imputation.prob.txt",sep="")
      prob.res <- read.delim(prob.file, stringsAsFactors=F)
      #All possible genotypes
      allele.pairs <- rownames(prob.res)
      
      #Initialize lists to store doses for this ancestry group & locus
      dosage <- list()
      call <- list()
      a <- 0
      for (allele in alleles){
        a <- a + 1
        #Across all possible genotypes, determine binary expanded dose for this allele
        #Initialize dose to 0 indicating genotype does not have current allele
        allele.pairs.dose <- rep(0, length(allele.pairs))
        #If genotype contains allele, change dose to 1
        allele.pairs.dose[grep(allele,allele.pairs)] <- 1
        if (model=="Additive"){
          #Bug? delim is / not .
          #If homozygote for allele, change dose to 2
          allele.pairs.dose[allele.pairs==paste(allele,allele,sep="/")] <- 2  
        }
        #For each sample, sum of probabilities of all genotypes w/o allele
        pr0 <- apply(prob.res,2,function(x){sum(x[allele.pairs.dose==0])})
        #For each sample, sum of probabilities of all genotypes with allele
        pr1 <- apply(prob.res,2,function(x){sum(x[allele.pairs.dose==1])})
        if (model=="Additive"){
          #For each sample, probability of homozygote genotype
          pr2 <- apply(prob.res,2,function(x){sum(x[allele.pairs.dose==2])})
          
          if (uncertainty){
            #List of binary-expanded doses for each allele
            dosage[[a]] <- pr2*2 + pr1*1 + pr0*0  
          } else {
            #Convert best-guess genotypes to discrete doses (0/1/2) for this allele
            dose.tmp <- apply(res1.locus,1,function(x){sum(x %in% allele)})
            #Assuming sample ID order same
            names(dose.tmp) <- names(pr0)
            #List of binary-expanded doses for each allele based on best-guess genotypes
            dosage[[a]] <- dose.tmp
          }
          #?Should be list of probabilities of the binary-expanded dose for each allele
          #Instead, single value for each allele - if intentional, what representing?
          call[[a]] <- max(pr2,pr1,pr0)
        } else {
          if (uncertainty){
            #List of binary-expanded doses for each allele
            dosage[[a]] <- pr1*1 + pr0*0 
          } else {
            #Convert best-guess genotypes to discrete doses (0/1) for this allele
            dose.tmp <- apply(res1.locus,1,function(x){sum(any(x %in% allele))})
            #Assuming sample ID order same
            names(dose.tmp) <- names(pr0)
            #List of binary-expanded doses for each allele based on best-guess genotypes
            dosage[[a]] <- dose.tmp 
          }          
          #?Should be list of probabilities of the binary-expanded dose for each allele
          #Instead, single value for each allele - if intentional, what representing?
          call[[a]] <- max(pr1,pr0)
        }
      }
      #Collapse list of binary-expanded doses for each allele to matrix of alleles (rows) x samples (Columns)
      dosage <- do.call(rbind, dosage)
      #Check each sample's by-allele doses sum to 2
      if (uncertainty) {
        apply(dosage,2,function(x) {
          #Due to precision of small numbers, will not sum to exactly 2
          if (sum(x) < 1.99) stop("Sum of doses across all alleles is not 2 for each sample in locus ", hla.id, " at least one is ", sum(x))
        })
      } else {
        apply(dosage,2,function(x) {
          if (sum(x) != 2) stop("Sum of doses across all alleles is not 2 for each sample in locus ", hla.id, " at least one is ", sum(x))
        })
      }
      #Transpose binary-expanded doses for each allele to dataframe of samples (rows) x alleles (columns)
      dosage <- as.data.frame(t(dosage))
      #Name columns with allele
      colnames(dosage) <- paste("HLA-", hla.id,"*",alleles,sep="")
      #List of binary-expanded dose dataframes for each locus
      dose.locus[[n]] <- dosage

      #Do same collapse, tranposition and allele naming of call (whatever value this is intended to be - still unclear from above)
      call <- do.call(rbind, call)
      call <- as.data.frame(t(call))
      colnames(call) <- paste("HLA-", hla.id,"*",alleles,sep="")
      call.locus[[n]] <- call      
      
    }
    #Combine list of binary-expanded dose dataframes for each locus into single dataframe of samples (rows) x alleles (columns)
    dose.locus <- do.call(cbind, dose.locus)
    #Same for call (whatever value this is intended to be - still unclear from above)
    call.locus <- do.call(cbind, call.locus)
    
    #Get sample IDs
    sid <- rownames(dose.locus)
    #Convert sample IDs to FID->IID to emulate minimac output
    sid <- paste(sid,sid,sep="->")
    #Initialize a dataframe with first column as sample ID and second as "DOSE" to emulate minimac dose file
    sid <- data.frame(SAMPLEID=sid, DOSE="DOSE")
    #Append binary-expanded doses
    dose.locus <- cbind(sid, dose.locus)
    
    #Calculate metrics required to emulate minimac info file
    info.locus <- info(dose.locus)

    #Don't think this is correct - unclear what AvgCall in minimac INFO file represents
    #info.locus$AvgCall <- as.numeric(call.locus)
    #From MaCH documentation, most likely related to Quality (average of posterior probabilities for most likely genotype over all samples)
    #But calculating this value from minimac prob files does not match AvgCall - in minimac, when low Rsq, AvgCall tends to match Freq1
    #Value not used in downstream analyses so leaving blank
    info.locus$AvgCall <- '-'
    
    Freq1 <- info.locus$Freq1
    if (model=="Additive"){
      #Need to divide dose by 2 before calculating variance else will get Rsq values > 1
      info.locus$Rsq <- as.numeric(apply(dose.locus[,grep("HLA",names(dose.locus))]/2,2,var))/(2*Freq1*(1-Freq1))  
    } else {
      #This calculation does not seem correct - may not be relevant for dominant model
      #Note, multiplying Freq1 by 2 since info function assumes additive and will return half of Freq1
      #But shouldn't the adjustment be made directly to info.locus$Freq1 since that is also output in the info file
      info.locus$Rsq <- as.numeric(apply(dose.locus[,grep("HLA",names(dose.locus))],2,var))/(2*Freq1*(1-2*Freq1))  
    }
    #Set remaining info columns to blank
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
    
    #List of minimac dose dataframe and call dataframe for each ancestry group
    out.dose[[i]] <- dose.locus
    out.call[[i]] <- call.locus
  }
  #Combine dose and call dataframes for each ancestry group
  out.dose <- do.call(rbind, out.dose)
  out.call <- do.call(rbind, out.call)
  #?
  out.call1 <- apply(out.call,2,max, na.rm=T)
  
  #Calculate metrics required to emulate minimac info file
  info.all <- info(out.dose)
  #See notes above in by-ancestry group info calculation for concerns over correctness
  #info.all$AvgCall <- as.numeric(out.call1)
  info.all$AvgCall <- '-'
  
  #See notes above in by-ancestry calculations
  Freq1 <- info.all$Freq1
  if (uncertainty){
    if (model=="Additive"){
      info.all$Rsq <- as.numeric(apply(out.dose[,grep("HLA",names(out.dose))]/2,2,var))/(2*Freq1*(1-Freq1))  
    } else {
      info.all$Rsq <- as.numeric(apply(out.dose[,grep("HLA",names(out.dose))],2,var))/(2*Freq1*(1-2*Freq1))
    }
    #Now that fixed, should not need this ceiling
    #info.all$Rsq <- pmin(info.all$Rsq, 1)
  } else {
    #Why NA? Possibly because this calculation is invalid when converting best-guess to discrete dosages (0/1/2) instead of using quantitative doses
    #Would '-' be more appropriate to match other columns' blank values?
    info.all$Rsq <- NA
  }
  #Set all remaining info columns to blank
  info.all$Genotyped <- info.all$LooRsq  <- info.all$EmpR <- info.all$EmpRsq <- info.all$Dose1 <- info.all$Dose2 <- "-"
  
  #Assemble relevant out file paths
  if (uncertainty){
    info.file1 <- paste("./Results_ImputedHLAAlleles_Converted/Imputed_HLAalleles_AllSubjects_",model,".info",sep="")
    dose.file1 <- paste("./Results_ImputedHLAAlleles_Converted/Imputed_HLAalleles_AllSubjects_",model,".dose",sep="")
  } else {
    info.file1 <- paste("./Results_ImputedHLAAlleles_Converted/Imputed_HLAalleles_AllSubjects_",model,"_BestGuessed.info",sep="")
    dose.file1 <- paste("./Results_ImputedHLAAlleles_Converted/Imputed_HLAalleles_AllSubjects_",model,"_BestGuessed.dose",sep="")
  }
  #Write results
  write.table(info.all, file = info.file1, quote = F,col.names = T,na="", row.names = F, sep = "\t" )  
  write.table(out.dose, file = dose.file1, quote = F,col.names = F,na="", row.names = F, sep = "\t" )
  invisible(list(dose=out.dose, info=info.all))
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


#Convert probabilities of genotypes to additive doses (range 0:2)
#Uses best-guess genotypes as starting point
#dosage2(pgx=basename(pgx0), race=races, model="Additive", uncertainty=T)
dosage3()

#Excluding - can easily derive from quantitative doses as do for minimac output where no best-guess output
#Convert best-guess genotypes to additive doses (0/1/2)
#dosage2(pgx=basename(pgx0), race=races, model="Additive", uncertainty=F)


#Excluding dominant mode due to concerns over correctness of Freq1 and Rsq calculations
#Can easily derive carrier status from additive doses as do for minimac output where no dominant output
#Convert best-guess genotypes to dominant doses (0/1)
#dosage2(pgx=basename(pgx0), race=races, model="Dominant", uncertainty=F)
#Convert probabilities of genotypes to dominant doses (range 0:1)
#dosage2(pgx=basename(pgx0), race=races, model="Dominant", uncertainty=T)
