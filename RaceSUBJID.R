######################################################################################################
#
# RaceSUBJID.R
# Author: Judong Shen (judong.x.shen@gsk.com)
#
######################################################################################################

race.subjid <- function(race.file){
  race.d <- read.table(race.file,na.strings = c("",".",NA),check.name=FALSE,as.is=T, head=T, sep="")
  races <- unique(race.d$Ethnicity)
  for (race in races){
    race.sub <- subset(race.d, Ethnicity==race)
    race.sub <- data.frame(V1=race.sub$SUBJID,V2=race.sub$SUBJID)
    #out.file <- paste(race, "SUBJID.txt",sep="_")
    out.file <- paste(race, ".txt",sep="")
    out.file <- paste("./ProcessedData",out.file,sep="/")
    write.table(race.sub, out.file, quote = F, col.names = F, row.names = F, sep = "\t")
  }
}


myargs = commandArgs(TRUE)
#check all arguments specified
if (length(myargs) < 1) {
  cat("Invalid arguments, should be: \"--args [ancestry file] \"\n")
  q()
}
race.file <- myargs[1]
race.subjid(race.file) 