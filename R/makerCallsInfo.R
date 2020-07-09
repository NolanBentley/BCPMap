makeCallsInfo<-function(calls){
    callsInfo<-data.frame(stringsAsFactors = F,
                          Markers=rownames(calls))
    callsInfo$Chr<-gsub("_.*","",callsInfo$Markers)
    callsInfo$Pos<-as.numeric(gsub(".*_","",callsInfo$Markers))
    callsInfo<-cbind(callsInfo,callsList$Alleles)
    return(callsInfo)
}
