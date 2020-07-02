genoFreq<-function(boundCallList){
    alleles<-c("A","C","G","T","-")
    determined<-matrix(NA,nrow(boundCallList$Calls),ncol(boundCallList$Calls),dimnames = dimnames(boundCallList$Calls))
    callNChar <-nchar(boundCallList$Calls)
    determined[callNChar==2]<-match(boundCallList$Calls[callNChar==2],paste0(c(alleles,"?"),"/"))*11

    hetMatch   <-sapply(paste0(alleles,"/"),paste0,alleles)
    hetMatchNum<-sapply((1:length(alleles))*10,function(x,y)x+y,1:length(alleles))
    hetMatch   <-hetMatch   [c(which(lower.tri(hetMatch   )),which(upper.tri(hetMatch   )))]
    hetMatchNum<-hetMatchNum[c(which(lower.tri(hetMatchNum)),which(lower.tri(hetMatchNum)))]
    determined[callNChar==3]<-hetMatchNum[match(boundCallList$Calls[callNChar==3],hetMatch)]

    unlist(apply(head(determined,1000),2,table))
    cl<-min(makeCluster(36),detectCores-1)
    fullGenoFreq<-t(pbapply(cbind(matrix(c(11*(1:6),hetMatchNum),nrow(determined),length(c(11*(1:6),hetMatchNum)),byrow = T),
                                  determined),
                            MARGIN = 1,
                            FUN = table,
                            cl = cl)-1)
    stopCluster(cl)
    colnames(fullGenoFreq)<-c(alleles,"?",hetMatch)[order(c(11*(1:6),unique(hetMatchNum)))]
    return(fullGenoFreq)
}
