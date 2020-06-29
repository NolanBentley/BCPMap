readCallsAndCounts<-function(rawPath,
                             rowStart=3,
                             colStart=5,
                             rowsKept=0,
                             colEnd=0,
                             type="impute"){
    #Open and format calls and counts files
    #Make rowsKept NA where undesired
    no_col_calls <- max(count.fields(callStrPath, sep = ","))
    headerNames<-gsub("\"$","",gsub("^\"","",unlist(strsplit(readLines(callStrPath,n = 1),","))))
    headerSupplemental<-paste0("x",1:no_col_calls)
    headerNames<-c(headerNames,headerSupplemental[(1:length(headerSupplemental))>length(headerNames)])
    currCsv<-read.csv(rawPath,header = F,col.names = headerNames,colClasses = "character",fill = T,stringsAsFactors = F)
    dimnames(currCsv)<-list(currCsv[,1],headerNames)
    if(colEnd==0){colEnd=length(headerNames)}
    if(all(rowsKept==0)){
        currRows<-rowStart:nrow(currCsv)
    }else{
        if(!all(rownames(currCsv)==names(rowsKept))){stop("mismatched rows")}
        currRows<-which(currCsv[,1]==rowsKept)
    }
    currCsv<-currCsv[currRows,]
    currMarkers<-currCsv[,1]
    currChr    <-currCsv[,2]
    currPos    <-as.numeric(currCsv[,3])
    currRef    <-currCsv[,4]
    currCsv    <-currCsv[,colStart:colEnd]
    #If any numbers in top 10 rows, imputes as count files
    if(type=="impute"){
        if(mean(is.na(suppressWarnings(apply(currCsv[1:10,],c(1,2),as.numeric))))==1){
            type="calls"
        }else{
            type="counts"
        }
    }
    if(type=="calls"){
        return(list(calls=currCsv,
                    callsInfo=cbind(Markers=currMarkers,Chr=currChr,Pos=currPos,Ref=currRef)))
    }
    if(type=="counts"){
        leftInd<-seq(1,ncol(currCsv)-1,by = 2)
        return(list(leftCnt= currCsv[,leftInd]  ,
                    rightCnt=currCsv[,leftInd+1],
                    callsInfo=cbind(Markers=currMarkers,Chr=currChr,Pos=currPos,Ref=currRef)))
    }
}
