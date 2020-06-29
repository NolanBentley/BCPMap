#Functions meant to subset call and count files to common marker sets
SubsetCallsAndCountsFiles<-function(fileDf){
    #Cut off the first column of each file and create a matrix for joining them
    if(!is.data.frame(fileDf)){stop("fileDf must be dataframe")}
    if(!(nrow(fileDf)>1&ncol(fileDf)>1)){stop("fileDf must be dataframe with more than 1 row and column where the columns")}

    markerList<-sapply(fileDf[,1],cutFunction)
    markerDf<-data.frame(ds1=1:length(markerList[[1]]))
    dimnames(markerDf)<-list(markerList[[1]],fileDf[1,1])
    markerDf[,2:length(markerList)]<-lapply(markerList[2:length(markerList)],match,x = markerList[[1]])
    markerDf<-markerDf[rowSums(is.na(markerDf))==0,]

    #Read and subset the files based on common markers
    for(rowIter in 1:nrow(fileDf)){
        for(colIter in 1:ncol(fileDf)){
            print(paste0("starting ",rowIter,"_",colIter,"..."))
            currFile<-fileDf[rowIter,colIter]
            outFile<-paste0(gsub("\\.csv","",currFile),"_subset.csv")
            currLines<-readLines(currFile)[markerDf[,rowIter]]
            writeLines(currLines,con = outFile)
        }
    }
}

cutFunction<-function(path,columnNum=1,delim=","){
    if(!file.exists(path)){stop("File doesn't exist")}
    if(!(!is.na(as.numeric(columnNum))&length(columnNum)==1)){stop("Only 1 column can be extracted in current implementation")}
    output<-system(paste0("cut -f ",columnNum," -d",delim," ",path),intern = T)
    if(substr(output[1],1,1)=="\""){
        output<-gsub("^\"","",output)#Remove bracketing quotes around quoted strings
        output<-gsub("\"$","",output)#Remove bracketing quotes around quoted strings
    }
    return(output)
}
