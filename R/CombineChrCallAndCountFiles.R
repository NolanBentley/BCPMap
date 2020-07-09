joinCallAndCountFiles<-function(callString,countString,numPartsExpected,
                                combinedCallFile ="../GitExamples/CallsCombined.csv",
                                combinedCountFile="../GitExamples/CountsCombined.csv"){
    #Check inputs
    if(!dir.exists(dirname(callString))){stop(dirname(callString)," not found")}
    if(!dir.exists(dirname(countString))){stop(dirname(countString)," not found")}
    if(!dir.exists(dirname(combinedCallFile))){stop(dirname(combinedCallFile)," not found")}
    if(!dir.exists(dirname(combinedCountFile))){stop(dirname(combinedCountFile)," not found")}
    if(!numPartsExpected>0){stop("numPartsExpected needs to be numeric and > 0")}

    #Grep the basename(string) in the path resulting from dirname(string)
    callFileList <-list.files(path = dirname(callString ),pattern = basename(callString ),full.names = T)
    countFileList<-list.files(path = dirname(countString),pattern = basename(countString),full.names = T)
    if(!(length(callFileList)==numPartsExpected&length(countFileList)==numPartsExpected)){
        cat(paste0("\nFile list not expected length:\n",paste0(basename(c(callFileList,countFileList)),collapse = "\n")))
        stop("File list not expected length")
    }

    #Submit bash to combine files beginning with the first ind
    for(i in 1:numPartsExpected){
        file_calls <-callFileList [i]
        file_counts<-countFileList[i]
        if(i==1){
            system(paste0("cp ",file_calls," ",combinedCallFile))
            system(paste0("cp ",file_counts," ",combinedCountFile),wait = T)
        }else{
            wcCalls <-as.numeric(gsub(" .*","",system(paste0("wc -l ",file_calls ),intern = T)))
            wcCounts<-as.numeric(gsub(" .*","",system(paste0("wc -l ",file_counts),intern = T)))
            if(wcCalls!=wcCounts){stop("Call and count files not same length")}
            system(paste0("tail -n ",wcCalls -2," ",file_calls ," >> ",combinedCallFile))
            system(paste0("tail -n ",wcCounts-2," ",file_counts," >> ",combinedCountFile),wait = T)
        }
        print(paste0(basename(file_calls),": iter ",i," of ",numPartsExpected," catted."))
    }
    return(c(combinedCallFile,combinedCountFile))
}
