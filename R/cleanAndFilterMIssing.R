cleanAndFilterMissing<-function(combinedCallFile,combinedCountFile,
                            maxMissProp = 0.2,headerRow = 2,missString="missing"){
    if(file.exists(paste0(gsub("\\.csv$","",combinedCountFile),"_stripped.csv"))){
        stop("Stripped count file alread exists")
    }
    #Read header row
    row2<-system(paste0("head -n ",headerRow," ",combinedCallFile," | tail -n 1"),intern = T)
    row2<-unlist(strsplit(row2,split = ","))
    numGenos<-length(grep("^G",row2))
    #Cut out missing call counts
    missingCount<-system(paste0("cut -d, -f ",which(row2==missString)," ",combinedCallFile),intern = T)
    missingCount[1:2]<-NA
    missingCount<-as.numeric(missingCount)
    hist(missingCount,100);abline(v = ceiling(maxMissProp*numGenos),col="red")
    #Open and subset calls file
    body<-system(paste0("cut -d, -f -",which(row2=="missing")-1," ",combinedCallFile),intern = T)
    body<-body[c(1,2,which(missingCount<=ceiling(maxMissProp*numGenos)))]
    write(paste0(gsub("\\.csv$","",combinedCallFile),"_stripped.csv"), x = body, sep = ",")
    #Open and subset counts file
    body<-readLines(combinedCountFile)
    body<-body[c(1,2,which(missingCount<=ceiling(maxMissProp*numGenos)))]
    write(paste0(gsub("\\.csv$","",combinedCountFile),"_stripped.csv"), x = body, sep = ",")
}
