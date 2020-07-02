#Prints each element of input list as seperate files to a common directory for easier reconstruction later
require(rlist)
saveLargeList<-function(listToPrint,
                        baseDir,
                        baseName,
                        elementDelim="__",
                        fileType="rds",
                        elementsInd=NULL,
                        minPad=2){
    if(any(names(listToPrint)==""|duplicated(names(listToPrint)))){
        stop("List names need to be unique strings with nchar>=1")
    }
    if(is.null(elementsInd)){elementsInd<-1:length(listToPrint)}
    dir.create(file.path(baseDir,baseName),showWarnings = F,recursive = T)
    for(i in elementsInd){
        currFilename<-file.path(baseDir,baseName,
                                paste0(baseName,".",
                                       gsub(" ","0",format(i,width = max(c(minPad,nchar(max(elementsInd)))))),
                                       elementDelim,names(listToPrint)[i],".",fileType))
        write(1,currFilename)
        if(file.exists(currFilename)){
            list.save(listToPrint[[i]],file = currFilename)
        }else{
            stop(paste0("Filename failed to create: ",currFilename))
        }
        print(paste0(which(i==elementsInd)," of ",length(elementsInd)," saved @ ",Sys.time(),": ",currFilename))
    }
}
loadLargeList<-function(baseDir,
                        baseName,
                        elementDelim="__",
                        fileType="yaml"){
    fileVector<-list.files(file.path(baseDir,baseName),pattern = paste0(".*",elementDelim,".*",fileType,"$"),full.names = T)
    for(i in 1:length(fileVector)){
        if(i==1){
            outputList<-list()
        }
        outputList[[i]]<-list.load(fileVector[i])
        names(outputList)[i]<-gsub(paste0("\\.",fileType,"$"),"",gsub(paste0(".*",elementDelim),"",fileVector[i]))
        print(paste0(i," of ",length(fileVector)," loaded @ ",Sys.time(),": ",fileVector[i]))
    }
    return(outputList)
}
