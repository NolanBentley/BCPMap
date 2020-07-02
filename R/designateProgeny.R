designateProgeny<-function(sampleData,callColnames,progenyMappedFile,useSubsetOfGenotypes){
    if(!useSubsetOfGenotypes){
        if(any(callColnames!=sampleData$joinedName)){
            stop("non matching names between calls colnames and sampledata$joinedName detected")
        }
        return(rep(T,nrow(sampleData)))
    }
    #Load progeny file
    if(!file.exists(progenyMappedFile)){stop(paste0("Progeny file ",progenyMappedFile," doesn't exist."))}
    samplesDesignated<-readLines(progenyMappedFile)
    if(samplesDesignated[1]=="names"){samplesDesignated<-samplesDesignated[2:length(samplesDesignated)]}
    #Check objects for expected matching
    if(any(callColnames!=sampleData$names)){
        stop("non matching names between calls colnames and sampledata$names detected")
    }
    if(any(!samplesDesignated%in%sampleData$names)){
        stop("non matching names between progeny calls and sampledata$joinedName detected")
    }
    if(any(duplicated(samplesDesignated))){
        stop("Duplicate names detected")
    }
    #Return logical vector of if name in progeny file
    mapped<-sampleData$joinedName%in%samplesDesignated
    return(mapped)
}
