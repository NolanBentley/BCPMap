callsListJoiner<-function(callsListAdded,
                          callsListOrigin,
                          callsListNames=c("added","origin"),
                          fill=F,
                          boundElements=c("Alleles","MajorAlleleCount","MinorAlleleCount","IsMultiAllele","Calls")
){
    if(!all(names(callsListAdded)==names(callsListOrigin)[!names(callsListOrigin)%in%c("SampleData","callsInfo")])){
        stop("irggularity between names of callsLists")
    }
    sapply(inputCheckVector,exists)
    numGenoOrigin<-ncol(callsListOrigin$Calls)
    numGenoAdded<-ncol(callsListAdded$Calls)
    namesOrigin<-colnames(callsListOrigin$Calls)
    namesAdded<-colnames(callsListAdded$Calls)
    #Create a sample data list if needed and fill in dataset information
    if(exists("callsListOrigin[[\"SampleData\"]]")){
        if(callsListNames[1]=="added"){
            addName<-paste0("added",length(unique(callsListOrigin$SampleData$datasetNames)))
        }else{addName<-callsListNames[1]}
        callsListOrigin$SampleData<-rbind(callsListOrigin$SampleData,
                                          data.frame(stringsAsFactors = F,
                                                     names=namesAdded,
                                                     datasetNames=rep(addName,numGenoAdded)))
    }else{
        if(callsListNames[1]=="added"){
            addName<-"added1"
        }else{addName<-callsListNames[1]}
        callsListOrigin$SampleData<-data.frame(stringsAsFactors = F,
                                               names=c(namesOrigin,namesAdded),
                                               datasetNames=c(rep(callsListNames[2],numGenoOrigin),
                                                              rep(addName,numGenoAdded)))
    }
    if(any(duplicated(callsListOrigin$SampleData$names))){
        print("Duplicate names need attention")
        warning(paste0("Duplicated genotypes saved\n:",
                       paste0(callsListOrigin$SampleData$names[duplicated(callsListOrigin$SampleData$names)],
                              collapse = "\n")))
    }

    if(fill==T){
        #Adds missing data to calls object where snps are absent currently
        print("Fill not yet implemented")
    }else{
        #Make sure all rownames match before combining
        rowTest<-mapply(function(x,y){if(length(x)!=length(y)){F}else{all(x==y)}},
                        x=lapply(callsListOrigin[boundElements],rownames),
                        y=lapply(callsListAdded [boundElements],rownames))
        if(any(!rowTest)){stop(paste0(paste0(boundElements[!rowTest],collapse = " & ")," don't match."))}
        callsListOrigin[boundElements]<-mapply(function(x,y){cbind(x,y)},
                                               x=callsListOrigin[boundElements],
                                               y=callsListAdded[boundElements])
    }
}
