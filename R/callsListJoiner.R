callsListJoiner<-function(callsListOrigin,
                          callsListAdded,
                          callsListNames=c(origin="origin",added="added"),
                          fill=F,
                          boundElements=c("Alleles","MajorAlleleCount","MinorAlleleCount",
                                          "LeftCalls","RightCalls","LeftCounts","RightCounts",
                                          "CombCount","IsMultiAllele","Calls")){
    ####Check formatting of inputs
    #Check if names of list elements are the same excluding SampleData and callsInfo
    if(!all(names(callsListAdded)==names(callsListOrigin)[!names(callsListOrigin)%in%c("SampleData","callsInfo")])){
        stop("irregularity between names of callsLists")
    }
    #Check if the elements to be combined specifically are in both
    if(!all(c(boundElements%in%names(callsListAdded),
              boundElements%in%names(callsListOrigin)))){
        stop(paste0("callsLists need dataframe elements; '",paste0(boundElements,collapse = "', '"),"'"))
    }
    #Since callsListOrigin can have already been bound previously, convert it to a matrix before combineing
    callsListOrigin$IsMultiAllele<-as.matrix(callsListOrigin$IsMultiAllele)
    callsListAdded$IsMultiAllele <-as.matrix(callsListAdded$IsMultiAllele)
    #Check class of bound elements.
    if(!all(c(sapply(callsListAdded[boundElements],class)=="matrix",
              sapply(callsListAdded[boundElements],class)=="matrix"))){
        stop(paste0("callsLists need matrix objects in elements; '",paste0(boundElements,collapse = "', '"),"'"))
    }
    if(!all(unique(c("origin","added"))%in%names(callsListNames))){
        stop("Calls list names should be a named vector of unique dataset names")
    }

    #Set reference objects
    numGenoOrigin<-ncol(callsListOrigin$Calls)
    numGenoAdded<-ncol(callsListAdded$Calls)
    namesOrigin<-colnames(callsListOrigin$Calls)
    namesAdded<-colnames(callsListAdded$Calls)

    #Create a sample data list if needed and fill in dataset information
    if("SampleData"%in%names(callsListOrigin)){
        if(callsListNames["added"]=="added"){
            addName<-paste0("added",length(unique(callsListOrigin$SampleData$datasetNames)))
        }else{addName<-callsListNames["added"]}
        callsListOrigin$SampleData<-rbind(callsListOrigin$SampleData,
                                          data.frame(stringsAsFactors = F,
                                                     names=namesAdded,
                                                     datasetNames=rep(addName,numGenoAdded),
                                                     duplicated=NA,
                                                     joinedName=NA))
    }else{
        if(callsListNames["added"]=="added"){
            addName<-"added1"
        }else{addName<-callsListNames["added"]}
        callsListOrigin$SampleData<-data.frame(stringsAsFactors = F,
                                               names=c(namesOrigin,namesAdded),
                                               datasetNames=c(rep(callsListNames["origin"],numGenoOrigin),
                                                              rep(addName,numGenoAdded)),
                                               duplicated=NA,
                                               joinedName=NA)
    }

    #Check for and correct duplicate names
    callsListOrigin$SampleData$duplicated<-duplicated(callsListOrigin$SampleData$names)
    callsListOrigin$SampleData$joinedName<-callsListOrigin$SampleData$names
    namesToCorrect<-callsListOrigin$SampleData$joinedName[callsListOrigin$SampleData$duplicated]

    callsListOrigin$SampleData$joinedName[callsListOrigin$SampleData$names%in%namesToCorrect]<-
        paste0(callsListOrigin$SampleData$joinedName[callsListOrigin$SampleData$names%in%namesToCorrect],
               "...",
               callsListOrigin$SampleData$datasetNames[callsListOrigin$SampleData$names%in%namesToCorrect])

    #Stop if catting dataset didnt make unique
    if(any(duplicated(callsListOrigin$SampleData$joinedName))){
        duplicatedNamesLogic<-callsListOrigin$SampleData$joinedName%in%callsListOrigin$SampleData$joinedName[duplicated(callsListOrigin$SampleData$joinedName)]
        stop(paste0("Duplicated genotypes remianing after catting dataset name. Joining failed. Rename offending genotypes.:\n",
                    paste0(callsListOrigin$SampleData$joinedName[duplicatedNamesLogic],
                           collapse = "\n")))
    }

    #Warn if any duplicate names found and describe new name
    if(any(duplicated(callsListOrigin$SampleData$names))){
        duplicatedNamesLogic<-callsListOrigin$SampleData$names%in%callsListOrigin$SampleData$names[duplicated(callsListOrigin$SampleData$names)]
        warning(paste0("Duplicated genotypes saved and analysis continued:\n",
                       paste0(callsListOrigin$SampleData$joinedName[duplicatedNamesLogic],collapse = "\n")))
    }

    #Join the lists
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

    return(callsListOrigin)
}

joinDfChecker<-function(fileDf){
    if(!all(file.exists(fileDf))){stop("A subset file was not found")}
    markerNames<-apply(fileDf,c(1,2),cutFunction)
    if(!all(c(all(markerNames[,1,1]==markerNames[,2,1]),all(markerNames[,1,1]==markerNames[,3,1]),
              all(markerNames[,1,1]==markerNames[,1,2]),all(markerNames[,1,1]==markerNames[,2,2]),
              all(markerNames[,1,1]==markerNames[,3,2])))){stop("Subsets don't match")}
    print("Subset files found and marker names matched.")
}

joinCallsLists<-function(callsListList,datasetNames){
    for(i in 1:(length(datasetNames)-1)){
        if(i==1){
            callsList<-callsListList[[1]]
        }
        callsList<-callsListJoiner(callsListOrigin = callsList,
                                   callsListAdded  = callsListList[[i+1]],
                                   callsListNames  = c(origin=datasetNames[i],added=datasetNames[i+1]))
        print(paste0(i," of ",(length(datasetNames)-1)," loaded @ ",Sys.time()))
    }
}

checkJoinedListNames<-function(joinedNames,callsListList){
    return(all(joinedNames==unlist(lapply(callsListList,function(x){colnames(x$Calls)}))))
}
