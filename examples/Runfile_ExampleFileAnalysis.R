#Initial processing of files to create default files
#Create combined calls file of diversity, simulated, and LxO progeny
#Begins with calls and counts files from the Klein-lab
# Begin merging of files based on reference sequence position based on SNPs in first column of filterCutFile
# filterCutFile is a csv file with the unique SNP names in column one and row initParam$markerRow and lasting until the end of the file



####Load packages
if(T){
    require(devtools)
    require(viridisLite)
    require(parallel)
    require(snow)
    require(dplyr)
    require(rlist)
    require(ggplot2)
    require(reshape)
    load_all()
}


####Set parameter values in list initParam
if(T){
    initParam<-list()
    initParam$markerRow<-3
    initParam$storageDir<-"../GitExamples"
    initParam$fileDf<-data.frame(stringsAsFactors = F,
        calls=c("../GitExamples/LxOMapCombinedDivSubset_calls.csv",#This should be the most limiting file to save memory
                "../GitExamples/LxOSimulatedShotgun_4G_130bp12x25min_Vs_OaxMainV1_3-3-15_20190716_182609_sam_calls.csv",
                "../GitExamples/PecanCRRDiv_844G_PstI_Vs_Carya_illinoinensis_var_87MX3.mainGenome_3-3-15_20190513_144932_sam_calls.csv"))
    initParam$fileDf$counts<-gsub("calls\\.csv","counts.csv",initParam$fileDf$calls)
    initParam$datasetNames<-c("LxO","Sim","Div")
    initParam$progenyMappedFile<-"../GitExamples/SampleNamesMapped.txt"
    initParam$paternalAnalysis<-"./examples/_LakotaFrameworkMarkerSelection.R"
    initParam$chrGrepString<-"Chr"
    initParam$minAbsPC1Framework<-2
    initParam$P1GP1<-"Mahan-CSV-6-6-s2"
    initParam$P1GP2<-"CRR-Mat-Major-CSV-19-5-s5+s6+s7+s8"
}


####Sourced example customized script files
#Add source of maternalFramworkMarkers function that creats a vector of markers to use in maternal framework map.
if(T){
    source(initParam$paternalAnalysis)
}


####Tracking progress and turning off completed sections
if(TRUE){
    progLst<-list()
    progLst$subsetCompleted<-T
    progLst$useSubsetFiles<-T
    progLst$joinCallsListsCompleted<-T
    progLst$loadCallsListFromJoinedFile<-T
    progLst$useSubsetOfGenotypes<-T
    progLst$callsCoded<-F
    progLst$markerStatisticsCreated<-F
    progLst$loadStatisticsFromFile<-F
    progLst$PutativeMaternalDetected<-F
    progLst$PutativeMaternalMarkersClustered<-F
    progLst$MaternalFrameworkCentimorgansCalculated<-F
}


####Collect list of obList of objects to not remove during bug testing, clean using: eval(parse(text=cleanList$cleanEnv))
if(T){
    cleanList<-list()
    cleanList$futureObjectsToNotDelete<-c("cleanList","callsList")
    cleanList$baseObjectList<-sort(unique(c(ls(),cleanList$futureObjectsToNotDelete)))
    cleanList$cleanEnv<-"rm(list = newObj(environment(),cleanList$baseObjectList));load_all()"#eval(parse(text=cleanList$cleanEnv)
}



####Subset input files to common loci
if(!progLst$subsetCompleted){
    #Check for file paths
    if(!all(file.exists(initParam$fileDf$calls))){stop("Call files not found")}
    if(!all(file.exists(initParam$fileDf$counts))){stop("Count files corresponding to the names of the call files not found")}

    #Produces subset versions of the input call and count files with a combined marker set
    SubsetCallsAndCountsFiles(initParam$fileDf)
}



####Change input to subset files
if(progLst$useSubsetFiles){
    initParam$fileDf<-apply(initParam$fileDf,c(1,2),function(x){paste0(gsub("\\.csv","",x),"_subset.csv")})
}



####Check if calls and counts are joinable then join
if(!progLst$joinCallsListsCompleted){
    #Check markers are in common
    joinDfChecker(initParam$fileDf)

    #Import Call and Count files
    callsListList<-multiCallListImporter(initParam$fileDf,initParam$datasetNames)

    #Join the call and count files
    callsList<-joinCallsLists(callsListList,initParam$datasetNames)

    #Check if correctly joined and remove the remaining objects if good
    namesMatch<-colnames(callsList$Calls)==unlist(lapply(callsListList,function(x){colnames(x$Calls)}))
    if(all(namesMatch)){rm("callsListList");rm("namesMatch")
    }else{stop("Joining incomplete as indicated by a check of call names")}

    #Save an image state
    save.image(file.path(initParam$storageDir,"image0010_unformattedJoined.R"))

    #Write joined callsList
    saveLargeList(listToPrint = callsList,baseDir = initParam$storageDir,baseName = "joinedCallsList",elementDelim = "__",fileType = "rds")
}



####Load callsList from file if needed
if(progLst$loadCallsListFromJoinedFile){
    if(exists("callsList")){
        warning("callList exists, skipping load from file")
    }else{
        callsList<-loadLargeList(baseDir = initParam$storageDir,baseName = "joinedCallsList",fileType = "rds")
    }
}



####Initiate callsInfo if needed
if(!exists("callsList$CallsInfo")){
    callsList$callsInfo<-makeCallsInfo(callsList$Calls)
    callsList$SampleData$mapped<-T
}



####Load list of progeny names to be used in mapping
if(progLst$useSubsetOfGenotypes){
    callsList$SampleData$mapped<-designateProgeny(sampleData = callsList$SampleData,
                                                  callColnames = colnames(callsList$Calls),
                                                  progenyMappedFile = initParam$progenyMappedFile,
                                                  useSubsetOfGenotypes = progLst$useSubsetOfGenotypes)
}



####Code mapped calls as 0, 1, 2, or NA if either homozygous major, heterozygous, homozygous minor, or missing
if(!progLst$callsCoded){
    callsList$CodedCalls<-codeCalls(Calls = callsList$Calls,
                                    majAlleles = callsList$Alleles[,1],
                                    minAlleles = callsList$Alleles[,2])
}



########Begin analyzing the maternal markers
#Custom selection of potentially maternal informative markers
#Script should produce a list with an element named "selected" that gives a logical vector indicated the selected markers
if(!progLst$PutativeMaternalDetected){
    maternalSelection<-maternalFrameworkMarkers(callsList,
                                                histogramFilename =   file.path(initParam$storageDir,"MaternalFrameworkHist.png"),
                                                scatterplotFilename = file.path(initParam$storageDir,"MaternalFrameworkScatter.png"),
                                                minHet = 0.05,closestProximity = 150,
                                                maxHet = 0.75,zRange = c(10,qnorm(0.999)),
                                                minPercDiff.HomoCallFreqs = 0.4)
    callsList$callsInfo$PutativeMatTestcross<-maternalSelection$selected
    addColumns<-colnames(maternalSelection$info)[!colnames(maternalSelection$info)%in%colnames(callsList$callsInfo)]
    callsList$callsInfo<-cbind(callsList$callsInfo,maternalSelection$info[,c("Markers",addColumns)])
    summarizeSelection(groupings = callsList$callsInfo$Chr,
                       positions = callsList$callsInfo$Pos,
                       selected  = callsList$callsInfo$PutativeMatTestcross,
                       excludeGroups = "scaffold")
}


#Cluster maternal framework markers
if(!progLst$PutativeMaternalMarkersClustered){
    callsList$callsInfo<-phaseCodedCalls(callsList,parent = "LakFramework",
                                         rowLogic  = callsList$callsInfo$PutativeMatTestcross,
                                         minAbsPC1 = initParam$minAbsPC1Framework,
                                         storageDir= initParam$storageDir)
    callsList$callsInfo$FrameworkPC1<-callsList$callsInfo$PC1
    callsList$callsInfo$FrameworkPhase<-callsList$callsInfo$Phase
}

#Filter and calculate cM positions for maternal framework markers
if(!progLst$MaternalFrameworkCentimorgansCalculated){
    #R/QTL Code framework markers
    callsList$phasedFramework<-adjustCodedCallsByPhase(codes=callsList$CodedCalls,
                                                       phases=callsList$callsInfo$FrameworkPhase)
    colnames(callsList$Alleles)<-apply(expand.grid(c("Major","Minor"),unique(callsList$SampleData$datasetNames)),1,paste0,collapse=".")
    callsList$Alleles<-as.data.frame(callsList$Alleles,stringsAsFactors = F)
    callsList$Alleles$Phase1<-NA
    callsList$Alleles$Phase1[which(callsList$callsInfo$FrameworkPhase==1)]<-callsList$callsInfo$MajorAllele[which(callsList$callsInfo$FrameworkPhase==1)]
    callsList$Alleles$Phase1[which(callsList$callsInfo$FrameworkPhase==2)]<-callsList$callsInfo$MinorAllele[which(callsList$callsInfo$FrameworkPhase==2)]
    callsList$Alleles$Phase2<-NA
    callsList$Alleles$Phase2[which(callsList$callsInfo$FrameworkPhase==1)]<-callsList$callsInfo$MinorAllele[which(callsList$callsInfo$FrameworkPhase==1)]
    callsList$Alleles$Phase2[which(callsList$callsInfo$FrameworkPhase==2)]<-callsList$callsInfo$MajorAllele[which(callsList$callsInfo$FrameworkPhase==2)]




    uniChr<-unique(callsList$callsInfo$Chr)
    uniPhases<-unique(callsList$callsInfo$FrameworkPhase)
    for(i in 1:length(uniChr)){
        callsList$SampleData$ChrPhase<-NA
        callsList$SampleData$ChrPhase<-colMeans(callsList$phasedFramework[which(callsList$callsInfo$Phase==1&callsList$callsInfo$Chr==uniChr[i]),],na.rm=T)
        colnames(callsList$SampleData)[colnames(callsList$SampleData)=="ChrPhase"]<-paste0(uniChr[i],"_1")
        callsList$SampleData$ChrPhase<-NA
        callsList$SampleData$ChrPhase<-colMeans(callsList$phasedFramework[which(callsList$callsInfo$Phase==2&callsList$callsInfo$Chr==uniChr[i]),],na.rm=T)
        colnames(callsList$SampleData)[colnames(callsList$SampleData)=="ChrPhase"]<-paste0(uniChr[i],"_2")
    }
    install.packages("heatmaply")
    install.packages("htmlwidgets")
    require("heatmaply")
    require("htmlwidgets")
    hm<-callsList$SampleData[,6:(ncol(callsList$SampleData)-2)]

    rownames(hm)<-paste0(callsList$SampleData$datasetNames,".",callsList$SampleData$names)
    plotOrder<-order(callsList$SampleData$datasetNames,rowMeans(abs(0.5-hm)))

    sizingPolicy(padding = 0, browser.fill = F)

    hm1<-heatmaply(hm[plotOrder,],
                  Colv = F,Rowv = F,
                  height = 12*nrow(callsList$SampleData))
    saveWidget(hm1,file = "Test1.html")

    plot(1:32,
        (as.numeric(hm[initParam$P1GP1==callsList$SampleData$names,])-
             as.numeric(hm[initParam$P1GP2==callsList$SampleData$names,]))/
            (as.numeric(hm[initParam$P1GP1==callsList$SampleData$names,])+
                 as.numeric(hm[initParam$P1GP2==callsList$SampleData$names,]))
        )
    dev.off()
}



