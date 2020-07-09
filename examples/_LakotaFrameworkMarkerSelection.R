maternalFrameworkMarkers<-function(callsList,
                                   histogramFilename="./maternalFrameworkMarkerHistogram.png",
                                   scatterplotFilename="./maternalFrameworkMarkerScatterplot.png",
                                   minHet=0.05,
                                   maxHet=0.75,
                                   minPercDiff.HomoCallFreqs=0.4,
                                   b = 3,
                                   zRange = c(qnorm(0.999),20),
                                   zLength = 20,
                                   closestProximity=150,
                                   seqGrepStr="Chr"
                                   ){
    ####Select markers to be mapped
    parentalInfo<-list(name             = c("Lakota-CSHQ-12-6-s2","LakotaaltHap-s1","LakotamainGenome-s1",
                                            "87MX3-2-11-CSPB-1-30-s2+s3","87MX3altHap-s1","87MX3mainGenome-s1"),
                       type             = c("GBS","Sim","Sim" ,"GBS","Sim","Sim" ),
                       parent           = c("Mat","Mat","Mat" ,"Pat","Pat","Pat" ),
                       reference        = c("GBS" ,"Alt","Main","GBS" ,"Alt","Main"),
                       limit.Called     = c(T    ,T    ,T     ,T    ,F    ,T     ),
                       limit.Het        = c(T    ,F    ,F     ,T    ,F    ,F     ),
                       limit.SimCombHet = c(T    ,T    ,T     ,T    ,F    ,F     ),
                       limit.RDMinHomo  = c(6    ,11   ,11    ,6    ,11   ,11    ),
                       limit.RDMaxHomo  = c(Inf  ,12   ,12    ,Inf  ,12   ,12    ),
                       limit.RDMinHet   = c(2    ,0    ,0     ,2    ,0    ,0     ),
                       limit.RDMaxHet   = c(Inf  ,0    ,0     ,Inf  ,0    ,0     ))
    parentalInfo$SamCol<-match(parentalInfo$name,callsList$SampleData$names)
    parentalInfo<-lapply(parentalInfo,function(x,currNames){names(x)<-currNames;return(x)},
                         currNames=paste0(parentalInfo$parent,parentalInfo$reference))
    callsList$parentalInfo<-parentalInfo

    ####Find putative Lakota sites
    #Check Lakota against SIM and GBS
    callsList$callsInfo$IsMatHet<-
        (callsList$CodedCalls[,parentalInfo$SamCol["MatMain"]]==0&
             callsList$CodedCalls[,parentalInfo$SamCol["MatAlt" ]]==2&
             callsList$CodedCalls[,parentalInfo$SamCol["MatGBS" ]]==1)|
        (callsList$CodedCalls[,parentalInfo$SamCol["MatMain"]]==2&
             callsList$CodedCalls[,parentalInfo$SamCol["MatAlt" ]]==0&
             callsList$CodedCalls[,parentalInfo$SamCol["MatGBS" ]]==1)
    callsList$callsInfo$IsMatHet[is.na(callsList$callsInfo$IsMatHet)]<-F
    #Check Oaxaca to confirm not het but allow for PatAlt to be NA
    callsList$callsInfo$IsOaxNotHet<-
        (callsList$CodedCalls[,parentalInfo$SamCol["PatMain"]]==2&
             callsList$CodedCalls[,parentalInfo$SamCol["PatAlt" ]]==2&
             callsList$CodedCalls[,parentalInfo$SamCol["PatGBS" ]]==2)|
        (callsList$CodedCalls[,parentalInfo$SamCol["PatMain"]]==0&
             callsList$CodedCalls[,parentalInfo$SamCol["PatAlt" ]]==0&
             callsList$CodedCalls[,parentalInfo$SamCol["PatGBS" ]]==0)|
        (callsList$CodedCalls[,parentalInfo$SamCol["PatMain"]]==0&
             is.na(callsList$CodedCalls[,parentalInfo$SamCol["PatAlt" ]])&
             callsList$CodedCalls[,parentalInfo$SamCol["PatGBS" ]]==0)
    callsList$callsInfo$IsOaxNotHet[is.na(callsList$callsInfo$IsOaxNotHet)]<-F
    #Are Lakota sims read depths as expected
    callsList$callsInfo$MatSimRDLogic<-
        callsList$MajorAlleleCount[,parentalInfo$SamCol["MatMain"]]>=parentalInfo$limit.RDMinHomo["MatMain"]&
        callsList$MajorAlleleCount[,parentalInfo$SamCol["MatAlt" ]]<=parentalInfo$limit.RDMaxHomo["MatAlt" ]&
        callsList$MajorAlleleCount[,parentalInfo$SamCol["PatMain"]]<=parentalInfo$limit.RDMaxHomo["PatMain"]
    #Are Paternal and Maternal read depths within normal distributions
    dimnames(callsList$CombCount)<-dimnames(callsList$Calls)
    callsList$parCombCnt<-callsList$CombCount[,parentalInfo$SamCol[c("MatGBS","PatGBS")]]

    #Filter based on relative read depths
    callsList$callsInfo$PutativeMatTestcross<-
        callsList$callsInfo$IsMatHet&
        callsList$callsInfo$IsOaxNotHet&
        callsList$callsInfo$MatSimRDLogic
    callsList$callsInfo$PutativeMatTestcross<-
        normalCountCalculation(combCnt = callsList$parCombCnt,plotName = histogramFilename,b=b,zRange = zRange,zLength = zLength,
                               originalTargetCalls = callsList$callsInfo$PutativeMatTestcross)

    #Filter based on Het and HomoMin frequency in progeny
    callsList$callsInfo$MappedHet    <-rowMeans(callsList$CodedCalls[,callsList$SampleData$mapped]==1,na.rm=T)
    callsList$callsInfo$MappedHomoMaj<-rowMeans(callsList$CodedCalls[,callsList$SampleData$mapped]==0,na.rm=T)
    callsList$callsInfo$MappedHomoMin<-rowMeans(callsList$CodedCalls[,callsList$SampleData$mapped]==2,na.rm=T)

    callsList$callsInfo$AbsHomoPercDiff<-
        abs((callsList$callsInfo$MappedHomoMaj-callsList$callsInfo$MappedHomoMin))/
        (callsList$callsInfo$MappedHomoMaj+callsList$callsInfo$MappedHomoMin)

    callsList$callsInfo$MinAbsHomoPercDiff<-callsList$callsInfo$AbsHomoPercDiff>=minPercDiff.HomoCallFreqs
    callsList$callsInfo$MinHet<-callsList$callsInfo$MappedHet>=minHet
    callsList$callsInfo$MaxHet<-callsList$callsInfo$MappedHet<=maxHet

    callsList$callsInfo$CallStatisticFiltered<-
        callsList$callsInfo$MinAbsHomoPercDiff&
        callsList$callsInfo$MinHet&
        callsList$callsInfo$MaxHet

    print(paste0("An additional ",
                 sum(!callsList$callsInfo$CallStatisticFiltered[callsList$callsInfo$PutativeMatTestcross]),
                 " markers filtered due to extreme marker statistics.")
    )

    plottedDf<-data.frame(x=callsList$callsInfo$AbsHomoPercDiff,
                          y=callsList$callsInfo$MappedHet,
                          chr=callsList$callsInfo$Chr,
                          selected=callsList$callsInfo$PutativeMatTestcross,
                          filtered=callsList$callsInfo$CallStatisticFiltered
                          )
    plottedDf<-plottedDf[!is.na(plottedDf$x),]
    plottedDf<-plottedDf[order(plottedDf$selected),]

    scatterP1<-ggplot(plottedDf,aes(x,y,shape=filtered))+
        geom_rect(data = data.frame(xmin=(-1),
                                    xmax=minPercDiff.HomoCallFreqs,
                                    ymin=(-1),
                                    ymax=2),
                  aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),inherit.aes = F,alpha=0.2)+
        geom_rect(data = data.frame(xmin=minPercDiff.HomoCallFreqs,
                                    xmax=2,
                                    ymin=(-1),
                                    ymax=minHet),
                  aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),inherit.aes = F,alpha=0.2)+
        geom_rect(data = data.frame(xmin=minPercDiff.HomoCallFreqs,
                                    xmax=2,
                                    ymin=maxHet,
                                    ymax=2),
                  aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),inherit.aes = F,alpha=0.2)+
        geom_bin2d(binwidth=0.01)+
        scale_fill_viridis_c(trans="log10")+
        geom_vline(xintercept = minPercDiff.HomoCallFreqs,color="red")+
        geom_hline(yintercept = c(minHet,maxHet),color="red")+
        geom_point(data = plottedDf[plottedDf$selected,],alpha=0.85,color="red",fill="transparent")+
        scale_shape(solid = F)+
        labs(x="abs(HomoMajFreq-HomoMinFreq)/(HomoMajFreq+HomoMinFreq)",
            y="HetFreq",fill="Count",color="Selected")+
        theme_bw()+
        coord_cartesian(xlim=c(0,1),ylim=c(0,1))
    ggsave(scatterP1,width = 8,height = 8,dpi = 600,units = "in",filename = scatterplotFilename)


    ####Cluster based on physical proximity
    #Sort so that smallest HomoMinorFreqs are first considered
    callsList$callsInfo$MappedChr<-grepl(seqGrepStr,callsList$callsInfo$Chr)
    currSelected<-which(
        callsList$callsInfo$PutativeMatTestcross&
            callsList$callsInfo$CallStatisticFiltered&
            callsList$callsInfo$MappedChr
    )
    names(currSelected)<-callsList$callsInfo$Markers[currSelected]
    currSelected<-currSelected[order(callsList$callsInfo$MappedHomoMin[currSelected])]

    #Bin by proximity and subset selection
    callsList$callsInfo$PosBin<-paste0(callsList$callsInfo$Chr,"|",floor((callsList$callsInfo$Pos[currSelected])/closestProximity))
    #Iterate by offsets to space bins
    posBinX<-callsList$callsInfo$PosBin
    for(distMult in seq(0,closestProximity,by = 1)){
        posBinX[currSelected]<-paste0(callsList$callsInfo$Chr[currSelected],"|",floor((callsList$callsInfo$Pos[currSelected]+distMult)/closestProximity))
        currSelected<-currSelected[!duplicated(posBinX[currSelected])]
        posBinX[currSelected]<-paste0(callsList$callsInfo$Chr[currSelected],"|",floor((callsList$callsInfo$Pos[currSelected]-distMult)/closestProximity))
        currSelected<-currSelected[!duplicated(posBinX[currSelected])]
    }
    callsList$callsInfo$SelectedInBin<-F
    callsList$callsInfo$SelectedInBin[currSelected]<-T
    print(paste0("An additional ",
                 sum(!callsList$callsInfo$SelectedInBin[which(callsList$callsInfo$PutativeMatTestcross&callsList$callsInfo$CallStatisticFiltered)]),
                 " markers filtered to space them at least ",
                 closestProximity,
                 " bp apart."))

    return(list(selected=callsList$callsInfo$SelectedInBin,
                info=callsList$callsInfo))
}

