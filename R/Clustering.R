phaseCodedCalls<-function(callsList,
                           parent,
                           rowLogic,
                           minAbsPC1=1,
                           chrGrepString="Chr",
                          storageDir="./",
                          imageWidthIn=8,
                          imageHeightIn=8,
                          imageDPI=600){
    #Identify clustered markers and columns
    currRows<-rowLogic&grepl(chrGrepString,callsList$callsInfo$Chr)
    currCols<-callsList$SampleData$mapped
    callsList$callsInfo$PC1<-NA

    #Use PCA to cluster markers
    callsList$callsInfo$PC1[currRows]<-pcaPhasing(callsList$CodedCalls[currRows,currCols],
                                                  callsList$callsInfo$Chr[currRows])
    #Analyze PC1 to call phase based on the abs minimum PC1 value
    callsList$callsInfo$Phase<-(callsList$callsInfo$PC1<0)*1+1
    callsList$callsInfo$Phase[which(abs(callsList$callsInfo$PC1)<minAbsPC1)]<-0
    #Visualize and save PCA
    pPCA<-pcaVis(callsList$callsInfo,currRows,minPC1=minAbsPC1)
    ggsave(filename = file.path(storageDir,paste0("PCA_",parent,".png")),
           plot = pPCA,
           width = imageWidthIn,
           height = imageHeightIn,
           dpi = imageDPI,
           units = "in")

    return(callsList$callsInfo)
}

pcaPhasing<-function(codedCalls,groupings,imputeHomoMin=T,numDimToReturn=1){
    if(imputeHomoMin){codedCalls[codedCalls==2]<-1}
    codedCalls[is.na(codedCalls)]<-0.5
    outputPCA<-NULL
    uniGroup<-sort(unique(groupings))
    for(i in 1:length(uniGroup)){
        currCodes<-codedCalls[groupings==uniGroup[i],]
        currPCA<-prcomp(rbind(currCodes,1-currCodes),scale. = T)
        currPCA<-currPCA$x[1:nrow(currCodes),]
        outputPCA<-rbind(outputPCA,currPCA[,1:(numDimToReturn+1)])
    }
    return(outputPCA[,1:numDimToReturn])
}

pcaVis<-function(callsInfo,currRows,xDivisor=1000000,xUnit="Mbp",yUnit="PC1",colUnit="Phase",minPC1=1,
                 xycgVector=c("Pos","PC1","Phase","Chr")){
    df1<-data.frame(x        = callsInfo[,xycgVector[1]],
                    y        = callsInfo[,xycgVector[2]],
                    color    = callsInfo[,xycgVector[3]],
                    grouping = callsInfo[,xycgVector[4]])
    df1<-df1[currRows,]
    df1$x<-df1$x/xDivisor
    p1<-ggplot(df1,aes(x,y,color=as.factor(color)))+
        geom_hline(yintercept = c(-minPC1,minPC1),color="yellow")+
        geom_point()+
        scale_color_viridis_d(direction = -1)+
        facet_wrap(vars(grouping))+
        theme_bw()+
        labs(x=xUnit,y=yUnit,color=colUnit)
    return(p1)
}

