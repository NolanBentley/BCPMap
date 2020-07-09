summarizeSelection<-function(groupings,positions,selected,excludeGroups="."){
    uniGroups<-sort(unique(groupings))
    for(i in 1:length(uniGroups)){
        if(i==1){
            selectionSummary<-data.frame(stringsAsFactors = F)
        }
        currRows<-selected&uniGroups[i]==groupings
        currDf<-data.frame(NumSelected=sum(currRows),
                           MeanDiff=mean(diff(positions[currRows])),
                           SdDiff=sd(diff(positions[currRows])),
                           MadDiff=mad(diff(positions[currRows]))
        )
        probs<-c(0,0.05,0.25,0.5,0.75,0.95,1)
        currQuants<-matrix(quantile(diff(positions[currRows]),probs = probs),1,
                           dimnames = list(uniGroups[i],paste0(probs*100,"%")))
        currDf<-cbind(currDf,currQuants)
        selectionSummary<-rbind(selectionSummary,currDf)
    }
    numCols<-unlist(lapply(selectionSummary[1,],is.numeric))
    totalLine<-as.data.frame(matrix(NA,1,ncol(selectionSummary),dimnames = list("Average",colnames(selectionSummary))))
    if(sum(rownames(selectionSummary)%in%excludeGroups)>0){
        warning(paste0("Excluding ",paste0(excludeGroups,collapse = " & ")," from average row."))
    }
    totalLine[numCols]<-colMeans(selectionSummary[!rownames(selectionSummary)%in%excludeGroups,numCols])
    selectionSummary<-rbind(selectionSummary,totalLine)
    selectionSummary[,numCols]<-apply(selectionSummary[,numCols],2,format,big.mark=",")
    print(paste0(sum(selected)," markers selected:"))
    return(selectionSummary)
}
