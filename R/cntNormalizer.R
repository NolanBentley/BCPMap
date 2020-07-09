normalCountCalculation<-function(combCnt,
                                 originalTargetCalls,
                                 zRange=c(qnorm(0.99),20),
                                 zLength=10,
                                 plotToggle=T,
                                 plotName="./cntNorm.png",
                                 b=2){
    #Since this script is commonly nested, print its use
    print(paste0("Count normalizer used on ",ncol(combCnt)," samples."))

    #Formula of zScore cutoffs
    x<-1:zLength
    h<-min(x)
    r<-max(zRange)-min(zRange)
    k<-(sqrt(r^2*b)-max(zRange))*(-1)
    a<-(zLength-h)^2/(r^2-((min(zRange)-k)^2/b))
    y<-sqrt((r^2-(x-h)^2/a)*b)+k
    limitZScore<-y

    #Create a range of ZScore limits that lowers slowly then faster
    totalInput<-sum(originalTargetCalls)
    for(i in 1:length(limitZScore)){
        if(i==1){targetCalls<-originalTargetCalls}
        zMult<-limitZScore[i]
        currMeds<-apply(combCnt[targetCalls,],2,median,na.rm=T)
        currMads<-apply(combCnt[targetCalls,],2,mad,na.rm=T)
        medZCnt<-t((t(combCnt)-currMeds)/currMads)
        currQuants<-apply(medZCnt[targetCalls,],2,quantile,probs=c(pnorm(-zMult),pnorm(zMult)))
        currQuantsLogic<-medZCnt>=(-zMult)&medZCnt<=zMult
        targetCalls<-rowMeans(currQuantsLogic)==1&targetCalls
        numRemoved<-sum(targetCalls==F&originalTargetCalls)
        plotMain<-paste0("Iter ",i,", z-score limit= +-",round(zMult,6),": ",
                         numRemoved," of ",totalInput," removed (",round(numRemoved/sum(totalInput)*100,2),"%)")
        print(plotMain)
    }

    #Visualize results
    if(plotToggle){
        p1<-ggplot(melt(medZCnt),aes(x=value,fill=X2))+geom_histogram(binwidth = 0.25)+
            geom_vline(data = melt(currQuants),aes(xintercept=value,color=X2))+
            coord_cartesian(xlim=c(max(min(-limitZScore),min(medZCnt)),max(limitZScore)))+
            labs(title=paste0(plotMain,"\n","lines at ",
                              paste0(round(c(pnorm(-zMult),pnorm(zMult)),4),collapse = " & "),
                              " quantiles. Final cutoff at leftmost rugplot."))+
            geom_rug(data = data.frame(x=limitZScore),
                     mapping = aes(x),color=viridis(zLength),inherit.aes = F)+
            theme_bw()
        ggsave(plotName,p1)
    }
    return(targetCalls)
}
