#Initial processing of files to create default files
#Create combined calls file of diversity, simulated, and LxO progeny

#Version 4.2
##Example text (https://samtools.github.io/hts-specs/VCFv4.2.pdf)
'
##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM	POS	    ID	        REF	ALT	    QUAL    FILTER	INFO	                            FORMAT	    NA00001	        NA00002	        NA00003
20	    14370	rs6054257	G	A	    29	    PASS	NS=3;DP=14;AF=0.5;DB;H2	            GT:GQ:DP:HQ	0|0:48:1:51,51	1|0:48:8:51,51	1/1:43:5:.,.
20	    17330	.	        T	A	    3	    q10	    NS=3;DP=11;AF=0.017	                GT:GQ:DP:HQ	0|0:49:3:58,50	0|1:3:5:65,3	0/0:41:3
20	    1110696	rs6040355	A	G,T	    67	    PASS	NS=2;DP=10;AF=0.333,0.667;AA=T;DB	GT:GQ:DP:HQ	1|2:21:6:23,27	2|1:2:0:18,2	2/2:35:4
20	    1230237	.	        T	.	    47	    PASS	NS=3;DP=13;AA=T	                    GT:GQ:DP:HQ	0|0:54:7:56,60	0|0:48:4:51,51	0/0:61:2
20	    1234567	microsat1	GTC	G,GTCT	50	    PASS	NS=3;DP=9;AA=G	                    GT:GQ:DP	0/1:35:4	    0/2:17:2    	1/1:40:3
'
# This example shows (in order):
## a good simple SNP,
## a possible SNP that has been filtered out because its quality is below 10,
## a site at which two alternate alleles are called, with one of them (T) being ancestral (possibly a reference sequencing error),
## a site that is called monomorphic reference (i.e. with no alternate alleles), and
## a microsatellite with two alternative alleles, one a deletion of 2 bases (TC), and the other an insertion of one base (T).
# Genotype data are given for three samples,
## two of which are phased and the third unphased,
## with per sample genotype quality,
## depth and haplotype qualities (the latter only for the phased samples) given as well as the genotypes.
## The microsatellite calls are unphased.

######### This is completely customized for processing my files
# Begin merging of files based on reference sequence position based on SNPs in first column of filterCutFile
# filterCutFile is a csv file with the unique SNP names in column one and row $markerRow and lasting until the end of the file

require(devtools)
require(parallel)
require(snow)
require(dplyr)
require(rlist)

load_all()
markerRow<-3
storageDir<-"../GitExamples"
fileDf<-data.frame(calls=c("../GitExamples/LxOMapCombinedDivSubset_calls.csv",#This should be the most limiting file to save memory
                           "../GitExamples/LxOSimulatedShotgun_4G_130bp12x25min_Vs_OaxMainV1_3-3-15_20190716_182609_sam_calls.csv",
                           "../GitExamples/PecanCRRDiv_844G_PstI_Vs_Carya_illinoinensis_var_87MX3.mainGenome_3-3-15_20190513_144932_sam_calls.csv"),
                   stringsAsFactors = F)
fileDf$counts<-gsub("calls\\.csv","counts.csv",fileDf$calls)
datasetNames<-c("LxO","Sim","Div")
progenyMappedFile<-"../GitExamples/SampleNamesMapped.txt"


####Tracking progress and turning off completed sections
subsetCompleted<-T
useSubsetFiles<-T
joinCallsListsCompleted<-T
loadCallsListFromJoinedFile<-T
useSubsetOfGenotypes<-T
callsCoded<-F
markerStatisticsCreated<-F
loadStatisticsFromFile<-F

####Collect list of obList of objects to not remove during bug testing, clean using: eval(parse(text=cleanEnv))
futureObjectsToNotDelete<-c("baseObjectList","baseEnv","cleanEnv","callsList")
baseObjectList<-sort(unique(c(ls(),futureObjectsToNotDelete)))
cleanEnv<-"rm(list = newObj(environment(),baseObjectList));load_all()"

####Subset input files to common loci
if(!subsetCompleted){
    #Check for file paths
    if(!all(file.exists(fileDf$calls))){stop("Call files not found")}
    if(!all(file.exists(fileDf$counts))){stop("Count files corresponding to the names of the call files not found")}

    #Produces subset versions of the input call and count files with a combined marker set
    SubsetCallsAndCountsFiles(fileDf)
}

####Change input to subset files
if(useSubsetFiles){
    fileDf<-apply(fileDf,c(1,2),function(x){paste0(gsub("\\.csv","",x),"_subset.csv")})
}

####Check if calls and counts are joinable then join
if(!joinCallsListsCompleted){
    #Check markers are in common
    joinDfChecker(fileDf)

    #Import Call and Count files
    callsListList<-multiCallListImporter(fileDf,datasetNames)

    #Join the call and count files
    callsList<-joinCallsLists(callsListList,datasetNames)

    #Check if correctly joined and remove the remaining objects if good
    namesMatch<-colnames(callsList$Calls)==unlist(lapply(callsListList,function(x){colnames(x$Calls)}))
    if(all(namesMatch)){rm("callsListList");rm("namesMatch")
    }else{stop("Joining incomplete as indicated by a check of call names")}

    #Save an image state
    save.image(file.path(storageDir,"image0010_unformattedJoined.R"))

    #Write joined callsList
    saveLargeList(listToPrint = callsList,baseDir = storageDir,baseName = "joinedCallsList",elementDelim = "__",fileType = "rds")
}

####Load callsList from file if needed
if(!exists("callsList")){
    if(loadCallsListFromJoinedFile){
        callsList<-loadLargeList(baseDir = storageDir,baseName = "joinedCallsList",fileType = "rds")
    }else if(!exists("callsList")){
        stop("No callsList object")
    }
}else if(loadCallsListFromJoinedFile){
    warning("callList exists, skipping load from file")
}

####Load list of progeny names to be used in mapping
callsList$SampleData$mapped<-designateProgeny(sampleData = callsList$SampleData,
                                              callColnames = colnames(callsList$Calls),
                                              progenyMappedFile = progenyMappedFile,
                                              useSubsetOfGenotypes = useSubsetOfGenotypes)

####Code mapped calls as 0, 1, 2, or NA if either homozygous major, heterozygous, homozygous minor, or missing
if(!callsCoded){
    callsList$CodedCalls<-codeCalls(Calls = callsList$Calls,
                                    majAlleles = callsList$Alleles[,1],
                                    minAlleles = callsList$Alleles[,2])
}

####Select markers to be mapped
parentalInfo<-data.frame(name             = c("Lakota-CSHQ-12-6-s2","LakotaaltHap-s1","LakotamainGenome-s1",
                                              "87MX3-2-11-CSPB-1-30-s2+s3","87MX3altHap-s1","87MX3mainGenome-s1"),
                         type             = c("GBS","Sim","Sim" ,"GBS","Sim","Sim" ),
                         parent           = c("Pat","Pat","Pat" ,"Mat","Mat","Mat" ),
                         reference        = c("NA" ,"Alt","Main","NA" ,"Alt","Main"),
                         limit.Called     = c(T    ,T    ,T     ,T    ,F    ,T     ),
                         limit.Het        = c(T    ,F    ,F     ,T    ,F    ,F     ),
                         limit.SimCombHet = c(T    ,T    ,T     ,T    ,F    ,F     ),
                         limit.RDMinHomo  = c(6    ,11   ,11    ,6    ,11   ,11    ),
                         limit.RDMaxHomo  = c(Inf  ,12   ,12    ,Inf  ,12   ,12    ),
                         limit.RDMinHet   = c(2    ,0    ,0     ,2    ,0    ,0     ),
                         limit.RDMaxHet   = c(Inf  ,0    ,0     ,Inf  ,0    ,0     ),
                         stringsAsFactors = F)
parentalInfo$SamCol<-match(parentalInfo$name,callsList$SampleData$names)

makeCallsInfo<-function(namedDf,gsubToIsolateSeq="_.*",gsubToIsolatePos=".*_"){
    callsInfo<-data.frame(Marker=rownames(namedDf),
                          Chr=gsub(gsubToIsolateSeq,"",rownames(namedDf)),
                          Pos=as.numeric(gsub(gsubToIsolatePos,"",rownames(namedDf))),
                          stringsAsFactors = F)
    callsInfo$RelaPos<-NA
    for(i in unique(callsInfo$Chr)){
        currRows<-callsInfo$Chr==i
        currPos<-callsInfo$Pos[currRows]
        callsInfo$RelaPos[currRows]<-(currPos-min(currPos))/(max(currPos)-min(currPos))
    }
    return(callsInfo)
}

plotHetSumWin<-function(isHet,callsInfo,win=500000){
    currHets<-cbind(Het=isHet,callsInfo)
    uniChr<-unique(currHets$Chr)
    currHets$WinHet<-NA
    currHets$WinHetNum<-NA
    for(j in 1:length(uniChr)){
        currcurrHets<-currHets[currHets$Chr==uniChr[j],]
        for(i in 1:nrow(currcurrHets)){
            lowerLimit<-currcurrHets$Pos[i]-win
            upperLimit<-currcurrHets$Pos[i]+win
            currWinHet<-currcurrHets$Het[currcurrHets$Pos>=lowerLimit&currcurrHets$Pos<=upperLimit]
            currcurrHets$WinHet[i]<-sum(currWinHet)
            currcurrHets$WinHetNum[i]<-length(currWinHet)
        }
        currHets[currHets$Chr==uniChr[j],c("WinHet","WinHetNum")]<-currcurrHets[,c("WinHet","WinHetNum")]
        print(j)
    }
    plot(currHets$RelaPos+as.numeric(gsub("Chr","",currHets$Chr)),currHets$WinHet)
}

stripSlashesFromCalls<-function(calls){
    alleles<-paste0(c("A","C","G","T","-"),"/")
    for(i in alleles){
        calls[calls==i]<-gsub("/","",i)
    }
}



checkSim<-function(callsList,parental,useCodes=T,homoMinRequired=F){

}

checkSim(hetMain = callsList$


applyParentalLimits<-function(callsList,limitControl){
    ####Applies a variety of filters to markers based on parental sequencing data and saves analysis to calls Info
    #Initiate callsInfo if needed
    if(!exists("callsList$CallsInfo")){callsList$callsInfo<-makeCallsInfo(callsList$Calls)}

    #Markers that are potentially maternally segregating
    currInfo<-parentalInfo[parentalInfo$parent=="Mat"&parentalInfo$limit.Het&parentalInfo$limit.Called,]
    callsList$callsInfo$MatRD<-
        matrix(callsList$CombCount [,currInfo$SamCol],ncol=nrow(currInfo))>=matrix(currInfo$limit.RDMinHet,nrow = nrow(callsList$CombCount),ncol = nrow(currInfo))&
        matrix(callsList$CombCount [,currInfo$SamCol],ncol=nrow(currInfo))<=matrix(currInfo$limit.RDMaxHet,nrow = nrow(callsList$CombCount),ncol = nrow(currInfo))&
        !is.na(matrix(callsList$CodedCalls[,currInfo$SamCol],ncol=nrow(currInfo)))
    currHets<-matrix(callsList$CodedCalls[,currInfo$SamCol],ncol=nrow(currInfo))==1
    currHets[,!currInfo$limit.Het]<-NA
    callsList$callsInfo$MatHet<-rowMeans(currHets,na.rm=T)==1



    currHets<-cbind(Het=currLogic[currEligble,],callsList$callsInfo[currEligble,])
    plotHetSumWin(currLogic[currEligble,],callsList$callsInfo[currEligble,])
}





####Create marker statistics
if(!markerStatisticsCreated){

    #Establish marker info sheet
    callsList$CallsInfo<-data.frame(stringsAsFactors = F,
                                    Markers=rownames(callsList$Calls))
    callsList$CallsInfo$Chr<-gsub("_.*","",callsList$CallsInfo$Markers)
    callsList$CallsInfo$Pos<-as.numeric(gsub(".*_","",callsList$CallsInfo$Markers))
    callsList$CallsInfo$cntInner99Logic<-callsList$cntInner99Logic
    callsList$CallsInfo<-cbind(callsList$CallsInfo,callsList$Alleles)



    #Start at high cutoff
    callsList$SampleData$cntMedian<-apply(callsList$CombCount,MARGIN = 2,median,na.rm=T)
    callsList$SampleData$cntMAD<-apply(callsList$CombCount,MARGIN = 2,mad,na.rm=T)
    callsList$SampleData$InCntAnalysis<-callsList$SampleData$cntMAD>1&callsList$SampleData$cntMedian>2&callsList$SampleData$mapped
    callsList$zCnt<-t(t(callsList$CombCount)-callsList$SampleData$cntMedian/callsList$SampleData$cntMAD)
    callsList$CallsInfo$zCntRow<-rowMeans(callsList$zCnt[,callsList$SampleData$InCntAnalysis],na.rm=T)
    zCntRowMedian<-median(callsList$CallsInfo$zCntRow,na.rm=T)
    zCntRowMad<-mad(callsList$CallsInfo$zCntRow,na.rm=T)
    callsList$CallsInfo$cntInner99.99Logic<-
        callsList$CallsInfo$zCntRow>=(zCntRowMedian-3.89*zCntRowMad)&
        callsList$CallsInfo$zCntRow<=(zCntRowMedian+3.89*zCntRowMad)
    hist(callsList$CallsInfo$zCntRow,100000,xlim=c(zCntRowMedian-3.89*2*zCntRowMad,zCntRowMedian+3.89*2*zCntRowMad))
    abline(v = c(zCntRowMedian-3.89*zCntRowMad,zCntRowMedian+3.89*zCntRowMad),col="red")

    #Redo with less biased starting parameters
    callsList$SampleData$cntMedian<-apply(callsList$CombCount[callsList$CallsInfo$cntInner99.99Logic,],MARGIN = 2,median,na.rm=T)
    callsList$SampleData$cntMAD   <-apply(callsList$CombCount[callsList$CallsInfo$cntInner99.99Logic,],MARGIN = 2,mad   ,na.rm=T)
    callsList$SampleData$InCntAnalysis<-callsList$SampleData$cntMAD>1&callsList$SampleData$cntMedian>2&callsList$SampleData$mapped
    callsList$zCnt<-t(t(callsList$CombCount)-callsList$SampleData$cntMedian/callsList$SampleData$cntMAD)
    callsList$CallsInfo$zCntRow<-rowMeans(callsList$zCnt[,callsList$SampleData$InCntAnalysis],na.rm=T)
    zCntRowMedian<-median(callsList$CallsInfo$zCntRow,na.rm=T)
    zCntRowMad<-mad(callsList$CallsInfo$zCntRow,na.rm=T)
    callsList$CallsInfo$cntInner99Logic<-
        callsList$CallsInfo$zCntRow>=(zCntRowMedian-2.58*zCntRowMad)&
        callsList$CallsInfo$zCntRow<=(zCntRowMedian+2.58*zCntRowMad)
    hist(callsList$CallsInfo$zCntRow[callsList$CallsInfo$cntInner99.99Logic],100000,xlim=c(zCntRowMedian-3.89*zCntRowMad,zCntRowMedian+3.89*zCntRowMad))

    callsList$CallsInfo$hetFreqProg


    #Redo at lower cutoffjhhhn bbbn
    cntMedian<-apply(callsList$CombCount[GbsDf$cntInner99.99Logic,],MARGIN = 2,median)
    cntMAD   <-apply(callsList$CombCount[GbsDf$cntInner99.99Logic,],MARGIN = 2,mad)
    GbsDf$zCnt<-t((t(callsList$CombCount)-cntMedian)/cntMAD)
    GbsDf$zCntRow<-rowMeans(GbsDf$zCnt)
    zCntRowMedian<-median(GbsDf$zCntRow[GbsDf$cntInner99.99Logic])
    zCntRowMad<-mad(GbsDf$zCntRow[GbsDf$cntInner99.99Logic])
    GbsDf$cntInner99Logic<-
        GbsDf$zCntRow>=(zCntRowMedian-2.58*zCntRowMad)&
        GbsDf$zCntRow<=(zCntRowMedian+2.58*zCntRowMad)

}
if(loadStatisticsFromFile){
}




