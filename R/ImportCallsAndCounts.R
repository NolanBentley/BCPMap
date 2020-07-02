require("stringi")
require("pbapply")
require("parallel")
require("doSNOW")
GBSImport<-function(file_calls,file_counts,cores=NULL){
    print("Initializing GBS Import...")

    #Parralellize
    if(is.null(cores)){cores<-detectCores()-1}
    cl <- makeCluster(cores)
    registerDoSNOW (cl)

    #Import Data
    no_col_calls <- max(count.fields(file_calls, sep = ","))
    #no_col_counts <- max(count.fields(file_counts, sep = ","))
    csv_calls<-read.csv(file_calls,stringsAsFactors = FALSE,header = FALSE,fill=TRUE,col.names=1:no_col_calls)
    csv_counts<-read.csv(file_counts,stringsAsFactors = FALSE,header = FALSE,fill=TRUE)

    #Call File Parameters
    Call_GenoRow = 1
    Call_SnpCol = 1
    Call_GenoColBeg = 5
    Call_GenoColEnd<-max(which(!""==csv_calls [Call_GenoRow,]))
    Call_SnpRowBeg = 3
    Call_SnpRowEnd <-max(which(!""==csv_calls [,Call_SnpCol]))

    #Count File Parameters
    Cnt_GenoRow = 1
    Cnt_SnpCol = 1
    Cnt_SnpRowEnd = 0
    Cnt_GenoColBeg = 5
    Cnt_GenoColEnd <-max(which(!""==csv_counts[Cnt_GenoRow ,]))
    Cnt_SnpRowBeg = 3
    Cnt_SnpRowEnd  <-max(which(!""==csv_counts[,Cnt_SnpCol ]))

    #Check that they match
    print("Testing Calls and Counts")
    print("Number of Genotypes:")
    if(!all(csv_calls [Call_GenoRow,Call_GenoColBeg:Call_GenoColEnd]==csv_counts[Cnt_GenoRow ,seq(Cnt_GenoColBeg,Cnt_GenoColEnd-1,by=2)])){stop("Geno Columns don't equal")
    }else{print(paste0(Call_GenoColEnd-Call_GenoColBeg+1," Genotypes Match"))}
    print("Number of SNPs:")
    if(!all(csv_calls [Call_SnpRowBeg:Call_SnpRowEnd,Call_SnpCol]==csv_counts[Cnt_SnpRowBeg:Cnt_SnpRowEnd,Cnt_SnpCol])){stop("SNP Rows don't equal")
    }else{print(paste0(Call_SnpRowEnd-Call_SnpRowBeg+1  ," SNPs Match"))}

    #Import Data
    Snps       <-csv_calls[Call_SnpRowBeg:Call_SnpRowEnd,Call_SnpCol]
    Genotypes  <-csv_calls[Call_GenoRow,Call_GenoColBeg:Call_GenoColEnd]
    Calls      <-csv_calls[Call_SnpRowBeg:Call_SnpRowEnd,Call_GenoColBeg:Call_GenoColEnd]
    Cnt_Geno   <-csv_counts[Cnt_GenoRow,Cnt_GenoColBeg:Cnt_GenoColEnd]
    Cnt_Snp    <-csv_counts[Cnt_SnpRowBeg:Cnt_SnpRowEnd,Cnt_SnpCol]
    Counts     <-csv_counts[Cnt_SnpRowBeg:Cnt_SnpRowEnd,Cnt_GenoColBeg:Cnt_GenoColEnd]
    rm(list=c("csv_calls","csv_counts"))

    #Remove Left from Right Counts
    Cnt_Geno   <-Cnt_Geno[1:ncol(Counts)%%2==1]
    LeftCounts <-Counts[,1:ncol(Counts)%%2==1]
    RightCounts<-Counts[,1:ncol(Counts)%%2==0]
    LeftCounts [LeftCounts =="?"]<-0
    RightCounts[RightCounts=="?"]<-0
    print("Converting left counts to numbers...")
    LeftCounts <-pbapply(LeftCounts ,MARGIN = 2, as.numeric)
    print("Converting right counts to numbers...")
    RightCounts<-pbapply(RightCounts,MARGIN = 2, as.numeric)
    rm(list=c("Counts"))

    #Recode Calls and seperate left from right
    Calls[Calls=="|"]<-"?"
    Calls[Calls=="--"]<-"?"
    print("Substringing Left Calls...")
    LeftCalls <-pbapply(Calls,MARGIN = 2,FUN = stri_sub,from=1,to=1)
    print("Substringing Right Calls...")
    RightCalls<-pbapply(Calls,MARGIN = 2,FUN = stri_sub,from=3,to=3)
    rm(list=c("Calls"))

    #Determine major and minor alleles per locus
    print("Counting alleles...")
    unique_calls<-c("A","T","C","G","-")
    Allele_counts<-matrix(nrow=nrow(LeftCalls),ncol=length(unique_calls))
    for(i in 1: length(unique_calls)){
        Allele_counts[,i]<-rowSums(LeftCalls==unique_calls[i])+rowSums(RightCalls==unique_calls[i])
    }
    print("Allele counts found.")
    print("Ranking alleles...")
    colnames(Allele_counts)<-unique_calls
    multiallele<-rowSums(Allele_counts>0)>2
    Allele_ranks<-t(pbapply(Allele_counts,1,rank,ties.method = "random"))
    print("Finding major alleles...")
    MajorAllele<-unique_calls[pbapply(Allele_ranks,1,match,x=length(unique_calls)  )]
    print("Finding minor alleles...")
    MinorAllele<-unique_calls[pbapply(Allele_ranks,1,match,x=length(unique_calls)-1)]
    print("Major and Minor Alleles per SNP identified")

    MajorCount<-matrix(0,ncol=ncol(LeftCalls),nrow = nrow(LeftCalls))
    MinorCount<-matrix(0,ncol=ncol(LeftCalls),nrow = nrow(LeftCalls))

    LeftMajCnt<-LeftCalls==MajorAllele
    LeftMinCnt<-LeftCalls==MinorAllele
    RightMajCnt<-RightCalls==MajorAllele
    RightMinCnt<-RightCalls==MinorAllele
    Calls=matrix(paste0(LeftCalls,"/",RightCalls),
                 nrow=nrow(LeftCalls),
                 ncol=ncol(LeftCalls),
                 dimnames = list(Snps,Genotypes)
    )

    MajorCount[LeftMajCnt]<-LeftCounts[LeftMajCnt]
    MinorCount[LeftMinCnt]<-LeftCounts[LeftMinCnt]
    rm(list=c("LeftMajCnt","LeftMinCnt"))

    MajorCount[RightMajCnt]<-RightCounts[RightMajCnt]
    MinorCount[RightMinCnt]<-RightCounts[RightMinCnt]
    rm(list=c("RightMajCnt","RightMinCnt"))
    print("RD per allele determined.")

    CombCount<-LeftCounts+RightCounts

    #Name the outputs
    dimnames(MajorCount)<-list(Snps,Genotypes)
    dimnames(MinorCount)<-list(Snps,Genotypes)
    Alleles<-cbind(MajorAllele,MinorAllele)
    rownames(Alleles)<-Snps
    print("Output data formatted.")

    #Data to Return
    print("Returning output...!")

    stopCluster(cl)

    return(list(Alleles=Alleles,
                MajorAlleleCount=MajorCount,
                MinorAlleleCount=MinorCount,
                CombCount=CombCount,
                IsMultiAllele=multiallele,
                Calls=Calls,
                LeftCalls=LeftCalls,
                RightCalls=RightCalls,
                LeftCounts=LeftCounts,
                RightCounts=RightCounts
                )
           )
}

GBSCallMajMinCoder<-function(calls,majAllele,minAllele){
    Het1<-paste0(majAllele,"/",minAllele)
    Het2<-paste0(minAllele,"/",majAllele)
    callsCoded<-matrix(numeric(),nrow=nrow(calls),ncol=ncol(calls),dimnames = dimnames(calls))
    callsCoded[calls==matrix(paste0(majAllele,"/"),nrow=nrow(calls),ncol=ncol(calls))]<-0
    callsCoded[calls==matrix(Het1,nrow=nrow(calls),ncol=ncol(calls))]<-1
    callsCoded[calls==matrix(Het2,nrow=nrow(calls),ncol=ncol(calls))]<-1
    callsCoded[calls==matrix(paste0(minAllele,"/"),nrow=nrow(calls),ncol=ncol(calls))]<-2
    return(callsCoded)
}

readCallsAndCounts<-function(rawPath,
                             rowStart=3,
                             colStart=5,
                             rowsKept=0,
                             colEnd=0,
                             type="impute"){
    #Open and format calls and counts files
    #Make rowsKept NA where undesired
    no_col_calls <- max(count.fields(rawPath, sep = ","))
    headerNames<-gsub("\"$","",gsub("^\"","",unlist(strsplit(readLines(rawPath,n = 1),","))))
    headerSupplemental<-paste0("x",1:no_col_calls)
    headerNames<-c(headerNames,headerSupplemental[(1:length(headerSupplemental))>length(headerNames)])
    currCsv<-read.csv(rawPath,header = F,col.names = headerNames,colClasses = "character",fill = T,stringsAsFactors = F)
    dimnames(currCsv)<-list(currCsv[,1],headerNames)
    if(colEnd==0){colEnd=length(headerNames)}
    if(all(rowsKept==0)){
        currRows<-rowStart:nrow(currCsv)
    }else{
        if(!all(rownames(currCsv)==names(rowsKept))){stop("mismatched rows")}
        currRows<-which(currCsv[,1]==rowsKept)
    }
    currCsv<-currCsv[currRows,]
    currMarkers<-currCsv[,1]
    currChr    <-currCsv[,2]
    currPos    <-as.numeric(currCsv[,3])
    currRef    <-currCsv[,4]
    currCsv    <-currCsv[,colStart:colEnd]
    #If any numbers in top 10 rows, imputes as count files
    if(type=="impute"){
        if(mean(is.na(suppressWarnings(apply(currCsv[1:10,],c(1,2),as.numeric))))==1){
            type="calls"
        }else{
            type="counts"
        }
    }
    if(type=="calls"){
        return(list(calls=currCsv,
                    callsInfo=cbind(Markers=currMarkers,Chr=currChr,Pos=currPos,Ref=currRef)))
    }
    if(type=="counts"){
        leftInd<-seq(1,ncol(currCsv)-1,by = 2)
        return(list(leftCnt= currCsv[,leftInd]  ,
                    rightCnt=currCsv[,leftInd+1],
                    callsInfo=cbind(Markers=currMarkers,Chr=currChr,Pos=currPos,Ref=currRef)))
    }
}

multiCallListImporter<-function(fileDf,datasetNames){
    if(length(datasetNames)!=nrow(fileDf)){stop("Dataset names not equal to number of datasets")}
    for(i in 1:length(datasetNames)){
        if(i==1){callsListList<-list()}
        callsListList[[i]]<-GBSImport(file_calls = fileDf[i,1],fileDf[i,2])
        if(i==length(datasetNames)){
            names(callsListList)<-datasetNames
        }
    }
    return(callsListList)
}
