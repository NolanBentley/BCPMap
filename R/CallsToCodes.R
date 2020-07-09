codeCalls<-function(Calls,majAlleles,minAlleles,homoCallAddedText="/"){
    #Call the progeny as homoMaj homoMin or Het
    #Set homoCallAddedText = "" if homozygous calls do not have a "/" appended
    Codes<-matrix(NA,nrow(Calls),ncol(Calls),dimnames=dimnames(Calls))
    Codes[Calls==matrix(paste0(majAlleles,homoCallAddedText),nrow(Calls),ncol(Calls))]<-0
    Codes[Calls==matrix(paste0(minAlleles,homoCallAddedText),nrow(Calls),ncol(Calls))]<-2
    Codes[Calls==matrix(paste0(majAlleles,"/",minAlleles   ),nrow(Calls),ncol(Calls))]<-1
    Codes[Calls==matrix(paste0(minAlleles,"/",majAlleles   ),nrow(Calls),ncol(Calls))]<-1
    return(Codes)
}

adjustCodedCallsByPhase<-function(codes,
                                  phases,
                                  homoMajNumber=0,
                                  hetNumber=1,
                                  homoMinNumber=1,
                                  missNumber=0.5){
    #Recode numbers as needed
    codesOrig<-codes
    codes[codesOrig==0]<-homoMajNumber
    codes[codesOrig==2]<-homoMinNumber
    codes[codesOrig==1]<-hetNumber
    #Swap phase on even testcross phases
    codes[which(phases==2|phases==4),]<-hetNumber-codes[which(phases==2|phases==4),]
    return(codes)
}
