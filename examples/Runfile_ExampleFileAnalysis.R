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
load_all()
require(parallel)
require(snow)
require(dplyr)
cl<-makeCluster(min(detectCores()-1,36))

markerRow<-3
fileDf<-data.frame(calls=c("../GitExamples/LxOMapCombinedDivSubset_calls.csv",#This should be the most limiting file to save memory
                           "../GitExamples/LxOSimulatedShotgun_4G_130bp12x25min_Vs_OaxMainV1_3-3-15_20190716_182609_sam_calls.csv",
                           "../GitExamples/PecanCRRDiv_844G_PstI_Vs_Carya_illinoinensis_var_87MX3.mainGenome_3-3-15_20190513_144932_sam_calls.csv"),
                   stringsAsFactors = F)
fileDf$counts<-gsub("calls\\.csv","counts.csv",fileDf$calls)
if(!all(file.exists(fileDf$calls))){stop("Call files not found")}
if(!all(file.exists(fileDf$counts))){stop("Count files corresponding to the names of the call files not found")}

#Produces subset versions of the input call and count files with a combined marker set
SubsetCallsAndCountsFiles(fileDf)

#Check subsets
fileDf_subset<-apply(fileDf,c(1,2),function(x){paste0(gsub("\\.csv","",x),"_subset.csv")})
if(!all(file.exists(fileDf_subset))){stop("A subset file was not found")}
markerNames<-parApply(cl,fileDf_subset,c(1,2),cutFunction)
if(!all(c(all(markerNames[,1,1]==markerNames[,2,1]),all(markerNames[,1,1]==markerNames[,3,1]),
          all(markerNames[,1,1]==markerNames[,1,2]),all(markerNames[,1,1]==markerNames[,2,2]),
          all(markerNames[,1,1]==markerNames[,3,2])))){stop("Subsets don't match")}

#Calculate call file statistics
csvList<-parApply(cl,fileDf_subset,c(1,2),readCallsAndCounts)
callsList<-list(calls=bind_cols(lapply(csvList,function(x){x[[1]]})),callsInfo=csvList[[1]][[2]])
colnames(callsList$calls)

stopCluster(cl)
callsListList$LxO<-GBSImport(file_calls = fileDf_subset[1,1],fileDf_subset[1,2])
callsListList$Sim<-GBSImport(file_calls = fileDf_subset[2,1],fileDf_subset[2,2])
callsListList$Div<-GBSImport(file_calls = fileDf_subset[3,1],fileDf_subset[3,2])

