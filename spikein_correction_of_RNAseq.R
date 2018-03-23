##########################################
# Spike-in Correction of RNAseq datasets #
##########################################

#Date
#18.03.19

#RNA-seq libraries from experimental S.cerevisiae have a consistent S.pombe spiked in to it. I need to normalise againset this to ensure that there is consistent amounts of total RNA. 
#The files are an output of ........., and exist as counts for each gene.

#Load viridis color palette for graphing
if (!require("viridis")) {
  install.packages("viridis")
  library(viridis)
}

#Load dplyr for maniulating matrices
if (!require("dplyr")) {
  install.packages("dplyr")
  library(dplyr)
}


# Load count data and create table of all samples
path <- '~/Google_Drive/Lab/Data/RNA/exp_counts/'

pathToFiles_A <- list(paste0(path, 'matt_49_featureCounts.txt'),
                      paste0(path, 'matt_50_featureCounts.txt'),
                      paste0(path, 'matt_51_featureCounts.txt'),
                      paste0(path, 'mrp66_featureCounts.txt'),
                      paste0(path, 'mrp67_featureCounts.txt'),
                      paste0(path, 'mrp68_featureCounts.txt'),
                      paste0(path, 'mrp69_featureCounts.txt'),
                      paste0(path, 'mrp70_featureCounts.txt'),
                      paste0(path, 'mrp71_featureCounts.txt'),
                      paste0(path, 'mrp72_featureCounts.txt'),
                      paste0(path, 'mrp73_featureCounts.txt'),
                      paste0(path, 'mrp74_featureCounts.txt'),
                      paste0(path, 'mrp75_featureCounts.txt'),
                      paste0(path, 'mrp76_featureCounts.txt'),
                      paste0(path, 'mrp77_featureCounts.txt'),
                      paste0(path, 'mrp78_featureCounts.txt'),
                      paste0(path, 'mrp79_featureCounts.txt'),
                      paste0(path, 'mrp80_featureCounts.txt'),
                      paste0(path, 'mrp81_featureCounts.txt')
)

path <- '~/Google_Drive/Lab/Data/RNA/YCG1/'
pathToFiles_B <- list(paste0(path, 'mrp82_L_featureCounts.txt'),
                      paste0(path, 'mrp83_L_featureCounts.txt'),
                      paste0(path, 'mrp84_L_featureCounts.txt'),
                      paste0(path, 'mrp85_L_featureCounts.txt'),
                      paste0(path, 'mrp86_L_featureCounts.txt'),
                      paste0(path, 'mrp87_L_featureCounts.txt'),
                      paste0(path, 'mrp88_L_featureCounts.txt'),
                      paste0(path, 'mrp89_L_featureCounts.txt')
)

pathToFiles<-append(pathToFiles_A, pathToFiles_B)

counts <- lapply(pathToFiles, read.table, header = T, stringsAsFactors=F)

#Assign the different data set names
sampleNames <- list('matt49', 'matt50', 'matt51', 'mrp66','mrp67', 'mrp68', 'mrp69', 'mrp70','mrp71', 'mrp72', 'mrp73', 'mrp74','mrp75', 'mrp76', 'mrp77', 'mrp78','mrp79', 'mrp80', 'mrp81', 'mrp82', 'mrp83', 'mrp84','mrp85', 'mrp86', 'mrp87', 'mrp88','mrp89')
for(i in 1:length(counts)){
  colnames(counts[[i]])[7] <- sampleNames[i]
}

#Merge my datasets together
counts_all <- Reduce(dplyr::inner_join, counts)
head(counts_all)

## Load count data and create table of all samples
path <- '~/Google_Drive/Lab/Data/RNA/spike_counts/'

pathToFiles_A <- list(paste0(path, 'matt_49_spike_featureCounts.txt'),
                      paste0(path, 'matt_50_spike_featureCounts.txt'),
                      paste0(path, 'matt_51_spike_featureCounts.txt'),
                      paste0(path, 'mrp66_spike_featureCounts.txt'),
                      paste0(path, 'mrp67_spike_featureCounts.txt'),
                      paste0(path, 'mrp68_spike_featureCounts.txt'),
                      paste0(path, 'mrp69_spike_featureCounts.txt'),
                      paste0(path, 'mrp70_spike_featureCounts.txt'),
                      paste0(path, 'mrp71_spike_featureCounts.txt'),
                      paste0(path, 'mrp72_spike_featureCounts.txt'),
                      paste0(path, 'mrp73_spike_featureCounts.txt'),
                      paste0(path, 'mrp74_spike_featureCounts.txt'),
                      paste0(path, 'mrp75_spike_featureCounts.txt'),
                      paste0(path, 'mrp76_spike_featureCounts.txt'),
                      paste0(path, 'mrp77_spike_featureCounts.txt'),
                      paste0(path, 'mrp78_spike_featureCounts.txt'),
                      paste0(path, 'mrp79_spike_featureCounts.txt'),
                      paste0(path, 'mrp80_spike_featureCounts.txt'),
                      paste0(path, 'mrp81_spike_featureCounts.txt')
)

path <- '~/Google_Drive/Lab/Data/RNA/YCG1/'
pathToFiles_B <- list(paste0(path, 'mrp82_L_spike_featureCounts.txt'),
                      paste0(path, 'mrp83_L_spike_featureCounts.txt'),
                      paste0(path, 'mrp84_L_spike_featureCounts.txt'),
                      paste0(path, 'mrp85_L_spike_featureCounts.txt'),
                      paste0(path, 'mrp86_L_spike_featureCounts.txt'),
                      paste0(path, 'mrp87_L_spike_featureCounts.txt'),
                      paste0(path, 'mrp88_L_spike_featureCounts.txt'),
                      paste0(path, 'mrp89_L_spike_featureCounts.txt')
)

pathToFiles<-append(pathToFiles_A, pathToFiles_B)


counts <- lapply(pathToFiles, read.table, header = T, stringsAsFactors=F)

#Assign the different data set names
sampleNames <- list('matt49', 'matt50', 'matt51', 'mrp66','mrp67', 'mrp68', 'mrp69', 'mrp70','mrp71', 'mrp72', 'mrp73', 'mrp74','mrp75', 'mrp76', 'mrp77', 'mrp78','mrp79', 'mrp80', 'mrp81', 'mrp82', 'mrp83', 'mrp84','mrp85', 'mrp86', 'mrp87', 'mrp88','mrp89')
for(i in 1:length(counts)){
  colnames(counts[[i]])[7] <- sampleNames[i]
}

#merge my datasets together
counts_all_spike <- Reduce(dplyr::inner_join, counts)
head(counts_all_spike)


### Calculate CPMs ###

#Use the existing file as a template for for inputting cpms
cpm_all_spike<-counts_all_spike
cpm_all<-counts_all


#Calculate CPM for every for each, step by step
for (i in 7:33){
  #Do it for spike
  total_counts<-sum(counts_all_spike[,i])  
  cpm_counts<-(counts_all_spike[,i]/total_counts)*10^6
  cpm_all_spike[,i]<-cpm_counts
  #Do it for experimental
  total_counts<-sum(counts_all[,i])  
  cpm_counts<-(counts_all[,i]/total_counts)*10^6
  cpm_all[,i]<-cpm_counts
  print(i)
}

#Calculate RPKM for every for each, step by step
rpkm_all_spike<-counts_all_spike
rpkm_all<-counts_all

for (i in 7:33){
  #Do it for spike
  total_counts<-sum(counts_all_spike[,i])  
  cpm_counts<-(counts_all_spike[,i]/(counts_all_spike[,6]/1000*total_counts/10^6))
  rpkm_all_spike[,i]<-cpm_counts
  #Do it for experimental
  total_counts<-sum(counts_all[,i])  
  cpm_counts<-(counts_all[,i]/(counts_all[,6]/1000*total_counts/10^6))
  rpkm_all[,i]<-cpm_counts
  print(i)
}

#Calculate TPM for every for each, step by step
tpm_all_spike<-counts_all_spike
tpm_all<-counts_all

for (i in 7:33){
  #Do it for spike
  total_counts<-sum(rpkm_all_spike[,i])  
  cpm_counts<-(rpkm_all_spike[,i]/(total_counts/10^6))
  tpm_all_spike[,i]<-cpm_counts
  #Do it for experimental
  total_counts<-sum(rpkm_all[,i])  
  cpm_counts<-(rpkm_all[,i]/(total_counts/10^6))
  tpm_all[,i]<-cpm_counts
  print(i)
}


# Brn1 - AAway experiments
rpkm_all_spike_brn1<-cbind(rpkm_all_spike[,1:6],rpkm_all_spike[,18],rpkm_all_spike[,19],rpkm_all_spike[,22],rpkm_all_spike[,24],rpkm_all_spike[,20],rpkm_all_spike[,21],rpkm_all_spike[,23],rpkm_all_spike[,25])
rpkm_all_brn1<-cbind(rpkm_all[,1:6],rpkm_all[,18],rpkm_all[,19],rpkm_all[,22],rpkm_all[,24],rpkm_all[,20],rpkm_all[,21],rpkm_all[,23],rpkm_all[,25])
tpm_all_spike_brn1<-cbind(tpm_all_spike[,1:6],tpm_all_spike[,18],tpm_all_spike[,19],tpm_all_spike[,22],tpm_all_spike[,24],tpm_all_spike[,20],tpm_all_spike[,21],tpm_all_spike[,23],tpm_all_spike[,25])
tpm_all_brn1<-cbind(tpm_all[,1:6],tpm_all[,18],tpm_all[,19],tpm_all[,22],tpm_all[,24],tpm_all[,20],tpm_all[,21],tpm_all[,23],tpm_all[,25])
cpm_all_spike_brn1<-cbind(cpm_all_spike[,1:6],cpm_all_spike[,18],cpm_all_spike[,19],cpm_all_spike[,22],cpm_all_spike[,24],cpm_all_spike[,20],cpm_all_spike[,21],cpm_all_spike[,23],cpm_all_spike[,25])
cpm_all_brn1<-cbind(cpm_all[,1:6],cpm_all[,18],cpm_all[,19],cpm_all[,22],cpm_all[,24],cpm_all[,20],cpm_all[,21],cpm_all[,23],cpm_all[,25])

header_brn1<-c(row.names(rpkm_all_spike_brn1)[1:6],'brn1_wt_norap_rep1','brn1_wt_rap_rep1','brn1_wt_rap_rep2','brn1_wt_rap_rep3','brn1_brn1_norap_rep1','brn1_brn1_rap_rep1','brn1_brn1_rap_rep2','brn1_brn1_rap_rep3')

colnames(rpkm_all_spike_brn1)<-header_brn1
colnames(rpkm_all_brn1)<-header_brn1
colnames(cpm_all_spike_brn1)<-header_brn1
colnames(cpm_all_brn1)<-header_brn1
colnames(tpm_all_spike_brn1)<-header_brn1
colnames(tpm_all_brn1)<-header_brn1


# ycs4-2 experiments
rpkm_all_spike_ycs4<-cbind(rpkm_all_spike[,1:6],rpkm_all_spike[,10],rpkm_all_spike[,14],rpkm_all_spike[,7],rpkm_all_spike[,11],rpkm_all_spike[,15],rpkm_all_spike[,8],rpkm_all_spike[,12],rpkm_all_spike[,16],rpkm_all_spike[,9],rpkm_all_spike[,13],rpkm_all_spike[,17])
rpkm_all_ycs4<-cbind(rpkm_all[,1:6],rpkm_all[,10],rpkm_all[,14],rpkm_all[,7],rpkm_all[,11],rpkm_all[,15],rpkm_all[,8],rpkm_all[,12],rpkm_all[,16],rpkm_all[,9],rpkm_all[,13],rpkm_all[,17])
tpm_all_spike_ycs4<-cbind(tpm_all_spike[,1:6],tpm_all_spike[,10],tpm_all_spike[,14],tpm_all_spike[,7],tpm_all_spike[,11],tpm_all_spike[,15],tpm_all_spike[,8],tpm_all_spike[,12],tpm_all_spike[,16],tpm_all_spike[,9],tpm_all_spike[,13],tpm_all_spike[,17])
tpm_all_ycs4<-cbind(tpm_all[,1:6],tpm_all[,10],tpm_all[,14],tpm_all[,7],tpm_all[,11],tpm_all[,15],tpm_all[,8],tpm_all[,12],tpm_all[,16],tpm_all[,9],tpm_all[,13],tpm_all[,17])
cpm_all_spike_ycs4<-cbind(cpm_all_spike[,1:6],cpm_all_spike[,10],cpm_all_spike[,14],cpm_all_spike[,7],cpm_all_spike[,11],cpm_all_spike[,15],cpm_all_spike[,8],cpm_all_spike[,12],cpm_all_spike[,16],cpm_all_spike[,9],cpm_all_spike[,13],cpm_all_spike[,17])
cpm_all_ycs4<-cbind(cpm_all[,1:6],cpm_all[,10],cpm_all[,14],cpm_all[,7],cpm_all[,11],cpm_all[,15],cpm_all[,8],cpm_all[,12],cpm_all[,16],cpm_all[,9],cpm_all[,13],cpm_all[,17])

header_ycs4<-c(row.names(rpkm_all_spike_ycs4)[1:6],'ycs4_wt_23C_rep1','ycs4_wt_23C_rep2','ycs4_wt_37C_rep1','ycs4_wt_37C_rep2','ycs4_wt_37C_rep3','ycs4_ycs4_23C_rep1','ycs4_ycs4_23C_rep2','ycs4_ycs4_23C_rep3','ycs4_ycs4_37C_rep1','ycs4_ycs4_37C_rep2','ycs4_ycs4_37C_rep3')

colnames(rpkm_all_spike_ycs4)<-header_ycs4
colnames(rpkm_all_ycs4)<-header_ycs4
colnames(cpm_all_spike_ycs4)<-header_ycs4
colnames(cpm_all_ycs4)<-header_ycs4
colnames(tpm_all_spike_ycs4)<-header_ycs4
colnames(tpm_all_ycs4)<-header_ycs4

# ycg1-2 experiments
rpkm_all_spike_ycg1<-cbind(rpkm_all_spike[,1:6],rpkm_all_spike[,26],rpkm_all_spike[,30],rpkm_all_spike[,27],rpkm_all_spike[,31],rpkm_all_spike[,28],rpkm_all_spike[,32],rpkm_all_spike[,29],rpkm_all_spike[,33])
rpkm_all_ycg1<-cbind(rpkm_all[,1:6],rpkm_all[,26],rpkm_all[,30],rpkm_all[,27],rpkm_all[,31],rpkm_all[,28],rpkm_all[,32],rpkm_all[,29],rpkm_all[,33])
tpm_all_spike_ycg1<-cbind(tpm_all_spike[,1:6],tpm_all_spike[,26],tpm_all_spike[,30],tpm_all_spike[,27],tpm_all_spike[,31],tpm_all_spike[,28],tpm_all_spike[,32],tpm_all_spike[,29],tpm_all_spike[,33])
tpm_all_ycg1<-cbind(tpm_all[,1:6],tpm_all[,26],tpm_all[,30],tpm_all[,27],tpm_all[,31],tpm_all[,28],tpm_all[,32],tpm_all[,29],tpm_all[,33])
cpm_all_spike_ycg1<-cbind(cpm_all_spike[,1:6],cpm_all_spike[,26],cpm_all_spike[,30],cpm_all_spike[,27],cpm_all_spike[,31],cpm_all_spike[,28],cpm_all_spike[,32],cpm_all_spike[,29],cpm_all_spike[,33])
cpm_all_ycg1<-cbind(cpm_all[,1:6],cpm_all[,26],cpm_all[,30],cpm_all[,27],cpm_all[,31],cpm_all[,28],cpm_all[,32],cpm_all[,29],cpm_all[,33])

header_ycg1<-c(row.names(rpkm_all_spike_ycg1)[1:6],'ycg1_wt_23C_rep1','ycg1_wt_23C_rep2','ycg1_wt_37C_rep1','ycg1_wt_37C_rep2','ycg1_ycg1_23C_rep1','ycg1_ycg1_23C_rep2','ycg1_ycg1_37C_rep1','ycg1_ycg1_37C_rep2')

colnames(rpkm_all_spike_ycg1)<-header_ycg1
colnames(rpkm_all_ycg1)<-header_ycg1
colnames(cpm_all_spike_ycg1)<-header_ycg1
colnames(cpm_all_ycg1)<-header_ycg1
colnames(tpm_all_spike_ycg1)<-header_ycg1
colnames(tpm_all_ycg1)<-header_ycg1

#Give me expressed only those that are expressed

rpkm_all_spike_ycg1_exp<-rpkm_all_spike_ycg1[which(rpkm_all_spike_ycg1[,7]>1 & rpkm_all_spike_ycg1[,8]>1 & rpkm_all_spike_ycg1[,9]>1 & rpkm_all_spike_ycg1[,10]>1 & rpkm_all_spike_ycg1[,11]>1 & rpkm_all_spike_ycg1[,12]>1 & rpkm_all_spike_ycg1[,13]>1 & rpkm_all_spike_ycg1[,14]>1),]
rpkm_all_ycg1_exp<-rpkm_all_ycg1[which(rpkm_all_ycg1[,7]>1 & rpkm_all_ycg1[,8]>1 & rpkm_all_ycg1[,9]>1 & rpkm_all_ycg1[,10]>1 & rpkm_all_ycg1[,11]>1 & rpkm_all_ycg1[,12]>1 & rpkm_all_ycg1[,13]>1 & rpkm_all_ycg1[,14]>1),]
tpm_all_spike_ycg1_exp<-tpm_all_spike_ycg1[which(tpm_all_spike_ycg1[,7]>1 & tpm_all_spike_ycg1[,8]>1 & tpm_all_spike_ycg1[,9]>1 & tpm_all_spike_ycg1[,10]>1 & tpm_all_spike_ycg1[,11]>1 & tpm_all_spike_ycg1[,12]>1 & tpm_all_spike_ycg1[,13]>1 & tpm_all_spike_ycg1[,14]>1),]
tpm_all_ycg1_exp<-tpm_all_ycg1[which(tpm_all_ycg1[,7]>1 & tpm_all_ycg1[,8]>1 & tpm_all_ycg1[,9]>1 & tpm_all_ycg1[,10]>1 & tpm_all_ycg1[,11]>1 & tpm_all_ycg1[,12]>1 & tpm_all_ycg1[,13]>1 & tpm_all_ycg1[,14]>1),]
cpm_all_spike_ycg1_exp<-cpm_all_spike_ycg1[which(cpm_all_spike_ycg1[,7]>1 & cpm_all_spike_ycg1[,8]>1 & cpm_all_spike_ycg1[,9]>1 & cpm_all_spike_ycg1[,10]>1 & cpm_all_spike_ycg1[,11]>1 & cpm_all_spike_ycg1[,12]>1 & cpm_all_spike_ycg1[,13]>1 & cpm_all_spike_ycg1[,14]>1),]
cpm_all_ycg1_exp<-cpm_all_ycg1[which(cpm_all_ycg1[,7]>1 & cpm_all_ycg1[,8]>1 & cpm_all_ycg1[,9]>1 & cpm_all_ycg1[,10]>1 & cpm_all_ycg1[,11]>1 & cpm_all_ycg1[,12]>1 & cpm_all_ycg1[,13]>1 & cpm_all_ycg1[,14]>1),]

rpkm_all_spike_ycs4_exp<-rpkm_all_spike_ycs4[which(rpkm_all_spike_ycs4[,7]>1 & rpkm_all_spike_ycs4[,8]>1 & rpkm_all_spike_ycs4[,9]>1 & rpkm_all_spike_ycs4[,10]>1 & rpkm_all_spike_ycs4[,11]>1 & rpkm_all_spike_ycs4[,12]>1 & rpkm_all_spike_ycs4[,13]>1 & rpkm_all_spike_ycs4[,14]>1& rpkm_all_spike_ycs4[,15]>1& rpkm_all_spike_ycs4[,16]>1& rpkm_all_spike_ycs4[,17]>1),]
rpkm_all_ycs4_exp<-rpkm_all_ycs4[which(rpkm_all_ycs4[,7]>1 & rpkm_all_ycs4[,8]>1 & rpkm_all_ycs4[,9]>1 & rpkm_all_ycs4[,10]>1 & rpkm_all_ycs4[,11]>1 & rpkm_all_ycs4[,12]>1 & rpkm_all_ycs4[,13]>1 & rpkm_all_ycs4[,14]>1 & rpkm_all_ycs4[,15]>1 & rpkm_all_ycs4[,16]>1 & rpkm_all_ycs4[,17]>1),]
tpm_all_spike_ycs4_exp<-tpm_all_spike_ycs4[which(tpm_all_spike_ycs4[,7]>1 & tpm_all_spike_ycs4[,8]>1 & tpm_all_spike_ycs4[,9]>1 & tpm_all_spike_ycs4[,10]>1 & tpm_all_spike_ycs4[,11]>1 & tpm_all_spike_ycs4[,12]>1 & tpm_all_spike_ycs4[,13]>1 & tpm_all_spike_ycs4[,14]>1& tpm_all_spike_ycs4[,15]>1& tpm_all_spike_ycs4[,16]>1& tpm_all_spike_ycs4[,17]>1),]
tpm_all_ycs4_exp<-tpm_all_ycs4[which(tpm_all_ycs4[,7]>1 & tpm_all_ycs4[,8]>1 & tpm_all_ycs4[,9]>1 & tpm_all_ycs4[,10]>1 & tpm_all_ycs4[,11]>1 & tpm_all_ycs4[,12]>1 & tpm_all_ycs4[,13]>1 & tpm_all_ycs4[,14]>1& tpm_all_ycs4[,15]>1& tpm_all_ycs4[,16]>1& tpm_all_ycs4[,17]>1),]
cpm_all_spike_ycs4_exp<-cpm_all_spike_ycs4[which(cpm_all_spike_ycs4[,7]>1 & cpm_all_spike_ycs4[,8]>1 & cpm_all_spike_ycs4[,9]>1 & cpm_all_spike_ycs4[,10]>1 & cpm_all_spike_ycs4[,11]>1 & cpm_all_spike_ycs4[,12]>1 & cpm_all_spike_ycs4[,13]>1 & cpm_all_spike_ycs4[,14]>1& cpm_all_spike_ycs4[,15]>1 & cpm_all_spike_ycs4[,16]>1 & cpm_all_spike_ycs4[,17]>1),]
cpm_all_ycs4_exp<-cpm_all_ycs4[which(cpm_all_ycs4[,7]>1 & cpm_all_ycs4[,8]>1 & cpm_all_ycs4[,9]>1 & cpm_all_ycs4[,10]>1 & cpm_all_ycs4[,11]>1 & cpm_all_ycs4[,12]>1 & cpm_all_ycs4[,13]>1 & cpm_all_ycs4[,14]>1& cpm_all_ycs4[,15]>1& cpm_all_ycs4[,16]>1& cpm_all_ycs4[,17]>1),]

rpkm_all_spike_brn1_exp<-rpkm_all_spike_brn1[which(rpkm_all_spike_brn1[,7]>1 & rpkm_all_spike_brn1[,8]>1 & rpkm_all_spike_brn1[,9]>1 & rpkm_all_spike_brn1[,10]>1 & rpkm_all_spike_brn1[,11]>1 & rpkm_all_spike_brn1[,12]>1 & rpkm_all_spike_brn1[,13]>1 & rpkm_all_spike_brn1[,14]>1),]
rpkm_all_brn1_exp<-rpkm_all_brn1[which(rpkm_all_brn1[,7]>1 & rpkm_all_brn1[,8]>1 & rpkm_all_brn1[,9]>1 & rpkm_all_brn1[,10]>1 & rpkm_all_brn1[,11]>1 & rpkm_all_brn1[,12]>1 & rpkm_all_brn1[,13]>1 & rpkm_all_brn1[,14]>1),]
tpm_all_spike_brn1_exp<-tpm_all_spike_brn1[which(tpm_all_spike_brn1[,7]>1 & tpm_all_spike_brn1[,8]>1 & tpm_all_spike_brn1[,9]>1 & tpm_all_spike_brn1[,10]>1 & tpm_all_spike_brn1[,11]>1 & tpm_all_spike_brn1[,12]>1 & tpm_all_spike_brn1[,13]>1 & tpm_all_spike_brn1[,14]>1),]
tpm_all_brn1_exp<-tpm_all_brn1[which(tpm_all_brn1[,7]>1 & tpm_all_brn1[,8]>1 & tpm_all_brn1[,9]>1 & tpm_all_brn1[,10]>1 & tpm_all_brn1[,11]>1 & tpm_all_brn1[,12]>1 & tpm_all_brn1[,13]>1 & tpm_all_brn1[,14]>1),]
cpm_all_spike_brn1_exp<-cpm_all_spike_brn1[which(cpm_all_spike_brn1[,7]>1 & cpm_all_spike_brn1[,8]>1 & cpm_all_spike_brn1[,9]>1 & cpm_all_spike_brn1[,10]>1 & cpm_all_spike_brn1[,11]>1 & cpm_all_spike_brn1[,12]>1 & cpm_all_spike_brn1[,13]>1 & cpm_all_spike_brn1[,14]>1),]
cpm_all_brn1_exp<-cpm_all_brn1[which(cpm_all_brn1[,7]>1 & cpm_all_brn1[,8]>1 & cpm_all_brn1[,9]>1 & cpm_all_brn1[,10]>1 & cpm_all_brn1[,11]>1 & cpm_all_brn1[,12]>1 & cpm_all_brn1[,13]>1 & cpm_all_brn1[,14]>1),]

#Expressed under all experimetnal conditions

rpkm_all_spike_brn1_exp_all<-rpkm_all_spike_brn1[which(rpkm_all_spike_brn1[,7]>1 & rpkm_all_spike_brn1[,8]>1 & rpkm_all_spike_brn1[,9]>1 & rpkm_all_spike_brn1[,10]>1 & rpkm_all_spike_brn1[,11]>1 & rpkm_all_spike_brn1[,12]>1 & rpkm_all_spike_brn1[,13]>1 & rpkm_all_spike_brn1[,14]>1&rpkm_all_spike_ycs4[,7]>1 & rpkm_all_spike_ycs4[,8]>1 & rpkm_all_spike_ycs4[,9]>1 & rpkm_all_spike_ycs4[,10]>1 & rpkm_all_spike_ycs4[,11]>1 & rpkm_all_spike_ycs4[,12]>1 & rpkm_all_spike_ycs4[,13]>1 & rpkm_all_spike_ycs4[,14]>1& rpkm_all_spike_ycs4[,15]>1 & rpkm_all_spike_ycs4[,16]>1 & rpkm_all_spike_ycs4[,17]>1&rpkm_all_spike_ycg1[,7]>1 & rpkm_all_spike_ycg1[,8]>1 & rpkm_all_spike_ycg1[,9]>1 & rpkm_all_spike_ycg1[,10]>1 & rpkm_all_spike_ycg1[,11]>1 & rpkm_all_spike_ycg1[,12]>1 & rpkm_all_spike_ycg1[,13]>1 & rpkm_all_spike_ycg1[,14]>1),]
rpkm_all_brn1_exp_all<-rpkm_all_brn1[which(rpkm_all_brn1[,7]>1 & rpkm_all_brn1[,8]>1 & rpkm_all_brn1[,9]>1 & rpkm_all_brn1[,10]>1 & rpkm_all_brn1[,11]>1 & rpkm_all_brn1[,12]>1 & rpkm_all_brn1[,13]>1 & rpkm_all_brn1[,14]>1&rpkm_all_ycs4[,7]>1 & rpkm_all_ycs4[,8]>1 & rpkm_all_ycs4[,9]>1 & rpkm_all_ycs4[,10]>1 & rpkm_all_ycs4[,11]>1 & rpkm_all_ycs4[,12]>1 & rpkm_all_ycs4[,13]>1 & rpkm_all_ycs4[,14]>1& rpkm_all_ycs4[,15]>1& rpkm_all_ycs4[,16]>1& rpkm_all_ycs4[,17]>1&rpkm_all_ycg1[,7]>1 & rpkm_all_ycg1[,8]>1 & rpkm_all_ycg1[,9]>1 & rpkm_all_ycg1[,10]>1 & rpkm_all_ycg1[,11]>1 & rpkm_all_ycg1[,12]>1 & rpkm_all_ycg1[,13]>1 & rpkm_all_ycg1[,14]>1),]
tpm_all_spike_brn1_exp_all<-tpm_all_spike_brn1[which(tpm_all_spike_brn1[,7]>1 & tpm_all_spike_brn1[,8]>1 & tpm_all_spike_brn1[,9]>1 & tpm_all_spike_brn1[,10]>1 & tpm_all_spike_brn1[,11]>1 & tpm_all_spike_brn1[,12]>1 & tpm_all_spike_brn1[,13]>1 & tpm_all_spike_brn1[,14]>1&tpm_all_spike_ycs4[,7]>1 & tpm_all_spike_ycs4[,8]>1 & tpm_all_spike_ycs4[,9]>1 & tpm_all_spike_ycs4[,10]>1 & tpm_all_spike_ycs4[,11]>1 & tpm_all_spike_ycs4[,12]>1 & tpm_all_spike_ycs4[,13]>1 & tpm_all_spike_ycs4[,14]>1& tpm_all_spike_ycs4[,15]>1& tpm_all_spike_ycs4[,16]>1& tpm_all_spike_ycs4[,17]>1&tpm_all_spike_ycg1[,7]>1 & tpm_all_spike_ycg1[,8]>1 & tpm_all_spike_ycg1[,9]>1 & tpm_all_spike_ycg1[,10]>1 & tpm_all_spike_ycg1[,11]>1 & tpm_all_spike_ycg1[,12]>1 & tpm_all_spike_ycg1[,13]>1 & tpm_all_spike_ycg1[,14]>1),]
tpm_all_brn1_exp_all<-tpm_all_brn1[which(tpm_all_brn1[,7]>1 & tpm_all_brn1[,8]>1 & tpm_all_brn1[,9]>1 & tpm_all_brn1[,10]>1 & tpm_all_brn1[,11]>1 & tpm_all_brn1[,12]>1 & tpm_all_brn1[,13]>1 & tpm_all_brn1[,14]>1&tpm_all_ycs4[,7]>1 & tpm_all_ycs4[,8]>1 & tpm_all_ycs4[,9]>1 & tpm_all_ycs4[,10]>1 & tpm_all_ycs4[,11]>1 & tpm_all_ycs4[,12]>1 & tpm_all_ycs4[,13]>1 & tpm_all_ycs4[,15]>1& tpm_all_ycs4[,16]>1& tpm_all_ycs4[,14]>1& tpm_all_ycs4[,17]>1&tpm_all_ycg1[,7]>1 & tpm_all_ycg1[,8]>1 & tpm_all_ycg1[,9]>1 & tpm_all_ycg1[,10]>1 & tpm_all_ycg1[,11]>1 & tpm_all_ycg1[,12]>1 & tpm_all_ycg1[,13]>1 & tpm_all_ycg1[,14]>1),]
cpm_all_spike_brn1_exp_all<-cpm_all_spike_brn1[which(cpm_all_spike_brn1[,7]>1 & cpm_all_spike_brn1[,8]>1 & cpm_all_spike_brn1[,9]>1 & cpm_all_spike_brn1[,10]>1 & cpm_all_spike_brn1[,11]>1 & cpm_all_spike_brn1[,12]>1 & cpm_all_spike_brn1[,13]>1 & cpm_all_spike_brn1[,14]>1&cpm_all_spike_ycs4[,7]>1 & cpm_all_spike_ycs4[,8]>1 & cpm_all_spike_ycs4[,9]>1 & cpm_all_spike_ycs4[,10]>1 & cpm_all_spike_ycs4[,11]>1 & cpm_all_spike_ycs4[,12]>1 & cpm_all_spike_ycs4[,13]>1 & cpm_all_spike_ycs4[,14]>1& cpm_all_spike_ycs4[,15]>1& cpm_all_spike_ycs4[,16]>1& cpm_all_spike_ycs4[,17]>1&cpm_all_spike_ycg1[,7]>1 & cpm_all_spike_ycg1[,8]>1 & cpm_all_spike_ycg1[,9]>1 & cpm_all_spike_ycg1[,10]>1 & cpm_all_spike_ycg1[,11]>1 & cpm_all_spike_ycg1[,12]>1 & cpm_all_spike_ycg1[,13]>1 & cpm_all_spike_ycg1[,14]>1),]
cpm_all_brn1_exp_all<-cpm_all_brn1[which(cpm_all_brn1[,7]>1 & cpm_all_brn1[,8]>1 & cpm_all_brn1[,9]>1 & cpm_all_brn1[,10]>1 & cpm_all_brn1[,11]>1 & cpm_all_brn1[,12]>1 & cpm_all_brn1[,13]>1 & cpm_all_brn1[,14]>1&cpm_all_ycs4[,7]>1 & cpm_all_ycs4[,8]>1 & cpm_all_ycs4[,9]>1 & cpm_all_ycs4[,10]>1 & cpm_all_ycs4[,11]>1 & cpm_all_ycs4[,12]>1 & cpm_all_ycs4[,13]>1 & cpm_all_ycs4[,14]>1& cpm_all_ycs4[,15]>1& cpm_all_ycs4[,16]>1& cpm_all_ycs4[,17]>1&cpm_all_ycg1[,7]>1 & cpm_all_ycg1[,8]>1 & cpm_all_ycg1[,9]>1 & cpm_all_ycg1[,10]>1 & cpm_all_ycg1[,11]>1 & cpm_all_ycg1[,12]>1 & cpm_all_ycg1[,13]>1 & cpm_all_ycg1[,14]>1),]

rpkm_all_spike_ycg1_exp_all<-rpkm_all_spike_ycg1[which(rpkm_all_spike_brn1[,7]>1 & rpkm_all_spike_brn1[,8]>1 & rpkm_all_spike_brn1[,9]>1 & rpkm_all_spike_brn1[,10]>1 & rpkm_all_spike_brn1[,11]>1 & rpkm_all_spike_brn1[,12]>1 & rpkm_all_spike_brn1[,13]>1 & rpkm_all_spike_brn1[,14]>1&rpkm_all_spike_ycs4[,7]>1 & rpkm_all_spike_ycs4[,8]>1 & rpkm_all_spike_ycs4[,9]>1 & rpkm_all_spike_ycs4[,10]>1 & rpkm_all_spike_ycs4[,11]>1 & rpkm_all_spike_ycs4[,12]>1 & rpkm_all_spike_ycs4[,13]>1 & rpkm_all_spike_ycs4[,14]>1& rpkm_all_spike_ycs4[,15]>1 & rpkm_all_spike_ycs4[,16]>1 & rpkm_all_spike_ycs4[,17]>1&rpkm_all_spike_ycg1[,7]>1 & rpkm_all_spike_ycg1[,8]>1 & rpkm_all_spike_ycg1[,9]>1 & rpkm_all_spike_ycg1[,10]>1 & rpkm_all_spike_ycg1[,11]>1 & rpkm_all_spike_ycg1[,12]>1 & rpkm_all_spike_ycg1[,13]>1 & rpkm_all_spike_ycg1[,14]>1),]
rpkm_all_ycg1_exp_all<-rpkm_all_ycg1[which(rpkm_all_brn1[,7]>1 & rpkm_all_brn1[,8]>1 & rpkm_all_brn1[,9]>1 & rpkm_all_brn1[,10]>1 & rpkm_all_brn1[,11]>1 & rpkm_all_brn1[,12]>1 & rpkm_all_brn1[,13]>1 & rpkm_all_brn1[,14]>1&rpkm_all_ycs4[,7]>1 & rpkm_all_ycs4[,8]>1 & rpkm_all_ycs4[,9]>1 & rpkm_all_ycs4[,10]>1 & rpkm_all_ycs4[,11]>1 & rpkm_all_ycs4[,12]>1 & rpkm_all_ycs4[,13]>1 & rpkm_all_ycs4[,14]>1& rpkm_all_ycs4[,15]>1& rpkm_all_ycs4[,16]>1& rpkm_all_ycs4[,17]>1&rpkm_all_ycg1[,7]>1 & rpkm_all_ycg1[,8]>1 & rpkm_all_ycg1[,9]>1 & rpkm_all_ycg1[,10]>1 & rpkm_all_ycg1[,11]>1 & rpkm_all_ycg1[,12]>1 & rpkm_all_ycg1[,13]>1 & rpkm_all_ycg1[,14]>1),]
tpm_all_spike_ycg1_exp_all<-tpm_all_spike_ycg1[which(tpm_all_spike_brn1[,7]>1 & tpm_all_spike_brn1[,8]>1 & tpm_all_spike_brn1[,9]>1 & tpm_all_spike_brn1[,10]>1 & tpm_all_spike_brn1[,11]>1 & tpm_all_spike_brn1[,12]>1 & tpm_all_spike_brn1[,13]>1 & tpm_all_spike_brn1[,14]>1&tpm_all_spike_ycs4[,7]>1 & tpm_all_spike_ycs4[,8]>1 & tpm_all_spike_ycs4[,9]>1 & tpm_all_spike_ycs4[,10]>1 & tpm_all_spike_ycs4[,11]>1 & tpm_all_spike_ycs4[,12]>1 & tpm_all_spike_ycs4[,13]>1 & tpm_all_spike_ycs4[,14]>1& tpm_all_spike_ycs4[,15]>1& tpm_all_spike_ycs4[,16]>1& tpm_all_spike_ycs4[,17]>1&tpm_all_spike_ycg1[,7]>1 & tpm_all_spike_ycg1[,8]>1 & tpm_all_spike_ycg1[,9]>1 & tpm_all_spike_ycg1[,10]>1 & tpm_all_spike_ycg1[,11]>1 & tpm_all_spike_ycg1[,12]>1 & tpm_all_spike_ycg1[,13]>1 & tpm_all_spike_ycg1[,14]>1),]
tpm_all_ycg1_exp_all<-tpm_all_ycg1[which(tpm_all_brn1[,7]>1 & tpm_all_brn1[,8]>1 & tpm_all_brn1[,9]>1 & tpm_all_brn1[,10]>1 & tpm_all_brn1[,11]>1 & tpm_all_brn1[,12]>1 & tpm_all_brn1[,13]>1 & tpm_all_brn1[,14]>1&tpm_all_ycs4[,7]>1 & tpm_all_ycs4[,8]>1 & tpm_all_ycs4[,9]>1 & tpm_all_ycs4[,10]>1 & tpm_all_ycs4[,11]>1 & tpm_all_ycs4[,12]>1 & tpm_all_ycs4[,13]>1 & tpm_all_ycs4[,14]>1& tpm_all_ycs4[,15]>1& tpm_all_ycs4[,16]>1& tpm_all_ycs4[,14]>1&tpm_all_ycg1[,7]>1 & tpm_all_ycg1[,8]>1 & tpm_all_ycg1[,9]>1 & tpm_all_ycg1[,10]>1 & tpm_all_ycg1[,11]>1 & tpm_all_ycg1[,12]>1 & tpm_all_ycg1[,13]>1 & tpm_all_ycg1[,14]>1),]
cpm_all_spike_ycg1_exp_all<-cpm_all_spike_ycg1[which(cpm_all_spike_brn1[,7]>1 & cpm_all_spike_brn1[,8]>1 & cpm_all_spike_brn1[,9]>1 & cpm_all_spike_brn1[,10]>1 & cpm_all_spike_brn1[,11]>1 & cpm_all_spike_brn1[,12]>1 & cpm_all_spike_brn1[,13]>1 & cpm_all_spike_brn1[,14]>1&cpm_all_spike_ycs4[,7]>1 & cpm_all_spike_ycs4[,8]>1 & cpm_all_spike_ycs4[,9]>1 & cpm_all_spike_ycs4[,10]>1 & cpm_all_spike_ycs4[,11]>1 & cpm_all_spike_ycs4[,12]>1 & cpm_all_spike_ycs4[,13]>1 & cpm_all_spike_ycs4[,14]>1& cpm_all_spike_ycs4[,15]>1& cpm_all_spike_ycs4[,16]>1& cpm_all_spike_ycs4[,17]>1&cpm_all_spike_ycg1[,7]>1 & cpm_all_spike_ycg1[,8]>1 & cpm_all_spike_ycg1[,9]>1 & cpm_all_spike_ycg1[,10]>1 & cpm_all_spike_ycg1[,11]>1 & cpm_all_spike_ycg1[,12]>1 & cpm_all_spike_ycg1[,13]>1 & cpm_all_spike_ycg1[,14]>1),]
cpm_all_ycg1_exp_all<-cpm_all_ycg1[which(cpm_all_brn1[,7]>1 & cpm_all_brn1[,8]>1 & cpm_all_brn1[,9]>1 & cpm_all_brn1[,10]>1 & cpm_all_brn1[,11]>1 & cpm_all_brn1[,12]>1 & cpm_all_brn1[,13]>1 & cpm_all_brn1[,14]>1&cpm_all_ycs4[,7]>1 & cpm_all_ycs4[,8]>1 & cpm_all_ycs4[,9]>1 & cpm_all_ycs4[,10]>1 & cpm_all_ycs4[,11]>1 & cpm_all_ycs4[,12]>1 & cpm_all_ycs4[,13]>1 & cpm_all_ycs4[,14]>1& cpm_all_ycs4[,15]>1& cpm_all_ycs4[,16]>1& cpm_all_ycs4[,17]>1&cpm_all_ycg1[,7]>1 & cpm_all_ycg1[,8]>1 & cpm_all_ycg1[,9]>1 & cpm_all_ycg1[,10]>1 & cpm_all_ycg1[,11]>1 & cpm_all_ycg1[,12]>1 & cpm_all_ycg1[,13]>1 & cpm_all_ycg1[,14]>1),]

rpkm_all_spike_ycs4_exp_all<-rpkm_all_spike_ycs4[which(rpkm_all_spike_brn1[,7]>1 & rpkm_all_spike_brn1[,8]>1 & rpkm_all_spike_brn1[,9]>1 & rpkm_all_spike_brn1[,10]>1 & rpkm_all_spike_brn1[,11]>1 & rpkm_all_spike_brn1[,12]>1 & rpkm_all_spike_brn1[,13]>1 & rpkm_all_spike_brn1[,14]>1&rpkm_all_spike_ycs4[,7]>1 & rpkm_all_spike_ycs4[,8]>1 & rpkm_all_spike_ycs4[,9]>1 & rpkm_all_spike_ycs4[,10]>1 & rpkm_all_spike_ycs4[,11]>1 & rpkm_all_spike_ycs4[,12]>1 & rpkm_all_spike_ycs4[,13]>1 & rpkm_all_spike_ycs4[,14]>1& rpkm_all_spike_ycs4[,15]>1 & rpkm_all_spike_ycs4[,16]>1 & rpkm_all_spike_ycs4[,17]>1&rpkm_all_spike_ycg1[,7]>1 & rpkm_all_spike_ycg1[,8]>1 & rpkm_all_spike_ycg1[,9]>1 & rpkm_all_spike_ycg1[,10]>1 & rpkm_all_spike_ycg1[,11]>1 & rpkm_all_spike_ycg1[,12]>1 & rpkm_all_spike_ycg1[,13]>1 & rpkm_all_spike_ycg1[,14]>1),]
rpkm_all_ycs4_exp_all<-rpkm_all_ycs4[which(rpkm_all_brn1[,7]>1 & rpkm_all_brn1[,8]>1 & rpkm_all_brn1[,9]>1 & rpkm_all_brn1[,10]>1 & rpkm_all_brn1[,11]>1 & rpkm_all_brn1[,12]>1 & rpkm_all_brn1[,13]>1 & rpkm_all_brn1[,14]>1&rpkm_all_ycs4[,7]>1 & rpkm_all_ycs4[,8]>1 & rpkm_all_ycs4[,9]>1 & rpkm_all_ycs4[,10]>1 & rpkm_all_ycs4[,11]>1 & rpkm_all_ycs4[,12]>1 & rpkm_all_ycs4[,13]>1 & rpkm_all_ycs4[,14]& rpkm_all_ycs4[,15]>1& rpkm_all_ycs4[,16]>1& rpkm_all_ycs4[,17]>1>1&rpkm_all_ycg1[,7]>1 & rpkm_all_ycg1[,8]>1 & rpkm_all_ycg1[,9]>1 & rpkm_all_ycg1[,10]>1 & rpkm_all_ycg1[,11]>1 & rpkm_all_ycg1[,12]>1 & rpkm_all_ycg1[,13]>1 & rpkm_all_ycg1[,14]>1),]
tpm_all_spike_ycs4_exp_all<-tpm_all_spike_ycs4[which(tpm_all_spike_brn1[,7]>1 & tpm_all_spike_brn1[,8]>1 & tpm_all_spike_brn1[,9]>1 & tpm_all_spike_brn1[,10]>1 & tpm_all_spike_brn1[,11]>1 & tpm_all_spike_brn1[,12]>1 & tpm_all_spike_brn1[,13]>1 & tpm_all_spike_brn1[,14]>1&tpm_all_spike_ycs4[,7]>1 & tpm_all_spike_ycs4[,8]>1 & tpm_all_spike_ycs4[,9]>1 & tpm_all_spike_ycs4[,10]>1 & tpm_all_spike_ycs4[,11]>1 & tpm_all_spike_ycs4[,12]>1 & tpm_all_spike_ycs4[,13]>1 & tpm_all_spike_ycs4[,14]>1& tpm_all_spike_ycs4[,15]>1& tpm_all_spike_ycs4[,16]>1& tpm_all_spike_ycs4[,17]>1&tpm_all_spike_ycg1[,7]>1 & tpm_all_spike_ycg1[,8]>1 & tpm_all_spike_ycg1[,9]>1 & tpm_all_spike_ycg1[,10]>1 & tpm_all_spike_ycg1[,11]>1 & tpm_all_spike_ycg1[,12]>1 & tpm_all_spike_ycg1[,13]>1 & tpm_all_spike_ycg1[,14]>1),]
tpm_all_ycs4_exp_all<-tpm_all_ycs4[which(tpm_all_brn1[,7]>1 & tpm_all_brn1[,8]>1 & tpm_all_brn1[,9]>1 & tpm_all_brn1[,10]>1 & tpm_all_brn1[,11]>1 & tpm_all_brn1[,12]>1 & tpm_all_brn1[,13]>1 & tpm_all_brn1[,14]>1&tpm_all_ycs4[,7]>1 & tpm_all_ycs4[,8]>1 & tpm_all_ycs4[,9]>1 & tpm_all_ycs4[,10]>1 & tpm_all_ycs4[,11]>1 & tpm_all_ycs4[,12]>1 & tpm_all_ycs4[,13]>1 & tpm_all_ycs4[,14]>1& tpm_all_ycs4[,15]>1& tpm_all_ycs4[,16]>1& tpm_all_ycs4[,14]>1&tpm_all_ycg1[,7]>1 & tpm_all_ycg1[,8]>1 & tpm_all_ycg1[,9]>1 & tpm_all_ycg1[,10]>1 & tpm_all_ycg1[,11]>1 & tpm_all_ycg1[,12]>1 & tpm_all_ycg1[,13]>1 & tpm_all_ycg1[,14]>1),]
cpm_all_spike_ycs4_exp_all<-cpm_all_spike_ycs4[which(cpm_all_spike_brn1[,7]>1 & cpm_all_spike_brn1[,8]>1 & cpm_all_spike_brn1[,9]>1 & cpm_all_spike_brn1[,10]>1 & cpm_all_spike_brn1[,11]>1 & cpm_all_spike_brn1[,12]>1 & cpm_all_spike_brn1[,13]>1 & cpm_all_spike_brn1[,14]>1&cpm_all_spike_ycs4[,7]>1 & cpm_all_spike_ycs4[,8]>1 & cpm_all_spike_ycs4[,9]>1 & cpm_all_spike_ycs4[,10]>1 & cpm_all_spike_ycs4[,11]>1 & cpm_all_spike_ycs4[,12]>1 & cpm_all_spike_ycs4[,13]>1 & cpm_all_spike_ycs4[,14]>1& cpm_all_spike_ycs4[,15]>1& cpm_all_spike_ycs4[,16]>1& cpm_all_spike_ycs4[,17]>1&cpm_all_spike_ycg1[,7]>1 & cpm_all_spike_ycg1[,8]>1 & cpm_all_spike_ycg1[,9]>1 & cpm_all_spike_ycg1[,10]>1 & cpm_all_spike_ycg1[,11]>1 & cpm_all_spike_ycg1[,12]>1 & cpm_all_spike_ycg1[,13]>1 & cpm_all_spike_ycg1[,14]>1),]
cpm_all_ycs4_exp_all<-cpm_all_ycs4[which(cpm_all_brn1[,7]>1 & cpm_all_brn1[,8]>1 & cpm_all_brn1[,9]>1 & cpm_all_brn1[,10]>1 & cpm_all_brn1[,11]>1 & cpm_all_brn1[,12]>1 & cpm_all_brn1[,13]>1 & cpm_all_brn1[,14]>1&cpm_all_ycs4[,7]>1 & cpm_all_ycs4[,8]>1 & cpm_all_ycs4[,9]>1 & cpm_all_ycs4[,10]>1 & cpm_all_ycs4[,11]>1 & cpm_all_ycs4[,12]>1 & cpm_all_ycs4[,13]>1 & cpm_all_ycs4[,14]>1& cpm_all_ycs4[,15]>1& cpm_all_ycs4[,16]>1& cpm_all_ycs4[,17]>1&cpm_all_ycg1[,7]>1 & cpm_all_ycg1[,8]>1 & cpm_all_ycg1[,9]>1 & cpm_all_ycg1[,10]>1 & cpm_all_ycg1[,11]>1 & cpm_all_ycg1[,12]>1 & cpm_all_ycg1[,13]>1 & cpm_all_ycg1[,14]>1),]

##Create means
cpm_all_brn1_exp_all_mean<-cbind(cpm_all_brn1_exp_all,cpm_all_brn1_exp_all[,7],apply(cpm_all_brn1_exp_all[,8:10],1,mean),cpm_all_brn1_exp_all[,11], apply(cpm_all_brn1_exp_all[,12:14],1,mean))
colnames(cpm_all_brn1_exp_all_mean)<-c(colnames(cpm_all_brn1_exp_all), 'brn1_wt_norap_mean', 'brn1_wt_rap_mean', 'brn1_brn1_norap_mean', 'brn1_brn1_rap_mean')
tpm_all_brn1_exp_all_mean<-cbind(tpm_all_brn1_exp_all,tpm_all_brn1_exp_all[,7],apply(tpm_all_brn1_exp_all[,8:10],1,mean),tpm_all_brn1_exp_all[,11], apply(tpm_all_brn1_exp_all[,12:14],1,mean))
colnames(tpm_all_brn1_exp_all_mean)<-c(colnames(tpm_all_brn1_exp_all), 'brn1_wt_norap_mean', 'brn1_wt_rap_mean', 'brn1_brn1_norap_mean', 'brn1_brn1_rap_mean')
rpkm_all_brn1_exp_all_mean<-cbind(rpkm_all_brn1_exp_all,rpkm_all_brn1_exp_all[,7],apply(rpkm_all_brn1_exp_all[,8:10],1,mean),rpkm_all_brn1_exp_all[,11], apply(rpkm_all_brn1_exp_all[,12:14],1,mean))
colnames(rpkm_all_brn1_exp_all_mean)<-c(colnames(rpkm_all_brn1_exp_all), 'brn1_wt_norap_mean', 'brn1_wt_rap_mean', 'brn1_brn1_norap_mean', 'brn1_brn1_rap_mean')

cpm_all_ycg1_exp_all_mean<-cbind(cpm_all_ycg1_exp_all,apply(cpm_all_ycg1_exp_all[,7:8],1,mean),apply(cpm_all_ycg1_exp_all[,9:10],1,mean),apply(cpm_all_ycg1_exp_all[,11:12],1,mean), apply(cpm_all_ycg1_exp_all[,13:14],1,mean))
colnames(cpm_all_ycg1_exp_all_mean)<-c(colnames(cpm_all_ycg1_exp_all), 'ycg1_wt_norap_mean', 'ycg1_wt_rap_mean', 'ycg1_ycg1_norap_mean', 'ycg1_ycg1_rap_mean')
tpm_all_ycg1_exp_all_mean<-cbind(tpm_all_ycg1_exp_all,apply(tpm_all_ycg1_exp_all[,7:8],1,mean),apply(tpm_all_ycg1_exp_all[,9:10],1,mean),apply(tpm_all_ycg1_exp_all[,11:12],1,mean), apply(tpm_all_ycg1_exp_all[,13:14],1,mean))
colnames(tpm_all_ycg1_exp_all_mean)<-c(colnames(tpm_all_ycg1_exp_all), 'ycg1_wt_norap_mean', 'ycg1_wt_rap_mean', 'ycg1_ycg1_norap_mean', 'ycg1_ycg1_rap_mean')
rpkm_all_ycg1_exp_all_mean<-cbind(rpkm_all_ycg1_exp_all,apply(rpkm_all_ycg1_exp_all[,7:8],1,mean),apply(rpkm_all_ycg1_exp_all[,9:10],1,mean),apply(rpkm_all_ycg1_exp_all[,11:12],1,mean), apply(rpkm_all_ycg1_exp_all[,13:14],1,mean))
colnames(rpkm_all_ycg1_exp_all_mean)<-c(colnames(rpkm_all_ycg1_exp_all), 'ycg1_wt_norap_mean', 'ycg1_wt_rap_mean', 'ycg1_ycg1_norap_mean', 'ycg1_ycg1_rap_mean')

cpm_all_ycs4_exp_all_mean<-cbind(cpm_all_ycs4_exp_all,apply(cpm_all_ycs4_exp_all[,7:8],1,mean),apply(cpm_all_ycs4_exp_all[,9:11],1,mean),apply(cpm_all_ycs4_exp_all[,12:14],1,mean), apply(cpm_all_ycs4_exp_all[,15:17],1,mean))
colnames(cpm_all_ycs4_exp_all_mean)<-c(colnames(cpm_all_ycs4_exp_all), 'ycs4_wt_norap_mean', 'ycs4_wt_rap_mean', 'ycs4_ycs4_norap_mean', 'ycs4_ycs4_rap_mean')
tpm_all_ycs4_exp_all_mean<-cbind(tpm_all_ycs4_exp_all,apply(tpm_all_ycs4_exp_all[,7:8],1,mean),apply(tpm_all_ycs4_exp_all[,9:11],1,mean),apply(tpm_all_ycs4_exp_all[,12:14],1,mean), apply(tpm_all_ycs4_exp_all[,15:17],1,mean))
colnames(tpm_all_ycs4_exp_all_mean)<-c(colnames(tpm_all_ycs4_exp_all), 'ycs4_wt_norap_mean', 'ycs4_wt_rap_mean', 'ycs4_ycs4_norap_mean', 'ycs4_ycs4_rap_mean')
rpkm_all_ycs4_exp_all_mean<-cbind(rpkm_all_ycs4_exp_all,apply(rpkm_all_ycs4_exp_all[,7:8],1,mean),apply(rpkm_all_ycs4_exp_all[,9:11],1,mean),apply(rpkm_all_ycs4_exp_all[,12:14],1,mean), apply(rpkm_all_ycs4_exp_all[,15:17],1,mean))
colnames(rpkm_all_ycs4_exp_all_mean)<-c(colnames(rpkm_all_ycs4_exp_all), 'ycs4_wt_norap_mean', 'ycs4_wt_rap_mean', 'ycs4_ycs4_norap_mean', 'ycs4_ycs4_rap_mean')

 
cpm_all_brn1_exp_mean<-cbind(cpm_all_brn1_exp,cpm_all_brn1_exp[,7],apply(cpm_all_brn1_exp[,8:10],1,mean),cpm_all_brn1_exp[,11], apply(cpm_all_brn1_exp[,12:14],1,mean))
colnames(cpm_all_brn1_exp_mean)<-c(colnames(cpm_all_brn1_exp), 'brn1_wt_norap_mean', 'brn1_wt_rap_mean', 'brn1_brn1_norap_mean', 'brn1_brn1_rap_mean')
tpm_all_brn1_exp_mean<-cbind(tpm_all_brn1_exp,tpm_all_brn1_exp[,7],apply(tpm_all_brn1_exp[,8:10],1,mean),tpm_all_brn1_exp[,11], apply(tpm_all_brn1_exp[,12:14],1,mean))
colnames(tpm_all_brn1_exp_mean)<-c(colnames(tpm_all_brn1_exp), 'brn1_wt_norap_mean', 'brn1_wt_rap_mean', 'brn1_brn1_norap_mean', 'brn1_brn1_rap_mean')
rpkm_all_brn1_exp_mean<-cbind(rpkm_all_brn1_exp,rpkm_all_brn1_exp[,7],apply(rpkm_all_brn1_exp[,8:10],1,mean),rpkm_all_brn1_exp[,11], apply(rpkm_all_brn1_exp[,12:14],1,mean))
colnames(rpkm_all_brn1_exp_mean)<-c(colnames(rpkm_all_brn1_exp), 'brn1_wt_norap_mean', 'brn1_wt_rap_mean', 'brn1_brn1_norap_mean', 'brn1_brn1_rap_mean')

cpm_all_ycg1_exp_mean<-cbind(cpm_all_ycg1_exp,apply(cpm_all_ycg1_exp[,7:8],1,mean),apply(cpm_all_ycg1_exp[,9:10],1,mean),apply(cpm_all_ycg1_exp[,11:12],1,mean), apply(cpm_all_ycg1_exp[,13:14],1,mean))
colnames(cpm_all_ycg1_exp_mean)<-c(colnames(cpm_all_ycg1_exp), 'ycg1_wt_norap_mean', 'ycg1_wt_rap_mean', 'ycg1_ycg1_norap_mean', 'ycg1_ycg1_rap_mean')
tpm_all_ycg1_exp_mean<-cbind(tpm_all_ycg1_exp,apply(tpm_all_ycg1_exp[,7:8],1,mean),apply(tpm_all_ycg1_exp[,9:10],1,mean),apply(tpm_all_ycg1_exp[,11:12],1,mean), apply(tpm_all_ycg1_exp[,13:14],1,mean))
colnames(tpm_all_ycg1_exp_mean)<-c(colnames(tpm_all_ycg1_exp), 'ycg1_wt_norap_mean', 'ycg1_wt_rap_mean', 'ycg1_ycg1_norap_mean', 'ycg1_ycg1_rap_mean')
rpkm_all_ycg1_exp_mean<-cbind(rpkm_all_ycg1_exp,apply(rpkm_all_ycg1_exp[,7:8],1,mean),apply(rpkm_all_ycg1_exp[,9:10],1,mean),apply(rpkm_all_ycg1_exp[,11:12],1,mean), apply(rpkm_all_ycg1_exp[,13:14],1,mean))
colnames(rpkm_all_ycg1_exp_mean)<-c(colnames(rpkm_all_ycg1_exp), 'ycg1_wt_norap_mean', 'ycg1_wt_rap_mean', 'ycg1_ycg1_norap_mean', 'ycg1_ycg1_rap_mean')

cpm_all_ycs4_exp_mean<-cbind(cpm_all_ycs4_exp,apply(cpm_all_ycs4_exp[,7:8],1,mean),apply(cpm_all_ycs4_exp[,9:11],1,mean),apply(cpm_all_ycs4_exp[,12:14],1,mean), apply(cpm_all_ycs4_exp[,15:17],1,mean))
colnames(cpm_all_ycs4_exp_mean)<-c(colnames(cpm_all_ycs4_exp), 'ycs4_wt_norap_mean', 'ycs4_wt_rap_mean', 'ycs4_ycs4_norap_mean', 'ycs4_ycs4_rap_mean')
tpm_all_ycs4_exp_mean<-cbind(tpm_all_ycs4_exp,apply(tpm_all_ycs4_exp[,7:8],1,mean),apply(tpm_all_ycs4_exp[,9:11],1,mean),apply(tpm_all_ycs4_exp[,12:14],1,mean), apply(tpm_all_ycs4_exp[,15:17],1,mean))
colnames(tpm_all_ycs4_exp_mean)<-c(colnames(tpm_all_ycs4_exp), 'ycs4_wt_norap_mean', 'ycs4_wt_rap_mean', 'ycs4_ycs4_norap_mean', 'ycs4_ycs4_rap_mean')
rpkm_all_ycs4_exp_mean<-cbind(rpkm_all_ycs4_exp,apply(rpkm_all_ycs4_exp[,7:8],1,mean),apply(rpkm_all_ycs4_exp[,9:11],1,mean),apply(rpkm_all_ycs4_exp[,12:14],1,mean), apply(rpkm_all_ycs4_exp[,15:17],1,mean))
colnames(rpkm_all_ycs4_exp_mean)<-c(colnames(rpkm_all_ycs4_exp), 'ycs4_wt_norap_mean', 'ycs4_wt_rap_mean', 'ycs4_ycs4_norap_mean', 'ycs4_ycs4_rap_mean')


cpm_all_brn1_mean<-cbind(cpm_all_brn1,cpm_all_brn1[,7],apply(cpm_all_brn1[,8:10],1,mean),cpm_all_brn1[,11], apply(cpm_all_brn1[,12:14],1,mean))
colnames(cpm_all_brn1_mean)<-c(colnames(cpm_all_brn1), 'brn1_wt_norap_mean', 'brn1_wt_rap_mean', 'brn1_brn1_norap_mean', 'brn1_brn1_rap_mean')
tpm_all_brn1_mean<-cbind(tpm_all_brn1,tpm_all_brn1[,7],apply(tpm_all_brn1[,8:10],1,mean),tpm_all_brn1[,11], apply(tpm_all_brn1[,12:14],1,mean))
colnames(tpm_all_brn1_mean)<-c(colnames(tpm_all_brn1), 'brn1_wt_norap_mean', 'brn1_wt_rap_mean', 'brn1_brn1_norap_mean', 'brn1_brn1_rap_mean')
rpkm_all_brn1_mean<-cbind(rpkm_all_brn1,rpkm_all_brn1[,7],apply(rpkm_all_brn1[,8:10],1,mean),rpkm_all_brn1[,11], apply(rpkm_all_brn1[,12:14],1,mean))
colnames(rpkm_all_brn1_mean)<-c(colnames(rpkm_all_brn1), 'brn1_wt_norap_mean', 'brn1_wt_rap_mean', 'brn1_brn1_norap_mean', 'brn1_brn1_rap_mean')

cpm_all_ycg1_mean<-cbind(cpm_all_ycg1,apply(cpm_all_ycg1[,7:8],1,mean),apply(cpm_all_ycg1[,9:10],1,mean),apply(cpm_all_ycg1[,11:12],1,mean), apply(cpm_all_ycg1[,13:14],1,mean))
colnames(cpm_all_ycg1_mean)<-c(colnames(cpm_all_ycg1), 'ycg1_wt_norap_mean', 'ycg1_wt_rap_mean', 'ycg1_ycg1_norap_mean', 'ycg1_ycg1_rap_mean')
tpm_all_ycg1_mean<-cbind(tpm_all_ycg1,apply(tpm_all_ycg1[,7:8],1,mean),apply(tpm_all_ycg1[,9:10],1,mean),apply(tpm_all_ycg1[,11:12],1,mean), apply(tpm_all_ycg1[,13:14],1,mean))
colnames(tpm_all_ycg1_mean)<-c(colnames(tpm_all_ycg1), 'ycg1_wt_norap_mean', 'ycg1_wt_rap_mean', 'ycg1_ycg1_norap_mean', 'ycg1_ycg1_rap_mean')
rpkm_all_ycg1_mean<-cbind(rpkm_all_ycg1,apply(rpkm_all_ycg1[,7:8],1,mean),apply(rpkm_all_ycg1[,9:10],1,mean),apply(rpkm_all_ycg1[,11:12],1,mean), apply(rpkm_all_ycg1[,13:14],1,mean))
colnames(rpkm_all_ycg1_mean)<-c(colnames(rpkm_all_ycg1), 'ycg1_wt_norap_mean', 'ycg1_wt_rap_mean', 'ycg1_ycg1_norap_mean', 'ycg1_ycg1_rap_mean')

cpm_all_ycs4_mean<-cbind(cpm_all_ycs4,apply(cpm_all_ycs4[,7:8],1,mean),apply(cpm_all_ycs4[,9:11],1,mean),apply(cpm_all_ycs4[,12:14],1,mean), apply(cpm_all_ycs4[,15:17],1,mean))
colnames(cpm_all_ycs4_mean)<-c(colnames(cpm_all_ycs4), 'ycs4_wt_norap_mean', 'ycs4_wt_rap_mean', 'ycs4_ycs4_norap_mean', 'ycs4_ycs4_rap_mean')
tpm_all_ycs4_mean<-cbind(tpm_all_ycs4,apply(tpm_all_ycs4[,7:8],1,mean),apply(tpm_all_ycs4[,9:11],1,mean),apply(tpm_all_ycs4[,12:14],1,mean), apply(tpm_all_ycs4[,15:17],1,mean))
colnames(tpm_all_ycs4_mean)<-c(colnames(tpm_all_ycs4), 'ycs4_wt_norap_mean', 'ycs4_wt_rap_mean', 'ycs4_ycs4_norap_mean', 'ycs4_ycs4_rap_mean')
rpkm_all_ycs4_mean<-cbind(rpkm_all_ycs4,apply(rpkm_all_ycs4[,7:8],1,mean),apply(rpkm_all_ycs4[,9:11],1,mean),apply(rpkm_all_ycs4[,12:14],1,mean), apply(rpkm_all_ycs4[,15:17],1,mean))
colnames(rpkm_all_ycs4_mean)<-c(colnames(rpkm_all_ycs4), 'ycs4_wt_norap_mean', 'ycs4_wt_rap_mean', 'ycs4_ycs4_norap_mean', 'ycs4_ycs4_rap_mean')


cpm_all_spike_brn1_exp_all_mean<-cbind(cpm_all_spike_brn1_exp_all,cpm_all_spike_brn1_exp_all[,7],apply(cpm_all_spike_brn1_exp_all[,8:10],1,mean),cpm_all_spike_brn1_exp_all[,11], apply(cpm_all_spike_brn1_exp_all[,12:14],1,mean))
colnames(cpm_all_spike_brn1_exp_all_mean)<-c(colnames(cpm_all_spike_brn1_exp_all), 'brn1_wt_norap_mean', 'brn1_wt_rap_mean', 'brn1_brn1_norap_mean', 'brn1_brn1_rap_mean')
tpm_all_spike_brn1_exp_all_mean<-cbind(tpm_all_spike_brn1_exp_all,tpm_all_spike_brn1_exp_all[,7],apply(tpm_all_spike_brn1_exp_all[,8:10],1,mean),tpm_all_spike_brn1_exp_all[,11], apply(tpm_all_spike_brn1_exp_all[,12:14],1,mean))
colnames(tpm_all_spike_brn1_exp_all_mean)<-c(colnames(tpm_all_spike_brn1_exp_all), 'brn1_wt_norap_mean', 'brn1_wt_rap_mean', 'brn1_brn1_norap_mean', 'brn1_brn1_rap_mean')
rpkm_all_spike_brn1_exp_all_mean<-cbind(rpkm_all_spike_brn1_exp_all,rpkm_all_spike_brn1_exp_all[,7],apply(rpkm_all_spike_brn1_exp_all[,8:10],1,mean),rpkm_all_spike_brn1_exp_all[,11], apply(rpkm_all_spike_brn1_exp_all[,12:14],1,mean))
colnames(rpkm_all_spike_brn1_exp_all_mean)<-c(colnames(rpkm_all_spike_brn1_exp_all), 'brn1_wt_norap_mean', 'brn1_wt_rap_mean', 'brn1_brn1_norap_mean', 'brn1_brn1_rap_mean')

cpm_all_spike_ycg1_exp_all_mean<-cbind(cpm_all_spike_ycg1_exp_all,apply(cpm_all_spike_ycg1_exp_all[,7:8],1,mean),apply(cpm_all_spike_ycg1_exp_all[,9:10],1,mean),apply(cpm_all_spike_ycg1_exp_all[,11:12],1,mean), apply(cpm_all_spike_ycg1_exp_all[,13:14],1,mean))
colnames(cpm_all_spike_ycg1_exp_all_mean)<-c(colnames(cpm_all_spike_ycg1_exp_all), 'ycg1_wt_norap_mean', 'ycg1_wt_rap_mean', 'ycg1_ycg1_norap_mean', 'ycg1_ycg1_rap_mean')
tpm_all_spike_ycg1_exp_all_mean<-cbind(tpm_all_spike_ycg1_exp_all,apply(tpm_all_spike_ycg1_exp_all[,7:8],1,mean),apply(tpm_all_spike_ycg1_exp_all[,9:10],1,mean),apply(tpm_all_spike_ycg1_exp_all[,11:12],1,mean), apply(tpm_all_spike_ycg1_exp_all[,13:14],1,mean))
colnames(tpm_all_spike_ycg1_exp_all_mean)<-c(colnames(tpm_all_spike_ycg1_exp_all), 'ycg1_wt_norap_mean', 'ycg1_wt_rap_mean', 'ycg1_ycg1_norap_mean', 'ycg1_ycg1_rap_mean')
rpkm_all_spike_ycg1_exp_all_mean<-cbind(rpkm_all_spike_ycg1_exp_all,apply(rpkm_all_spike_ycg1_exp_all[,7:8],1,mean),apply(rpkm_all_spike_ycg1_exp_all[,9:10],1,mean),apply(rpkm_all_spike_ycg1_exp_all[,11:12],1,mean), apply(rpkm_all_spike_ycg1_exp_all[,13:14],1,mean))
colnames(rpkm_all_spike_ycg1_exp_all_mean)<-c(colnames(rpkm_all_spike_ycg1_exp_all), 'ycg1_wt_norap_mean', 'ycg1_wt_rap_mean', 'ycg1_ycg1_norap_mean', 'ycg1_ycg1_rap_mean')

cpm_all_spike_ycs4_exp_all_mean<-cbind(cpm_all_spike_ycs4_exp_all,apply(cpm_all_spike_ycs4_exp_all[,7:8],1,mean),apply(cpm_all_spike_ycs4_exp_all[,9:11],1,mean),apply(cpm_all_spike_ycs4_exp_all[,12:14],1,mean), apply(cpm_all_spike_ycs4_exp_all[,15:17],1,mean))
colnames(cpm_all_spike_ycs4_exp_all_mean)<-c(colnames(cpm_all_spike_ycs4_exp_all), 'ycs4_wt_norap_mean', 'ycs4_wt_rap_mean', 'ycs4_ycs4_norap_mean', 'ycs4_ycs4_rap_mean')
tpm_all_spike_ycs4_exp_all_mean<-cbind(tpm_all_spike_ycs4_exp_all,apply(tpm_all_spike_ycs4_exp_all[,7:8],1,mean),apply(tpm_all_spike_ycs4_exp_all[,9:11],1,mean),apply(tpm_all_spike_ycs4_exp_all[,12:14],1,mean), apply(tpm_all_spike_ycs4_exp_all[,15:17],1,mean))
colnames(tpm_all_spike_ycs4_exp_all_mean)<-c(colnames(tpm_all_spike_ycs4_exp_all), 'ycs4_wt_norap_mean', 'ycs4_wt_rap_mean', 'ycs4_ycs4_norap_mean', 'ycs4_ycs4_rap_mean')
rpkm_all_spike_ycs4_exp_all_mean<-cbind(rpkm_all_spike_ycs4_exp_all,apply(rpkm_all_spike_ycs4_exp_all[,7:8],1,mean),apply(rpkm_all_spike_ycs4_exp_all[,9:11],1,mean),apply(rpkm_all_spike_ycs4_exp_all[,12:14],1,mean), apply(rpkm_all_spike_ycs4_exp_all[,15:17],1,mean))
colnames(rpkm_all_spike_ycs4_exp_all_mean)<-c(colnames(rpkm_all_spike_ycs4_exp_all), 'ycs4_wt_norap_mean', 'ycs4_wt_rap_mean', 'ycs4_ycs4_norap_mean', 'ycs4_ycs4_rap_mean')


cpm_all_spike_brn1_exp_mean<-cbind(cpm_all_spike_brn1_exp,cpm_all_spike_brn1_exp[,7],apply(cpm_all_spike_brn1_exp[,8:10],1,mean),cpm_all_spike_brn1_exp[,11], apply(cpm_all_spike_brn1_exp[,12:14],1,mean))
colnames(cpm_all_spike_brn1_exp_mean)<-c(colnames(cpm_all_spike_brn1_exp), 'brn1_wt_norap_mean', 'brn1_wt_rap_mean', 'brn1_brn1_norap_mean', 'brn1_brn1_rap_mean')
tpm_all_spike_brn1_exp_mean<-cbind(tpm_all_spike_brn1_exp,tpm_all_spike_brn1_exp[,7],apply(tpm_all_spike_brn1_exp[,8:10],1,mean),tpm_all_spike_brn1_exp[,11], apply(tpm_all_spike_brn1_exp[,12:14],1,mean))
colnames(tpm_all_spike_brn1_exp_mean)<-c(colnames(tpm_all_spike_brn1_exp), 'brn1_wt_norap_mean', 'brn1_wt_rap_mean', 'brn1_brn1_norap_mean', 'brn1_brn1_rap_mean')
rpkm_all_spike_brn1_exp_mean<-cbind(rpkm_all_spike_brn1_exp,rpkm_all_spike_brn1_exp[,7],apply(rpkm_all_spike_brn1_exp[,8:10],1,mean),rpkm_all_spike_brn1_exp[,11], apply(rpkm_all_spike_brn1_exp[,12:14],1,mean))
colnames(rpkm_all_spike_brn1_exp_mean)<-c(colnames(rpkm_all_spike_brn1_exp), 'brn1_wt_norap_mean', 'brn1_wt_rap_mean', 'brn1_brn1_norap_mean', 'brn1_brn1_rap_mean')

cpm_all_spike_ycg1_exp_mean<-cbind(cpm_all_spike_ycg1_exp,apply(cpm_all_spike_ycg1_exp[,7:8],1,mean),apply(cpm_all_spike_ycg1_exp[,9:10],1,mean),apply(cpm_all_spike_ycg1_exp[,11:12],1,mean), apply(cpm_all_spike_ycg1_exp[,13:14],1,mean))
colnames(cpm_all_spike_ycg1_exp_mean)<-c(colnames(cpm_all_spike_ycg1_exp), 'ycg1_wt_norap_mean', 'ycg1_wt_rap_mean', 'ycg1_ycg1_norap_mean', 'ycg1_ycg1_rap_mean')
tpm_all_spike_ycg1_exp_mean<-cbind(tpm_all_spike_ycg1_exp,apply(tpm_all_spike_ycg1_exp[,7:8],1,mean),apply(tpm_all_spike_ycg1_exp[,9:10],1,mean),apply(tpm_all_spike_ycg1_exp[,11:12],1,mean), apply(tpm_all_spike_ycg1_exp[,13:14],1,mean))
colnames(tpm_all_spike_ycg1_exp_mean)<-c(colnames(tpm_all_spike_ycg1_exp), 'ycg1_wt_norap_mean', 'ycg1_wt_rap_mean', 'ycg1_ycg1_norap_mean', 'ycg1_ycg1_rap_mean')
rpkm_all_spike_ycg1_exp_mean<-cbind(rpkm_all_spike_ycg1_exp,apply(rpkm_all_spike_ycg1_exp[,7:8],1,mean),apply(rpkm_all_spike_ycg1_exp[,9:10],1,mean),apply(rpkm_all_spike_ycg1_exp[,11:12],1,mean), apply(rpkm_all_spike_ycg1_exp[,13:14],1,mean))
colnames(rpkm_all_spike_ycg1_exp_mean)<-c(colnames(rpkm_all_spike_ycg1_exp), 'ycg1_wt_norap_mean', 'ycg1_wt_rap_mean', 'ycg1_ycg1_norap_mean', 'ycg1_ycg1_rap_mean')

cpm_all_spike_ycs4_exp_mean<-cbind(cpm_all_spike_ycs4_exp,apply(cpm_all_spike_ycs4_exp[,7:8],1,mean),apply(cpm_all_spike_ycs4_exp[,9:11],1,mean),apply(cpm_all_spike_ycs4_exp[,12:14],1,mean), apply(cpm_all_spike_ycs4_exp[,15:17],1,mean))
colnames(cpm_all_spike_ycs4_exp_mean)<-c(colnames(cpm_all_spike_ycs4_exp), 'ycs4_wt_norap_mean', 'ycs4_wt_rap_mean', 'ycs4_ycs4_norap_mean', 'ycs4_ycs4_rap_mean')
tpm_all_spike_ycs4_exp_mean<-cbind(tpm_all_spike_ycs4_exp,apply(tpm_all_spike_ycs4_exp[,7:8],1,mean),apply(tpm_all_spike_ycs4_exp[,9:11],1,mean),apply(tpm_all_spike_ycs4_exp[,12:14],1,mean), apply(tpm_all_spike_ycs4_exp[,15:17],1,mean))
colnames(tpm_all_spike_ycs4_exp_mean)<-c(colnames(tpm_all_spike_ycs4_exp), 'ycs4_wt_norap_mean', 'ycs4_wt_rap_mean', 'ycs4_ycs4_norap_mean', 'ycs4_ycs4_rap_mean')
rpkm_all_spike_ycs4_exp_mean<-cbind(rpkm_all_spike_ycs4_exp,apply(rpkm_all_spike_ycs4_exp[,7:8],1,mean),apply(rpkm_all_spike_ycs4_exp[,9:11],1,mean),apply(rpkm_all_spike_ycs4_exp[,12:14],1,mean), apply(rpkm_all_spike_ycs4_exp[,15:17],1,mean))
colnames(rpkm_all_spike_ycs4_exp_mean)<-c(colnames(rpkm_all_spike_ycs4_exp), 'ycs4_wt_norap_mean', 'ycs4_wt_rap_mean', 'ycs4_ycs4_norap_mean', 'ycs4_ycs4_rap_mean')


cpm_all_spike_brn1_mean<-cbind(cpm_all_spike_brn1,cpm_all_spike_brn1[,7],apply(cpm_all_spike_brn1[,8:10],1,mean),cpm_all_spike_brn1[,11], apply(cpm_all_spike_brn1[,12:14],1,mean))
colnames(cpm_all_spike_brn1_mean)<-c(colnames(cpm_all_spike_brn1), 'brn1_wt_norap_mean', 'brn1_wt_rap_mean', 'brn1_brn1_norap_mean', 'brn1_brn1_rap_mean')
tpm_all_spike_brn1_mean<-cbind(tpm_all_spike_brn1,tpm_all_spike_brn1[,7],apply(tpm_all_spike_brn1[,8:10],1,mean),tpm_all_spike_brn1[,11], apply(tpm_all_spike_brn1[,12:14],1,mean))
colnames(tpm_all_spike_brn1_mean)<-c(colnames(tpm_all_spike_brn1), 'brn1_wt_norap_mean', 'brn1_wt_rap_mean', 'brn1_brn1_norap_mean', 'brn1_brn1_rap_mean')
rpkm_all_spike_brn1_mean<-cbind(rpkm_all_spike_brn1,rpkm_all_spike_brn1[,7],apply(rpkm_all_spike_brn1[,8:10],1,mean),rpkm_all_spike_brn1[,11], apply(rpkm_all_spike_brn1[,12:14],1,mean))
colnames(rpkm_all_spike_brn1_mean)<-c(colnames(rpkm_all_spike_brn1), 'brn1_wt_norap_mean', 'brn1_wt_rap_mean', 'brn1_brn1_norap_mean', 'brn1_brn1_rap_mean')

cpm_all_spike_ycg1_mean<-cbind(cpm_all_spike_ycg1,apply(cpm_all_spike_ycg1[,7:8],1,mean),apply(cpm_all_spike_ycg1[,9:10],1,mean),apply(cpm_all_spike_ycg1[,11:12],1,mean), apply(cpm_all_spike_ycg1[,13:14],1,mean))
colnames(cpm_all_spike_ycg1_mean)<-c(colnames(cpm_all_spike_ycg1), 'ycg1_wt_norap_mean', 'ycg1_wt_rap_mean', 'ycg1_ycg1_norap_mean', 'ycg1_ycg1_rap_mean')
tpm_all_spike_ycg1_mean<-cbind(tpm_all_spike_ycg1,apply(tpm_all_spike_ycg1[,7:8],1,mean),apply(tpm_all_spike_ycg1[,9:10],1,mean),apply(tpm_all_spike_ycg1[,11:12],1,mean), apply(tpm_all_spike_ycg1[,13:14],1,mean))
colnames(tpm_all_spike_ycg1_mean)<-c(colnames(tpm_all_spike_ycg1), 'ycg1_wt_norap_mean', 'ycg1_wt_rap_mean', 'ycg1_ycg1_norap_mean', 'ycg1_ycg1_rap_mean')
rpkm_all_spike_ycg1_mean<-cbind(rpkm_all_spike_ycg1,apply(rpkm_all_spike_ycg1[,7:8],1,mean),apply(rpkm_all_spike_ycg1[,9:10],1,mean),apply(rpkm_all_spike_ycg1[,11:12],1,mean), apply(rpkm_all_spike_ycg1[,13:14],1,mean))
colnames(rpkm_all_spike_ycg1_mean)<-c(colnames(rpkm_all_spike_ycg1), 'ycg1_wt_norap_mean', 'ycg1_wt_rap_mean', 'ycg1_ycg1_norap_mean', 'ycg1_ycg1_rap_mean')

cpm_all_spike_ycs4_mean<-cbind(cpm_all_spike_ycs4,apply(cpm_all_spike_ycs4[,7:8],1,mean),apply(cpm_all_spike_ycs4[,9:11],1,mean),apply(cpm_all_spike_ycs4[,12:14],1,mean), apply(cpm_all_spike_ycs4[,15:17],1,mean))
colnames(cpm_all_spike_ycs4_mean)<-c(colnames(cpm_all_spike_ycs4), 'ycs4_wt_norap_mean', 'ycs4_wt_rap_mean', 'ycs4_ycs4_norap_mean', 'ycs4_ycs4_rap_mean')
tpm_all_spike_ycs4_mean<-cbind(tpm_all_spike_ycs4,apply(tpm_all_spike_ycs4[,7:8],1,mean),apply(tpm_all_spike_ycs4[,9:11],1,mean),apply(tpm_all_spike_ycs4[,12:14],1,mean), apply(tpm_all_spike_ycs4[,15:17],1,mean))
colnames(tpm_all_spike_ycs4_mean)<-c(colnames(tpm_all_spike_ycs4), 'ycs4_wt_norap_mean', 'ycs4_wt_rap_mean', 'ycs4_ycs4_norap_mean', 'ycs4_ycs4_rap_mean')
rpkm_all_spike_ycs4_mean<-cbind(rpkm_all_spike_ycs4,apply(rpkm_all_spike_ycs4[,7:8],1,mean),apply(rpkm_all_spike_ycs4[,9:11],1,mean),apply(rpkm_all_spike_ycs4[,12:14],1,mean), apply(rpkm_all_spike_ycs4[,15:17],1,mean))
colnames(rpkm_all_spike_ycs4_mean)<-c(colnames(rpkm_all_spike_ycs4), 'ycs4_wt_norap_mean', 'ycs4_wt_rap_mean', 'ycs4_ycs4_norap_mean', 'ycs4_ycs4_rap_mean')


#Export all the raw datasets.
path='~/Google_Drive/Lab/Data_Analysis/RNA/all_datasets/'
write.csv(rpkm_all_spike_ycg1_mean, paste0(path, '2018_03_22_rpkm_spike_ycg12.csv'))
write.csv(rpkm_all_ycg1_mean, paste0(path, '2018_03_22_rpkm_exp_ycg12.csv'))
write.csv(tpm_all_spike_ycg1_mean, paste0(path, '2018_03_22_tpm_spike_ycg12.csv'))
write.csv(tpm_all_ycg1_mean, paste0(path, '2018_03_22_tpm_exp_ycg12.csv'))
write.csv(cpm_all_spike_ycg1_mean, paste0(path, '2018_03_22_cpm_spike_ycg12.csv'))
write.csv(cpm_all_ycg1_mean, paste0(path, '2018_03_22_cpm_exp_ycg12.csv'))

write.csv(rpkm_all_spike_ycs4_mean, paste0(path, '2018_03_22_rpkm_spike_ycs42.csv'))
write.csv(rpkm_all_ycs4_mean, paste0(path, '2018_03_22_rpkm_exp_ycs42.csv'))
write.csv(tpm_all_spike_ycs4_mean, paste0(path, '2018_03_22_tpm_spike_ycs42.csv'))
write.csv(tpm_all_ycs4_mean, paste0(path, '2018_03_22_tpm_exp_ycs42.csv'))
write.csv(cpm_all_spike_ycs4_mean, paste0(path, '2018_03_22_cpm_spike_ycs42.csv'))
write.csv(cpm_all_ycs4_mean, paste0(path, '2018_03_22_cpm_exp_ycs42.csv'))

write.csv(rpkm_all_spike_brn1_mean, paste0(path, '2018_03_22_rpkm_spike_brn1aaway.csv'))
write.csv(rpkm_all_brn1_mean, paste0(path, '2018_03_22_rpkm_exp_brn1aaway.csv'))
write.csv(tpm_all_spike_brn1_mean, paste0(path, '2018_03_22_tpm_spike_brn1aaway.csv'))
write.csv(tpm_all_brn1_mean, paste0(path, '2018_03_22_tpm_exp_brn1aaway.csv'))
write.csv(cpm_all_spike_brn1_mean, paste0(path, '2018_03_22_cpm_spike_brn1aaway.csv'))
write.csv(cpm_all_brn1_mean, paste0(path, '2018_03_22_cpm_exp_brn1aaway.csv'))

#Export all the expressed only (per experiment) datasets.
write.csv(rpkm_all_spike_ycg1_exp_mean, paste0(path, '2018_03_22_rpkm_spike_ycg12_expressed.csv'))
write.csv(rpkm_all_ycg1_exp_mean, paste0(path, '2018_03_22_rpkm_exp_ycg12_expressed.csv'))
write.csv(tpm_all_spike_ycg1_exp_mean, paste0(path, '2018_03_22_tpm_spike_ycg12_expressed.csv'))
write.csv(tpm_all_ycg1_exp_mean, paste0(path, '2018_03_22_tpm_exp_ycg12_expressed.csv'))
write.csv(cpm_all_spike_ycg1_exp_mean, paste0(path, '2018_03_22_cpm_spike_ycg12_expressed.csv'))
write.csv(cpm_all_ycg1_exp_mean, paste0(path, '2018_03_22_cpm_exp_ycg12_expressed.csv'))

write.csv(rpkm_all_spike_ycs4_exp_mean, paste0(path, '2018_03_22_rpkm_spike_ycs42_expressed.csv'))
write.csv(rpkm_all_ycs4_exp_mean, paste0(path, '2018_03_22_rpkm_exp_ycs42_expressed.csv'))
write.csv(tpm_all_spike_ycs4_exp_mean, paste0(path, '2018_03_22_tpm_spike_ycs42_expressed.csv'))
write.csv(tpm_all_ycs4_exp_mean, paste0(path, '2018_03_22_tpm_exp_ycs42_expressed.csv'))
write.csv(cpm_all_spike_ycs4_exp_mean, paste0(path, '2018_03_22_cpm_spike_ycs42_expressed.csv'))
write.csv(cpm_all_ycs4_exp_mean, paste0(path, '2018_03_22_cpm_exp_ycs42_expressed.csv'))

write.csv(rpkm_all_spike_brn1_exp_mean, paste0(path, '2018_03_22_rpkm_spike_brn1aaway_expressed.csv'))
write.csv(rpkm_all_brn1_exp_mean, paste0(path, '2018_03_22_rpkm_exp_brn1aaway_expressed.csv'))
write.csv(tpm_all_spike_brn1_exp_mean, paste0(path, '2018_03_22_tpm_spike_brn1aaway_expressed.csv'))
write.csv(tpm_all_brn1_exp_mean, paste0(path, '2018_03_22_tpm_exp_brn1aaway_expressed.csv'))
write.csv(cpm_all_spike_brn1_exp_mean, paste0(path, '2018_03_22_cpm_spike_brn1aaway_expressed.csv'))
write.csv(cpm_all_brn1_exp_mean, paste0(path, '2018_03_22_cpm_exp_brn1aaway_expressed.csv'))

#Export all the expressed in all datasets only.
write.csv(rpkm_all_spike_ycg1_exp_all_mean, paste0(path, '2018_03_22_rpkm_spike_ycg12_expressed.in.all.csv'))
write.csv(rpkm_all_ycg1_exp_all_mean, paste0(path, '2018_03_22_rpkm_exp_ycg12_expressed.in.all.csv'))
write.csv(tpm_all_spike_ycg1_exp_all_mean, paste0(path, '2018_03_22_tpm_spike_ycg12_expressed.in.all.csv'))
write.csv(tpm_all_ycg1_exp_all_mean, paste0(path, '2018_03_22_tpm_exp_ycg12_expressed.in.all.csv'))
write.csv(cpm_all_spike_ycg1_exp_all_mean, paste0(path, '2018_03_22_cpm_spike_ycg12_expressed.in.all.csv'))
write.csv(cpm_all_ycg1_exp_all_mean, paste0(path, '2018_03_22_cpm_exp_ycg12_expressed.in.all.csv'))

write.csv(rpkm_all_spike_ycs4_exp_all_mean, paste0(path, '2018_03_22_rpkm_spike_ycs42_expressed.in.all.csv'))
write.csv(rpkm_all_ycs4_exp_all_mean, paste0(path, '2018_03_22_rpkm_exp_ycs42_expressed.in.all.csv'))
write.csv(tpm_all_spike_ycs4_exp_all_mean, paste0(path, '2018_03_22_tpm_spike_ycs42_expressed.in.all.csv'))
write.csv(tpm_all_ycs4_exp_all_mean, paste0(path, '2018_03_22_tpm_exp_ycs42_expressed.in.all.csv'))
write.csv(cpm_all_spike_ycs4_exp_all_mean, paste0(path, '2018_03_22_cpm_spike_ycs42_expressed.in.all.csv'))
write.csv(cpm_all_ycs4_exp_all_mean, paste0(path, '2018_03_22_cpm_exp_ycs42_expressed.in.all.csv'))

write.csv(rpkm_all_spike_brn1_exp_all_mean, paste0(path, '2018_03_22_rpkm_spike_brn1aaway_expressed.in.all.csv'))
write.csv(rpkm_all_brn1_exp_all_mean, paste0(path, '2018_03_22_rpkm_exp_brn1aaway_expressed.in.all.csv'))
write.csv(tpm_all_spike_brn1_exp_all_mean, paste0(path, '2018_03_22_tpm_spike_brn1aaway_expressed.in.all.csv'))
write.csv(tpm_all_brn1_exp_all_mean, paste0(path, '2018_03_22_tpm_exp_brn1aaway_expressed.in.all.csv'))
write.csv(cpm_all_spike_brn1_exp_all_mean, paste0(path, '2018_03_22_cpm_spike_brn1aaway_expressed.in.all.csv'))
write.csv(cpm_all_brn1_exp_all_mean, paste0(path, '2018_03_22_cpm_exp_brn1aaway_expressed.in.all.csv'))


#Run correction on all. Seperate and means

par(mfrow=c(2,3))
boxplot(cpm_all_brn1_exp_all_mean[,15:18],outline=F)
boxplot(cpm_all_ycg1_exp_all_mean[,15:18],outline=F)
boxplot(cpm_all_ycs4_exp_all_mean[,18:21],outline=F)
boxplot(cpm_all_spike_brn1_exp_all_mean[,15:18],outline=F)
boxplot(cpm_all_spike_ycg1_exp_all_mean[,15:18],outline=F)
boxplot(cpm_all_spike_ycs4_exp_all_mean[,18:21],outline=F)

boxplot(cpm_all_brn1_exp_mean[,15:18],outline=F)
boxplot(cpm_all_ycg1_exp_mean[,15:18],outline=F)
boxplot(cpm_all_ycs4_exp_mean[,18:21],outline=F)
boxplot(cpm_all_spike_brn1_exp_mean[,15:18],outline=F)
boxplot(cpm_all_spike_ycg1_exp_mean[,15:18],outline=F)
boxplot(cpm_all_spike_ycs4_exp_mean[,18:21],outline=F)

boxplot(cpm_all_brn1_mean[,15:18],outline=F)
boxplot(cpm_all_ycg1_mean[,15:18],outline=F)
boxplot(cpm_all_ycs4_mean[,18:21],outline=F)
boxplot(cpm_all_spike_brn1_mean[,15:18],outline=F)
boxplot(cpm_all_spike_ycg1_mean[,15:18],outline=F)
boxplot(cpm_all_spike_ycs4_mean[,18:21],outline=F)


#Variables for correction

#cpm_all_ycs4_exp_mean
cpm_all_ycs4_exp_mean_corrected_log2<-cpm_all_ycs4_exp_mean
bottom_mean<-length(colnames(cpm_all_ycs4_mean))-3
top_mean<-length(colnames(cpm_all_ycs4_mean))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in 8:(bottom_mean-1)) {
  print(i)
experiment_name<-colnames(cpm_all_ycs4_exp_mean)[i]
plot(log2(cpm_all_ycs4_exp_mean[,7]), log2(cpm_all_ycs4_exp_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
#linear model
r <-lm(log2(cpm_all_ycs4_exp_mean[,i]) ~ log2(cpm_all_ycs4_exp_mean[,7]))
abline(r, col="red")
abline(0,1, col="gray")
legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))

#(2) Correct the spike-in data
plot(log2(cpm_all_ycs4_exp_mean[,7]), (log2(cpm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
# corrected linear model
rc <-lm((log2(cpm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_ycs4_exp_mean[,7]))
abline(rc, col="red")
abline(0,1, col="gray")

#(3) experimental-experimental comparison
plot(log2(cpm_all_ycs4_exp_mean[,7]), log2(cpm_all_ycs4_exp_mean[,i]), pch=20, main="S. cerevisiae experimental", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
abline(0,1, col="gray")
abline(r, col="red")
r <-lm(log2(cpm_all_ycs4_exp_mean[,i]) ~ log2(cpm_all_ycs4_exp_mean[,7]))

#(4) Correct the spike-in data
plot(log2(cpm_all_ycs4_exp_mean[,7]), (log2(cpm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
# corrected linear model
rc <-lm((log2(cpm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_ycs4_exp_mean[,7]))
abline(rc, col="red")
abline(0,1, col="gray")


#(5) corrected matrix
cpm_all_ycs4_exp_mean_corrected_log2[,i]<-((log2(cpm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
cpm_all_ycs4_exp_mean_corrected_log2[,7]<-(log2(cpm_all_ycs4_exp_mean[,7]))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in (bottom_mean+1):top_mean) {
  print(i)
  experiment_name<-colnames(cpm_all_ycs4_exp_mean)[i]
  plot(log2(cpm_all_ycs4_exp_mean[,(bottom_mean)]), log2(cpm_all_ycs4_exp_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(cpm_all_ycs4_exp_mean[,i]) ~ log2(cpm_all_ycs4_exp_mean[,(bottom_mean)]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(cpm_all_ycs4_exp_mean[,(bottom_mean)]), (log2(cpm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_ycs4_exp_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) experimental-experimental comparison
  plot(log2(cpm_all_ycs4_exp_mean[,(bottom_mean)]), log2(cpm_all_ycs4_exp_mean[,i]), pch=20, main="S. cerevisiae experimental", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(cpm_all_ycs4_exp_mean[,i]) ~ log2(cpm_all_ycs4_exp_mean[,(bottom_mean)]))
  
  #(4) Correct the spike-in data
  plot(log2(cpm_all_ycs4_exp_mean[,(bottom_mean)]), (log2(cpm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_ycs4_exp_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  cpm_all_ycs4_exp_mean_corrected_log2[,i]<-((log2(cpm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
cpm_all_ycs4_exp_mean_corrected_log2[,bottom_mean]<-(log2(cpm_all_ycs4_exp_mean[,bottom_mean]))

cpm_all_ycs4_exp_mean_corrected<-cpm_all_ycs4_exp_mean_corrected_log2
cpm_all_ycs4_exp_mean_corrected[,7:top_mean]<-2^(cpm_all_ycs4_exp_mean_corrected_log2[,7:top_mean])


#tpm_all_ycs4_exp_mean
tpm_all_ycs4_exp_mean_corrected_log2<-tpm_all_ycs4_exp_mean
bottom_mean<-length(colnames(tpm_all_ycs4_mean))-3
top_mean<-length(colnames(tpm_all_ycs4_mean))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in 8:(bottom_mean-1)) {
  print(i)
  experiment_name<-colnames(tpm_all_ycs4_exp_mean)[i]
  plot(log2(tpm_all_ycs4_exp_mean[,7]), log2(tpm_all_ycs4_exp_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(tpm_all_ycs4_exp_mean[,i]) ~ log2(tpm_all_ycs4_exp_mean[,7]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(tpm_all_ycs4_exp_mean[,7]), (log2(tpm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_ycs4_exp_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) experimental-experimental comparison
  plot(log2(tpm_all_ycs4_exp_mean[,7]), log2(tpm_all_ycs4_exp_mean[,i]), pch=20, main="S. cerevisiae experimental", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(tpm_all_ycs4_exp_mean[,i]) ~ log2(tpm_all_ycs4_exp_mean[,7]))
  
  #(4) Correct the spike-in data
  plot(log2(tpm_all_ycs4_exp_mean[,7]), (log2(tpm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_ycs4_exp_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  tpm_all_ycs4_exp_mean_corrected_log2[,i]<-((log2(tpm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
tpm_all_ycs4_exp_mean_corrected_log2[,7]<-(log2(tpm_all_ycs4_exp_mean[,7]))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in (bottom_mean+1):top_mean) {
  print(i)
  experiment_name<-colnames(tpm_all_ycs4_exp_mean)[i]
  plot(log2(tpm_all_ycs4_exp_mean[,(bottom_mean)]), log2(tpm_all_ycs4_exp_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(tpm_all_ycs4_exp_mean[,i]) ~ log2(tpm_all_ycs4_exp_mean[,(bottom_mean)]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(tpm_all_ycs4_exp_mean[,(bottom_mean)]), (log2(tpm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_ycs4_exp_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) experimental-experimental comparison
  plot(log2(tpm_all_ycs4_exp_mean[,(bottom_mean)]), log2(tpm_all_ycs4_exp_mean[,i]), pch=20, main="S. cerevisiae experimental", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(tpm_all_ycs4_exp_mean[,i]) ~ log2(tpm_all_ycs4_exp_mean[,(bottom_mean)]))
  
  #(4) Correct the spike-in data
  plot(log2(tpm_all_ycs4_exp_mean[,(bottom_mean)]), (log2(tpm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_ycs4_exp_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  tpm_all_ycs4_exp_mean_corrected_log2[,i]<-((log2(tpm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
tpm_all_ycs4_exp_mean_corrected_log2[,bottom_mean]<-(log2(tpm_all_ycs4_exp_mean[,bottom_mean]))

tpm_all_ycs4_exp_mean_corrected<-tpm_all_ycs4_exp_mean_corrected_log2
tpm_all_ycs4_exp_mean_corrected[,7:top_mean]<-2^(tpm_all_ycs4_exp_mean_corrected_log2[,7:top_mean])



#rpkm_all_ycs4_exp_mean
rpkm_all_ycs4_exp_mean_corrected_log2<-rpkm_all_ycs4_exp_mean
bottom_mean<-length(colnames(rpkm_all_ycs4_mean))-3
top_mean<-length(colnames(rpkm_all_ycs4_mean))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in 8:(bottom_mean-1)) {
  print(i)
  experiment_name<-colnames(rpkm_all_ycs4_exp_mean)[i]
  plot(log2(rpkm_all_ycs4_exp_mean[,7]), log2(rpkm_all_ycs4_exp_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(rpkm_all_ycs4_exp_mean[,i]) ~ log2(rpkm_all_ycs4_exp_mean[,7]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(rpkm_all_ycs4_exp_mean[,7]), (log2(rpkm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_ycs4_exp_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) experimental-experimental comparison
  plot(log2(rpkm_all_ycs4_exp_mean[,7]), log2(rpkm_all_ycs4_exp_mean[,i]), pch=20, main="S. cerevisiae experimental", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(rpkm_all_ycs4_exp_mean[,i]) ~ log2(rpkm_all_ycs4_exp_mean[,7]))
  
  #(4) Correct the spike-in data
  plot(log2(rpkm_all_ycs4_exp_mean[,7]), (log2(rpkm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_ycs4_exp_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  rpkm_all_ycs4_exp_mean_corrected_log2[,i]<-((log2(rpkm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
rpkm_all_ycs4_exp_mean_corrected_log2[,7]<-(log2(rpkm_all_ycs4_exp_mean[,7]))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in (bottom_mean+1):top_mean) {
  print(i)
  experiment_name<-colnames(rpkm_all_ycs4_exp_mean)[i]
  plot(log2(rpkm_all_ycs4_exp_mean[,(bottom_mean)]), log2(rpkm_all_ycs4_exp_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(rpkm_all_ycs4_exp_mean[,i]) ~ log2(rpkm_all_ycs4_exp_mean[,(bottom_mean)]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(rpkm_all_ycs4_exp_mean[,(bottom_mean)]), (log2(rpkm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_ycs4_exp_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) experimental-experimental comparison
  plot(log2(rpkm_all_ycs4_exp_mean[,(bottom_mean)]), log2(rpkm_all_ycs4_exp_mean[,i]), pch=20, main="S. cerevisiae experimental", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(rpkm_all_ycs4_exp_mean[,i]) ~ log2(rpkm_all_ycs4_exp_mean[,(bottom_mean)]))
  
  #(4) Correct the spike-in data
  plot(log2(rpkm_all_ycs4_exp_mean[,(bottom_mean)]), (log2(rpkm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_ycs4_exp_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  rpkm_all_ycs4_exp_mean_corrected_log2[,i]<-((log2(rpkm_all_ycs4_exp_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
rpkm_all_ycs4_exp_mean_corrected_log2[,bottom_mean]<-(log2(rpkm_all_ycs4_exp_mean[,bottom_mean]))

rpkm_all_ycs4_exp_mean_corrected<-rpkm_all_ycs4_exp_mean_corrected_log2
rpkm_all_ycs4_exp_mean_corrected[,7:top_mean]<-2^(rpkm_all_ycs4_exp_mean_corrected_log2[,7:top_mean])



#cpm_all_ycg1_exp_mean
cpm_all_ycg1_exp_mean_corrected_log2<-cpm_all_ycg1_exp_mean
bottom_mean<-length(colnames(cpm_all_ycg1_mean))-3
top_mean<-length(colnames(cpm_all_ycg1_mean))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in 8:(bottom_mean-1)) {
  print(i)
  experiment_name<-colnames(cpm_all_ycg1_exp_mean)[i]
  plot(log2(cpm_all_ycg1_exp_mean[,7]), log2(cpm_all_ycg1_exp_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(cpm_all_ycg1_exp_mean[,i]) ~ log2(cpm_all_ycg1_exp_mean[,7]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(cpm_all_ycg1_exp_mean[,7]), (log2(cpm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_ycg1_exp_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) experimental-experimental comparison
  plot(log2(cpm_all_ycg1_exp_mean[,7]), log2(cpm_all_ycg1_exp_mean[,i]), pch=20, main="S. cerevisiae experimental", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(cpm_all_ycg1_exp_mean[,i]) ~ log2(cpm_all_ycg1_exp_mean[,7]))
  
  #(4) Correct the spike-in data
  plot(log2(cpm_all_ycg1_exp_mean[,7]), (log2(cpm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_ycg1_exp_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  cpm_all_ycg1_exp_mean_corrected_log2[,i]<-((log2(cpm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
cpm_all_ycg1_exp_mean_corrected_log2[,7]<-(log2(cpm_all_ycg1_exp_mean[,7]))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in (bottom_mean+1):top_mean) {
  print(i)
  experiment_name<-colnames(cpm_all_ycg1_exp_mean)[i]
  plot(log2(cpm_all_ycg1_exp_mean[,(bottom_mean)]), log2(cpm_all_ycg1_exp_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(cpm_all_ycg1_exp_mean[,i]) ~ log2(cpm_all_ycg1_exp_mean[,(bottom_mean)]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(cpm_all_ycg1_exp_mean[,(bottom_mean)]), (log2(cpm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_ycg1_exp_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) experimental-experimental comparison
  plot(log2(cpm_all_ycg1_exp_mean[,(bottom_mean)]), log2(cpm_all_ycg1_exp_mean[,i]), pch=20, main="S. cerevisiae experimental", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(cpm_all_ycg1_exp_mean[,i]) ~ log2(cpm_all_ycg1_exp_mean[,(bottom_mean)]))
  
  #(4) Correct the spike-in data
  plot(log2(cpm_all_ycg1_exp_mean[,(bottom_mean)]), (log2(cpm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_ycg1_exp_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  cpm_all_ycg1_exp_mean_corrected_log2[,i]<-((log2(cpm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
cpm_all_ycg1_exp_mean_corrected_log2[,bottom_mean]<-(log2(cpm_all_ycg1_exp_mean[,bottom_mean]))

cpm_all_ycg1_exp_mean_corrected<-cpm_all_ycg1_exp_mean_corrected_log2
cpm_all_ycg1_exp_mean_corrected[,7:top_mean]<-2^(cpm_all_ycg1_exp_mean_corrected_log2[,7:top_mean])


#tpm_all_ycg1_exp_mean
tpm_all_ycg1_exp_mean_corrected_log2<-tpm_all_ycg1_exp_mean
bottom_mean<-length(colnames(tpm_all_ycg1_mean))-3
top_mean<-length(colnames(tpm_all_ycg1_mean))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in 8:(bottom_mean-1)) {
  print(i)
  experiment_name<-colnames(tpm_all_ycg1_exp_mean)[i]
  plot(log2(tpm_all_ycg1_exp_mean[,7]), log2(tpm_all_ycg1_exp_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(tpm_all_ycg1_exp_mean[,i]) ~ log2(tpm_all_ycg1_exp_mean[,7]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(tpm_all_ycg1_exp_mean[,7]), (log2(tpm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_ycg1_exp_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) experimental-experimental comparison
  plot(log2(tpm_all_ycg1_exp_mean[,7]), log2(tpm_all_ycg1_exp_mean[,i]), pch=20, main="S. cerevisiae experimental", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(tpm_all_ycg1_exp_mean[,i]) ~ log2(tpm_all_ycg1_exp_mean[,7]))
  
  #(4) Correct the spike-in data
  plot(log2(tpm_all_ycg1_exp_mean[,7]), (log2(tpm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_ycg1_exp_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  tpm_all_ycg1_exp_mean_corrected_log2[,i]<-((log2(tpm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
tpm_all_ycg1_exp_mean_corrected_log2[,7]<-(log2(tpm_all_ycg1_exp_mean[,7]))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in (bottom_mean+1):top_mean) {
  print(i)
  experiment_name<-colnames(tpm_all_ycg1_exp_mean)[i]
  plot(log2(tpm_all_ycg1_exp_mean[,(bottom_mean)]), log2(tpm_all_ycg1_exp_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(tpm_all_ycg1_exp_mean[,i]) ~ log2(tpm_all_ycg1_exp_mean[,(bottom_mean)]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(tpm_all_ycg1_exp_mean[,(bottom_mean)]), (log2(tpm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_ycg1_exp_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) experimental-experimental comparison
  plot(log2(tpm_all_ycg1_exp_mean[,(bottom_mean)]), log2(tpm_all_ycg1_exp_mean[,i]), pch=20, main="S. cerevisiae experimental", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(tpm_all_ycg1_exp_mean[,i]) ~ log2(tpm_all_ycg1_exp_mean[,(bottom_mean)]))
  
  #(4) Correct the spike-in data
  plot(log2(tpm_all_ycg1_exp_mean[,(bottom_mean)]), (log2(tpm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_ycg1_exp_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  tpm_all_ycg1_exp_mean_corrected_log2[,i]<-((log2(tpm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
tpm_all_ycg1_exp_mean_corrected_log2[,bottom_mean]<-(log2(tpm_all_ycg1_exp_mean[,bottom_mean]))

tpm_all_ycg1_exp_mean_corrected<-tpm_all_ycg1_exp_mean_corrected_log2
tpm_all_ycg1_exp_mean_corrected[,7:top_mean]<-2^(tpm_all_ycg1_exp_mean_corrected_log2[,7:top_mean])



#rpkm_all_ycg1_exp_mean
rpkm_all_ycg1_exp_mean_corrected_log2<-rpkm_all_ycg1_exp_mean
bottom_mean<-length(colnames(rpkm_all_ycg1_mean))-3
top_mean<-length(colnames(rpkm_all_ycg1_mean))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in 8:(bottom_mean-1)) {
  print(i)
  experiment_name<-colnames(rpkm_all_ycg1_exp_mean)[i]
  plot(log2(rpkm_all_ycg1_exp_mean[,7]), log2(rpkm_all_ycg1_exp_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(rpkm_all_ycg1_exp_mean[,i]) ~ log2(rpkm_all_ycg1_exp_mean[,7]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(rpkm_all_ycg1_exp_mean[,7]), (log2(rpkm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_ycg1_exp_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) experimental-experimental comparison
  plot(log2(rpkm_all_ycg1_exp_mean[,7]), log2(rpkm_all_ycg1_exp_mean[,i]), pch=20, main="S. cerevisiae experimental", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(rpkm_all_ycg1_exp_mean[,i]) ~ log2(rpkm_all_ycg1_exp_mean[,7]))
  
  #(4) Correct the spike-in data
  plot(log2(rpkm_all_ycg1_exp_mean[,7]), (log2(rpkm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_ycg1_exp_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  rpkm_all_ycg1_exp_mean_corrected_log2[,i]<-((log2(rpkm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
rpkm_all_ycg1_exp_mean_corrected_log2[,7]<-(log2(rpkm_all_ycg1_exp_mean[,7]))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in (bottom_mean+1):top_mean) {
  print(i)
  experiment_name<-colnames(rpkm_all_ycg1_exp_mean)[i]
  plot(log2(rpkm_all_ycg1_exp_mean[,(bottom_mean)]), log2(rpkm_all_ycg1_exp_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(rpkm_all_ycg1_exp_mean[,i]) ~ log2(rpkm_all_ycg1_exp_mean[,(bottom_mean)]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(rpkm_all_ycg1_exp_mean[,(bottom_mean)]), (log2(rpkm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_ycg1_exp_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) experimental-experimental comparison
  plot(log2(rpkm_all_ycg1_exp_mean[,(bottom_mean)]), log2(rpkm_all_ycg1_exp_mean[,i]), pch=20, main="S. cerevisiae experimental", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(rpkm_all_ycg1_exp_mean[,i]) ~ log2(rpkm_all_ycg1_exp_mean[,(bottom_mean)]))
  
  #(4) Correct the spike-in data
  plot(log2(rpkm_all_ycg1_exp_mean[,(bottom_mean)]), (log2(rpkm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_ycg1_exp_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  rpkm_all_ycg1_exp_mean_corrected_log2[,i]<-((log2(rpkm_all_ycg1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
rpkm_all_ycg1_exp_mean_corrected_log2[,bottom_mean]<-(log2(rpkm_all_ycg1_exp_mean[,bottom_mean]))

rpkm_all_ycg1_exp_mean_corrected<-rpkm_all_ycg1_exp_mean_corrected_log2
rpkm_all_ycg1_exp_mean_corrected[,7:top_mean]<-2^(rpkm_all_ycg1_exp_mean_corrected_log2[,7:top_mean])


#cpm_all_brn1_exp_mean
cpm_all_brn1_exp_mean_corrected_log2<-cpm_all_brn1_exp_mean
bottom_mean<-length(colnames(cpm_all_brn1_mean))-3
top_mean<-length(colnames(cpm_all_brn1_mean))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in 8:(bottom_mean-1)) {
  print(i)
  experiment_name<-colnames(cpm_all_brn1_exp_mean)[i]
  plot(log2(cpm_all_brn1_exp_mean[,7]), log2(cpm_all_brn1_exp_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(cpm_all_brn1_exp_mean[,i]) ~ log2(cpm_all_brn1_exp_mean[,7]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(cpm_all_brn1_exp_mean[,7]), (log2(cpm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_brn1_exp_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) experimental-experimental comparison
  plot(log2(cpm_all_brn1_exp_mean[,7]), log2(cpm_all_brn1_exp_mean[,i]), pch=20, main="S. cerevisiae experimental", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(cpm_all_brn1_exp_mean[,i]) ~ log2(cpm_all_brn1_exp_mean[,7]))
  
  #(4) Correct the spike-in data
  plot(log2(cpm_all_brn1_exp_mean[,7]), (log2(cpm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_brn1_exp_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  cpm_all_brn1_exp_mean_corrected_log2[,i]<-((log2(cpm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
cpm_all_brn1_exp_mean_corrected_log2[,7]<-(log2(cpm_all_brn1_exp_mean[,7]))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in (bottom_mean+1):top_mean) {
  print(i)
  experiment_name<-colnames(cpm_all_brn1_exp_mean)[i]
  plot(log2(cpm_all_brn1_exp_mean[,(bottom_mean)]), log2(cpm_all_brn1_exp_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(cpm_all_brn1_exp_mean[,i]) ~ log2(cpm_all_brn1_exp_mean[,(bottom_mean)]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(cpm_all_brn1_exp_mean[,(bottom_mean)]), (log2(cpm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_brn1_exp_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) experimental-experimental comparison
  plot(log2(cpm_all_brn1_exp_mean[,(bottom_mean)]), log2(cpm_all_brn1_exp_mean[,i]), pch=20, main="S. cerevisiae experimental", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(cpm_all_brn1_exp_mean[,i]) ~ log2(cpm_all_brn1_exp_mean[,(bottom_mean)]))
  
  #(4) Correct the spike-in data
  plot(log2(cpm_all_brn1_exp_mean[,(bottom_mean)]), (log2(cpm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',experiment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_brn1_exp_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  cpm_all_brn1_exp_mean_corrected_log2[,i]<-((log2(cpm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
cpm_all_brn1_exp_mean_corrected_log2[,bottom_mean]<-(log2(cpm_all_brn1_exp_mean[,bottom_mean]))

cpm_all_brn1_exp_mean_corrected<-cpm_all_brn1_exp_mean_corrected_log2
cpm_all_brn1_exp_mean_corrected[,7:top_mean]<-2^(cpm_all_brn1_exp_mean_corrected_log2[,7:top_mean])


#tpm_all_brn1_exp_mean
tpm_all_brn1_exp_mean_corrected_log2<-tpm_all_brn1_exp_mean
bottom_mean<-length(colnames(tpm_all_brn1_mean))-3
top_mean<-length(colnames(tpm_all_brn1_mean))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in 8:(bottom_mean-1)) {
  print(i)
  experiment_name<-colnames(tpm_all_brn1_exp_mean)[i]
  plot(log2(tpm_all_brn1_exp_mean[,7]), log2(tpm_all_brn1_exp_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(tpm_all_brn1_exp_mean[,i]) ~ log2(tpm_all_brn1_exp_mean[,7]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(tpm_all_brn1_exp_mean[,7]), (log2(tpm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_brn1_exp_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) experimental-experimental comparison
  plot(log2(tpm_all_brn1_exp_mean[,7]), log2(tpm_all_brn1_exp_mean[,i]), pch=20, main="S. cerevisiae experimental", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(tpm_all_brn1_exp_mean[,i]) ~ log2(tpm_all_brn1_exp_mean[,7]))
  
  #(4) Correct the spike-in data
  plot(log2(tpm_all_brn1_exp_mean[,7]), (log2(tpm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_brn1_exp_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  tpm_all_brn1_exp_mean_corrected_log2[,i]<-((log2(tpm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
tpm_all_brn1_exp_mean_corrected_log2[,7]<-(log2(tpm_all_brn1_exp_mean[,7]))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in (bottom_mean+1):top_mean) {
  print(i)
  experiment_name<-colnames(tpm_all_brn1_exp_mean)[i]
  plot(log2(tpm_all_brn1_exp_mean[,(bottom_mean)]), log2(tpm_all_brn1_exp_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(tpm_all_brn1_exp_mean[,i]) ~ log2(tpm_all_brn1_exp_mean[,(bottom_mean)]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(tpm_all_brn1_exp_mean[,(bottom_mean)]), (log2(tpm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_brn1_exp_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) experimental-experimental comparison
  plot(log2(tpm_all_brn1_exp_mean[,(bottom_mean)]), log2(tpm_all_brn1_exp_mean[,i]), pch=20, main="S. cerevisiae experimental", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(tpm_all_brn1_exp_mean[,i]) ~ log2(tpm_all_brn1_exp_mean[,(bottom_mean)]))
  
  #(4) Correct the spike-in data
  plot(log2(tpm_all_brn1_exp_mean[,(bottom_mean)]), (log2(tpm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',experiment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_brn1_exp_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  tpm_all_brn1_exp_mean_corrected_log2[,i]<-((log2(tpm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
tpm_all_brn1_exp_mean_corrected_log2[,bottom_mean]<-(log2(tpm_all_brn1_exp_mean[,bottom_mean]))

tpm_all_brn1_exp_mean_corrected<-tpm_all_brn1_exp_mean_corrected_log2
tpm_all_brn1_exp_mean_corrected[,7:top_mean]<-2^(tpm_all_brn1_exp_mean_corrected_log2[,7:top_mean])



#rpkm_all_brn1_exp_mean
rpkm_all_brn1_exp_mean_corrected_log2<-rpkm_all_brn1_exp_mean
bottom_mean<-length(colnames(rpkm_all_brn1_mean))-3
top_mean<-length(colnames(rpkm_all_brn1_mean))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in 8:(bottom_mean-1)) {
  print(i)
  experiment_name<-colnames(rpkm_all_brn1_exp_mean)[i]
  plot(log2(rpkm_all_brn1_exp_mean[,7]), log2(rpkm_all_brn1_exp_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(rpkm_all_brn1_exp_mean[,i]) ~ log2(rpkm_all_brn1_exp_mean[,7]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(rpkm_all_brn1_exp_mean[,7]), (log2(rpkm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_brn1_exp_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) experimental-experimental comparison
  plot(log2(rpkm_all_brn1_exp_mean[,7]), log2(rpkm_all_brn1_exp_mean[,i]), pch=20, main="S. cerevisiae experimental", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(rpkm_all_brn1_exp_mean[,i]) ~ log2(rpkm_all_brn1_exp_mean[,7]))
  
  #(4) Correct the spike-in data
  plot(log2(rpkm_all_brn1_exp_mean[,7]), (log2(rpkm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_brn1_exp_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  rpkm_all_brn1_exp_mean_corrected_log2[,i]<-((log2(rpkm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
rpkm_all_brn1_exp_mean_corrected_log2[,7]<-(log2(rpkm_all_brn1_exp_mean[,7]))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in (bottom_mean+1):top_mean) {
  print(i)
  experiment_name<-colnames(rpkm_all_brn1_exp_mean)[i]
  plot(log2(rpkm_all_brn1_exp_mean[,(bottom_mean)]), log2(rpkm_all_brn1_exp_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(rpkm_all_brn1_exp_mean[,i]) ~ log2(rpkm_all_brn1_exp_mean[,(bottom_mean)]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(rpkm_all_brn1_exp_mean[,(bottom_mean)]), (log2(rpkm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_brn1_exp_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) experimental-experimental comparison
  plot(log2(rpkm_all_brn1_exp_mean[,(bottom_mean)]), log2(rpkm_all_brn1_exp_mean[,i]), pch=20, main="S. cerevisiae experimental", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(rpkm_all_brn1_exp_mean[,i]) ~ log2(rpkm_all_brn1_exp_mean[,(bottom_mean)]))
  
  #(4) Correct the spike-in data
  plot(log2(rpkm_all_brn1_exp_mean[,(bottom_mean)]), (log2(rpkm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',experiment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_brn1_exp_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  rpkm_all_brn1_exp_mean_corrected_log2[,i]<-((log2(rpkm_all_brn1_exp_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
rpkm_all_brn1_exp_mean_corrected_log2[,bottom_mean]<-(log2(rpkm_all_brn1_exp_mean[,bottom_mean]))

rpkm_all_brn1_exp_mean_corrected<-rpkm_all_brn1_exp_mean_corrected_log2
rpkm_all_brn1_exp_mean_corrected[,7:top_mean]<-2^(rpkm_all_brn1_exp_mean_corrected_log2[,7:top_mean])




#cpm_all_ycs4_exp_all_mean
cpm_all_ycs4_exp_all_mean_corrected_log2<-cpm_all_ycs4_exp_all_mean
bottom_mean<-length(colnames(cpm_all_ycs4_mean))-3
top_mean<-length(colnames(cpm_all_ycs4_mean))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in 8:(bottom_mean-1)) {
  print(i)
  exp_alleriment_name<-colnames(cpm_all_ycs4_exp_all_mean)[i]
  plot(log2(cpm_all_ycs4_exp_all_mean[,7]), log2(cpm_all_ycs4_exp_all_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(cpm_all_ycs4_exp_all_mean[,i]) ~ log2(cpm_all_ycs4_exp_all_mean[,7]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(cpm_all_ycs4_exp_all_mean[,7]), (log2(cpm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_ycs4_exp_all_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) exp_allerimental-exp_allerimental comparison
  plot(log2(cpm_all_ycs4_exp_all_mean[,7]), log2(cpm_all_ycs4_exp_all_mean[,i]), pch=20, main="S. cerevisiae exp_allerimental", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(cpm_all_ycs4_exp_all_mean[,i]) ~ log2(cpm_all_ycs4_exp_all_mean[,7]))
  
  #(4) Correct the spike-in data
  plot(log2(cpm_all_ycs4_exp_all_mean[,7]), (log2(cpm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_ycs4_exp_all_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  cpm_all_ycs4_exp_all_mean_corrected_log2[,i]<-((log2(cpm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
cpm_all_ycs4_exp_all_mean_corrected_log2[,7]<-(log2(cpm_all_ycs4_exp_all_mean[,7]))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in (bottom_mean+1):top_mean) {
  print(i)
  exp_alleriment_name<-colnames(cpm_all_ycs4_exp_all_mean)[i]
  plot(log2(cpm_all_ycs4_exp_all_mean[,(bottom_mean)]), log2(cpm_all_ycs4_exp_all_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(cpm_all_ycs4_exp_all_mean[,i]) ~ log2(cpm_all_ycs4_exp_all_mean[,(bottom_mean)]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(cpm_all_ycs4_exp_all_mean[,(bottom_mean)]), (log2(cpm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_ycs4_exp_all_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) exp_allerimental-exp_allerimental comparison
  plot(log2(cpm_all_ycs4_exp_all_mean[,(bottom_mean)]), log2(cpm_all_ycs4_exp_all_mean[,i]), pch=20, main="S. cerevisiae exp_allerimental", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(cpm_all_ycs4_exp_all_mean[,i]) ~ log2(cpm_all_ycs4_exp_all_mean[,(bottom_mean)]))
  
  #(4) Correct the spike-in data
  plot(log2(cpm_all_ycs4_exp_all_mean[,(bottom_mean)]), (log2(cpm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_ycs4_exp_all_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  cpm_all_ycs4_exp_all_mean_corrected_log2[,i]<-((log2(cpm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
cpm_all_ycs4_exp_all_mean_corrected_log2[,bottom_mean]<-(log2(cpm_all_ycs4_exp_all_mean[,bottom_mean]))

cpm_all_ycs4_exp_all_mean_corrected<-cpm_all_ycs4_exp_all_mean_corrected_log2
cpm_all_ycs4_exp_all_mean_corrected[,7:top_mean]<-2^(cpm_all_ycs4_exp_all_mean_corrected_log2[,7:top_mean])


#tpm_all_ycs4_exp_all_mean
tpm_all_ycs4_exp_all_mean_corrected_log2<-tpm_all_ycs4_exp_all_mean
bottom_mean<-length(colnames(tpm_all_ycs4_mean))-3
top_mean<-length(colnames(tpm_all_ycs4_mean))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in 8:(bottom_mean-1)) {
  print(i)
  exp_alleriment_name<-colnames(tpm_all_ycs4_exp_all_mean)[i]
  plot(log2(tpm_all_ycs4_exp_all_mean[,7]), log2(tpm_all_ycs4_exp_all_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(tpm_all_ycs4_exp_all_mean[,i]) ~ log2(tpm_all_ycs4_exp_all_mean[,7]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(tpm_all_ycs4_exp_all_mean[,7]), (log2(tpm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_ycs4_exp_all_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) exp_allerimental-exp_allerimental comparison
  plot(log2(tpm_all_ycs4_exp_all_mean[,7]), log2(tpm_all_ycs4_exp_all_mean[,i]), pch=20, main="S. cerevisiae exp_allerimental", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(tpm_all_ycs4_exp_all_mean[,i]) ~ log2(tpm_all_ycs4_exp_all_mean[,7]))
  
  #(4) Correct the spike-in data
  plot(log2(tpm_all_ycs4_exp_all_mean[,7]), (log2(tpm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_ycs4_exp_all_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  tpm_all_ycs4_exp_all_mean_corrected_log2[,i]<-((log2(tpm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
tpm_all_ycs4_exp_all_mean_corrected_log2[,7]<-(log2(tpm_all_ycs4_exp_all_mean[,7]))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in (bottom_mean+1):top_mean) {
  print(i)
  exp_alleriment_name<-colnames(tpm_all_ycs4_exp_all_mean)[i]
  plot(log2(tpm_all_ycs4_exp_all_mean[,(bottom_mean)]), log2(tpm_all_ycs4_exp_all_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(tpm_all_ycs4_exp_all_mean[,i]) ~ log2(tpm_all_ycs4_exp_all_mean[,(bottom_mean)]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(tpm_all_ycs4_exp_all_mean[,(bottom_mean)]), (log2(tpm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_ycs4_exp_all_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) exp_allerimental-exp_allerimental comparison
  plot(log2(tpm_all_ycs4_exp_all_mean[,(bottom_mean)]), log2(tpm_all_ycs4_exp_all_mean[,i]), pch=20, main="S. cerevisiae exp_allerimental", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(tpm_all_ycs4_exp_all_mean[,i]) ~ log2(tpm_all_ycs4_exp_all_mean[,(bottom_mean)]))
  
  #(4) Correct the spike-in data
  plot(log2(tpm_all_ycs4_exp_all_mean[,(bottom_mean)]), (log2(tpm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_ycs4_exp_all_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  tpm_all_ycs4_exp_all_mean_corrected_log2[,i]<-((log2(tpm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
tpm_all_ycs4_exp_all_mean_corrected_log2[,bottom_mean]<-(log2(tpm_all_ycs4_exp_all_mean[,bottom_mean]))

tpm_all_ycs4_exp_all_mean_corrected<-tpm_all_ycs4_exp_all_mean_corrected_log2
tpm_all_ycs4_exp_all_mean_corrected[,7:top_mean]<-2^(tpm_all_ycs4_exp_all_mean_corrected_log2[,7:top_mean])



#rpkm_all_ycs4_exp_all_mean
rpkm_all_ycs4_exp_all_mean_corrected_log2<-rpkm_all_ycs4_exp_all_mean
bottom_mean<-length(colnames(rpkm_all_ycs4_mean))-3
top_mean<-length(colnames(rpkm_all_ycs4_mean))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in 8:(bottom_mean-1)) {
  print(i)
  exp_alleriment_name<-colnames(rpkm_all_ycs4_exp_all_mean)[i]
  plot(log2(rpkm_all_ycs4_exp_all_mean[,7]), log2(rpkm_all_ycs4_exp_all_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(rpkm_all_ycs4_exp_all_mean[,i]) ~ log2(rpkm_all_ycs4_exp_all_mean[,7]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(rpkm_all_ycs4_exp_all_mean[,7]), (log2(rpkm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_ycs4_exp_all_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) exp_allerimental-exp_allerimental comparison
  plot(log2(rpkm_all_ycs4_exp_all_mean[,7]), log2(rpkm_all_ycs4_exp_all_mean[,i]), pch=20, main="S. cerevisiae exp_allerimental", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(rpkm_all_ycs4_exp_all_mean[,i]) ~ log2(rpkm_all_ycs4_exp_all_mean[,7]))
  
  #(4) Correct the spike-in data
  plot(log2(rpkm_all_ycs4_exp_all_mean[,7]), (log2(rpkm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_ycs4_exp_all_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  rpkm_all_ycs4_exp_all_mean_corrected_log2[,i]<-((log2(rpkm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
rpkm_all_ycs4_exp_all_mean_corrected_log2[,7]<-(log2(rpkm_all_ycs4_exp_all_mean[,7]))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in (bottom_mean+1):top_mean) {
  print(i)
  exp_alleriment_name<-colnames(rpkm_all_ycs4_exp_all_mean)[i]
  plot(log2(rpkm_all_ycs4_exp_all_mean[,(bottom_mean)]), log2(rpkm_all_ycs4_exp_all_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(rpkm_all_ycs4_exp_all_mean[,i]) ~ log2(rpkm_all_ycs4_exp_all_mean[,(bottom_mean)]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(rpkm_all_ycs4_exp_all_mean[,(bottom_mean)]), (log2(rpkm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_ycs4_exp_all_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) exp_allerimental-exp_allerimental comparison
  plot(log2(rpkm_all_ycs4_exp_all_mean[,(bottom_mean)]), log2(rpkm_all_ycs4_exp_all_mean[,i]), pch=20, main="S. cerevisiae exp_allerimental", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(rpkm_all_ycs4_exp_all_mean[,i]) ~ log2(rpkm_all_ycs4_exp_all_mean[,(bottom_mean)]))
  
  #(4) Correct the spike-in data
  plot(log2(rpkm_all_ycs4_exp_all_mean[,(bottom_mean)]), (log2(rpkm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_ycs4_exp_all_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  rpkm_all_ycs4_exp_all_mean_corrected_log2[,i]<-((log2(rpkm_all_ycs4_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
rpkm_all_ycs4_exp_all_mean_corrected_log2[,bottom_mean]<-(log2(rpkm_all_ycs4_exp_all_mean[,bottom_mean]))

rpkm_all_ycs4_exp_all_mean_corrected<-rpkm_all_ycs4_exp_all_mean_corrected_log2
rpkm_all_ycs4_exp_all_mean_corrected[,7:top_mean]<-2^(rpkm_all_ycs4_exp_all_mean_corrected_log2[,7:top_mean])



#cpm_all_ycg1_exp_all_mean
cpm_all_ycg1_exp_all_mean_corrected_log2<-cpm_all_ycg1_exp_all_mean
bottom_mean<-length(colnames(cpm_all_ycg1_mean))-3
top_mean<-length(colnames(cpm_all_ycg1_mean))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in 8:(bottom_mean-1)) {
  print(i)
  exp_alleriment_name<-colnames(cpm_all_ycg1_exp_all_mean)[i]
  plot(log2(cpm_all_ycg1_exp_all_mean[,7]), log2(cpm_all_ycg1_exp_all_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(cpm_all_ycg1_exp_all_mean[,i]) ~ log2(cpm_all_ycg1_exp_all_mean[,7]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(cpm_all_ycg1_exp_all_mean[,7]), (log2(cpm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_ycg1_exp_all_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) exp_allerimental-exp_allerimental comparison
  plot(log2(cpm_all_ycg1_exp_all_mean[,7]), log2(cpm_all_ycg1_exp_all_mean[,i]), pch=20, main="S. cerevisiae exp_allerimental", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(cpm_all_ycg1_exp_all_mean[,i]) ~ log2(cpm_all_ycg1_exp_all_mean[,7]))
  
  #(4) Correct the spike-in data
  plot(log2(cpm_all_ycg1_exp_all_mean[,7]), (log2(cpm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_ycg1_exp_all_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  cpm_all_ycg1_exp_all_mean_corrected_log2[,i]<-((log2(cpm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
cpm_all_ycg1_exp_all_mean_corrected_log2[,7]<-(log2(cpm_all_ycg1_exp_all_mean[,7]))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in (bottom_mean+1):top_mean) {
  print(i)
  exp_alleriment_name<-colnames(cpm_all_ycg1_exp_all_mean)[i]
  plot(log2(cpm_all_ycg1_exp_all_mean[,(bottom_mean)]), log2(cpm_all_ycg1_exp_all_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(cpm_all_ycg1_exp_all_mean[,i]) ~ log2(cpm_all_ycg1_exp_all_mean[,(bottom_mean)]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(cpm_all_ycg1_exp_all_mean[,(bottom_mean)]), (log2(cpm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_ycg1_exp_all_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) exp_allerimental-exp_allerimental comparison
  plot(log2(cpm_all_ycg1_exp_all_mean[,(bottom_mean)]), log2(cpm_all_ycg1_exp_all_mean[,i]), pch=20, main="S. cerevisiae exp_allerimental", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(cpm_all_ycg1_exp_all_mean[,i]) ~ log2(cpm_all_ycg1_exp_all_mean[,(bottom_mean)]))
  
  #(4) Correct the spike-in data
  plot(log2(cpm_all_ycg1_exp_all_mean[,(bottom_mean)]), (log2(cpm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_ycg1_exp_all_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  cpm_all_ycg1_exp_all_mean_corrected_log2[,i]<-((log2(cpm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
cpm_all_ycg1_exp_all_mean_corrected_log2[,bottom_mean]<-(log2(cpm_all_ycg1_exp_all_mean[,bottom_mean]))

cpm_all_ycg1_exp_all_mean_corrected<-cpm_all_ycg1_exp_all_mean_corrected_log2
cpm_all_ycg1_exp_all_mean_corrected[,7:top_mean]<-2^(cpm_all_ycg1_exp_all_mean_corrected_log2[,7:top_mean])


#tpm_all_ycg1_exp_all_mean
tpm_all_ycg1_exp_all_mean_corrected_log2<-tpm_all_ycg1_exp_all_mean
bottom_mean<-length(colnames(tpm_all_ycg1_mean))-3
top_mean<-length(colnames(tpm_all_ycg1_mean))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in 8:(bottom_mean-1)) {
  print(i)
  exp_alleriment_name<-colnames(tpm_all_ycg1_exp_all_mean)[i]
  plot(log2(tpm_all_ycg1_exp_all_mean[,7]), log2(tpm_all_ycg1_exp_all_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(tpm_all_ycg1_exp_all_mean[,i]) ~ log2(tpm_all_ycg1_exp_all_mean[,7]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(tpm_all_ycg1_exp_all_mean[,7]), (log2(tpm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_ycg1_exp_all_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) exp_allerimental-exp_allerimental comparison
  plot(log2(tpm_all_ycg1_exp_all_mean[,7]), log2(tpm_all_ycg1_exp_all_mean[,i]), pch=20, main="S. cerevisiae exp_allerimental", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(tpm_all_ycg1_exp_all_mean[,i]) ~ log2(tpm_all_ycg1_exp_all_mean[,7]))
  
  #(4) Correct the spike-in data
  plot(log2(tpm_all_ycg1_exp_all_mean[,7]), (log2(tpm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_ycg1_exp_all_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  tpm_all_ycg1_exp_all_mean_corrected_log2[,i]<-((log2(tpm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
tpm_all_ycg1_exp_all_mean_corrected_log2[,7]<-(log2(tpm_all_ycg1_exp_all_mean[,7]))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in (bottom_mean+1):top_mean) {
  print(i)
  exp_alleriment_name<-colnames(tpm_all_ycg1_exp_all_mean)[i]
  plot(log2(tpm_all_ycg1_exp_all_mean[,(bottom_mean)]), log2(tpm_all_ycg1_exp_all_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(tpm_all_ycg1_exp_all_mean[,i]) ~ log2(tpm_all_ycg1_exp_all_mean[,(bottom_mean)]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(tpm_all_ycg1_exp_all_mean[,(bottom_mean)]), (log2(tpm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_ycg1_exp_all_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) exp_allerimental-exp_allerimental comparison
  plot(log2(tpm_all_ycg1_exp_all_mean[,(bottom_mean)]), log2(tpm_all_ycg1_exp_all_mean[,i]), pch=20, main="S. cerevisiae exp_allerimental", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(tpm_all_ycg1_exp_all_mean[,i]) ~ log2(tpm_all_ycg1_exp_all_mean[,(bottom_mean)]))
  
  #(4) Correct the spike-in data
  plot(log2(tpm_all_ycg1_exp_all_mean[,(bottom_mean)]), (log2(tpm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_ycg1_exp_all_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  tpm_all_ycg1_exp_all_mean_corrected_log2[,i]<-((log2(tpm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
tpm_all_ycg1_exp_all_mean_corrected_log2[,bottom_mean]<-(log2(tpm_all_ycg1_exp_all_mean[,bottom_mean]))

tpm_all_ycg1_exp_all_mean_corrected<-tpm_all_ycg1_exp_all_mean_corrected_log2
tpm_all_ycg1_exp_all_mean_corrected[,7:top_mean]<-2^(tpm_all_ycg1_exp_all_mean_corrected_log2[,7:top_mean])



#rpkm_all_ycg1_exp_all_mean
rpkm_all_ycg1_exp_all_mean_corrected_log2<-rpkm_all_ycg1_exp_all_mean
bottom_mean<-length(colnames(rpkm_all_ycg1_mean))-3
top_mean<-length(colnames(rpkm_all_ycg1_mean))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in 8:(bottom_mean-1)) {
  print(i)
  exp_alleriment_name<-colnames(rpkm_all_ycg1_exp_all_mean)[i]
  plot(log2(rpkm_all_ycg1_exp_all_mean[,7]), log2(rpkm_all_ycg1_exp_all_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(rpkm_all_ycg1_exp_all_mean[,i]) ~ log2(rpkm_all_ycg1_exp_all_mean[,7]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(rpkm_all_ycg1_exp_all_mean[,7]), (log2(rpkm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_ycg1_exp_all_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) exp_allerimental-exp_allerimental comparison
  plot(log2(rpkm_all_ycg1_exp_all_mean[,7]), log2(rpkm_all_ycg1_exp_all_mean[,i]), pch=20, main="S. cerevisiae exp_allerimental", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(rpkm_all_ycg1_exp_all_mean[,i]) ~ log2(rpkm_all_ycg1_exp_all_mean[,7]))
  
  #(4) Correct the spike-in data
  plot(log2(rpkm_all_ycg1_exp_all_mean[,7]), (log2(rpkm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_ycg1_exp_all_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  rpkm_all_ycg1_exp_all_mean_corrected_log2[,i]<-((log2(rpkm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
rpkm_all_ycg1_exp_all_mean_corrected_log2[,7]<-(log2(rpkm_all_ycg1_exp_all_mean[,7]))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in (bottom_mean+1):top_mean) {
  print(i)
  exp_alleriment_name<-colnames(rpkm_all_ycg1_exp_all_mean)[i]
  plot(log2(rpkm_all_ycg1_exp_all_mean[,(bottom_mean)]), log2(rpkm_all_ycg1_exp_all_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(rpkm_all_ycg1_exp_all_mean[,i]) ~ log2(rpkm_all_ycg1_exp_all_mean[,(bottom_mean)]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(rpkm_all_ycg1_exp_all_mean[,(bottom_mean)]), (log2(rpkm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_ycg1_exp_all_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) exp_allerimental-exp_allerimental comparison
  plot(log2(rpkm_all_ycg1_exp_all_mean[,(bottom_mean)]), log2(rpkm_all_ycg1_exp_all_mean[,i]), pch=20, main="S. cerevisiae exp_allerimental", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(rpkm_all_ycg1_exp_all_mean[,i]) ~ log2(rpkm_all_ycg1_exp_all_mean[,(bottom_mean)]))
  
  #(4) Correct the spike-in data
  plot(log2(rpkm_all_ycg1_exp_all_mean[,(bottom_mean)]), (log2(rpkm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_ycg1_exp_all_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  rpkm_all_ycg1_exp_all_mean_corrected_log2[,i]<-((log2(rpkm_all_ycg1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
rpkm_all_ycg1_exp_all_mean_corrected_log2[,bottom_mean]<-(log2(rpkm_all_ycg1_exp_all_mean[,bottom_mean]))

rpkm_all_ycg1_exp_all_mean_corrected<-rpkm_all_ycg1_exp_all_mean_corrected_log2
rpkm_all_ycg1_exp_all_mean_corrected[,7:top_mean]<-2^(rpkm_all_ycg1_exp_all_mean_corrected_log2[,7:top_mean])


#cpm_all_brn1_exp_all_mean
cpm_all_brn1_exp_all_mean_corrected_log2<-cpm_all_brn1_exp_all_mean
bottom_mean<-length(colnames(cpm_all_brn1_mean))-3
top_mean<-length(colnames(cpm_all_brn1_mean))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in 8:(bottom_mean-1)) {
  print(i)
  exp_alleriment_name<-colnames(cpm_all_brn1_exp_all_mean)[i]
  plot(log2(cpm_all_brn1_exp_all_mean[,7]), log2(cpm_all_brn1_exp_all_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(cpm_all_brn1_exp_all_mean[,i]) ~ log2(cpm_all_brn1_exp_all_mean[,7]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(cpm_all_brn1_exp_all_mean[,7]), (log2(cpm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_brn1_exp_all_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) exp_allerimental-exp_allerimental comparison
  plot(log2(cpm_all_brn1_exp_all_mean[,7]), log2(cpm_all_brn1_exp_all_mean[,i]), pch=20, main="S. cerevisiae exp_allerimental", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(cpm_all_brn1_exp_all_mean[,i]) ~ log2(cpm_all_brn1_exp_all_mean[,7]))
  
  #(4) Correct the spike-in data
  plot(log2(cpm_all_brn1_exp_all_mean[,7]), (log2(cpm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_brn1_exp_all_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  cpm_all_brn1_exp_all_mean_corrected_log2[,i]<-((log2(cpm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
cpm_all_brn1_exp_all_mean_corrected_log2[,7]<-(log2(cpm_all_brn1_exp_all_mean[,7]))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in (bottom_mean+1):top_mean) {
  print(i)
  exp_alleriment_name<-colnames(cpm_all_brn1_exp_all_mean)[i]
  plot(log2(cpm_all_brn1_exp_all_mean[,(bottom_mean)]), log2(cpm_all_brn1_exp_all_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(cpm_all_brn1_exp_all_mean[,i]) ~ log2(cpm_all_brn1_exp_all_mean[,(bottom_mean)]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(cpm_all_brn1_exp_all_mean[,(bottom_mean)]), (log2(cpm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_brn1_exp_all_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) exp_allerimental-exp_allerimental comparison
  plot(log2(cpm_all_brn1_exp_all_mean[,(bottom_mean)]), log2(cpm_all_brn1_exp_all_mean[,i]), pch=20, main="S. cerevisiae exp_allerimental", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,' CPM)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(cpm_all_brn1_exp_all_mean[,i]) ~ log2(cpm_all_brn1_exp_all_mean[,(bottom_mean)]))
  
  #(4) Correct the spike-in data
  plot(log2(cpm_all_brn1_exp_all_mean[,(bottom_mean)]), (log2(cpm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt CPM)', ylab=paste0('log2(',exp_alleriment_name,'Corrected CPM)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(cpm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(cpm_all_brn1_exp_all_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  cpm_all_brn1_exp_all_mean_corrected_log2[,i]<-((log2(cpm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
cpm_all_brn1_exp_all_mean_corrected_log2[,bottom_mean]<-(log2(cpm_all_brn1_exp_all_mean[,bottom_mean]))

cpm_all_brn1_exp_all_mean_corrected<-cpm_all_brn1_exp_all_mean_corrected_log2
cpm_all_brn1_exp_all_mean_corrected[,7:top_mean]<-2^(cpm_all_brn1_exp_all_mean_corrected_log2[,7:top_mean])


#tpm_all_brn1_exp_all_mean
tpm_all_brn1_exp_all_mean_corrected_log2<-tpm_all_brn1_exp_all_mean
bottom_mean<-length(colnames(tpm_all_brn1_mean))-3
top_mean<-length(colnames(tpm_all_brn1_mean))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in 8:(bottom_mean-1)) {
  print(i)
  exp_alleriment_name<-colnames(tpm_all_brn1_exp_all_mean)[i]
  plot(log2(tpm_all_brn1_exp_all_mean[,7]), log2(tpm_all_brn1_exp_all_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(tpm_all_brn1_exp_all_mean[,i]) ~ log2(tpm_all_brn1_exp_all_mean[,7]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(tpm_all_brn1_exp_all_mean[,7]), (log2(tpm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_brn1_exp_all_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) exp_allerimental-exp_allerimental comparison
  plot(log2(tpm_all_brn1_exp_all_mean[,7]), log2(tpm_all_brn1_exp_all_mean[,i]), pch=20, main="S. cerevisiae exp_allerimental", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(tpm_all_brn1_exp_all_mean[,i]) ~ log2(tpm_all_brn1_exp_all_mean[,7]))
  
  #(4) Correct the spike-in data
  plot(log2(tpm_all_brn1_exp_all_mean[,7]), (log2(tpm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_brn1_exp_all_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  tpm_all_brn1_exp_all_mean_corrected_log2[,i]<-((log2(tpm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
tpm_all_brn1_exp_all_mean_corrected_log2[,7]<-(log2(tpm_all_brn1_exp_all_mean[,7]))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in (bottom_mean+1):top_mean) {
  print(i)
  exp_alleriment_name<-colnames(tpm_all_brn1_exp_all_mean)[i]
  plot(log2(tpm_all_brn1_exp_all_mean[,(bottom_mean)]), log2(tpm_all_brn1_exp_all_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(tpm_all_brn1_exp_all_mean[,i]) ~ log2(tpm_all_brn1_exp_all_mean[,(bottom_mean)]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(tpm_all_brn1_exp_all_mean[,(bottom_mean)]), (log2(tpm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_brn1_exp_all_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) exp_allerimental-exp_allerimental comparison
  plot(log2(tpm_all_brn1_exp_all_mean[,(bottom_mean)]), log2(tpm_all_brn1_exp_all_mean[,i]), pch=20, main="S. cerevisiae exp_allerimental", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,' tpm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(tpm_all_brn1_exp_all_mean[,i]) ~ log2(tpm_all_brn1_exp_all_mean[,(bottom_mean)]))
  
  #(4) Correct the spike-in data
  plot(log2(tpm_all_brn1_exp_all_mean[,(bottom_mean)]), (log2(tpm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt tpm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected tpm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(tpm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(tpm_all_brn1_exp_all_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  tpm_all_brn1_exp_all_mean_corrected_log2[,i]<-((log2(tpm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
tpm_all_brn1_exp_all_mean_corrected_log2[,bottom_mean]<-(log2(tpm_all_brn1_exp_all_mean[,bottom_mean]))

tpm_all_brn1_exp_all_mean_corrected<-tpm_all_brn1_exp_all_mean_corrected_log2
tpm_all_brn1_exp_all_mean_corrected[,7:top_mean]<-2^(tpm_all_brn1_exp_all_mean_corrected_log2[,7:top_mean])



#rpkm_all_brn1_exp_all_mean
rpkm_all_brn1_exp_all_mean_corrected_log2<-rpkm_all_brn1_exp_all_mean
bottom_mean<-length(colnames(rpkm_all_brn1_mean))-3
top_mean<-length(colnames(rpkm_all_brn1_mean))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in 8:(bottom_mean-1)) {
  print(i)
  exp_alleriment_name<-colnames(rpkm_all_brn1_exp_all_mean)[i]
  plot(log2(rpkm_all_brn1_exp_all_mean[,7]), log2(rpkm_all_brn1_exp_all_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(rpkm_all_brn1_exp_all_mean[,i]) ~ log2(rpkm_all_brn1_exp_all_mean[,7]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(rpkm_all_brn1_exp_all_mean[,7]), (log2(rpkm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_brn1_exp_all_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) exp_allerimental-exp_allerimental comparison
  plot(log2(rpkm_all_brn1_exp_all_mean[,7]), log2(rpkm_all_brn1_exp_all_mean[,i]), pch=20, main="S. cerevisiae exp_allerimental", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(rpkm_all_brn1_exp_all_mean[,i]) ~ log2(rpkm_all_brn1_exp_all_mean[,7]))
  
  #(4) Correct the spike-in data
  plot(log2(rpkm_all_brn1_exp_all_mean[,7]), (log2(rpkm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_brn1_exp_all_mean[,7]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  rpkm_all_brn1_exp_all_mean_corrected_log2[,i]<-((log2(rpkm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
rpkm_all_brn1_exp_all_mean_corrected_log2[,7]<-(log2(rpkm_all_brn1_exp_all_mean[,7]))

#(1) spike-spike comparison
par(mfrow=c(2,2))
for (i in (bottom_mean+1):top_mean) {
  print(i)
  exp_alleriment_name<-colnames(rpkm_all_brn1_exp_all_mean)[i]
  plot(log2(rpkm_all_brn1_exp_all_mean[,(bottom_mean)]), log2(rpkm_all_brn1_exp_all_mean[,i]), pch=20, main="S. pombe spike-in", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  #linear model
  r <-lm(log2(rpkm_all_brn1_exp_all_mean[,i]) ~ log2(rpkm_all_brn1_exp_all_mean[,(bottom_mean)]))
  abline(r, col="red")
  abline(0,1, col="gray")
  legend("topleft", paste("y=",round(r$coefficients[1],2), " + ",round(r$coefficients[2],2), "*x", sep=""))
  
  #(2) Correct the spike-in data
  plot(log2(rpkm_all_brn1_exp_all_mean[,(bottom_mean)]), (log2(rpkm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_brn1_exp_all_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  #(3) exp_allerimental-exp_allerimental comparison
  plot(log2(rpkm_all_brn1_exp_all_mean[,(bottom_mean)]), log2(rpkm_all_brn1_exp_all_mean[,i]), pch=20, main="S. cerevisiae exp_allerimental", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,' rpkm)'), xlim=c(0,15), ylim=c(0,15))
  abline(0,1, col="gray")
  abline(r, col="red")
  r <-lm(log2(rpkm_all_brn1_exp_all_mean[,i]) ~ log2(rpkm_all_brn1_exp_all_mean[,(bottom_mean)]))
  
  #(4) Correct the spike-in data
  plot(log2(rpkm_all_brn1_exp_all_mean[,(bottom_mean)]), (log2(rpkm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2], pch=20, main="S. pombe spike-in after correction", xlab='log2(wt rpkm)', ylab=paste0('log2(',exp_alleriment_name,'Corrected rpkm)'), xlim=c(0,15), ylim=c(0,15))
  # corrected linear model
  rc <-lm((log2(rpkm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2] ~ log2(rpkm_all_brn1_exp_all_mean[,(bottom_mean)]))
  abline(rc, col="red")
  abline(0,1, col="gray")
  
  
  #(5) corrected matrix
  rpkm_all_brn1_exp_all_mean_corrected_log2[,i]<-((log2(rpkm_all_brn1_exp_all_mean[,i])-r$coefficients[1])/r$coefficients[2])
}
rpkm_all_brn1_exp_all_mean_corrected_log2[,bottom_mean]<-(log2(rpkm_all_brn1_exp_all_mean[,bottom_mean]))

rpkm_all_brn1_exp_all_mean_corrected<-rpkm_all_brn1_exp_all_mean_corrected_log2
rpkm_all_brn1_exp_all_mean_corrected[,7:top_mean]<-2^(rpkm_all_brn1_exp_all_mean_corrected_log2[,7:top_mean])


