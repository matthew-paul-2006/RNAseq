###################################
# tRNA expression in RNA-seq data #
###################################

#Date
#18.03.15

#An analysis of tRNA expression in S. cerevisiae, compared to spike-in control S. pombe. Read in featurecounts data files of each dataset 
#following spike-in pipeline analysis. RNA seq data was produced with polyA depletion so should not select for tRNAs unless they are being 
#degraded. It may give insight, however, into their transcripiton. 

#Load in dplyr to help manipulate matrices
if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}

#Load in rlist to help manipulate lists
if (!require("rlist")) {
  install.packages("rlist", dependencies = TRUE)
  library(rlist)
}

#Load in gplots for graphing
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}

#Load viridis color palette for graphing
if (!require("viridis")) {
  install.packages("viridis")
  library(viridis)
}


### Read in tRNA count data ###

# Load count data and create table of all samples
path <- '~/Google_Drive/Lab/Data_Analysis/RNA/tRNAs/'

pathToFiles <- list(paste0(path, 'matt49_tRNA_featureCounts.txt'),
                    paste0(path, 'matt50_tRNA_featureCounts.txt'),
                    paste0(path, 'matt51_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp66_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp67_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp68_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp69_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp70_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp71_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp72_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp73_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp74_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp75_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp76_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp77_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp78_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp79_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp80_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp81_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp82_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp83_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp84_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp85_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp86_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp87_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp88_tRNA_featureCounts.txt'),
                    paste0(path, 'mrp89_tRNA_featureCounts.txt')
)

counts <- lapply(pathToFiles, read.table, header = T, stringsAsFactors=F)

#Assign the diffrent data set names
sampleNames <- list('matt49', 'matt50', 'matt51', 'mrp66','mrp67', 'mrp68', 'mrp69', 'mrp70','mrp71', 'mrp72', 'mrp73', 'mrp74','mrp75', 'mrp76', 'mrp77', 'mrp78','mrp79', 'mrp80', 'mrp81', 'mrp82', 'mrp83', 'mrp84','mrp85', 'mrp86', 'mrp87', 'mrp88','mrp89')
for(i in 1:length(counts)){
  colnames(counts[[i]])[7] <- sampleNames[i]
}

#merge my datasets together
counts_tRNA <- Reduce(dplyr::inner_join, counts)
head(counts_tRNA)

#Heatmap of raw counts
heatmap.2(data.matrix(counts_tRNA[,7:33]),dendrogram="none", Rowv=NA, Colv=NA,trace="none",density.info="none", col=viridis(100), main='Heatmap of tRNA counts by \n experimental sample')

### Read in tRNA spike count data ###

## Load count data and create table of all samples
path <- '~/Google_Drive/Lab/Data_Analysis/RNA/tRNAs/'

pathToFiles <- list(paste0(path, 'matt49_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'matt50_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'matt51_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp66_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp67_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp68_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp69_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp70_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp71_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp72_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp73_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp74_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp75_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp76_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp77_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp78_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp79_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp80_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp81_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp82_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp83_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp84_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp85_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp86_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp87_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp88_tRNA_spike_featureCounts.txt'),
                    paste0(path, 'mrp89_tRNA_spike_featureCounts.txt')
)

counts <- lapply(pathToFiles, read.table, header = T, stringsAsFactors=F)

#Assign the diffrent data set names
sampleNames <- list('matt49', 'matt50', 'matt51', 'mrp66','mrp67', 'mrp68', 'mrp69', 'mrp70','mrp71', 'mrp72', 'mrp73', 'mrp74','mrp75', 'mrp76', 'mrp77', 'mrp78','mrp79', 'mrp80', 'mrp81', 'mrp82', 'mrp83', 'mrp84','mrp85', 'mrp86', 'mrp87', 'mrp88','mrp89')
for(i in 1:length(counts)){
  colnames(counts[[i]])[7] <- sampleNames[i]
}

#merge my datasets together
counts_tRNA_spike <- Reduce(dplyr::inner_join, counts)
head(counts_tRNA_spike)

#Heatmap of raw counts
heatmap.2(data.matrix(counts_tRNA_spike[,7:33]),dendrogram="none", Rowv=NA, Colv=NA,trace="none",density.info="none", col=viridis(100), main='Heatmap of tRNA_spike counts \n by experimental sample')

### Read in CDS counts data ###
# tRNA expression is low. To ba able to efficiently normalise the data we will need to normalise
# using the whole genome

## Load count data and create table of all samples
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

#Assign the diffrent data set names
sampleNames <- list('matt49', 'matt50', 'matt51', 'mrp66','mrp67', 'mrp68', 'mrp69', 'mrp70','mrp71', 'mrp72', 'mrp73', 'mrp74','mrp75', 'mrp76', 'mrp77', 'mrp78','mrp79', 'mrp80', 'mrp81', 'mrp82', 'mrp83', 'mrp84','mrp85', 'mrp86', 'mrp87', 'mrp88','mrp89')
for(i in 1:length(counts)){
  colnames(counts[[i]])[7] <- sampleNames[i]
}

#merge my datasets together
counts_all <- Reduce(dplyr::inner_join, counts)
head(counts_all)

#Heatmap of raw counts
#heatmap.2(data.matrix(counts_all[,7:33]),dendrogram="none", Rowv=NA, Colv=NA,trace="none",density.info="none", col=viridis(100), main='Heatmap of counts \n by experimental sample')

### Read in CDS counts spike data ###

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

#Assign the diffrent data set names
sampleNames <- list('matt49', 'matt50', 'matt51', 'mrp66','mrp67', 'mrp68', 'mrp69', 'mrp70','mrp71', 'mrp72', 'mrp73', 'mrp74','mrp75', 'mrp76', 'mrp77', 'mrp78','mrp79', 'mrp80', 'mrp81', 'mrp82', 'mrp83', 'mrp84','mrp85', 'mrp86', 'mrp87', 'mrp88','mrp89')
for(i in 1:length(counts)){
  colnames(counts[[i]])[7] <- sampleNames[i]
}

#merge my datasets together
counts_all_spike <- Reduce(dplyr::inner_join, counts)
head(counts_all_spike)

#Heatmap of raw counts
#heatmap.2(data.matrix(counts_all_spike[,7:33]),dendrogram="none", Rowv=NA, Colv=NA,trace="none",density.info="none", col=viridis(100), main='Heatmap of spike counts \n by experimental sample')

### Merge tRNA and CDS datasets

spike_both_counts<-rbind(counts_all_spike, counts_tRNA_spike)
both_counts<-rbind(counts_all, counts_tRNA)


### Calculate CPMs ###
  
#Use the existing file as a template for for inputting cpms
cpm_all_spike<-spike_both_counts
cpm_all<-both_counts


#Calculate CPM for every for each, step by step
for (i in 7:33){
  #Do it for spike
  total_counts<-sum(spike_both_counts[,i])  
  cpm_counts<-(spike_both_counts[,i]/total_counts)*10^6
  cpm_all_spike[,i]<-cpm_counts
  #Do it for experimental
  total_counts<-sum(both_counts[,i])  
  cpm_counts<-(both_counts[,i]/total_counts)*10^6
  cpm_all[,i]<-cpm_counts
  print(i)
}


### Which are expressed ###

#Expressed in all conditions
expressed_all<-subset(cpm_all, matt49 >= 1 & matt50 >= 1 & matt51 >= 1 & mrp66 >= 1 & mrp67 >= 1 & mrp68 >= 1 & mrp69 >= 1 & mrp70 >= 1 & mrp71 >= 1 & mrp72 >= 1 & mrp73 >= 1 & mrp74 >= 1 &mrp75 >= 1 & mrp76 >= 1 & mrp77 >= 1 & mrp78 >= 1 & mrp79 >= 1 & mrp80 >= 1 & mrp81 >= 1 & mrp82 >= 1 & mrp83 >= 1 & mrp84 >= 1 & mrp85 >= 1 & mrp86 >= 1 & mrp87 >= 1 & mrp88 >= 1 & mrp89 >= 1)
expressed_spike<-subset(cpm_all_spike, matt49 >= 1 & matt50 >= 1 & matt51 >= 1 & mrp66 >= 1 & mrp67 >= 1 & mrp68 >= 1 & mrp69 >= 1 & mrp70 >= 1 & mrp71 >= 1 & mrp72 >= 1 & mrp73 >= 1 & mrp74 >= 1 &mrp75 >= 1 & mrp76 >= 1 & mrp77 >= 1 & mrp78 >= 1 & mrp79 >= 1 & mrp80 >= 1 & mrp81 >= 1 & mrp82 >= 1 & mrp83 >= 1 & mrp84 >= 1 & mrp85 >= 1 & mrp86 >= 1 & mrp87 >= 1 & mrp88 >= 1 & mrp89 >= 1)
#Expressed in one condition
expressed_once_all<-subset(cpm_all, matt49 != 0 | matt50 != 0 | matt51 != 0 | mrp66 != 0 | mrp67 != 0 | mrp68 != 0 | mrp69 != 0 | mrp70 != 0 | mrp71 != 0 | mrp72 != 0 | mrp73 != 0 | mrp74 != 0 |mrp75 != 0 | mrp76 != 0 | mrp77 != 0 | mrp78 != 0 | mrp79 != 0 | mrp80 != 0 | mrp81 != 0 | mrp82 != 0 | mrp83 != 0 | mrp84 != 0 | mrp85 != 0 | mrp86 != 0 | mrp87 != 0 | mrp88 != 0 | mrp89 != 0)
expressed_once_spike<-subset(cpm_all_spike, matt49 != 0 | matt50 != 0 | matt51 != 0 | mrp66 != 0 | mrp67 != 0 | mrp68 != 0 | mrp69 != 0 | mrp70 != 0 | mrp71 != 0 | mrp72 != 0 | mrp73 != 0 | mrp74 != 0 |mrp75 != 0 | mrp76 != 0 | mrp77 != 0 | mrp78 != 0 | mrp79 != 0 | mrp80 != 0 | mrp81 != 0 | mrp82 != 0 | mrp83 != 0 | mrp84 != 0 | mrp85 != 0 | mrp86 != 0 | mrp87 != 0 | mrp88 != 0 | mrp89 != 0)

### Graph the tRNAs ###

#Extract tRNAs from CPM lists
tRNAs<-cpm_all[grep('t',cpm_all[,1]),]
tRNAs_spike<-cpm_all_spike[grep('TRNA',cpm_all_spike[,1]),]

#Straight plot of all tRNAs
par(mar=c(6,4,4,2))
pdf('~/Google_Drive/Lab/Data_Analysis/RNA/tRNAs/tRNAexpressionpersample.pdf')
boxplot(tRNAs[,7:33], las=2, outline =F, main='tRNA expression per sample')
dev.off()
pdf('~/Google_Drive/Lab/Data_Analysis/RNA/tRNAs/tRNAexpressionpersample_spike.pdf')
boxplot(tRNAs_spike[,7:33], las=2, outline =F, main='tRNA expression per sample')
dev.off()

#Limited to tRNAs where there is expression in at least one condition
tRNAs_once<-expressed_once_all[grep('t',expressed_once_all[,1]),]
tRNAs_spike_once<-expressed_once_spike[grep('TRNA',expressed_once_spike[,1]),]

#Plot of all tRNAs where there is expression in at least one condition
par(mar=c(6,4,4,2))
pdf('~/Google_Drive/Lab/Data_Analysis/RNA/tRNAs/tRNAexpressionpersample_oneepxressed.pdf')
boxplot(tRNAs_once[,7:33], las=2, outline =F, main='tRNA expression in tRNAs with at  \n least one condition above 0 CPM per sample')
dev.off()
pdf('~/Google_Drive/Lab/Data_Analysis/RNA/tRNAs/SpiketRNAexpressionpersample_oneepxressed.pdf')
boxplot(tRNAs_spike_once[,7:33], las=2, outline =F, main='Spike tRNA expression in tRNAs with at  \n least one condition above 0 CPM per sample')
dev.off()

#Take means of each tRNA where there is expression in at least one condition
ycs4_wt_23<-apply(cbind(tRNAs_once[,10],tRNAs_once[,14]),1,mean)
ycs4_wt_37<-apply(cbind(tRNAs_once[,7],tRNAs_once[,11],tRNAs_once[,15]),1,mean)
ycs4_ycs42_23<-apply(cbind(tRNAs_once[,8],tRNAs_once[,12],tRNAs_once[,16]),1,mean)
ycs4_ycs42_37<-apply(cbind(tRNAs_once[,9],tRNAs_once[,13],tRNAs_once[,17]),1,mean)

ycg1_wt_23<-apply(cbind(tRNAs_once[,26],tRNAs_once[,30]),1,mean)
ycg1_wt_37<-apply(cbind(tRNAs_once[,27],tRNAs_once[,31]),1,mean)
ycg1_ycg12_23<-apply(cbind(tRNAs_once[,28],tRNAs_once[,32]),1,mean)
ycg1_ycg12_37<-apply(cbind(tRNAs_once[,29],tRNAs_once[,33]),1,mean)

brn1_wt_23<-apply(cbind(tRNAs_once[,18]),1,mean)
brn1_wt_37<-apply(cbind(tRNAs_once[,19],tRNAs_once[,22],tRNAs_once[,24]),1,mean)
brn1_aaway_23<-apply(cbind(tRNAs_once[,20]),1,mean)
brn1_aaway_37<-apply(cbind(tRNAs_once[,21],tRNAs_once[,23],tRNAs_once[,25]),1,mean)

#Plot means of each tRNA where there is expression in at least one condition
par(mar=c(9,4,4,2))
pdf('~/Google_Drive/Lab/Data_Analysis/RNA/tRNAs/tRNAexpression_oneepxressed.pdf')
boxplot(cbind(ycs4_wt_23,ycs4_wt_37,ycs4_ycs42_23,ycs4_ycs42_37,ycg1_wt_23, ycg1_wt_37,ycg1_ycg12_23,ycg1_ycg12_37,brn1_wt_23,brn1_wt_37,brn1_aaway_23,brn1_aaway_37  ),las=2, main='tRNA expression in tRNAs with at  \n least one condition above 0 CPM',ylab='CPMs', outline =F)
dev.off()

#Take means of each spike in tRNA where there is expression in at least one condition
ycs4_wt_23spike<-apply(cbind(tRNAs_spike_once[,10],tRNAs_spike_once[,14]),1,mean)
ycs4_wt_37spike<-apply(cbind(tRNAs_spike_once[,7],tRNAs_spike_once[,11],tRNAs_spike_once[,15]),1,mean)
ycs4_ycs42_23spike<-apply(cbind(tRNAs_spike_once[,8],tRNAs_spike_once[,12],tRNAs_spike_once[,16]),1,mean)
ycs4_ycs42_37spike<-apply(cbind(tRNAs_spike_once[,9],tRNAs_spike_once[,13],tRNAs_spike_once[,17]),1,mean)

ycg1_wt_23spike<-apply(cbind(tRNAs_spike_once[,26],tRNAs_spike_once[,30]),1,mean)
ycg1_wt_37spike<-apply(cbind(tRNAs_spike_once[,27],tRNAs_spike_once[,31]),1,mean)
ycg1_ycg12_23spike<-apply(cbind(tRNAs_spike_once[,28],tRNAs_spike_once[,32]),1,mean)
ycg1_ycg12_37spike<-apply(cbind(tRNAs_spike_once[,29],tRNAs_spike_once[,33]),1,mean)

brn1_wt_23spike<-apply(cbind(tRNAs_spike_once[,18]),1,mean)
brn1_wt_37spike<-apply(cbind(tRNAs_spike_once[,19],tRNAs_spike_once[,22],tRNAs_spike_once[,24]),1,mean)
brn1_aaway_23spike<-apply(cbind(tRNAs_spike_once[,20]),1,mean)
brn1_aaway_37spike<-apply(cbind(tRNAs_spike_once[,21],tRNAs_spike_once[,23],tRNAs_spike_once[,25]),1,mean)

#Plot means of each spike in tRNA where there is expression in at least one condition
par(mar=c(9,4,4,2))
pdf('~/Google_Drive/Lab/Data_Analysis/RNA/tRNAs/SpiketRNAexpression_oneepxressed.pdf')
boxplot(cbind(ycs4_wt_23spike,ycs4_wt_37spike,ycs4_ycs42_23spike,ycs4_ycs42_37spike,ycg1_wt_23spike, ycg1_wt_37spike,ycg1_ycg12_23spike,ycg1_ycg12_37spike,brn1_wt_23spike,brn1_wt_37spike,brn1_aaway_23spike,brn1_aaway_37spike  ),las=2, main='Spike in tRNA expression in tRNAs with at  \n least one condition above 0 CPM',ylab='CPMs', outline =F)
dev.off()

