# RNAseq #
A collection of analyses of S. cerevisiae RNAseq data

## spike-in RNA mapping ##
A slurm script that maps S.cerevisieae RNA datasets with S. pombe spike-ins. Will produce sam files to throw up on IGV and featurecounts files for further anlysis for reads aligned to either S. cerevisiae genome or to S. pombe genome. Based on regular RNA analysis script from Luis Silva.  

## tRNA expression ##
An Rscript analysis of tRNA expression in S. cerevisiae, compared to spike-in control S. pombe. Read in featurecounts data files of each dataset following spike-in pipeline analysis. RNA seq data was produced with polyA depletion so should not select for tRNAs unless they are being degraded. It may give insight, however, into their transcripiton. 


