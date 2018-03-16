#!/bin/sh
#
#BATCH --verbose
#SBATCH --job-name=RNA
#SBATCH --output=/scratch/mrp420/reports/slurm_RNA_%j.out
#SBATCH --error=/scratch/mrp420/reports/slurm_RNA_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=40GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mrp420@nyu.edu

#####
#To do!!!
#Automate to take whole dataset and provide correction according to spike-in. 
#Take correction and feed into edgeR for calling differentially expressed genes
####

### spikein_RNA_mapping.sh ###
#A slurm script that maps S.cerevisieae RNA datasets with S. pombe spike-ins. Will produce sam files to throw up on IGV and featurecounts 
#files for further anlysis for reads aligned to either S. cerevisiae genome or to S. pombe genome. Based on regular RNA analysis script 
#from Luis Silva

#-----------------------#
# How to run script #
#-----------------------#

# Run as: sbatch --export arg1='fastq file goes here' ~/yeast/scripts/Fq2ReadCounts_SacCer3_Spomspik.sh


# Example: sbatch --export arg1='/scratch/mrp420/fastq_RNA/C8CR7ACXX_l04n01_mrp75.35100000057e9a.fastq' ~/yeast/scripts/Fq2ReadCounts_SacCer3_Spomspike.sh

#Start time and ID reporting
echo 'Start Time'
echo $(date +%x_%r)

echo "
JobID is $PBS_JOBID
"

#-----------------------#
# Generate variables    #
#-----------------------#

#FASTQ
TAG1=$arg1

#SEQ ID
string=${TAG1}
IFS=_ read -a a <<< "$string"
var1=${a[0]}_${a[1]}
var2=${string#${var1}_}
TAG6=$( echo ${var2}|awk -F '[_.]' '{print $2}')
TAG3=$(echo ${TAG6}'_tRNA')

#Report TAGs
echo ${TAG1} ${TAG3}

#Generate directories
mkdir /scratch/mrp420/RNA
mkdir /scratch/mrp420/RNA/${TAG3}
mkdir /scratch/mrp420/RNA/${TAG3}/spike
TAG2=/scratch/mrp420/RNA/${TAG3}/${TAG3}
TAG4=/scratch/mrp420/RNA/${TAG3}/spike

#--------------------------------------------#
# Map reads to reference genome with Tophat2 #
#--------------------------------------------#
cd /scratch/mrp420/fastq/
echo "
Map reads with Tophat2...

"$(date +%r)

module load tophat/intel/2.1.1
module load bowtie2/intel/2.2.9
module load samtools/intel/1.3.1

#Map cerevisiae reads (first filter out pombe mappable reads)
echo "
Experimental:
"$(date +%r)

tophat2 --no-novel-juncs --library-type fr-firststrand -o ${TAG2}_Spombe_prefilter/ -G /home/mrp420/yeast/genomes/schizosaccharomyces_pombe.chr.gff3 /home/mrp420/yeast/genomes/Spombe_RNA/pombe ${TAG1}

tophat2 --no-novel-juncs --library-type fr-firststrand -o ${TAG2}_Scer_mapped/ -G /home/mrp420/yeast/genomes/saccharomyces_cerevisiae.gff /home/mrp420/yeast/genomes/SacCer_RNA/cerev ${TAG2}_Spombe_prefilter/unmapped.bam

mv ${TAG2}_Scer_mapped/accepted_hits.bam ${TAG2}_mapped.bam 

#Map pombe reads (first filter out cerev mappable reads)
echo "
Spike-in:
"$(date +%r)

tophat2 --no-novel-juncs --library-type fr-firststrand -o ${TAG2}_Scer_prefilter/ -G /home/mrp420/yeast/genomes/saccharomyces_cerevisiae.gff /home/mrp420/yeast/genomes/SacCer_RNA/cerev ${TAG1}

tophat2 --no-novel-juncs --library-type fr-firststrand -o ${TAG2}_Spom_mapped/ -G /home/mrp420/yeast/genomes/schizosaccharomyces_pombe.chr.gff3 /home/mrp420/yeast/genomes/Spombe_RNA/pombe ${TAG2}_Scer_prefilter/unmapped.bam

mv ${TAG2}_Spom_mapped/accepted_hits.bam ${TAG2}_spike_mapped.bam

#-------------------#
# Filter Duplicates #
#-------------------#

#echo "
#Remove duplicates
#"$(date +%r)
#
#module load picard/2.8.2
#

samtools sort -o ${TAG2}_mapped_s.sam ${TAG2}_mapped.bam 
samtools sort -o ${TAG2}_spike_mapped_s.sam ${TAG2}_spike_mapped.bam 

#java -jar $PICARD_JAR MarkDuplicates REMOVE_DUPLICATES=true I=${TAG2}_mapped_s.sam O=${TAG2}_nodups.sam M=${TAG2}_dups_metrics.txt
#java -jar $PICARD_JAR MarkDuplicates REMOVE_DUPLICATES=true I=spike/${TAG2}_spike_mapped_s.sam O=spike/${TAG2}_spike_nodups.sam M=spike/${TAG2}_spike_dups_metrics.txt

#------------------------------------------------#
# Sort and index SAM/BAM files for IGV and HTSeq #
#------------------------------------------------#

echo "
Sort and index sam/bam files for IGV and  HTSeq... 
"$(date +%r)

#Experimental
samtools view -bS -o ${TAG2}_mapped.bam ${TAG2}_mapped_s.sam     # for IGV
samtools sort -o ${TAG2}_s.bam ${TAG2}_mapped.bam        # for IGV
samtools index ${TAG2}_s.bam                     # for IGV

#Spike-in
samtools view -bS -o ${TAG2}_spike_mapped.bam ${TAG2}_spike_mapped_s.sam   # for IGV
samtools sort -o ${TAG2}_spike_s.bam ${TAG2}_spike_mapped.bam       # for IGV
samtools index ${TAG2}_spike_s.bam                               # for IGV

#Experimental no dupliactes
#samtools view -bS -o ${TAG2}_nodups.bam ${TAG2}_nodups.sam     # for IGV
#samtools sort -o ${TAG2}_nodups_s.bam ${TAG2}_nodups.bam       # for IGV
#samtools index ${TAG2}_nodups_s.bam                     # for IGV

#samtools sort -o ${TAG2}_nodups_sn.sam -n ${TAG2}_nodups.sam   # for htseq-count

#Spike-in no duplicates
#samtools view -bS -o spike/${TAG2}_nodups_spike.bam spike/${TAG2}_spike_nodups.sam     # for IGV
#samtools sort -o spike/${TAG2}_nodups_spike_s.bam spike/${TAG2}_nodups_spike.bam       # for IGV
#samtools index spike/${TAG2}_nodups_spike_s.bam                               # for IGV

#samtools sort -o spike/${TAG2}_nodups_spike_sn.sam -n spike/${TAG2}_spike_nodups.sam   # for htseq-count

#-----------------#
# Read reporting  #
#-----------------#

echo "
Report read information into report file
"$(date +%r)
echo $(date +%r)
echo "Total reads" ${TAG2} > ${TAG2}_report
#wc -l  ${TAG2}.sam >> ${TAG2}_report
total=$( wc -l /scratch/mrp420/fastq/${TAG1}| awk -F ' ' '{print $1}') 
echo $((total / 4)) >> ${TAG2}_report
echo "Total mapped reads" ${TAG2} >> ${TAG2}_report
wc -l ${TAG2}_mapped_s.sam >> ${TAG2}_report
#echo "Total mapped reads - no duplicates" ${TAG2} >> ${TAG2}_report
#wc -l ${TAG2}_nodups.sam >> ${TAG2}_report
echo "Total mapped reads - spike" ${TAG2} >> ${TAG2}_report
wc -l ${TAG2}_spike_mapped_s.sam >> ${TAG2}_report
#echo "Total mapped reads - spike - no duplicates" ${TAG2} >> ${TAG2}_report
#wc -l spike/${TAG2}_spike_nodups.sam >> ${TAG2}_report

#---------------------------------#
# Count reads with feature counts #
#---------------------------------#

echo "
Count reads with feature counts... 
"$(date +%r)

module load subread/intel/1.5.1

featureCounts -s 2 -t CDS -g Name -a /home/mrp420/yeast/genomes/saccharomyces_cerevisiae.gff -o ${TAG2}_featureCounts.txt ${TAG2}_mapped.bam

featureCounts -s 2 -t CDS -g ID -a /home/mrp420/yeast/genomes/schizosaccharomyces_pombe.chr.gff3 -o ${TAG2}_spike_featureCounts.txt ${TAG2}_spike_mapped.bam

#----------------------------#
# Remove intermediate files  #
#----------------------------#
echo "
Remove intermediate files
"$(date +%r)
#rm *.sam
#rm spike/*.sam
#rm *unmapped*
#rm spike/*unmapped*

echo "
Done!
"
#End time
echo 'End Time'
echo $(date +%x_%r)

exit 0;

