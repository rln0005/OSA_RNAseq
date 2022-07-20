#!/bin/sh

################### Used to map cleam/trimmed RNAseq reads to the canine reference sequence
################### Suggested run parameters: medium, 6 cores, 10 hrs, 16gb

module load hisat2
module load stringtie/1.3.3
module load python/2.7.1
module load gcc/9.3.0
module load samtools
module load bcftools/1.2
module load gffread/
module load gffcompare/

#Set stack size to unlimited 
ulimit -s
#Turn echo on so all commands are echoed in the output log
set -x

#Define and create directories
DATADIR=/scratch/OSAseq/Trimmed3
REFDIR=/scratch/OSAseq/RefSeq
OUTDIR=/scratch/OSAseq/Mapped3
COUNTSDIR=/scratch/OSAseq/ballgown3
RESULTSDIR=/home/aubrln001/OSAseq/Counts3

#make the directories
mkdir -p $OUTDIR
mkdir -p $COUNTSDIR
mkdir -p $RESULTSDIR

#define the reference genome
REF=K9_ref

######## Map and count the data using HiSat2 and Stringtie
#move into the data directory
cd $DATADIR

##create list of fastq files to map
ls | grep ".fastq.gz" | cut -d "_" -f 1 | sort | uniq > list

#move to directory for mapping
cd $OUTDIR

#copy the list of unique IDs from the original files to map
cp $DATADIR/list .

while read i
do
hisat2 -p 7 --dta --phred33 \
	-x "$REFDIR"/"$REF"_index \
	-1 "$DATADIR"/"$i"_1_paired.fastq.gz -2 "$DATADIR"/"$i"_2_paired.fastq.gz \
	-S "$i".sam

#view: convert the sam file into a bam file -bS: bam is the binary format corresponding to the SAM text format
#sort: convert the bam file into a sorted bam file
samtools view -@ 6 -bS ${i}.sam > ${i}.bam
samtools sort -@ 6 "$i".bam "$i"_sorted

#index the bam and get stats
samtools flagstat "$i"_sorted.bam > "$i"_stats.txt

mkdir "$COUNTSDIR"/"$i"

count reads that are mapped to each gene, exon, transcript model
stringtie -p 7 -e -B -G "$REFDIR"/"$REF".gtf -o "$COUNTSDIR"/"$i"/"$i".gtf -l "$i" "$OUTDIR"/"$i"_sorted.bam

done<list


#copy results to home directory
cp *.txt $RESULTSDIR
#move up 2 directories
cd ..

#Use PrepDE.py to convert files in ballgown folder to a count matrix
python /scratch/OSAseq/PrepDE.py /scratch/OSAseq/ballgown

#copy results to home directory
cp *.csv $RESULTSDIR
