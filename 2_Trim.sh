#!/bin/sh


####################### Used to trim the RNAseq data so only good quality is used for analysis

#source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load fastqc/0.10.1
module load trimmomatic
module load gnu_parallel/201612222
module load gnu_parallel

#Define variables for directories
DATADIR=/scratch/OSAseq/RawData
WORKDIR=/scratch/OSAseq/Trimmed3
OUTDIR=/home/aubrln001/OSAseq/Trimmed3

#Make directories
##mkdir -p $WORKDIR
#mkdir -p $OUTDIR

#Change to data directory
#cd $DATADIR

#Make list of file names to trim
#ls | grep ".fastq.gz" | cut -d "_" -f 1 | sort | uniq > list

#Move list to working directory
#mv list $WORKDIR

#Change to working directory
cd $WORKDIR

#Run trimmomatic on list
while read i
do
java -jar /mnt/beegfs/home/aubmxa/.conda/envs/BioInfo_Tools/share/trimmomatic-0.39-1/trimmomatic.jar PE -threads 7 -phred33 \
	$DATADIR/"$i"_1.fastq.gz $DATADIR/"$i"_2.fastq.gz \
	"$i"_1_paired.fastq.gz "$i"_1_unpaired.fastq.gz \
	"$i"_2_paired.fastq.gz "$i"_2_unpaired.fastq.gz \
	ILLUMINACLIP:AdaptersToTrim_All.fa:2:35:10 HEADCROP:1 MINLEN:36
done<list

#Assess quality on cleaned data with fastqc in parallel
ls *_1_paired.fastq.gz | parallel -j+0 --eta 'fastqc {}'
ls *_2_paired.fastq.gz | parallel -j+0 --eta 'fastqc {}'

#Copy results to home directory for safe keeping
cp *fastqc $OUTDIR

