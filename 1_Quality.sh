#!/bin/sh

################## Used to assess initial quality of RNAseq data from canine osteosarcoma RNA sequencing data

source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load fastqc/0.10.1
module load gnu_parallel/201612222
module load gnu_parallel

#Define variables for directories
DATADIR=/scratch/OSAseq
DATACOPY=/scratch/OSAseq/Raw_Copy
OUTDIR=/home/aubrln001/OSAseq/Quality_Raw

#Make the output directory for results
mkdir -p $OUTDIR
mkdir -p $DATACOPY

cd $DATACOPY
cp $DATADIR/SL*.fastq.gz .

ls *_1.fastq.gz | time parallel -j+0 --eta 'fastqc {}'
ls *_2.fastq.gz | time parallel -j+0 --eta 'fastqc {}'



#Run fastqc to assess quality - output will be a results folder + zip folder + html link to open & assess quality
#fastqc $DATADIR/*.fastq.gz

#Copy the results to home directory for safe keeing
cp *fastqc* $OUTDIR


###### You will then need to copy over the html files over to home desktop to visualize in a web browser to determine appropriate trimming parameters


