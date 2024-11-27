###Pipeline created by Shamsed Mahmud at 10th April 2023###

#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -N rm-pipeline
#$ -pe smp 8
#$ -l h_vmem=8G
#$ -q all.q

export LC_ALL=C
export MALLOC_ARENA_MAX=4

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 a" >&2
  exit 1
fi

### ENVIRONMENT ###
picard=$DIR_picard/picard.jar
refgenome=/data/gv_h/smahmud/proj_machine/GSE56293/hg19d/star_hg19_grch37_50
gene=/data/gv_h/smahmud/proj_machine/GSE56293/hg19d/Homo_sapiens.GRCh37.87.gtf
genome=/data/gv_h/smahmud/proj_machine/GSE56293/hg19d/star_hg19_grch37_50/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa
refname=/data/gv_h/smahmud/proj_machine/GSE56293/hg19d/star_hg19_grch37_50/Homo_sapiens.GRCh37.87.ercc.rsem

np=8
sn=${1}

## the start ###
echo "RS pipeline version 2023.03.17 by Shamsed Mahmud"; echo "Sample ID: "${sn}; echo "Processing starts on"; date; echo "==========";

####################
#########1##########
####################
cd ../

### 1.0 trim ###

report="1.0 trim"
echo ${report} "start on"; date; echo "====="
mkdir -p ./alignment/1-0-trim/

trim_galore -o ./alignment/1-0-trim \
./fastq/${sn}.fastq

mv ./alignment/1-0-trim/${sn}_trimmed.fq ./alignment/1-0-trim/${sn}.trimmed.fastq

echo "====="; echo ${report} "end on"; date

### 1.1 fastqc ###
report="1.1 fastqc"
echo ${report} "start on"; date; echo "====="
mkdir -p ./alignment/1-1-fastqc

fastqc -o ./alignment/1-1-fastqc ./alignment/1-0-trim/${sn}.trimmed.fastq

echo "====="; echo ${report} "end on"; date

### 1.2 fastq_screen ###
report="1.2 fastq_screen"
echo ${report} "start on"; date; echo "====="
mkdir -p ./alignment/1-2-fastq_screen

fastq_screen \
  --aligner bowtie2 \
  --conf /data/gv_h/smahmud/apps/FastQ-Screen-0.15.2/fastq_screen.conf \
  --outdir ./alignment/1-2-fastq_screen \
  ./alignment/1-0-trim/${sn}.trimmed.fastq

echo "====="; echo ${report} "end on"; date

### 1.3 alignment to human ref ###
report="1.3 alignment to human ref"
echo ${report} "start on"; date; echo "====="
mkdir -p ./alignment/1-3-human/
STAR \
  --runThreadN ${np} \
  --genomeDir ${refgenome} \
  --readFilesIn ./fastq/${sn}.fastq \
  --outFileNamePrefix ./alignment/1-3-human/${sn} \
  --outTmpDir ./alignment/1-3-human/${sn} \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode TranscriptomeSAM

## cleanup STAR transcriptome alignment ##
mkdir -p tmp/
java17 -Xmx20g -jar $picard AddOrReplaceReadGroups \
  VALIDATION_STRINGENCY=SILENT \
  TMP_DIR=`pwd`/tmp \
  I=./alignment/1-3-human/${sn}Aligned.toTranscriptome.out.bam \
  O=./alignment/1-3-human/${sn}.transcriptome.rd.bam \
  RGID=${sn} RGLB=${sn} RGPL=illumina RGPU=${sn} RGSM=${sn}

java17 -Xmx20g -jar $picard SortSam \
  VALIDATION_STRINGENCY=SILENT \
  TMP_DIR=`pwd`/tmp \
  I=./alignment/1-3-human/${sn}.transcriptome.rd.bam \
  O=./alignment/1-3-human/${sn}.transcriptome.sort.bam \
  SO=coordinate \
  CREATE_INDEX=true

# Remve duplications #
java17 -Xmx20g -jar $picard MarkDuplicates \
  TMP_DIR=`pwd`/tmp \
  I=./alignment/1-3-human/${sn}.transcriptome.sort.bam \
  O=./alignment/1-3-human/${sn}.transcriptome.star.bam \
  METRICS_FILE=./alignment/1-3-human/${sn}.transcriptome.dedup_metrics.txt \
  CREATE_INDEX=true \
  REMOVE_DUPLICATES=true

rm ./alignment/1-3-human/${sn}Aligned.toTranscriptome.out.bam \
  ./alignment/1-3-human/${sn}.transcriptome.rd.bam \
  ./alignment/1-3-human/${sn}.transcriptome.sort.ba?

#############
##### 2 #####
#############
report="2.1 convert-sam-for-rsem"
echo ${report} "start on"; date; echo "====="
convert-sam-for-rsem \
  -p ${np} \
  --memory-per-thread 6G \
  ./alignment/1-3-human/${sn}.transcriptome.star.bam \
  ./alignment/1-3-human/${sn}.transcriptome.rsem

fastqc -o ./alignment/1-3-human/ ./alignment/1-3-human/${sn}.transcriptome.rsem.bam
echo "====="; echo ${report} "end on"; date

## RSEM ##
mkdir -p ./alignment/2-rsem/${sn}
report="2.2 RSEM"
echo ${report} "start on"; date; echo "====="
rsem-calculate-expression \
  --no-bam-output \
  --alignments \
  --bam \
  -p ${np} \
  ./alignment/1-3-human/${sn}.transcriptome.rsem.bam \
  ${refname} \
  ./alignment/2-rsem/${sn}
echo "====="; echo ${report} "end on"; date

#####################
##### 3 summary #####
#####################
report="3 summary"
echo ${report} "start on"; date; echo "====="

### 3.1 Summary of the statistics ###

# col 1: sample name

# col 2: no. read pairs in r1 & r2
noreadpairs=$(grep "Number of input reads" alignment/1-3-human/${sn}Log.final.out | awk '{print $6}' | awk '{gsub(/,/,""); print $0}')

# col 3: % reads in r1 aligned to Human ref genome
# col 4: % reads in r2 aligned to Human ref genome
fractionreads_r1_hs=$(grep 'Human' alignment/1-2-fastq_screen/${sn}.trimmed_screen.txt | sed -n '1p' | awk '{print $6}')

# col 5: % reads in r1 aligned to Mouse ref genome
# col 6: % reads in r2 aligned to Mouse ref genome
fractionreads_r1_mm=$(grep 'Mouse' alignment/1-2-fastq_screen/${sn}.trimmed_screen.txt | sed -n '1p' | awk '{print $6}')


# col 7: no. read pairs uniquely aligned to human ref
noreadpairs_alignuniq=$(grep "Uniquely mapped reads number" alignment/1-3-human/${sn}Log.final.out | awk '{print $6}' | awk '{gsub(/,/,""); print $0}')

# col 8: no. read pairs multiplely (not "too many") aligned to human ref
noreadpairs_alignmulti=$(grep "Number of reads mapped to multiple loci" alignment/1-3-human/${sn}Log.final.out | awk '{print $9}' | awk '{gsub(/,/,""); print $0}')

# col 9: fraction of PCR duplications identified (%)
pcrpcrdupremoved=$(awk -v sn=${sn} '$1==sn' alignment/1-3-human/${sn}.transcriptome.dedup_metrics.txt | awk '{print $9 * 100}')

# write summary
file=processing_summary_gse56293.csv
if test -f "$file"; then
  echo "$file exists."
else
  echo "$file does not exist; create one."
  echo "SampleID, \
noreadpairs, \
fractionreads_r1_hs, \
fractionreads_r1_mm, \
noreadpairs_alignuniq,noreadpairs_alignmulti, \
pcrpcrdupremoved" \
> processing_summary_gse56293.csv
fi

echo ${sn},\
${noreadpairs},\
${fractionreads_r1_hs}, \
${fractionreads_r1_mm}, \
${noreadpairs_alignuniq},${noreadpairs_alignmulti},\
${pcrpcrdupremoved}\
>>processing_summary_gse56293.csv

### 5.2 clean ups ###
# collect useful output files: nothing need to be done now.
# remove intermediate files: nothing need to be done now.
echo "====="; echo ${report} "end on"; date

## the end ###
echo "=========="; echo "Processing ends on"; date; echo "Sample ID: "${sn}