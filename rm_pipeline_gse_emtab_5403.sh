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
refgenome=/data/gv_h/smahmud/proj_machine/all_files/hg38d/star-indices
gene=/data/gv_h/smahmud/proj_machine/all_files/hg38d/Homo_sapiens.GRCh38.109.gtf
genome=/data/gv_h/smahmud/proj_machine/all_files/hg38d/star-indices/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa
refname=/data/gv_h/smahmud/proj_machine/all_files/hg38d/star-indices/Homo_sapiens.GRCh38.dna_sm.rsem

np=8
sn=${1}
readends=$(seq 1 2)


## the start ###
echo "RS pipeline version 2023.03.27 by Shamsed Mahmud"; echo "Sample ID: "${sn}; echo "Processing starts on"; date; echo "==========";

#############
##### 1 #####
#############
cd ../

### 1.0 trim ###
report="1.0 trim"
echo ${report} "start on"; date; echo "====="
mkdir -p ./alignment/1-0-trim/
trim_galore --paired -o ./alignment/1-0-trim \
./fastq/${sn}_1.fastq \
./fastq/${sn}_2.fastq
mv ./alignment/1-0-trim/${sn}_1_val_1.fq ./alignment/1-0-trim/${sn}_1.trimmed.fastq
mv ./alignment/1-0-trim/${sn}_2_val_2.fq ./alignment/1-0-trim/${sn}_2.trimmed.fastq
echo "====="; echo ${report} "end on"; date

### 1.1 fastqc ###
report="1.1 fastqc"
echo ${report} "start on"; date; echo "====="
mkdir -p ./alignment/1-1-fastqc
for i in ${readends}; do
fastqc -o ./alignment/1-1-fastqc ./alignment/1-0-trim/${sn}_${i}.trimmed.fastq
done
echo "====="; echo ${report} "end on"; date

### 1.2 fastq_screen ###
report="1.2 fastq_screen"
echo ${report} "start on"; date; echo "====="
mkdir -p ./alignment/1-2-fastq_screen
for i in ${readends}; do
fastq_screen \
  --aligner bowtie2 \
  --conf /data/gv_h/smahmud/apps/FastQ-Screen-0.15.2/fastq_screen.conf \
  --outdir ./alignment/1-2-fastq_screen \
  ./alignment/1-0-trim/${sn}_${i}.trimmed.fastq
done
echo "====="; echo ${report} "end on"; date

### 1.3 alignment to human ref ###
report="1.3 alignment to human ref"
echo ${report} "start on"; date; echo "====="
mkdir -p ./alignment/1-3-human/
STAR \
  --runThreadN ${np} \
  --genomeDir ${refgenome} \
  --readFilesIn ./alignment/1-0-trim/${sn}_1.trimmed.fastq ./alignment/1-0-trim/${sn}_2.trimmed.fastq \
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

echo "====="; echo ${report} "end on"; date

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

echo "====="; echo ${report} "end on"; date

## RSEM ##
mkdir -p ./alignment/2-rsem/${sn}
report="2.2 RSEM"
echo ${report} "start on"; date; echo "====="
rsem-calculate-expression \
  --paired-end \
  --alignments \
  --single-cell-prior \
  -p ${np} \
  --bam \
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
fractionreads_r1_hs=$(grep 'Human' alignment/1-2-fastq_screen/${sn}_1.trimmed_screen.txt | sed -n '1p' | awk '{print $6}')
fractionreads_r2_hs=$(grep 'Human' alignment/1-2-fastq_screen/${sn}_2.trimmed_screen.txt | sed -n '1p' | awk '{print $6}')

# col 5: % reads in r1 aligned to Mouse ref genome
# col 6: % reads in r2 aligned to Mouse ref genome
fractionreads_r1_mm=$(grep 'Mouse' alignment/1-2-fastq_screen/${sn}_1.trimmed_screen.txt | sed -n '1p' | awk '{print $6}')
fractionreads_r2_mm=$(grep 'Mouse' alignment/1-2-fastq_screen/${sn}_2.trimmed_screen.txt | sed -n '1p' | awk '{print $6}')

# col 7: % reads in r1 aligned to Yeast ref genome
# col 8: % reads in r2 aligned to Yeast ref genome
fractionreads_r1_ye=$(grep 'Yeast' alignment/1-2-fastq_screen/${sn}_1.trimmed_screen.txt | sed -n '1p' | awk '{print $6}')
fractionreads_r2_ye=$(grep 'Yeast' alignment/1-2-fastq_screen/${sn}_2.trimmed_screen.txt | sed -n '1p' | awk '{print $6}')

# col 9: % reads in r1 aligned to Ecoli ref genome
# col 10: % reads in r2 aligned to Ecoli ref genome
fractionreads_r1_ec=$(grep 'Ecoli' alignment/1-2-fastq_screen/${sn}_1.trimmed_screen.txt | sed -n '1p' | awk '{print $6}')
fractionreads_r2_ec=$(grep 'Ecoli' alignment/1-2-fastq_screen/${sn}_2.trimmed_screen.txt | sed -n '1p' | awk '{print $6}')

# col 11: % reads in r1 aligned to PhiX ref genome
# col 12: % reads in r2 aligned to PhiX ref genome
fractionreads_r1_px=$(grep 'PhiX' alignment/1-2-fastq_screen/${sn}_1.trimmed_screen.txt | sed -n '1p' | awk '{print $6}')
fractionreads_r2_px=$(grep 'PhiX' alignment/1-2-fastq_screen/${sn}_2.trimmed_screen.txt | sed -n '1p' | awk '{print $6}')

# col 13: % reads in r1 aligned to Lambda ref genome
# col 14: % reads in r2 aligned to Lambda ref genome
fractionreads_r1_lm=$(grep 'Lambda' alignment/1-2-fastq_screen/${sn}_1.trimmed_screen.txt | sed -n '1p' | awk '{print $6}')
fractionreads_r2_lm=$(grep 'Lambda' alignment/1-2-fastq_screen/${sn}_2.trimmed_screen.txt | sed -n '1p' | awk '{print $6}')

# col 15: no. read pairs uniquely aligned to human ref
noreadpairs_alignuniq=$(grep "Uniquely mapped reads number" alignment/1-3-human/${sn}Log.final.out | awk '{print $6}' | awk '{gsub(/,/,""); print $0}')

# col 16: no. read pairs multiplely (not "too many") aligned to human ref
noreadpairs_alignmulti=$(grep "Number of reads mapped to multiple loci" alignment/1-3-human/${sn}Log.final.out | awk '{print $9}' | awk '{gsub(/,/,""); print $0}')

# col 17: fraction of PCR duplications identified (%)
pcrpcrdupremoved=$(awk -v sn=${sn} '$1==sn' alignment/1-3-human/${sn}.transcriptome.dedup_metrics.txt | awk '{print $9 * 100}')

# unzip -p alignment/1-3-human/${sn}.transcriptome.rsem_fastqc.zip ${sn}.transcriptome.rsem_fastqc/fastqc_data.txt | less

# write summary
file=processing_summary_gseuth.csv
if test -f "$file"; then
  echo "$file exists."
else
  echo "$file does not exist; create one."
  echo "SampleID,\
noreadpairs,\
fractionreads_r1_hs,fractionreads_r2_hs,\
fractionreads_r1_mm,fractionreads_r2_mm,\
fractionreads_r1_ye,fractionreads_r2_ye,\
fractionreads_r1_ec,fractionreads_r2_ec,\
fractionreads_r1_px,fractionreads_r2_px,\
fractionreads_r1_lm,fractionreads_r2_lm,\
noreadpairs_alignuniq,noreadpairs_alignmulti,\
pcrpcrdupremoved" \
> processing_summary_gseuth.csv
fi

echo ${sn},\
${noreadpairs},\
${fractionreads_r1_hs},${fractionreads_r2_hs},\
${fractionreads_r1_mm},${fractionreads_r2_mm},\
${fractionreads_r1_ye},${fractionreads_r2_ye},\
${fractionreads_r1_ec},${fractionreads_r2_ec},\
${fractionreads_r1_px},${fractionreads_r2_px},\
${fractionreads_r1_lm},${fractionreads_r2_lm},\
${noreadpairs_alignuniq},${noreadpairs_alignmulti},\
${pcrpcrdupremoved}\
>>processing_summary_gseuth.csv

### 5.2 clean ups ###
# collect useful output files: nothing need to be done now.
# remove intermediate files: nothing need to be done now.
echo "====="; echo ${report} "end on"; date

## the end ###
echo "=========="; echo "Processing ends on"; date; echo "Sample ID: "${sn}

