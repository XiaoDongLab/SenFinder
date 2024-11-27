#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd
#$ -N rm-pipeline
#$ -pe smp 10
#$ -l h_vmem=20G
#$ -q all.q

Rscript mrmr_500_run_new.R
