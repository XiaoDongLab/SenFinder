# senFinder
A svm trained model for predicting senescence on bulk RNA seq data.

Version: 1.0.0

Updated date: 2024.12.06

Citation: Shamsed Mahmud, Chen Zheng, Fernando Santiago, Lei Zhang, Paul D. Robbins and Xiao Dong. A machine learning approach identifies cellular senescence on transcriptome data of human cells in vitro. in review at Geroscience.

#####
## Author and License

Author: Shamsed Mahmud

Email: mahmu052@umn.edu (S.M.)

Licensed under the GNU Affero General Public License version 3 or later.

#####
## Dependencies
• R 4.4.1

• dplyr 1.1.4

• tidyr 1.3.1

• ggplot2 3.5.1

• lattice 0.22.6

• caret 6.0.94

• e1071 1.7.14



#####
## Introduction

(A) senFinder is a trained model with svm algorithm. For training we used all the 114 samples from the paper mentioned. The samples are from Fibroblast with MRC5 and HFF cell lines, keratinocyte and Melanocyte. The stressor are irradiation, bleomycin treatment, rotenone treatment, and replication.  

#####
## Usage

### Step 1. Prepare input file

• senFinder is supposed to take two input files, one is the senFinder.rds trained model and the other one is the user input file. The input file should be formated as a csv as the following:

| Sample_id | RRM2B | CPXM1 | PHGDH | PIK3IP1 | YIU2 | WASF1 | MRPS9 | TCF12 | DNAI3 | CAMK2N2 | NSD1 | MRPL11 | F8 | RPL10A | ISY1 | LUCAT1 | PURPL | CAPN10_DT Lnc_PLOD2_2 | LINC00663 |
|-----------|-------|-------|-------|---------|------|-------|-------|-------|-------|---------|------|--------|----|--------|------|--------|-------|----------------------|-----------|
| A         | 30.72 | 1.59  | 27    | 5.75    | 5.47 | 4.25  | 10.24 | 44.57 | 4.69  | 0.06    | 14.85 | 22.76  | 2.21 | 802.64 | 17.77 | 1.35   | 3.4   | 0.65                 | 0.51      |
| E         | 29.95 | 9.07  | 113.32| 20.58   | 12.85| 14.96 | 26.69 | 64.89 | 1     | 0.26    | 13.83 | 57.81  | 1.67| 1323.08| 34.63 | 0.87   | 0.66  | 0.24                 | 0.24      |


The file should not have NA, Inf, use 0 instead through basic cleaning process. Here, each row is a sample and column comprises of  Sample_id and twenty  gene found from mRMR and incremental feature selcections. Use this twenty gene's TPM value.

### Step 2. Create directory and run model 
• Then create a directory in your unix environment (like proj_sen).
```Bash
# Create a directory proj-sen
mkdir proj_sen 
# move to proj_sen directory
cd proj_sen
```
• Dowload the senFinder.R, senFinder.rds,test_file.csv files and put it in your directory(proj_sen).

• Then run senFinder.R like below
```Bash
Rscript senFinder.R "./senFinder.rds" "./test_file.csv"
```
• This will generate a new file alias predicted_results.csv. It has two columns. The first coulumn is Sample_id and the second column is Senescence_status. 'Yes' means predicted as senescence and 'No' means predicted as non-senescece.

#####
## Release Notes

• v1.0.0, 2024.12.06, 1st version uploades.





