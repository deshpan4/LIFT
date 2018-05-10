# LIFT:
## LncRNA Identification and Function prediction Tool
This document provides technical description of the LIFT pipeline for identification, genomic annotation and function prediction of long non-coding RNAs (lncRNAs) based on time-series RNA-Seq data.
## Getting Started
These instructions will help you to run LIFT from linux environment.
The LIFT framework consists of two modules:
* LncRNA Identification Module (LIM)
* Function prediction Module (FPM)
LIM performs the following tasks:
* Feature Extraction
* Classification of coding and non-coding RNAs
* Optimization of the selected features
* Genomic annotation of the predicted lncRNAs
FPM performs following tasks:
* Generation of lncRNA-protein co-expression similarity (LPCS) matrix
* Generation of protein-protein interaction matrix
* Generation of protein-function matrix
* Function prediction of the lncRNAs
## Download
Use git clone:
```
git clone https://github.com/deshpan4/LIFT
```
## Installation
Download the folder using 'git clone' command. Make sure all the dependencies are installed before running the framework.
## Prerequisites
Following tools should be installed before executing this pipeline on the host system:
* NodeJS
* SAMtools
* Python2.7 and Python3
* Bowtie2
* [Index and annotation for A.thaliana](ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Arabidopsis_thaliana/Ensembl/TAIR10/Arabidopsis_thaliana_Ensembl_TAIR10.tar.gz)
* R version 3.3 or higher
### R packages
* BMRF (https://github.com/jwbargsten/bmrf)
* iRF (https://github.com/sumbose/iRF)
* AUC
### NodeJS installation
For installing NodeJS on Ubuntu-based linux distributions, follow the instructions as mentioned [here](https://nodejs.org/en/download/package-manager/#debian-and-ubuntu-based-linux-distributions).
For Enterprise linux and Fedora systems, follow the instructions as mentioned [here](https://nodejs.org/en/download/package-manager/#enterprise-linux-and-fedora)
For other systems, please follow the instructions mentioned [here](https://nodejs.org/en/download/package-manager/)
### Usage
LIFT consists of five main scripts:
* LIFT_extractFeatures.sh (Feature extraction)
* LIFT_LiRFFS.py (Feature optimization)
* LIFT_lncRNAPredict.py (Predict lncRNAs)
* LIFT_annotateLncRNAs.py (Annotate lncRNAs)
* LIFT function prediction scripts
## Feature Extraction
The script 'LIFT_extractFeatures.sh' extract features from FASTA file for identification of lncRNA transcripts. Following input files are needed for extraction of features:
* Target FASTA file
* Coding FASTA file
* Non-coding FASTA file
* Output filename
* Scripts directory

### Usage
```
bash LIFT_extractFeatures.sh -c <coding filename> -n <noncoding filename> -f <target filename for feature extraction in FASTA format> -scripts <path to scripts directory> -o <output filename>
```

## Feature Optimization
The script 'LIFT_LiRFFS.py' performs feature selection for selection of optimal features. Following input files and parameters are required:
* Training set feature matrix
* Validation set feature matrix
* Lower cutoff for lambda
* Upper cutoff for lambda
* Step size between Lambda values
* Tolerance level
* Output training set file with optimal features
* Output validation set file with optimal features

### Usage
```
python3 LIFT_LiRFFS.py --training <training feature set in CSV format> --test <test feature set in CSV format> --lambdaLower <input lower lambda value> --lambdaUpper <input upper lambda value> --lambdaStepSize <input lambda step size> --tolerance <input tolerance value> --outputTr <output training set filename with optimal features> --outputTe <output test set with optimal features>
```

## lncRNA Identification
The script 'LIFT_lncRNAPredict.py' predicts lncRNAs based on features extracted. The script requires balanced training set and test set features for accurate classification of lncRNAs. Therefore, following files are required:
* Training set feature matrix
* Test set feature matrix
* Test set FASTA file
* Output test set file with predicted features

### Usage
```
usage: LIFT_lncRNAPredict.py [-h] [-tr TRAINING] [-te TEST] [-tef TESTFASTA]
                             [-o OUTPUT] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -tr TRAINING, --training TRAINING
                        input training set feature matrix
  -te TEST, --test TEST
                        input test set feature matrix
  -tef TESTFASTA, --testfasta TESTFASTA
                        input test set FASTA file
  -o OUTPUT, --output OUTPUT
                        output test set prediction filename
```

## lncRNA genomic annotation (sub-classification)
The script 'LIFT_annotateLncRNAs.py' provides annotation and classification of lncRNAs based on their coordinates. It requires following files and parameters for annotation:
* Coding sequence and coordinates file having following columns in CSV format: sequence,chromosome,start,end,strand
* Non-coding sequence and coordinates file having following columns in CSV format: sequence,chromosome,start,end,strand
* Output filename

### Usage
```
python2.7 LIFT_annotateLncRNAs.py -c <coding CSV file> -n <non-coding lncRNA CSV file> -o <output filename1> -o1 <output filename2>
usage: annotateLncRNA_New1.py [-h] [-c CODING] [-n NONCODING] [-o OUTPUT] [-o1 OUTPUT1] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -c CODING, --coding CODING
                        input coding coordinates file in CSV format
  -n NONCODING, --noncoding NONCODING
                        input non-coding coordinates file in CSV format
  -o OUTPUT, --output OUTPUT
                        output filename1
  -o1 OUTPUT1, --output1 OUTPUT
                        output filename2                        
```

## lncRNA function prediction
The script 'LIFT_functionPrediction.sh' generates lncRNA-protein co-expression similarity (LPCS) matrix based on relative expression of FPKM values using RNA-Seq data. Then, it performs function prediction of lncRNAs using BMRF. This script requires following steps:

* FPKM values of lncRNA and mRNA transcript sequences from RNA-Seq data analysis
* Protein IDs associated with GO terms.
* Protein-protein interaction (PPI) matrix file in required format (example): 

```
Zm00001d002235	Zm00001d015999	0.818921668362157
Zm00001d002235	Zm00001d041988	0.818921668362157
Zm00001d002235	Zm00001d041988	0.818921668362157
Zm00001d002235	Zm00001d022551	0.945066124109868
Zm00001d002235	Zm00001d022551	0.945066124109868
```
where column 1 is first interacting protein, column 2 is second interacting protein and column 3 is strength of interaction. NOTE: header should not be present

* The format of FPKM file should be in following format:

```
genename	sample1_FPKM	sample2_FPKM	sample3_FPKM
Zm00001d023154	0.761938	0.192704	0.837553
Zm00001d023190	6.14269	2.57758	3.90534
Zm00001d023074	2.85049	1.73323	2.50522
```
NOTE: header should be present in the file

* Once the lncRNA and mRNA FPKM files are obtained, run the following script to compute relative expression:

```
usage:

Rscript normalizeLPExpression.R <lncRNA FPKM file> <mRNA FPKM file>
``` 

* Once the FPKM values are normalized, compute the LPCS matrix using the following script 'computeLPCSmatrix.py':

### Usage
```
usage: computeLPCSmatrix.py [-h] [-c CODING] [-n NONCODING] [-o OUTPUT] [-v]

optional arguments:
  -h, --help            show this help message and exit
  -c CODING, --coding CODING
                        input coding FPKM relative expression file in CSV
                        format
  -n NONCODING, --noncoding NONCODING
                        input non-coding FPKM relative expression file in CSV
                        format
  -o OUTPUT, --output OUTPUT
                        output filename of LPCS matrix
  -v
```
* Concatenate the LPCS matrix and PPI matrix files as a single file. 
* The GOTerms file should have the following format:

```
Zm00001d002266	GO:0055114	P
Zm00001d048583	GO:0005743	P
Zm00001d048583	GO:0005739	P
Zm00001d048585	GO:0055114	P
Zm00001d032905	GO:0003676	P
Zm00001d032905	GO:0008270	P
Zm00001d032905	GO:0000151	P
```
NOTE: header should not be present

* Run the following script to obtain function prediction:

### Usage
```
usage:

Rscript predictFunction.R <LPCS-PPI matrix file> <Protein IDs associated with GOterms> <Output filename>
```
