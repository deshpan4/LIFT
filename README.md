# LIFT:
## LncRNA Identification and Function prediction Tool
This document provides technical description of the LIFT pipeline for identification, genomic annotation and function prediction of long non-coding RNAs (lncRNAs) based on time-series RNA-Seq data.
## Getting Started
These instructions will help you to run LIFT from linux environment.
The LIFT framework consists of two modules:
* LncRNA Identification Module (LIM)
* Function prediction Module (FPM)
LIM performs following tasks:
* Feature Extraction
* Classification of coding and non-coding RNAs
* Optimization of selected features
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
* SAMtools
* Python2.7 and Python3
* Bowtie2
* [Index and annotation for A.thaliana](ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Arabidopsis_thaliana/Ensembl/TAIR10/Arabidopsis_thaliana_Ensembl_TAIR10.tar.gz)
* R version 3.3 or higher
### R packages
* BMRF (https://github.com/jwbargsten/bmrf)
* iRF (https://github.com/sumbose/iRF)
* AUC
## Usage
