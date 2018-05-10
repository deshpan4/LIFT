# LIFT:
## LncRNA Identification and Function prediction Tool
This document provides technical description of the LIFT pipeline for identification, genomic annotation and function prediction of long non-coding RNAs (lncRNAs) based on RNA-Seq data.
## Download
Use git clone:
```
git clone https://github.com/deshpan4/LIFT
```
## Installation
Download the folder using 'git clone https://github.com/deshpan4/LIFT.git' command. Make sure all the dependencies are installed before running the framework.
## Prerequisites
Following tools should be installed before executing this pipeline on the host system:
* NodeJS
* Python2.7 and Python3
* R version 3.3 or higher
* Bowtie2
* [Index and annotation for A.thaliana](ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Arabidopsis_thaliana/Ensembl/TAIR10/Arabidopsis_thaliana_Ensembl_TAIR10.tar.gz)
* SAMtools
### R packages
* BMRF (https://github.com/jwbargsten/bmrf)
* iRF (https://github.com/sumbose/iRF)
* AUC
### NodeJS installation
For installing NodeJS on Ubuntu-based linux distributions, follow the instructions as mentioned [here](https://nodejs.org/en/download/package-manager/#debian-and-ubuntu-based-linux-distributions).
For Enterprise linux and Fedora systems, follow the instructions as mentioned [here](https://nodejs.org/en/download/package-manager/#enterprise-linux-and-fedora)
For other systems, please follow the instructions mentioned [here](https://nodejs.org/en/download/package-manager/)
## Tutorial (User Manual)
The installation consists of complete tutorial "LIFT Tutorial.pdf" demonstrating usage of the scripts. The installation also contains example datasets in the **"data"** folder.
For usage instructions of the scripts and data, please refer **"LIFT Tutorial.pdf"** tutorial.
## Usage scripts
LIFT consists of nine main scripts:
* LIFT_extractFeatures_testSet.sh (Feature extraction of test set sequences)
* LIFT_extractFeatures_trainingSet.sh (Feature extraction of training set sequences)
* LIFT_lncRNAPredict.py (Prediction of lncRNA sequences from test set sequences)
* LIFT_LiRFFS.py (Identification of optimal features)
* LIFT_annotateLncRNAs.py (Sub-classification of lncRNA sequences)
* normalizeLPExpression.R (Normalization of expression values)
* computeLPCSmatrix.py (Computation of LPCS matrix)
* predictFunction.R (Prediction of lncRNA functions)
* appendFunction.R (Annotation of Gene Ontology functions and function type)
### Feature Extraction of test set sequences
The script 'LIFT_extractFeatures_testSet.sh' extract features from test set FASTA file for identification of lncRNA transcripts. Following input files are needed for extraction of features:
* Coding FASTA sequences
* Noncoding FASTA sequences
* Test set FASTA file
* Output filename
* Scripts directory

#### Usage
```
bash LIFT_extractFeatures_testSet.sh -testc <coding filename> -testl <noncoding filename> -fasta <target filename for feature extraction in FASTA format> -output <output filename> -scripts /LIFT/lib -b /LIFT/lib -cpat /LIFT/lib
```
NOTE: For detailed usage, please see the tutorial. 

### Feature extraction of training set sequences
The script 'LIFT_extractFeatures_trainingSet.sh' extract features from training set FASTA file for identification of lncRNA transcripts. Following input files are needed for extraction of features:
* Training coding FASTA sequences
* Training noncoding FASTA sequences
* Output filename
* Scripts directory

#### Usage
```
bash LIFT_extractFeatures_trainingSet.sh -coding <coding filename> -noncoding <noncoding filename> -output <output filename> -scripts /LIFT/lib -b /LIFT/lib -cpat /LIFT/lib
```
**NOTE: For detailed usage, please see the tutorial.** 

### Prediction of lncRNA sequences from test set sequences
The script 'LIFT_lncRNAPredict.py' predicts lncRNA sequences from the test set FASTA sequences using training set sequences. The script requires following input files:
* Training set matrix
* Test set matrix
* Test set FASTA sequences

#### Usage
```
python3 LIFT_lncRNAPredict.py -tr <Training set matrix file> -te <Test set matrix file> -tef <Test set FASTA sequences> -o <Output filename>
```
NOTE: For detailed usage, please see the tutorial. 

### Feature Optimization
The script 'LIFT_LiRFFS.py' performs feature selection for selection of optimal features. Following input files and parameters are required:
* Training set feature matrix (tr)
* Validation set feature matrix (te)
* Lower cutoff for lambda (lambdaL)
* Upper cutoff for lambda (lambdaU)
* Step size between Lambda values (lambdaS)
* Tolerance level (tol)
* output training set with minimal optimal features (otrMin)
* Output validation set with minimal optimal features (oteMin)
* output training set with maximal optimal features (otrMax)
* Output validation set with maximal optimal features (oteMax)

#### Usage
```
python3 LIFT_LiRFFS.py -tr <training feature set in CSV format> -te <test feature set in CSV format> -lambdaL <input lower lambda value> -lambdaU <input upper lambda value> -lambdaS <input lambda step size> -tol <input tolerance value> -otrMin <output training set filename with minimal optimal features> -oteMin <output validation set with minimal optimal features> -otrMax <output training set filename with maximal optimal features> -oteMax <output validation set with maximal optimal features>
```
**NOTE: For detailed usage, please see the tutorial.**

### lncRNA genomic annotation or sub-classification of lncRNA sequences
The script 'LIFT_annotateLncRNAs.py' provides annotation and classification of lncRNAs based on their coordinates. It requires following files and parameters for annotation:
* Coding sequence and coordinates file having following columns in CSV format: sequence,chromosome,start,end,strand
* Non-coding sequence and coordinates file having following columns in CSV format: sequence,chromosome,start,end,strand
* Output filename

#### Usage
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
NOTE: For detailed usage, please see the tutorial. 

### lncRNA function prediction
For prediction of lncRNA functions, FPKM expression values are required. The format of FPKM file should be in following format:

```
genename	sample1_FPKM	sample2_FPKM	sample3_FPKM
Zm00001d023154	0.761938	0.192704	0.837553
Zm00001d023190	6.14269	2.57758	3.90534
Zm00001d023074	2.85049	1.73323	2.50522
```
**NOTE: For detailed usage, please see the tutorial.**

* Once the lncRNA and mRNA FPKM files are obtained, run the following script to compute relative expression:

```
usage:

Rscript normalizeLPExpression.R <lncRNA FPKM file> <mRNA FPKM file>
``` 

* Once the FPKM values are normalized, compute the LPCS matrix using the following script 'computeLPCSmatrix.py':

#### Usage
```
usage: /usr/bin/python2.7 computeLPCSmatrix.py -c <CODING FPKM expression file> -n <NONCODING FPKM expression file> -o <OUTPUT file>

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

#### Usage
For prediction of lncRNA function, following command is required:
```
usage:

Rscript predictFunction.R <LPCS-PPI matrix file> <Protein IDs associated with GOterms> <Output filename>
```
### Annotation and filtering of lncRNA sequences based on probability value
The 'appendFunction.R' script is used for annotation and filtering the lncRNA sequences based on probability and Gene Ontology data. Following input files and parameters are required:
* Output prediction file
* Protein-coding GO annotation file
* lncRNA identifiers
* Cutoff value
* Output file
#### Usage
```
usage:

Rscript appendFunction.R <Predicted lncRNA function association file> <protein-coding Gene Ontology annotated file> <lncRNA gene identifiers/genenames> <Probability cutoff value> <Output filename>
```

**NOTE: For detailed usage, please see the tutorial.**
