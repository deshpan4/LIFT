import argparse
import os, sys, getopt
import pandas as pd
import numpy as np
import math
import csv

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--coding',help='input coding FPKM relative expression file in CSV format')
    parser.add_argument('-n', '--noncoding',help='input non-coding FPKM relative expression file in CSV format')
    parser.add_argument('-tP', '--thresholdPositive',help='input threshold value for obtaining positive correlations. (Default value = 0.9)')
    parser.add_argument('-tN', '--thresholdNegative',help='input threshold value for obtaining negative correlations. (Default value = -0.9)')
    parser.add_argument('-o', '--output',help='output filename of LPCS matrix')
    args = parser.parse_args()
	
d1=pd.read_csv(args.noncoding)
d2=pd.read_csv(args.coding)
d1len=d1.shape
d2len=d2.shape

if args.thresholdPositive == None:
	args.thresholdPositive = 0.9

if args.thresholdNegative == None:
	args.thresholdNegative = -0.9	

thresPos = int(args.thresholdPositive)
thresNeg = int(args.thresholdNegative)

dash="\t"
e1="\t"
res1=[]
for i in range(0, d1len[1]):
 if np.sum(np.isnan(d1.ix[:,i].values)) >= 7:
  continue
 elif np.sum(d1.ix[:,i].values==0) >= 7:
  continue
 for x in range(0, d2len[1]):
  if np.sum(np.isnan(d2.ix[:,x].values)) >= 7:
   continue
  elif np.sum(d2.ix[:,x].values==0) >= 7:
   continue
  j=np.corrcoef(d1.ix[:,i],d2.ix[:,x])[0,1]
  if j >= thresPos and j != 1:
   #print(str(list(d1.iloc[:,i:i+1]))+dash+str(list(d2.iloc[:,x:x+1]))+e1,j)
   res1.append(str(list(d1.iloc[:,i:i+1]))+dash+str(list(d2.iloc[:,x:x+1]))+e1+str(j))

dash="\t"
e1="\t"
res2=[]
for i in range(0, d1len[1]):
 if np.sum(np.isnan(d1.ix[:,i].values)) >= 7:
  continue
 elif np.sum(d1.ix[:,i].values==0) >= 7:
  continue
 for x in range(0, d2len[1]):
  if np.sum(np.isnan(d2.ix[:,x].values)) >= 7:
   continue
  elif np.sum(d2.ix[:,x].values==0) >= 7:
   continue
  j=np.corrcoef(d1.ix[:,i],d2.ix[:,x])[0,1]
  if j <= thresNeg and j != 1:
   #print(str(list(d1.iloc[:,i:i+1]))+dash+str(list(d2.iloc[:,x:x+1]))+e1,j)
   res2.append(str(list(d1.iloc[:,i:i+1]))+dash+str(list(d2.iloc[:,x:x+1]))+e1+str(j))

negres2=[]
for i in range(0,len(res2)):
	tmp1=res2[i].replace("[", "")
	tmp2=tmp1.replace("]", "")
	tmp3=tmp2.replace("'", "")
	negres2.append(tmp3)

posres2=[]
for i in range(0,len(res1)):
	tmp1=res1[i].replace("[", "")
	tmp2=tmp1.replace("]", "")
	tmp3=tmp2.replace("'", "")
	posres2.append(tmp3)	
	
twolists = posres2 + negres2
np.savetxt(args.output, twolists, delimiter="", fmt="%s")