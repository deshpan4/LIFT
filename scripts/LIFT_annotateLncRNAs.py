from __future__ import division
import argparse
import os, sys, getopt
import itertools
from itertools import islice
from collections import defaultdict
from collections import Counter
import numpy as np
import csv
import re
import math
import pandas as pd
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--coding',help='input coding coordinates file in CSV format')
    parser.add_argument('-n', '--noncoding',help='input non-coding coordinates file in CSV format')
    parser.add_argument('-l', '--min_length',help='Minimum length of target sequence for alignment (default=10)')
    parser.add_argument('-ov', '--min_overlap',help='Minimum overlap length (default=5)')
    parser.add_argument('-b', '--bidirectional_cutoff',help='Bidirectional cutoff length (default=1000)')
    parser.add_argument('-d', '--distance_threshold',help='Distance threshold value (default=50)')
    parser.add_argument('-o', '--output',help='output filename')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

if args.min_length == None:
	args.min_length = 10

if args.min_overlap == None:
	args.min_overlap = 5

if args.distance_threshold == None:
	args.distance_threshold = 50

if args.bidirectional_cutoff == None:
	args.bidirectional_cutoff = 1000	
	
cdf=pd.read_csv(args.coding)
ndf=pd.read_csv(args.noncoding)

wordsArr=[]
hexamerArray = []
codingArray = []
noncodingArray = []
header = []
alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

cwordsArr1 = cdf['seqs'].values
nwordsArr1 = ndf['seqs'].values

codingCoordinates=cdf['start'].values
ncodingCoordinates=ndf['start'].values

codingCoordinatesEnd=cdf['end'].values
ncodingCoordinatesEnd=ndf['end'].values

codingChromosome=cdf['chromosome'].values
ncodingChromosome=ndf['chromosome'].values

codingStrand=cdf['strand'].values
ncodingStrand=ndf['strand'].values

cwordsArr=[]
nwordsArr=[]

for i in range(0,len(cwordsArr1)):
	cwordsArr.append(list(cwordsArr1[i]))

for i in range(0,len(nwordsArr1)):
	nwordsArr.append(list(nwordsArr1[i]))

codingArray=[]
noncodingArray=[]

for i in range(0,len(cwordsArr1)):
	codingArray.append({'chromosome':codingChromosome[i],'start': codingCoordinates[i], 'end': codingCoordinatesEnd[i], 'strand':codingStrand[i], 'data': cwordsArr1[i]})

for i in range(0,len(nwordsArr1)):
	noncodingArray.append({'chromosome':ncodingChromosome[i],'start': ncodingCoordinates[i], 'end': ncodingCoordinatesEnd[i], 'strand':ncodingStrand[i], 'data': nwordsArr1[i]})

def complement1(seq):    
    for k,v in alt_map.iteritems():
        seq = seq.replace(k,v)
    bases = list(seq) 
    bases = ([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.iteritems():
        bases = bases.replace(v,k)
    return bases



def aggregateCodonsAminoAcids(sequenceArray):
	c1=[]
	c2=[]
	for i in range(0,len(sequenceArray)):
		c1=[]
		for j in range(0,len(sequenceArray[i])):
			if sequenceArray[i][j] == 'TTC' or sequenceArray[i][j] == 'TTT':
				codon='F'
				c1.append(codon)
			elif sequenceArray[i][j] == 'TTA' or sequenceArray[i][j] == 'TTG' or sequenceArray[i][j] == 'CTT' or sequenceArray[i][j] == 'CTC' or sequenceArray[i][j] == 'CTA' or sequenceArray[i][j] == 'CTG':
				codon='L'
				c1.append(codon)
			elif sequenceArray[i][j] == 'ATT' or sequenceArray[i][j] == 'ATC' or sequenceArray[i][j] == 'ATA':
				codon='I'
				c1.append(codon)
			elif sequenceArray[i][j] == 'ATG':
				codon='M'
				c1.append(codon)
			elif sequenceArray[i][j] == 'GTT' or sequenceArray[i][j] == 'GTC' or sequenceArray[i][j] == 'GTA' or sequenceArray[i][j] == 'GTG':
				codon='V'
				c1.append(codon)
			elif sequenceArray[i][j] == 'TCT' or sequenceArray[i][j] == 'TCC' or sequenceArray[i][j] == 'TCA' or sequenceArray[i][j] == 'TCG' or sequenceArray[i][j] == 'AGT' or sequenceArray[i][j] == 'AGC':
				codon='S'
				c1.append(codon)
			elif sequenceArray[i][j] == 'CCT' or sequenceArray[i][j] == 'CCC' or sequenceArray[i][j] == 'CCA' or sequenceArray[i][j] == 'CCG':
				codon='P'
				c1.append(codon)
			elif sequenceArray[i][j] == 'ACT' or sequenceArray[i][j] == 'ACC' or sequenceArray[i][j] == 'ACA' or sequenceArray[i][j] == 'ACG':
				codon='T'
				c1.append(codon)
			elif sequenceArray[i][j] == 'GCT' or sequenceArray[i][j] == 'GCC' or sequenceArray[i][j] == 'GCA' or sequenceArray[i][j] == 'GCG':
				codon='A'
				c1.append(codon)
			elif sequenceArray[i][j] == 'TAT' or sequenceArray[i][j] == 'TAC':
				codon='Y'
				c1.append(codon)
			elif sequenceArray[i][j] == 'TAA' or sequenceArray[i][j] == 'TAG' or sequenceArray[i][j] == 'TGA':
				codon='-'
				c1.append(codon)
			elif sequenceArray[i][j] == 'CAT' or sequenceArray[i][j] == 'CAC':
				codon='H'
				c1.append(codon)
			elif sequenceArray[i][j] == 'CAA' or sequenceArray[i][j] == 'CAG':
				codon='Q'
				c1.append(codon)
			elif sequenceArray[i][j] == 'AAT' or sequenceArray[i][j] == 'AAC':
				codon='N'
				c1.append(codon)
			elif sequenceArray[i][j] == 'AAA' or sequenceArray[i][j] == 'AAG':
				codon='K'
				c1.append(codon)
			elif sequenceArray[i][j] == 'GAT' or sequenceArray[i][j] == 'GAC':
				codon='D'
				c1.append(codon)
			elif sequenceArray[i][j] == 'GAA' or sequenceArray[i][j] == 'GAG':
				codon='E'
				c1.append(codon)
			elif sequenceArray[i][j] == 'TGT' or sequenceArray[i][j] == 'TGC':
				codon='C'
				c1.append(codon)
			elif sequenceArray[i][j] == 'TGG':
				codon='W'
				c1.append(codon)
			elif sequenceArray[i][j] == 'CGT' or sequenceArray[i][j] == 'CGC' or sequenceArray[i][j] == 'CGA' or sequenceArray[i][j] == 'CGG' or sequenceArray[i][j] == 'AGA' or sequenceArray[i][j] == 'AGG':
				codon='R'
				c1.append(codon)
			elif sequenceArray[i][j] == 'GGT' or sequenceArray[i][j] == 'GGC' or sequenceArray[i][j] == 'GGA' or sequenceArray[i][j] == 'GGG':
				codon='G'
				c1.append(codon)
		c2.append(c1)
	
	return(c2)

def extractORF(forwardFrameArray,sequenceArray):
	a2=[]
	for i in range(0,len(forwardFrameArray)):
		a1=[]
		for x in range(0,len(forwardFrameArray[i]),3):
			list1 = forwardFrameArray[i][x:x+3]
			list2 = ''.join(map(str,list1))
			a1.append(list2)
		a2.append(a1)
	array11=[]
	array11.append(aggregateCodonsAminoAcids(a2))
	orfArray1=[]
	orfArray2=[]
	orfArray1=[]
	testing1=array11[0]
	coor=[20001,21113]
	c2d2=[]
	for i in range(0,len(testing1)):
		test3=[]
		a=[]
		b=[]
		c=[]
		d=[]
		c2d21=[]
		x=0
		a = [m for m, n in enumerate(testing1[i]) if n == "-"]
		b = [m for m, n in enumerate(testing1[i]) if n == "M"]
		if not a:
			a.append(len(testing1[i]))
		for j in range(0,len(b)):
			while b[j] > a[x] and b[j] < a[len(a)-1]:
				x=x+1
			if len(c) == 0:
				c.append(a[x])
				d.append(b[j])
			elif b[j] > a[x]:
				c.append(len(testing1[i]))
				d.append(b[j])
				break
			elif a[x] == c[len(c)-1]:
				continue
			else:
				c.append(a[x])	
				d.append(b[j])
		c1=[z * 3 for z in c]
		d1=[z * 3 for z in d]
		c2=[z + coor[0] for z in c]
		d2=[z + coor[0] for z in d]
		for u in range(0,len(c2)):
			c2d21.append(d2[u])
		c2d2.append(c2d21)		
		for t in range(0,len(c)):
				test2=cwordsArr1[i][d1[t]:c1[t]]
				test3.append(test2)
		orfArray1.append(test3)

	return (orfArray1,c2d2)

c1orfArr=extractORF(cwordsArr,cwordsArr1)
n1orfArr=extractORF(nwordsArr,nwordsArr1)

corfArr=c1orfArr[0]
cCoor1=c1orfArr[1]
norfArr=n1orfArr[0]
nCoor1=n1orfArr[1]


def getAGGTs(data1,coord,strand,chromosome):
	exonArray=[]
	intronArray=[]
	for i in range(0,len(data1)):
		E=[]
		I=[]
		for j in range(0,len(data1[i])):
		    data=data1[i][j]
		    bias=coord[i][j]
		    strand1=strand[i]
		    chromosome1=chromosome[i]
		    if data[-2:]!= "AG" :
		    	data = data + "AG" 
		    while len(data)!=0:
		        temp_data=data[:data.find("GT")]
		        start=bias
		        end=bias+data.find("GT")-1
		        E.append({'chromosome':chromosome1,'start': start, 'end': end, 'strand':strand1, 'data': temp_data})
		        temp_data2 = data[data.find("GT"):data.find("AG", data.find("GT"))+2]
		        start2 = bias+data.find("GT")
		        end2 = bias+data.find("AG", data.find("GT"))+1
		        I.append({'chromosome':chromosome1,'start': start2, 'end': end2, 'strand':strand1, 'data': temp_data2})
		        bias+=data.find("AG", data.find("GT"))+2
		        data = data[data.find("AG", data.find("GT"))+2:]
		        if data.find("GT")== -1 :
		            E.append({'chromosome':chromosome1,'start': bias, 'end': bias+len(data)-1, 'strand':strand1,'data': data})
		            data=""
		exonArray.append(E)
		intronArray.append(I)
	return (exonArray,intronArray)

(codingE,codingI)=getAGGTs(corfArr,cCoor1,codingStrand,codingChromosome);
(noncodingE,noncodingI)=getAGGTs(norfArr,nCoor1,ncodingStrand,ncodingChromosome);


def getOverlaps_old(cE,cI,nE,nI,start,end,strand,chromosome):
	res1=[]
	res2=[]
	for i in cE:
		for j in i:
			if j.get('start') >= 20004 and j.get('end') <= 20111:
				res1.append({'Type':'Exon','sequence':cE.index(i),'start': j.get('start'), 'end': j.get('end'), 'strand':j.get('strand'), 'chromosome':j.get('chromosome')})

	for i in cI:
		for j in i:
			if j.get('start') >= 20004 and j.get('end') <= 20111:
				res2.append({'Type':'Intron','sequence':cI.index(i),'start': j.get('start'), 'end': j.get('end'), 'strand':j.get('strand'), 'chromosome':j.get('chromosome')})

	return (res1,res2)

annotateResult=[]
def getOverlaps(cE,cI,nE,nI,minLength,minOverlap,distThreshold,biLength):
	res1=[]
	res2=[]
	for i in nE:
		for j in i:
			if len(j.get('data'))>=minOverlap:
				for i2 in cE:
					for j2 in i2:
						if j.get('strand')== '+' and j2.get('strand') == '+' and j.get('chromosome')==j2.get('chromosome'):
							if ((j2.get('start') <= j.get('start') and j2.get('start') >= (j.get('start')-distThreshold))) or ((j2.get('end') >= j.get('end') and j2.get('end') <= (j.get('end')+distThreshold))) or ((j2.get('start') >= j.get('start')) and ((j2.get('end') <= j.get('end')))):
								data=getOverlap(j.get('data'),j2.get('data'),minLength)
								if data!=None:
									annotateResult.append("Sequence: "+str(nE.index(i)+1)+" Sense Overlap Exonic "+str(data))
	for i in nE:
		for j in i:
			if len(j.get('data'))>=minOverlap:
				for i2 in cI:
					for j2 in i2:
						if j.get('strand')== '+' and j2.get('strand') == '+' and j.get('chromosome')==j2.get('chromosome'):
							if (j2.get('start') <= j.get('start') and j2.get('start') >= (j.get('start')-distThreshold)) or ((j2.get('end') >= j.get('end') and j2.get('end') <= (j.get('end')+distThreshold))) or ((j2.get('start') >= j.get('start')) and ((j2.get('end') <= j.get('end')))):
								data=getOverlap(j.get('data'),j2.get('data'),minLength)
								if data!=None:
									annotateResult.append("Sequence: "+str(nE.index(i))+" Sense Overlap Intronic Actual "+str(data))
	for i in nE:
		for j in i:
			if len(j.get('data'))>=minOverlap:
				for i2 in cE:
					for j2 in i2:
						if j.get('strand')== '-' and j2.get('strand') == '+' and j.get('chromosome')==j2.get('chromosome'):
							if ((j2.get('start') <= j.get('start') and j2.get('start') >= (j.get('start')-distThreshold))) or ((j2.get('end') >= j.get('end') and j2.get('end') <= (j.get('end')+distThreshold))) or ((j2.get('start') >= j.get('start')) and ((j2.get('end') <= j.get('end')))):
								data=getOverlap(j.get('data'),complement1(j2.get('data')),minLength)
								if data!=None:
									annotateResult.append("Sequence: "+str(nE.index(i))+" AntiSense Overlap Exonic "+str(data))
	for i in nE:
		for j in i:
			if len(j.get('data'))>=minOverlap:
				for i2 in cI:
					for j2 in i2:
						if j.get('strand')== '-' and j2.get('strand') == '+' and j.get('chromosome')==j2.get('chromosome'):
							if (j2.get('start') <= j.get('start') and j2.get('start') >= (j.get('start')-distThreshold)) or ((j2.get('end') >= j.get('end') and j2.get('end') <= (j.get('end')+distThreshold))) or ((j2.get('start') >= j.get('start')) and ((j2.get('end') <= j.get('end')))):
								data=getOverlap(j.get('data'),complement1(j2.get('data')),minLength)
								if data!=None:
									annotateResult.append("Sequence: "+str(nE.index(i))+" AntiSense Overlap Intronic Actual "+str(data))

	ntmp1=[]
	ntmp2=[]
	for i in nI:
		ntmp={}
		ntmp1=[]
		for j in i:
			ntmp=j
		ntmp1.append(ntmp)
		ntmp2.append(ntmp1)

	ctmp1=[]
	ctmp2=[]
	for i in cE:
		ctmp={}
		ctmp1=[]
		for j in i:
			ctmp=j
			ctmp1.append(ctmp)
			break
		ctmp2.append(ctmp1)		

	for i in ntmp1:
		for j in ctmp1:
			if i.get('strand') == '-' and j.get('strand') == '+' and i.get('chromosome')==j.get('chromosome'):
				if i.get('end') < j.get('start') and i.get('end') > (j.get('start')-biLength):
					annotateResult.append("Sequence: "+str(nI.index(z))+" Bidirectional lncRNA")
					print("Sequences ",ntmp1.index(i),ctmp1.index(j))
					print("Bidirectional lncRNA")
				

	return (res1,res2)

def getOverlap(left,right,min_overlap):
	if len(left)<min_overlap or len(right)<min_overlap:
		return None
	if left == right:
		return "equal:"+left+" - "+right+" perc: 100"
	else:
		if right not in left: 
			tmp=right
			iter=0;
			while len(tmp)>min_overlap:
				tmp=tmp[1:]
				iter+=1
				if left.find(tmp)==0:
					return "left:"+left+" - "+right+" perc:"+str(len(right)-iter/len(right)*100.0)
			tmp=right
			iter=0
			while len(tmp)>min_overlap:
				tmp=tmp[:-1]
				iter+=1
				if left.find(tmp)!=-1 and left.find(tmp)+len(tmp)==len(left):
					return "right:"+left+" - "+right+" perc:"+str((len(right)-iter)/len(right)*100.0)
			return None
		else:
			return "middle:"+left+" - "+right+" perc: 100"

result=getOverlaps(codingE,codingI,noncodingE,noncodingI,args.min_length,args.min_overlap,args.distance_threshold,args.bidirectional_cutoff)

def getIntergenicOverlaps(codingArray,noncodingArray):
	res1=[]
	res2=[]
	for i in noncodingArray:
		if len(i.get('data'))==len(i.get('data')):
			for i2 in codingArray:
				if i.get('strand')==i2.get('strand') and i.get('chromosome')==i2.get('chromosome'):
					if (i2.get('start') >= i.get('start')) and (i2.get('end')) <= i.get('end'):
						annotateResult.append("Sequence: "+str(noncodingArray.index(i))+"Non-Intergenic")
					elif (i2.get('start') >= i.get('start')) and (i2.get('start')) <= i.get('end'):
						annotateResult.append("Sequence: "+str(noncodingArray.index(i))+"Non-Intergenic")
					elif (i2.get('start') <= i.get('start')) and (i2.get('end')) >= i.get('start'):
						annotateResult.append("Sequence: "+str(noncodingArray.index(i))+"Non-Intergenic")
					elif (i2.get('start') <= i.get('end')) and (i2.get('end')) >= i.get('end'):
						annotateResult.append("Sequence: "+str(noncodingArray.index(i))+"Non-Intergenic")
	return res1	

np.savetxt(args.output, annotateResult, delimiter="\t", fmt="%s")