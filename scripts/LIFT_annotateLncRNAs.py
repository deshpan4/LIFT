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
import swalign

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--coding',help='input coding coordinates file in CSV format')
    parser.add_argument('-n', '--noncoding',help='input non-coding coordinates file in CSV format')
    parser.add_argument('-l', '--min-length',help='Minimum length of target sequence for alignment (default=10)')
    parser.add_argument('-ov', '--min-overlap',help='Minimum overlap length (default=5)')
    parser.add_argument('-b', '--bidirectional-cutoff',help='Bidirectional cutoff length (default=1000)')
    parser.add_argument('-d', '--distance-threshold',help='Distance threshold value (default=50)')
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

def extractORF(forwardFrameArray,sequenceArray,coor):
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
	print(len(testing1))
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
		c2=[z + coor[i] for z in c]
		d2=[z + coor[i] for z in d]
		for u in range(0,len(c2)):
			c2d21.append(d2[u])
		c2d2.append(c2d21)		
		for t in range(0,len(c)):
				print(t)
				test2=cwordsArr1[i][d1[t]:c1[t]]
				test3.append(test2)
		orfArray1.append(test3)

	return (orfArray1,c2d2)

c1orfArr=extractORF(cwordsArr,cwordsArr1,codingCoordinates)
n1orfArr=extractORF(nwordsArr,nwordsArr1,ncodingCoordinates)

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

resultArr=[]
annotateResult=[]
def getOverlaps(cE,cI,nE,nI,minLength,minOverlap,distThreshold):
	res1=[]
	res2=[]
	match = 2
	mismatch = -1
	scoring = swalign.NucleotideScoringMatrix(match, mismatch)
	sw = swalign.LocalAlignment(scoring)
	print("distance threshold: ",distThreshold)
	print("minOverlap: ",minOverlap)
	for i in nE:
		for j in i:
			if len(j.get('data'))>=int(minOverlap):
				for i2 in cE:
					for j2 in i2:
						if j.get('strand')== '+' and j2.get('strand') == '+' and j.get('chromosome')==j2.get('chromosome'):
							if ((j.get('start') <= j2.get('start') and j.get('end') > j2.get('start'))) or ((j.get('start') >= j2.get('start') and j.get('start') < j2.get('end'))):
								alignment = sw.align(j.get('data'),j2.get('data'))
								print("Sense Overlap Exonic "," sequence: ",str(nE.index(i)+1)," lncRNA start: ",j.get('start')," lncRNA end: ",j.get('end')," mRNA start: ",j2.get('start')," mRNA end: ",j2.get('end'),", Percent identity: ",alignment.identity*100)
								annotateResult.append("Sequence "+str(nE.index(i)+1)+" Sense Overlap Exonic "+" lncRNA start: "+str(j.get('start'))+" lncRNA end: "+str(j.get('end'))+" mRNA start: "+str(j2.get('start'))+" mRNA end: "+str(j2.get('end'))+" Percent identity: "+str(alignment.identity*100))
								resultArr.append({'sequence':nE.index(i)+1,'gene_type':'Sense Overlapping Exonic'})
	for i in nE:
		for j in i:
			if len(j.get('data'))>=int(minOverlap):
				for i2 in cI:
					for j2 in i2:
						if j.get('strand')== '+' and j2.get('strand') == '+' and j.get('chromosome')==j2.get('chromosome'):
							if ((j.get('start') <= j2.get('start') and j.get('end') > j2.get('start'))) or ((j.get('start') >= j2.get('start') and j.get('start') < j2.get('end'))) or ((j.get('start') <= j2.get('start') and j.get('end') >= j2.get('end'))) or ((j.get('start') >= j2.get('start') and j.get('end') <= j2.get('end'))):
								alignment = sw.align(j.get('data'),j2.get('data'))
								print("Sense Overlap Intronic "," sequence: ",str(nE.index(i)+1)," lncRNA start: ",j.get('start')," lncRNA end: ",j.get('end')," mRNA start: ",j2.get('start')," mRNA end: ",j2.get('end'),", Percent identity: ",alignment.identity*100)
								annotateResult.append("Sequence "+str(nE.index(i)+1)+" Sense Overlap Intronic "+" lncRNA start: "+str(j.get('start'))+" lncRNA end: "+str(j.get('end'))+" mRNA start: "+str(j2.get('start'))+" mRNA end: "+str(j2.get('end'))+" Percent identity: "+str(alignment.identity*100))
								resultArr.append({'sequence':nE.index(i)+1,'gene_type':'Sense Overlapping Intronic'})
	for i in nE:
		for j in i:
			if len(j.get('data'))>=int(minOverlap):
				for i2 in cE:
					for j2 in i2:
						if j.get('strand')== '-' and j2.get('strand') == '+' and j.get('chromosome')==j2.get('chromosome'):
							if ((j.get('start') <= j2.get('start') and j.get('end') > j2.get('start'))) or ((j.get('start') >= j2.get('start') and j.get('start') < j2.get('end'))) or ((j.get('start') <= j2.get('start') and j.get('end') >= j2.get('end'))) or ((j.get('start') >= j2.get('start') and j.get('end') <= j2.get('end'))):
								alignment = sw.align(j.get('data'),j2.get('data'))
								print("AntiSense Overlap Exonic "," sequence: ",str(nE.index(i)+1)," lncRNA start: ",j.get('start')," lncRNA end: ",j.get('end')," mRNA start: ",j2.get('start')," mRNA end: ",j2.get('end'),", Percent identity: ",alignment.identity*100)
								annotateResult.append("Sequence "+str(nE.index(i)+1)+" Antisense Overlap Exonic "+" lncRNA start: "+str(j.get('start'))+" lncRNA end: "+str(j.get('end'))+" mRNA start: "+str(j2.get('start'))+" mRNA end: "+str(j2.get('end'))+" Percent identity: "+str(alignment.identity*100))
								resultArr.append({'sequence':nE.index(i)+1,'gene_type':'AntiSense Overlap Exonic'})
	for i in nE:
		for j in i:
			if len(j.get('data'))>=int(minOverlap):
				for i2 in cI:
					for j2 in i2:
						if j.get('strand')== '-' and j2.get('strand') == '+' and j.get('chromosome')==j2.get('chromosome'):
							if ((j.get('start') <= j2.get('start') and j.get('end') > j2.get('start'))) or ((j.get('start') >= j2.get('start') and j.get('start') < j2.get('end'))) or ((j.get('start') <= j2.get('start') and j.get('end') >= j2.get('end'))) or ((j.get('start') >= j2.get('start') and j.get('end') <= j2.get('end'))):
								alignment = sw.align(j.get('data'),j2.get('data'))
								print("AntiSense Overlap Intronic "," sequence: ",str(nE.index(i)+1)," lncRNA start: ",j.get('start')," lncRNA end: ",j.get('end')," mRNA start: ",j2.get('start')," mRNA end: ",j2.get('end'),", Percent identity: ",alignment.identity*100)
								annotateResult.append("Sequence "+str(nE.index(i)+1)+" Antisense Overlap Intronic "+" lncRNA start: "+str(j.get('start'))+" lncRNA end: "+str(j.get('end'))+" mRNA start: "+str(j2.get('start'))+" mRNA end: "+str(j2.get('end'))+" Percent identity: "+str(alignment.identity*100))
								resultArr.append({'sequence':nE.index(i)+1,'gene_type':'AntiSense Overlap Intronic'})
	
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

result=getOverlaps(codingE,codingI,noncodingE,noncodingI,args.min_length,args.min_overlap,args.distance_threshold)

def getBidirectionalRNA(codingArray,ncodingArray,codingCoor,ncodingCoor,cStrand,ncStrand,codingChr,ncodingChr):
	cseqList=[]
	for i in range(0,len(codingArray)):
		cseqList.append({'sequence':(i+1),'orfstart':codingCoor[i],'strand':cStrand[i],'chromosome':codingChr[i]})
	
	cseqList1=[]
	for i in range(0,len(cseqList)):
		if cseqList[i].get('orfstart') == []:
			continue
		cseqList1.append({'sequence':cseqList[i].get('sequence'),'orfstart':cseqList[i].get('orfstart')[0],'strand':cseqList[i].get('strand'),'chromosome':cseqList[i].get('chromosome')})
	
	
	nseqList=[]
	for i in range(0,len(ncodingArray)):
		nseqList.append({'sequence':(i+1),'orfstart':ncodingCoor[i],'strand':ncStrand[i],'chromosome':ncodingChr[i]})
	
	nseqList1=[]
	for i in range(0,len(nseqList)):
		if nseqList[i].get('orfstart') == []:
			continue
		nseqList1.append({'sequence':nseqList[i].get('sequence'),'orfstart':nseqList[i].get('orfstart')[0],'strand':nseqList[i].get('strand'),'chromosome':nseqList[i].get('chromosome')})
	
	
	for i in range(0,len(nseqList1)):
		for j in range(0,len(cseqList1)):
			if (nseqList1[i].get('strand') == '-' and cseqList1[j].get('strand') == '+') or (nseqList1[i].get('strand') == '+' and cseqList1[j].get('strand')) or (nseqList1[i].get('strand') == cseqList1[j].get('strand')):
				if nseqList1[i].get('orfstart') < cseqList1[j].get('orfstart') and nseqList1[i].get('orfstart') > (cseqList1[j].get('orfstart')-1000):
					print("Bidirectional promoter ","lncRNA sequence: ",nseqList1[i].get('sequence'),", mRNA sequence: ",cseqList1[j].get('sequence'))
					annotateResult.append("Sequence: "+str(nseqList1[i].get('sequence'))+" Bidirectional promoter")
					resultArr.append({'sequence':nseqList1[i].get('sequence'),'gene_type':'Bidirectional promoter'})
	return resultArr 				

bidirectionalRes=getBidirectionalRNA(corfArr,norfArr,cCoor1,nCoor1,codingStrand,ncodingStrand,codingChromosome,ncodingChromosome)


intergenicArr=[]
seqFullArr=range(1,len(nwordsArr1)+1)
resultArr1=[]

for i in range(0,len(resultArr)):
	resultArr1.append(resultArr[i].get('sequence'))

intArr1=filter(lambda x: x not in resultArr1, seqFullArr)

for i in range(0,len(intArr1)):
	annotateResult.append("sequence "+str(intArr1[i])+" Intergenic")

np.savetxt(args.output, annotateResult, delimiter="\t", fmt="%s")
