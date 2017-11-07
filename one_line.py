#!/usr/bin/python
import re

#Get the CDS sequence file for D. seculans 
#Originally,the sequence in separated by '\n'. Get rid of '\n' and put them in one line
#Add '>' in the end of the original data file. Otherwise the program will not return the last sequence.
sec_dic={}
sec_header_list=[]
sec_seq_list=[]
sequence=''
file_1=open("dmel-all-gene-r6.12.fasta",'r')
for dta in file_1:
	dta=dta.strip()
	for i in dta:
		if i=='>':
			sec_header_list.append(dta)
			if sequence:
				sec_seq_list.append(sequence)
				sequence=''
			break
		else:
			continue
	if all([k==k.upper() for k in dta]):
		sequence=sequence+dta
file_1.close()

#Generate a new file "dsec-all-CDS.fasta" from the same data.
file_2=open("dmel-all-gene-r6.12.fasta.line", 'w')
for n in range(0,len(sec_seq_list)):
	header=sec_header_list[n]
	seq=sec_seq_list[n]
	file_2.write(header+'\n'+seq+'\n')



