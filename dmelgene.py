#! /usr/bin/python
import re

f1=open("dmel-all-gene-r6.12.fasta.cordinates", 'w')

with open("dmel-all-gene-r6.12.fasta.header") as f:
	light=f.readlines()
	for bulb in light:
		geneid=re.search('(FBgn)(\S+)' ,bulb).group(0)
		chromarm=re.search('(?<=loc=)([^:]+)' ,bulb).group(0)
		loc=re.search('(?<=loc=)([^;]+)', bulb).group(0)
		if '(' in loc:
			coor=re.search('(?<=\()([^)]+)',loc).group(0)
			coordina=re.findall('\d+', coor)
		else:
			coor_=re.search('(?<=:)(\S+)', loc).group(0)
			coordina=re.findall('\d+', coor_)
		lower=coordina[0]
		upper=coordina[1]
		f1.write(geneid+'\t'+chromarm+'\t'+lower+'\t'+upper+'\n')
f1.close()