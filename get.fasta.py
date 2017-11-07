import glob,re

dic={}

path1=glob.glob("cds.*.tfa")

for files in path1:
	with open(files) as f:
		book=f.readlines()
		for page in book:
			if ">" in page:
				geneid=page.strip()
				geneid=re.search("\w+", geneid).group(0)
			else:
				seq=page.strip()
				dic[geneid]=seq


# path2=glob.glob("*.p1c2a3")

# for files in path2:
# 	i=1
# 	j=1
# 	with open(files) as f:
# 		spe=files.split(".")[0]
# 		anc=files.split(".")[2]
# 		cdrom=f.readlines()
# 		for scripts in cdrom:
# 			dup1=scripts.split()[0]
# 			dup2=scripts.split()[1]
# 			ancs=scripts.split()[2]
# 			newf1=spe+"."+ anc+ ".parent." + str(i)+".fasta"
# 			f1=open(newf1, "w")
# 			try:
# 				f1.write(">"+dup1+ "\n" + dic[dup1]+"\n"+">"+ancs+"\n"+dic[ancs])
# 				i+=1
# 			except KeyError:
# 				print dup1, newf1
# 			f1.close()
# 			newf2=spe+"."+ anc+ ".child." + str(j)+".fasta"
# 			f2=open(newf2, "w")
# 			try:
# 				f2.write(">"+dup2+ "\n" + dic[dup2]+"\n"+">"+ancs+"\n"+dic[ancs])
# 				j+=1
# 			except KeyError:
# 				print dup1, newf1
# 			f2.close()
			

		
# path2=glob.glob("*.singles")

# for files in path2:
# 	i=1
# 	with open(files) as f:
# 		cdrom=f.readlines()
# 		for scripts in cdrom:
# 			dup1=scripts.split()[0]
# 			dup2=scripts.split()[1]
# 			newf1=files+"."+str(i)+".fasta"
# 			f1=open(newf1, "w")
# 			try:
# 				f1.write(">"+dup1+ "\n" + dic[dup1]+"\n"+">"+ dup2+"\n"+dic[dup2])
# 				i+=1
# 			except KeyError:
# 				print dup1, newf1
			
			

#path2=glob.glob("*.p1c2a3")
# cpassign=glob.glob(".p1c2")
# childids=[]
# parentids=[]
# for files in cpassign:
# 	with open(files) as f:
# 		crest=f.readlines()
# 		for paste in crest:
# 			id1=paste.split()[0]
# 			id2=paste.split()[1]
# 			parentids.append(id1)
# 			childids.append(id2)


geneidfile=open("singles.geneids.bdi_osa_sbi_mac_ath.noclade.txt", "r")

cdrom=geneidfile.readlines()
cdrom=cdrom[1:]
i=1
for scripts in cdrom:
	#spe=scripts.split()[17]
	#dup1=scripts.split()[18]
	sbi=scripts.split()[19]
	osa=scripts.split()[20]
	bdi=scripts.split()[21]
	mac=scripts.split()[22]
	ath=scripts.split()[23]
	newf1="single." + str(i)+".fasta"
	newf2="single." + str(i)+".tree"
	# if dup1 in childids and dup2 in parentids:
	# 	thischild=dup1
	# 	thisparent=dup2
	# else:
	# 	thischild=dup2
	# 	thisparent=dup1
	f2=open(newf2, "w")
	#if spe=="sbi":
		#((((sbi1, sbi2),(bdi, osa)), mac), ath);
		#f2.write("(((("+thisparent+", "+thischild+" ), ("+anc1+ ", "+anc2+ ")),"+mac+"), "+ath+");" )
		#f2.write("\t"+"1"+"\n"+"(((("+thisparent+" #1 , "+thischild+" ), ("+anc1+ ", "+anc2+ ")),"+mac+"), "+ath+");" )
		#f2.write("\t"+"1"+"\n"+"(((("+thisparent+", "+thischild+" #1), ("+anc1+ ", "+anc2+ ")),"+mac+"), "+ath+");" ) ### change this line if the tree differs for parent and child
	#else:
		#(((( (dup1, dup2), ory/bdi),sbi), mac), ath);
		#f2.write("((((("+thisparent+", "+thischild+" ), " +anc1+ "), "+anc2+")," +mac+"), "+ath+");" )
		#f2.write("\t"+"1"+"\n"+"((((("+thisparent+" #1, "+thischild+" ), " +anc1+ "), "+anc2+")," +mac+"), "+ath+");" )  ### change this line if the tree differs for parent and child
		#f2.write("\t"+"1"+"\n"+"((((("+thisparent+", "+thischild+" #1), " +anc1+ "), "+anc2+")," +mac+"), "+ath+");" )
	#i+=1
	#print newf1
	f1=open(newf1, "w")
	try:
		f1.write(">"+sbi+ "\n" + dic[sbi]+"\n"+ ">"+osa+ "\n" + dic[osa]+"\n"+ ">"+bdi+"\n"+dic[bdi]+"\n" +">"+mac+"\n"+dic[mac]+"\n" +">"+ath+"\n"+dic[ath] +"\n" )
		f2.write("(((("+bdi+", "+osa+" ), "+sbi+ "), "+mac+"), "+ath+");" )
		i+=1
	except KeyError:
		print dup1
	f1.close()





"""
	newf1=spe+"."+ anc+ ".branch." + str(i)+".fasta"
	newf2=spe+"."+ anc+ ".branch." + str(i)+".tree"



for files in path2:
	i=1
	j=1
	with open(files) as f:
		spe=files.split(".")[0]
		anc=files.split(".")[2]
		cdrom=f.readlines()
		for scripts in cdrom:
			dup1=scripts.split()[0]
			dup2=scripts.split()[1]
			ancs=scripts.split()[2]
			newf1=spe+"."+ anc+ ".branch." + str(i)+".fasta"
			newf2=spe+"."+ anc+ ".branch." + str(i)+".tree"
			#f1=open(newf1, "w")
			#f2=open(newf2, "w")
			try:
				#f1.write(">"+dup1+ "\n" + dic[dup1]+"\n"+ ">"+dup2+ "\n" + dic[dup2]+"\n"+ ">"+ancs+"\n"+dic[ancs])
				#f2.write("\t"+"1"+"\n"+"(("+dup1+"  #1,"+dup2+" ) ,"+ancs+");")
				i+=1
			except KeyError:
				print dup1, newf1
			#f1.close()
			f2.close()
			#print "(("+dup1+" #1,"+dup2+" #2),"+ancs+")"


"""			