import glob, re

path=glob.glob("mel_pse.ca.dup.align/*.paml")

f1=open("mel_pse.ca.dup.paml.results", "w")

for files in path:
	with open(files) as f:
		bloom=f.readlines()
		for anchor in bloom:
			if "(FBgn" in anchor:
				genids=re.findall("FBgn\w+", anchor)
				print genids
			if "t=" in anchor:
				print anchor
				allinfor=anchor.split("=")
				dn=re.search("\d+\.\d+", allinfor[5]).group(0)
				ds=re.search("\d+\.\d+", allinfor[6]).group(0)
				dnds=re.search("\d+\.\d+", allinfor[4]).group(0)
				f1.write(genids[0]+"\t"+genids[1]+"\t"+dn+"\t"+ds+"\t"+dnds+"\n")

f1.close()