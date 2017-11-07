import glob

path=glob.glob("./branch.5.species.fasta.aligned/*_NT.fasta")

for files in path:
	with open(files) as f:
		ring=f.read()
		ring=ring.replace("!", "-")
		with open(files,"w") as f1:
			f1.write(ring)
			

