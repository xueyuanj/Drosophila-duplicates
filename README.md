# Drosophila-duplicates
2016-2017 graduate project. Study the role of natural selection on Drosophila duplicates
For sequence based analysis, first parse the genome file
--> split.fasta.py or one_line.py

Generate the fasta files 
--> get.fasta.py (also can generate tree file)

Run alignment with macse (coding sequence only), or muscle if they're non-coding sequences.
--> automator.py

Remove the ! in alignment
--> remove_ex.py

The first part is Ka/Ks estimated using PAML
Convert the fasta file to phylip format
--> phylip_generator.py, format_convert.py

To get Ka, Ks or Ka/Ks ratio, run PAML
--> codeml.ctl, autopaml.py

Parse the PAML results
--> get.dnds.py

Keep in mind that all the jobs can be run on cluster $qsub -A rua15_collab *
--> cluster.job.txt
