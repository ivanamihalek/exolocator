#!/usr/bin/python3


from config import Config
import os, subprocess
from random import uniform as rand
####################################################
def main():
	# we'll do this for human only, for now
	species = 'homo_sapiens'

	fasta_path = "/storage/databases/ensembl-{}/fasta/{}/pep".format(Config.release_number, species)
	in_fasta = "{}/canonical_peptides.fa".format(fasta_path)
	for dep in [fasta_path, in_fasta, Config.blastdbcmd, Config.muscle]:
		if not os.path.exists(dep):
			print(dep,"not found")
			exit()
	paralogue_file = "paralogues.txt"
	if  not os.path.exists(paralogue_file):
		print(paralogue_file, "not found")
		exit()

	inf = open(paralogue_file,"r")
	for line in inf:
		if rand(0,1)<0.9: continue
		ids = line.strip().split("\t")
		if len(ids)<3: continue
		aux = "test.ids"
		outf = open(aux,"w")
		outf.write("\n".join(ids)+"\n");
		outf.close()

		# extract the seqs
		musclein = "test.fa"
		# %f is sequence with fasta header
		cmd = "{} -db {} -dbtype prot -out {} -outfmt %f -entry_batch {} -logfile blastcmd2b.log".\
				format(Config.blastdbcmd, in_fasta, musclein,  aux)
		subprocess.call(["bash","-c", cmd])
		if not os.path.exists(musclein) or os.path.getsize(musclein)==0: continue

		muscleout = "test.afa"
		cmd = "{} -in {} -out {} -maxiters 3".format(Config.muscle, musclein, muscleout)
		subprocess.call(["bash","-c", cmd])
		exit()


#####################################################
if __name__=="__main__":
	main()
