# Python 2.7.6
# vsearch_msaout_global_align_muscle.v2.py
# Mathias Scharmann


# 24 Feb 2015
# v2 21 Dec 2015: randomises cluster IDs / cluster order to avoid uneven coverages among reference intervals (slows down freebayes SNP calling)

# example usage:
# python 

# vcftools!
"""
prepare input like this:
vsearch --fasta_width 0 --threads 10 --id 0.9 --cluster_fast uniq.fasta --msaout vsearch_experiment.txt
"""


import argparse
import os
import random
#import numpy
import math
import subprocess
import multiprocessing as mp

# diese Funktion checkt ob ein file existiert und stoppt das script wenn ein file nicht exisitiert
def extant_file(x):
	"""
	'Type' for argparse - checks that file exists but does not open.
	"""
	
	if not os.path.exists(x):
		print "Error: {0} does not exist".format(x)
		exit()
	x = str(x)
	return x

#######


def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("--vsearch_msa", required=True, type=extant_file,
		help="a msaout file made with vsearch clustering", metavar="FILE")
		
	args = parser.parse_args()
	
	return args

	
def parse_rbdiv (vsearch_msa):
	
	cluster_dict = {}
	with open(vsearch_msa, "r") as INFILE:
		clstr = 0
		jump = False
		for line in INFILE:
#			print line
			
			if len(line) > 1:
			
				if jump:
					continue
				if line.startswith(">*"):
					clstr += 1
					same_clust = True
					continue
				if same_clust:
					if line.startswith(">"):
						continue
					elif line.startswith(">consensus"):
						jump = True
						continue
					else:
						try:
							cluster_dict[clstr].append(line.strip("\n"))
						except KeyError:
							cluster_dict[clstr] = [ line.strip("\n") ]
#	print cluster_dict	
	return cluster_dict	

def call_muscle (cluster_dict):
	
	print "there are {0} clusters".format(len(cluster_dict.keys()))
	
	aligned_clusters = {}
	
	
	
	for cluster in sorted(cluster_dict.keys(), key=int):
		
		if len(cluster_dict[cluster]) > 200:
			print "discarded cluster with {0} sequences".format(len(cluster_dict[cluster]))
			continue
			
		print "\rrunning muscle on cluster \t", cluster,
		
		with open("for_muscle.fasta", "w") as OUTFILE:
			for idx, seq in enumerate(cluster_dict[cluster]):
				OUTFILE.write(">"+str(idx)+"\n"+seq+"\n")
	
		bash_command = "muscle -quiet -in for_muscle.fasta" # -quiet

		p = subprocess.Popen(bash_command, shell=True, stdout=subprocess.PIPE) #, stderr=subprocess.STDOUT)
		
		rec = ""
		
		while True:
			line = p.stdout.readline()
#			print line,
			if line == '' and p.poll() != None:
				break
			if line.startswith(">"):
				new_rec = True
				if len(rec) > 0:
#					print len(rec)
					try:
						aligned_clusters[cluster].append(rec)
					except KeyError:
						aligned_clusters[cluster] = [rec]

				rec = ""
				continue
			else:
				new_rec = False
			if not new_rec: 	
				rec += line.strip("\n")	
				
			
#		print aligned_clusters[cluster]
						
	return aligned_clusters
	
def make_cluster_consensus(cluster_dict):
	
	# make consensus sequence out of individual sequences in the binned_dict:
	consensus_dict = {}
	for cluster, seqs in cluster_dict.items():
		outseq = []
		print "working on cluster {0}".format(cluster)
#		print seqs[0]
		if len(seqs) > 0:
			for column_index in range(0, len(seqs[0])): # to loop over all columns
				aligned_column = "".join([ seq[column_index] for seq in seqs ])
	#			print aligned_column
				maxcount = 0
				for nuc in ["A","T","G","C"]:
					if aligned_column.count(nuc) > maxcount:
						maxcount = aligned_column.count(nuc)
						consensus = nuc
				
			
				outseq.append(consensus)
		outseqstring = "".join(outseq)
#		print outseqstring
		consensus_dict[cluster] = outseqstring

	return consensus_dict

def write_fasta(fasta_dict):
	
	outlines = [ ]
	for contig in sorted(fasta_dict.keys(), key=int):
		seq = fasta_dict[contig]
		outlines.append( ">"+str(contig)+"_L"+str(len(seq)) )
		outlines.append(seq)
	
	with open("vsearch_muscle.fasta", "w") as OUTFILE:
		OUTFILE.write("\n".join(outlines))

def MT_wrapper_muscle(cluster_dict, nthreads):

	print "there are {0} clusters".format(len(cluster_dict.keys()))
	
	results = {}
	
	
	pool = mp.Pool(nthreads) #use all available cores, otherwise specify the number you want as an argument
	for clust in sorted(cluster_dict.keys(), key=int):
		if len(cluster_dict[clust]) > 100:
			print "discarded cluster with {0} sequences".format(len(cluster_dict[clust]))
			continue
		else:
			results[clust] = pool.apply_async(call_muscle_MT, args=(cluster_dict[clust], clust))
	pool.close()
	pool.join()
	
	# Get process results from the output queue
	#print output
	aligned_clusters = {}
	for idx, result in results.items():
#		print idx, result
		aligned_clusters[idx] = result.get()
#		print result.get()
	
	return aligned_clusters

def call_muscle_MT (seqs, clust):
	
	print "\rrunning muscle on cluster \t", clust,
	
	aligned_seqs = []
		
	with open("for_muscle_{0}.fasta".format(clust), "w") as OUTFILE:
		for idx, seq in enumerate(seqs):
			OUTFILE.write(">"+str(idx)+"\n"+seq+"\n")
	
	bash_command = "muscle -quiet -in for_muscle_{0}.fasta".format(clust) # -quiet

	p = subprocess.Popen(bash_command, shell=True, stdout=subprocess.PIPE) #, stderr=subprocess.STDOUT)
	
	rec = ""
	
	while True:
		line = p.stdout.readline()
#			print line,
		if line == '' and p.poll() != None:
			break
		if line.startswith(">"):
			new_rec = True
			if len(rec) > 0:
#					print len(rec)

				aligned_seqs.append(rec)

			rec = ""
			continue
		else:
			new_rec = False
		if not new_rec: 	
			rec += line.strip("\n")	
			
		
#		print aligned_clusters[cluster]
	
	os.remove("for_muscle_{0}.fasta".format(clust))
						
	return aligned_seqs

def randomise_cluster_ids(cluster_dict):
	
	old_IDs = cluster_dict.keys()
	new_IDs = random.sample(old_IDs, len(old_IDs))
	
	out_dict = {}
	for old, new in zip(old_IDs, new_IDs):
		out_dict[new] = cluster_dict[old]
	
	return out_dict

######### MAIN ########
		
args = 	get_commandline_arguments ()

cluster_dict = parse_rbdiv (args.vsearch_msa)
#print cluster_dict

# the following step prevents that clusters are sorted by coverage: if not randomising the ID, it is the order of clusters in the vsearch out format, which orders them by number of sequences in a cluster. When splitting the reference for freebayes, the first few chunks will contain much much more reads than the later ones and run >> longer. This uneveness is avoided if cluster IDs are randomised.
cluster_dict = randomise_cluster_ids(cluster_dict)


#aligned_clusters = call_muscle (cluster_dict)
nthreads = 30
aligned_clusters = MT_wrapper_muscle(cluster_dict, nthreads)

consensus_clusters = make_cluster_consensus(aligned_clusters)

write_fasta(consensus_clusters)


print "Done!"









