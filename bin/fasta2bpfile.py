#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
from Bio.SeqIO import parse
import os

###########################################################################################################
## python3 ./fasta2bpfile.py chr_10.MSA.fasta ADR INT INV MAC 5000 0.2 4 0.00000001 0.0000000001 100000  ##
## argument 1: name of the fasta file (ex: chr_10.MSA.fasta)                                             ##
## argument 2: name for species A (ex: ADR)                                                              ##
## argument 3: name for species B (ex: INT)                                                              ##
## argument 4: name for species C (ex: INV)                                                              ##
## argument 5: name for species D (ex: MAC)                                                              ##
## argument 6: size of the bin in bp (ex: 5000)                                                          ##
## argument 7: max proportion of missing data in the alignement (ex: 0.2)                                ##
## argument 8: minimum number of individuals in a given species at a given locus. (ex: 4)                ##
## argument 9: mutation rate mu per nucleotide and per generation (ex: 0.000000001)                      ##
## argument 10: recombinaration rate (ex: 0.0000000001 --> keep it small for big bins, i.e, 0.05mu       ##
## argument 11: size of the reference population (Nref) arbitrary fixed for ABC (ex: 100000)             ##
##                                                                                                       ##
## produces 3 output files:                                                                              ##
##         . general_infos.txt : informations about each bins (position, length, number of SNPs, etc...) ##
##         . observed_data.ms : the observed sequenced data in ms's format                               ##
##         . bpfile : used with observed_data.ms by mscalc to compute statistics                         ##
###########################################################################################################

bases = ['A', 'T', 'G', 'C']
def fasta2ms(aln, region):
	nTot = len(aln['spA']) + len(aln['spB']) + len(aln['spC']) + len(aln['spD'])
	nPos = len(aln['spA'][0])
	nPos_valid = 0
	nPos_biallelic = 0
	
	positions = []	
	alignements = []
	
	if region=='coding':
		liste_positions = range(2+3, nPos-3, 3) # from the third position of the second codon, to the codon before the last one
	else:
		liste_positions = range(nPos)
	
	for pos in liste_positions:
		pos_tmp = []
		for ind in aln['spA']:
			nuc = ind[pos]
			if nuc in bases:
				pos_tmp.append(nuc)
		for ind in aln['spB']:
			nuc = ind[pos]
			if nuc in bases:
				pos_tmp.append(nuc)
		for ind in aln['spC']:
			nuc = ind[pos]
			if nuc in bases:
				pos_tmp.append(nuc)
		for ind in aln['spD']:
			nuc = ind[pos]
			if nuc in bases:
				pos_tmp.append(nuc)
		
		if len(pos_tmp)==nTot:
			nPos_valid += 1
			nElements = set(pos_tmp)
			if len(nElements) == 2: # if biallelic
				positions.append(nPos_valid/(1.0*nPos))
				nPos_biallelic += 1
				nElements = [ i for i in nElements ]
				tmp = [ '0' if i==nElements[0] else '1' for i in pos_tmp ]
				alignements.append(tmp)

	res = {}
	res['nTot'] = nTot
	res['nPos'] = nPos
	res['nPos_valid'] = nPos_valid
	res['nPos_biallelic'] = nPos_biallelic
	res['positions'] = positions
	res['nA']=len(aln['spA'])
	res['nB']=len(aln['spB'])
	res['nC']=len(aln['spC'])
	res['nD']=len(aln['spD'])
	
	if nPos_biallelic > 0:
		seq = '\n//\nsegsites: {0}\npositions: {1}\n'.format(nPos_biallelic, ' '.join([ str(i) for i in positions ]))
		for ind in range(nTot):
			for pos in range(nPos_biallelic):
				seq += '{0}'.format(alignements[pos][ind])
			seq += '\n'
	else:
		seq = '\n//\nsegsites: 0\n'
	res['ms'] = seq
	return(res)


for tmp in sys.argv:
	arg=tmp.split('=')
	if arg[0]=='fastafile':
		fastafile=arg[1]
	if arg[0]=='region':
		region=arg[1] # coding ; noncoding
	if arg[0]=='spA':
		spA=arg[1]
	if arg[0]=='spB':
		spB=arg[1]
	if arg[0]=='spC':
		spC=arg[1]
	if arg[0]=='spD':
		spD=arg[1]
	if arg[0]=='maxN': # max prop. of N tolerated in the alignement
		maxN=float(arg[1])
	if arg[0]=='nMin': # minimum number of individuals in a given species at a given locus
		nMin=int(arg[1])
	if arg[0]=='Lmin': # minimum number of valid positions to retain a locus
		Lmin=int(arg[1])
	if arg[0]=='mu':
		mu=float(arg[1]) # mutation rate per generation and per bp
	if arg[0]=='rec':
		rec=float(arg[1]) # recombination rate per generation and per bp
	if arg[0]=='Nref':
		Nref=int(arg[1]) # Nref : 100000
	if arg[0]=='datapath':
		datapath=arg[1] # Nref : 100000
	if arg[0]=='binpath':
		binpath=arg[1] # Nref : 100000
	if arg[0]=='simulationpath':
		simulationpath=arg[1] # Nref : 100000

#fastafile = sys.argv[1] # name of the fasta file
#spA = sys.argv[2] # name of species A
#spB = sys.argv[3]
#spC = sys.argv[4]
#spD = sys.argv[5]
#maxN = float(sys.argv[6]) # max prop. of N tolerated in the alignement
#nMin = int(sys.argv[7]) # minimum number of individuals in a given species at a given locus
#mu = float(sys.argv[8]) # mutation rate per generation and per bp
#rec = float(sys.argv[9]) # recombination rate per generation and per bp
#Nref = int(sys.argv[10]) # Nref : 100000


#### TEST
#fastafile='100_RNAseq.fasta'
#spA='E1'
#spB='W1'
#spC='W2'
#spD='W3'
#maxN=0.5
#nMin=4
#mu=0.000000001
#rec=0.00000001
#Nref=100000

###

data = {}
infile = parse(fastafile, 'fasta')
for i in infile:
	tmp = i.id.split('|')
	gene = tmp[0]
	species = tmp[1]
	individual = tmp[2]
	allele = tmp[3]
	
	if gene not in data:
		data[gene]={}
	
	if species not in data[gene]:
		data[gene][species]=[]
	
	data[gene][species].append(i.seq)

outfile_infos = open('{datapath}/general_infos.txt'.format(datapath=datapath), 'w')
outfile_infos.write('gene\tcomment\tL\tL_valid\tnS\n')

outfile_ms = open('{datapath}/observed_data.ms'.format(datapath=datapath), 'w')
outfile_ms.write('msnsam 4pop\n111 222 333\n\n')

bpfile_header = "#Nref:{0}\tmu:{1}\tr:{2}\t".format(Nref, mu, rec)
bpfile_L1=''# L
bpfile_L2=''# nA
bpfile_L3=''# nB
bpfile_L4=''# nC
bpfile_L5=''# nD
bpfile_L6=''# theta
bpfile_L7=''# rho
bpfile_L8=''# nSNPs 

nRetainedLoci = 0
for gene in data: # loop over genes
	bin_tmp = {}
	bin_tmp['spA'] = []
	bin_tmp['spB'] = []
	bin_tmp['spC'] = []
	bin_tmp['spD'] = []
	nA = 0
	for seq in data[gene][spA]:
		seq_tmp = seq
		missing_data = (seq_tmp.count('N') + seq_tmp.count('n') + seq_tmp.count('?') + seq_tmp.count('-') + seq_tmp.count('_')) / (1.0 * len(seq_tmp))
		if missing_data < maxN:
			nA += 1
			bin_tmp['spA'].append(seq_tmp)

	nB = 0
	for seq in data[gene][spB]:
		seq_tmp = seq
		missing_data = (seq_tmp.count('N') + seq_tmp.count('n') + seq_tmp.count('?') + seq_tmp.count('-') + seq_tmp.count('_')) / (1.0 * len(seq_tmp))
		if missing_data < maxN:
			nB += 1
			bin_tmp['spB'].append(seq_tmp)

	nC = 0
	for seq in data[gene][spC]:
		seq_tmp = seq
		missing_data = (seq_tmp.count('N') + seq_tmp.count('n') + seq_tmp.count('?') + seq_tmp.count('-') + seq_tmp.count('_')) / (1.0 * len(seq_tmp))
		if missing_data < maxN:
			nC += 1
			bin_tmp['spC'].append(seq_tmp)

	nD = 0
	for seq in data[gene][spD]:
		seq_tmp = seq
		missing_data = (seq_tmp.count('N') + seq_tmp.count('n') + seq_tmp.count('?') + seq_tmp.count('-') + seq_tmp.count('_')) / (1.0 * len(seq_tmp))
		if missing_data < maxN:
			nD += 1
			bin_tmp['spD'].append(seq_tmp)


	if nA>=nMin and nB>=nMin and nC>=nMin and nD>=nMin:
		res = fasta2ms(bin_tmp, region)
		if res['nPos_valid']>=Lmin:
			nRetainedLoci += 1
			outfile_infos.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(gene, 'ok', res['nPos'], res['nPos_valid'], res['nPos_biallelic']))
			outfile_ms.write(res['ms'])
			bpfile_L1 += '{0}\t'.format(res['nPos_valid'])
			bpfile_L2 += '{0}\t'.format(res['nA'])
			bpfile_L3 += '{0}\t'.format(res['nB'])
			bpfile_L4 += '{0}\t'.format(res['nC'])
			bpfile_L5 += '{0}\t'.format(res['nD'])
			bpfile_L6 += '{0}\t'.format(4*Nref*mu*res['nPos_valid'])
			bpfile_L7 += '{0}\t'.format(4*Nref*rec*res['nPos_valid'])
			bpfile_L8 += '{0}\t'.format(res['nPos_biallelic'])
		else:
			outfile_infos.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(gene, 'rejected', res['nPos'], res['nPos_valid'], res['nPos_biallelic']))
	else:
		outfile_infos.write('{0}\t{1}\tNA\tNA\tNA\n'.format(gene, 'rejected'))
outfile_ms.write('\n')
outfile_infos.close()
outfile_ms.close()

bpfile_header += 'nLoci:{nLoci}'.format(nLoci=nRetainedLoci)
bpfile = open('{datapath}/bpfile'.format(datapath=datapath), 'w')
bpfile.write(bpfile_header.strip() + '\n')
bpfile.write(bpfile_L1.strip() + '\n')
bpfile.write(bpfile_L2.strip() + '\n')
bpfile.write(bpfile_L3.strip() + '\n')
bpfile.write(bpfile_L4.strip() + '\n')
bpfile.write(bpfile_L5.strip() + '\n')
bpfile.write(bpfile_L6.strip() + '\n')
bpfile.write(bpfile_L7.strip() + '\n')
bpfile.write(bpfile_L8.strip() + '\n')

bpfile.close()

### get ABCstat.txt
commande = 'cat {datapath}/observed_data.ms | pypy {binpath}/mscalc_4pop.py datapath={datapath} simulationpath={simulationpath}'.format(binpath=binpath, datapath=datapath, simulationpath=simulationpath)
os.system(commande)


