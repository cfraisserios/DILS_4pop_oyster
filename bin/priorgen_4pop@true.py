#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from numpy.random import uniform
from numpy.random import binomial
from numpy.random import beta
from numpy.random import seed

#seed(100)

list_models = ['topo1', 'topo2', 'topo3']
list_migration = ['AB', 'AC', 'AD', 'BC', 'BD', 'CD']
list_mutation = ['mu', 'nSNPs']

help = "\t\033[1;31;40mTakes one model specifier, a number of multilocus simulations and six indicators of migration (one per possible pair of populations) as arguments:\033[0m\n\t"
help += "Accepted model specifiers (model=):\n\t\t"
help += "\n\t\t".join(list_models)
help += "\n\n"
help += "\tAccepted indicators of migration:\n\t\t"
help += "\n\t\t".join(list_migration)
help += "\n\n"
help += "\tAccepted mode of mutations (mutation=), can be setted to 'mu' or 'nSNPs':\n\t\t"
help += "\n\t\t".join(list_migration)
help += "\n\n"
help += "\t\033[1;32;40m#topo1 (A; ((B; C); D))\033[0m\n\tmsnsam tbs nSIMULATIONS -t tbs -r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 3 tbs -m 3 1 tbs -m 2 4 tbs -m 4 2 tbs -m 1 2 tbs -m 2 1 tbs -m 2 3 tbs -m 3 2 tbs -m 1 4 tbs -m 4 1 tbs -m 3 4 tbs -m 4 3 tbs -em tbs 2 3 0 -em tbs 3 2 0 -em tbs 3 1 0 -em tbs 1 3 0 -em tbs 3 4 0 -em tbs 4 3 0 -em tbs 2 4 0 -em tbs 4 2 0 -em tbs 1 2 0 -em tbs 2 1 0 -em tbs 4 1 0 -em tbs 1 4 0 -ej tbs 3 2 -en tbs 2 tbs -ej tbs 2 4 -en tbs 4 tbs -ej tbs 4 1 -eN tbs tbs\033[0m\n\n" # topo1 : checked and valid !
help += "\t\033[1;32;40m#topo2 (A; (B; (C; D)))\033[0m\n\tmsnsam tbs nSIMULATIONS -t tbs -r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 3 tbs -m 3 1 tbs -m 2 4 tbs -m 4 2 tbs -m 1 2 tbs -m 2 1 tbs -m 2 3 tbs -m 3 2 tbs -m 1 4 tbs -m 4 1 tbs -m 3 4 tbs -m 4 3 tbs -em tbs 3 4 0 -em tbs 4 3 0 -em tbs 3 2 0 -em tbs 2 3 0 -em tbs 3 1 0 -em tbs 1 3 0 -em tbs 2 4 0 -em tbs 4 2 0 -em tbs 2 1 0 -em tbs 1 2 0 -em tbs 4 1 0 -em tbs 1 4 0 -ej tbs 3 4 -en tbs 4 tbs -ej tbs 4 2 -en tbs 2 tbs -ej tbs 2 1 -eN tbs tbs\033[0m\n\n" # topo2 : checked and valid !
help += "\t\033[1;32;40m#topo3 ((A; B); (C; D))\033[0m\n\tmsnsam tbs nSIMULATIONS -t tbs -r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 1 4 tbs -m 4 1 tbs -m 2 3 tbs -m 3 2 tbs -m 2 4 tbs -m 4 2 tbs -m 3 4 tbs -m 4 3 tbs -em tbs 1 2 0 -em tbs 2 1 0 -em tbs 1 3 0 -em tbs 3 1 0 -em tbs 1 4 0 -em tbs 4 1 0 -em tbs 2 3 0 -em tbs 3 2 0 -em tbs 2 4 0 -em tbs 4 2 0 -em tbs 3 4 0 -em tbs 4 3 0 -ej tbs 2 1 -en tbs 1 tbs -ej tbs 4 3 -en tbs 3 tbs -ej tbs 3 1 -eN tbs tbs\033[0m\n\n" # topo3
help += "\t\033[1;32;40m#AB=1\033[0m\t\tSecondary contact between pop_A and pop_B, bidirectional A<->B\033[0m\n"
help += "\t\033[1;32;40m#AB=0\033[0m\t\tNo gene flow between pop_A and pop_B,  A | B\033[0m\n"
help += "\t\033[1;32;40m#BD=0\033[0m\t\tNo gene flow between pop_B and pop_D,  B | D\033[0m\n"
help += "\t\033[1;32;40m#Thus, for all pairs ({pairs}), set to 0 if you want isolation, or set to 1 if you want migration\033[0m\n\n".format(pairs="; ".join(list_migration))
help += "\t\033[1;32;40mExample: ./priorgen_4pop.py datapath=/home/project simulationpath=/home/project/topo1_iteration123 mutation=mu model=topo1 nMultilocus=1000 AB=0 AC=0 AD=1 BC=1 BD=0 CD=1 N=0 N=10 T=0 T=10 M=0.4 M=8 shape=1 shape=20\033[0m\n"
help += "\t\033[1;32;40mIsolation between A<->B, A<->C and B<->D, migration between A<->D, B<->C and C<->D \033[0m\n"
help += "\t\033[1;32;40mBe sure that the pathway to simulationpath (where the simulations will be written) already exists\033[0m\n"
help += "\t\033[1;32;40mThe arguments 'target' and 'posteriorFile' are optional. 'target' can take as values ['modelComp', 'posterior', 'opti1', 'opti2', 'opti3', 'opti4', 'opti5', 'gof'] to define the purpose of the simulations. 'posteriorFile' is important only if 'target' defines one the five cycles of optimization.\033[0m\n"


if len(sys.argv) != 22:
	print(help)
	sys.exit()

datapath=None
simulationpath=None
model=None
nMultilocus=None
AB=None
AC=None
AD=None
BC=None
BD=None
CD=None
target=None
posteriorFile=None
N_bound = []
T_bound = []
M_bound = []
shape_bound = []

for tmp in sys.argv:
	arg=tmp.split('=')
	if arg[0]=='datapath':
		datapath=arg[1]
	if arg[0]=='simulationpath':
		simulationpath=arg[1]
	if arg[0]=='mutation': # theta or nSNPs
		mutation=arg[1]
	if arg[0]=='model':
		model=arg[1]
		if model not in list_models:
			print("{model} is not an autorized model ({liste})".format(model=model, liste="; ".join(list_models)))
			sys.exit()
	if arg[0]=='nMultilocus':
		nMultilocus = int(arg[1])
	if arg[0]=='AB':
		AB = arg[1]
	if arg[0]=='AC':
		AC = arg[1]
	if arg[0]=='AD':
		AD = arg[1]
	if arg[0]=='BC':
		BC = arg[1]
	if arg[0]=='BD':
		BD = arg[1]
	if arg[0]=='CD':
		CD = arg[1]
	if arg[0]=='N':
		N_bound.append(float(arg[1]))
	if arg[0]=='T':
		T_bound.append(float(arg[1]))
	if arg[0]=='M':
		M_bound.append(float(arg[1]))
	if arg[0]=='shape':
		shape_bound.append(float(arg[1]))
	if arg[0]=='target':
		target=arg[1]
	if arg[0]=='posteriorFile':
		posteriorFile=arg[1]

if mutation==None:
	print("please specify a method to implement mutations: 'mu' (mutations conditionned on 4.N.mu) or 'nSNPs' (mutations conditionned on the number of nSNPs)\n")
	sys.exit()

if datapath==None:
	print("please specify a path to the bpfile (for instance: datapath=/home/project\n")
	sys.exit()

if simulationpath==None:
	print("please specify a path to write the priorfile (for instance: datapath=/home/project\n")
	sys.exit()

if model==None:
	print("ERROR: please specify a correct model (i.e model=topo1, model=topo2 or model=topo3)\n")
	sys.exit()
if nMultilocus==None:
	print("ERROR: please specify a number of simulations to perform (for instance: nMultilocus=100)\n")
	sys.exit()
if AB==None or AB not in ['0', '1']:
	print("ERROR: please specify if there is migration or not between pop. A and B (i.e: AB=0 or AB=1)\n")
	sys.exit()
if AC==None or AC not in ['0', '1']:
	print("ERROR: please specify if there is migration or not between pop. A and C (i.e: AC=0 or AC=1)\n")
	sys.exit()
if AD==None or AD not in ['0', '1']:
	print("ERROR: please specify if there is migration or not between pop. A and D (i.e: AD=0 or AD=1)\n")
	sys.exit()
if BC==None or BC not in ['0', '1']:
	print("ERROR: please specify if there is migration or not between pop. B and C (i.e: BC=0 or BC=1)\n")
	sys.exit()
if BD==None or BD not in ['0', '1']:
	print("ERROR: please specify if there is migration or not between pop. B and D (i.e: BD=0 or BD=1)\n")
	sys.exit()
if CD==None or CD not in ['0', '1']:
	print("ERROR: please specify if there is migration or not between pop. C and D (i.e: CD=0 or CD=1)\n")
	sys.exit()
if len(N_bound)!=2:
	print("ERROR: please specify one bound of the prior distribution for N in coalescent units (N=0 or N=10)\n")
	sys.exit()
if len(T_bound)!=2:
	print("ERROR: please specify one bound of the prior distribution for T in coalescent units (T=0 or T=10)\n")
	sys.exit()
if len(M_bound)!=2:
	print("ERROR: please specify one bound of the prior distribution for M (M=0 or M=10 with M=4.N.m)\n")
	sys.exit()
if len(shape_bound)!=2:
	print("ERROR: please specify one bound of the prior distribution for shape parameters of beta distributions (shape=1, shape=10, etc...)\n")
	sys.exit()

# get the posterior of a previous estimate to optimize, just in case
if target in ['opti1', 'opti2', 'opti3', 'opti4', 'opti5', 'gof']:
	posterior = {}
	infile = open(posteriorFile, 'r')
	tmp = infile.readline()
	for line in infile:
		line = line.strip().split('\t')
		if line[0] not in posterior:
			posterior[line[0]] = float(line[1])
	infile.close()
	
	if target == 'opti1':
		shapeOpti = 10
	if target == 'opti2':
		shapeOpti = 20
	if target == 'opti3':
		shapeOpti = 50
	if target == 'opti4':
		shapeOpti = 100
	if target == 'opti5':
		shapeOpti = 500

# time_space
time_space = 1

# read bpfile
infile = open("{datapath}/bpfile".format(datapath=datapath), "r")
tmp = infile.readline()

if mutation=='mu':
	L = [ float(i) for i in infile.readline().strip().split("\t") ]
else:
	L = [ int(i)+1 for i in infile.readline().strip().split("\t") ]
nsamA = [ int(i) for i in infile.readline().strip().split("\t") ]
nsamB = [ int(i) for i in infile.readline().strip().split("\t") ]
nsamC = [ int(i) for i in infile.readline().strip().split("\t") ]
nsamD = [ int(i) for i in infile.readline().strip().split("\t") ]
theta = [ float(i) for i in infile.readline().strip().split("\t") ]
rho = [ float(i) for i in infile.readline().strip().split("\t") ]
nSNPs = [ float(i) for i in infile.readline().strip().split("\t") ]
infile.close()

# number of loci
nLoci = len(L)

# sum of nsamA + nsamB
nsam_tot = [ nsamA[i] + nsamB[i] + nsamC[i] + nsamD[i] for i in range(nLoci) ]

# param multilocus: values that will be printed in priorfile.txt
## N = N_pop_i / Nref
if target in ['opti1', 'opti2', 'opti3', 'opti4', 'opti5', 'gof']:
	if target == 'gof':
		N1 = [ posterior['N1'] for i in range(nMultilocus) ]
		N2 = [ posterior['N2'] for i in range(nMultilocus) ]
		N3 = [ posterior['N3'] for i in range(nMultilocus) ]
		N4 = [ posterior['N4'] for i in range(nMultilocus) ]
		Na = [ posterior['Na'] for i in range(nMultilocus) ]
		if model == 'topo1':
			Na_23 = [ posterior['Na_23'] for i in range(nMultilocus) ]
			Na_24 = [ posterior['Na_24'] for i in range(nMultilocus) ]
		if model == 'topo2':
			Na_34 = [ posterior['Na_34'] for i in range(nMultilocus) ]
			Na_24 = [ posterior['Na_24'] for i in range(nMultilocus) ]
		if model == 'topo3':
			Na_12 = [ posterior['Na_12'] for i in range(nMultilocus) ]
			Na_34 = [ posterior['Na_34'] for i in range(nMultilocus) ]
		
		## factor of local reduction in Ne. Model of "background selection"
		shape_N_a = [ posterior['shape_N_a'] for i in range(nMultilocus) ]
		shape_N_b = [ posterior['shape_N_b'] for i in range(nMultilocus) ]
	else:
		N1 = [ posterior['N1'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
		N2 = [ posterior['N2'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
		N3 = [ posterior['N3'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
		N4 = [ posterior['N4'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
		Na = [ posterior['Na'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
		if model == 'topo1':
			Na_23 = [ posterior['Na_23'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			Na_24 = [ posterior['Na_24'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
		if model == 'topo2':
			Na_34 = [ posterior['Na_34'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			Na_24 = [ posterior['Na_24'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
		if model == 'topo3':
			Na_12 = [ posterior['Na_12'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			Na_34 = [ posterior['Na_34'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
		
		## factor of local reduction in Ne. Model of "background selection"
		shape_N_a = [ posterior['shape_N_a'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
		shape_N_b = [ posterior['shape_N_b'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
else:
	N1 = uniform(low = min(N_bound), high = max(N_bound), size = nMultilocus)
	N2 = uniform(low = min(N_bound), high = max(N_bound), size = nMultilocus)
	N3 = uniform(low = min(N_bound), high = max(N_bound), size = nMultilocus)
	N4 = uniform(low = min(N_bound), high = max(N_bound), size = nMultilocus)

	Na = uniform(low = min(N_bound), high = max(N_bound), size = nMultilocus) # N1N2-N3N4 

	if model == 'topo1':
		Na_23 = uniform(low = min(N_bound), high = max(N_bound), size = nMultilocus) # ancestor between N1 and N2
		Na_24 = uniform(low = min(N_bound), high = max(N_bound), size = nMultilocus) # ancestor between N3 and N4
	if model == 'topo2':
		Na_34 = uniform(low = min(N_bound), high = max(N_bound), size = nMultilocus) # ancestor between N1 and N2
		Na_24 = uniform(low = min(N_bound), high = max(N_bound), size = nMultilocus) # ancestor between N3 and N4
	if model == 'topo3':
		Na_12 = uniform(low = min(N_bound), high = max(N_bound), size = nMultilocus) # ancestor between N1 and N2
		Na_34 = uniform(low = min(N_bound), high = max(N_bound), size = nMultilocus) # ancestor between N1 and N2

	## factor of local reduction in Ne. Model of "background selection"
	shape_N_a = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
	shape_N_b = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)

## Migration rates
if target in ['opti1', 'opti2', 'opti3', 'opti4', 'opti5', 'gof']:
	if target == 'gof':
		if AB == '1':
			M12 = [ posterior['M12'] for i in range(nMultilocus) ]
			M21 = [ posterior['M21'] for i in range(nMultilocus) ]
			shape_M12_a = [ posterior['shape_M12_a'] for i in range(nMultilocus) ]
			shape_M12_b = [ posterior['shape_M12_b'] for i in range(nMultilocus) ]
			shape_M21_a = [ posterior['shape_M21_a'] for i in range(nMultilocus) ]
			shape_M21_b = [ posterior['shape_M21_b'] for i in range(nMultilocus) ]
		else:
			M12 = [0]*nMultilocus
			M21 = [0]*nMultilocus
			shape_M12_a = [0]*nMultilocus
			shape_M12_b = [0]*nMultilocus
			shape_M21_a = [0]*nMultilocus
			shape_M21_b = [0]*nMultilocus

		if AC == '1':
			M13 = [ posterior['M13'] for i in range(nMultilocus) ]
			M31 = [ posterior['M31'] for i in range(nMultilocus) ]
			shape_M13_a = [ posterior['shape_M13_a'] for i in range(nMultilocus) ]
			shape_M13_b = [ posterior['shape_M13_b'] for i in range(nMultilocus) ]
			shape_M31_a = [ posterior['shape_M31_a'] for i in range(nMultilocus) ]
			shape_M31_b = [ posterior['shape_M31_b'] for i in range(nMultilocus) ]
		else:
			M13 = [0]*nMultilocus
			M31 = [0]*nMultilocus
			shape_M13_a = [0]*nMultilocus
			shape_M13_b = [0]*nMultilocus
			shape_M31_a = [0]*nMultilocus
			shape_M31_b = [0]*nMultilocus

		if AD == '1':
			M14 = [ posterior['M14'] for i in range(nMultilocus) ]
			M41 = [ posterior['M41'] for i in range(nMultilocus) ]
			shape_M14_a = [ posterior['shape_M14_a'] for i in range(nMultilocus) ]
			shape_M14_b = [ posterior['shape_M14_b'] for i in range(nMultilocus) ]
			shape_M41_a = [ posterior['shape_M41_a'] for i in range(nMultilocus) ]
			shape_M41_b = [ posterior['shape_M41_b'] for i in range(nMultilocus) ]
		else:
			M14 = [0]*nMultilocus
			M41 = [0]*nMultilocus
			shape_M14_a = [0]*nMultilocus
			shape_M14_b = [0]*nMultilocus
			shape_M41_a = [0]*nMultilocus
			shape_M41_b = [0]*nMultilocus

		if BC == '1':
			M23 = [ posterior['M23'] for i in range(nMultilocus) ]
			M32 = [ posterior['M32'] for i in range(nMultilocus) ]
			shape_M23_a = [ posterior['shape_M23_a'] for i in range(nMultilocus) ]
			shape_M23_b = [ posterior['shape_M23_b'] for i in range(nMultilocus) ]
			shape_M32_a = [ posterior['shape_M32_a'] for i in range(nMultilocus) ]
			shape_M32_b = [ posterior['shape_M32_b'] for i in range(nMultilocus) ]
		else:
			M23 = [0]*nMultilocus
			M32 = [0]*nMultilocus
			shape_M32_a = [0]*nMultilocus
			shape_M32_b = [0]*nMultilocus
			shape_M23_a = [0]*nMultilocus
			shape_M23_b = [0]*nMultilocus

		if BD == '1':
			M24 = [ posterior['M24'] for i in range(nMultilocus) ]
			M42 = [ posterior['M42'] for i in range(nMultilocus) ]
			shape_M24_a = [ posterior['shape_M24_a'] for i in range(nMultilocus) ]
			shape_M24_b = [ posterior['shape_M24_b'] for i in range(nMultilocus) ]
			shape_M42_a = [ posterior['shape_M42_a'] for i in range(nMultilocus) ]
			shape_M42_b = [ posterior['shape_M42_b'] for i in range(nMultilocus) ]
		else:
			M24 = [0]*nMultilocus
			M42 = [0]*nMultilocus
			shape_M24_a = [0]*nMultilocus
			shape_M24_b = [0]*nMultilocus
			shape_M42_a = [0]*nMultilocus
			shape_M42_b = [0]*nMultilocus

		if CD == '1':
			M34 = [ posterior['M34'] for i in range(nMultilocus) ]
			M43 = [ posterior['M43'] for i in range(nMultilocus) ]
			shape_M34_a = [ posterior['shape_M34_a'] for i in range(nMultilocus) ]
			shape_M34_b = [ posterior['shape_M34_b'] for i in range(nMultilocus) ]
			shape_M43_a = [ posterior['shape_M43_a'] for i in range(nMultilocus) ]
			shape_M43_b = [ posterior['shape_M43_b'] for i in range(nMultilocus) ]
		else:
			M34 = [0]*nMultilocus
			M43 = [0]*nMultilocus
			shape_M34_a = [0]*nMultilocus
			shape_M34_b = [0]*nMultilocus
			shape_M43_a = [0]*nMultilocus
			shape_M43_b = [0]*nMultilocus
	
	else:
		if AB == '1':
			M12 = [ posterior['M12'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			M21 = [ posterior['M21'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M12_a = [ posterior['shape_M12_a'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M12_b = [ posterior['shape_M12_b'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M21_a = [ posterior['shape_M21_a'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M21_b = [ posterior['shape_M21_b'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
		else:
			M12 = [0]*nMultilocus
			M21 = [0]*nMultilocus
			shape_M12_a = [0]*nMultilocus
			shape_M12_b = [0]*nMultilocus
			shape_M21_a = [0]*nMultilocus
			shape_M21_b = [0]*nMultilocus

		if AC == '1':
			M13 = [ posterior['M13'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			M31 = [ posterior['M31'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M13_a = [ posterior['shape_M13_a'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M13_b = [ posterior['shape_M13_b'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M31_a = [ posterior['shape_M31_a'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M31_b = [ posterior['shape_M31_b'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
		else:
			M13 = [0]*nMultilocus
			M31 = [0]*nMultilocus
			shape_M13_a = [0]*nMultilocus
			shape_M13_b = [0]*nMultilocus
			shape_M31_a = [0]*nMultilocus
			shape_M31_b = [0]*nMultilocus

		if AD == '1':
			M14 = [ posterior['M14'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			M41 = [ posterior['M41'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M14_a = [ posterior['shape_M14_a'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M14_b = [ posterior['shape_M14_b'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M41_a = [ posterior['shape_M41_a'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M41_b = [ posterior['shape_M41_b'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
		else:
			M14 = [0]*nMultilocus
			M41 = [0]*nMultilocus
			shape_M14_a = [0]*nMultilocus
			shape_M14_b = [0]*nMultilocus
			shape_M41_a = [0]*nMultilocus
			shape_M41_b = [0]*nMultilocus

		if BC == '1':
			M23 = [ posterior['M23'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			M32 = [ posterior['M32'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M23_a = [ posterior['shape_M23_a'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M23_b = [ posterior['shape_M23_b'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M32_a = [ posterior['shape_M32_a'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M32_b = [ posterior['shape_M32_b'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
		else:
			M23 = [0]*nMultilocus
			M32 = [0]*nMultilocus
			shape_M32_a = [0]*nMultilocus
			shape_M32_b = [0]*nMultilocus
			shape_M23_a = [0]*nMultilocus
			shape_M23_b = [0]*nMultilocus

		if BD == '1':
			M24 = [ posterior['M24'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			M42 = [ posterior['M42'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M24_a = [ posterior['shape_M24_a'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M24_b = [ posterior['shape_M24_b'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M42_a = [ posterior['shape_M42_a'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M42_b = [ posterior['shape_M42_b'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
		else:
			M24 = [0]*nMultilocus
			M42 = [0]*nMultilocus
			shape_M24_a = [0]*nMultilocus
			shape_M24_b = [0]*nMultilocus
			shape_M42_a = [0]*nMultilocus
			shape_M42_b = [0]*nMultilocus

		if CD == '1':
			M34 = [ posterior['M34'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			M43 = [ posterior['M43'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M34_a = [ posterior['shape_M34_a'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M34_b = [ posterior['shape_M34_b'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M43_a = [ posterior['shape_M43_a'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			shape_M43_b = [ posterior['shape_M43_b'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
		else:
			M34 = [0]*nMultilocus
			M43 = [0]*nMultilocus
			shape_M34_a = [0]*nMultilocus
			shape_M34_b = [0]*nMultilocus
			shape_M43_a = [0]*nMultilocus
			shape_M43_b = [0]*nMultilocus
else:
	if AB == '1':
		M12 = uniform(low = min(M_bound), high = max(M_bound), size = nMultilocus)
		M21 = uniform(low = min(M_bound), high = max(M_bound), size = nMultilocus)
		shape_M12_a = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
		shape_M12_b = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
		shape_M21_a = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
		shape_M21_b = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
	else:
		M12 = [0]*nMultilocus
		M21 = [0]*nMultilocus
		shape_M12_a = [0]*nMultilocus
		shape_M12_b = [0]*nMultilocus
		shape_M21_a = [0]*nMultilocus
		shape_M21_b = [0]*nMultilocus

	if AC == '1':
		M13 = uniform(low = min(M_bound), high = max(M_bound), size = nMultilocus)
		M31 = uniform(low = min(M_bound), high = max(M_bound), size = nMultilocus)
		shape_M13_a = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
		shape_M13_b = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
		shape_M31_a = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
		shape_M31_b = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
	else:
		M13 = [0]*nMultilocus
		M31 = [0]*nMultilocus
		shape_M13_a = [0]*nMultilocus
		shape_M13_b = [0]*nMultilocus
		shape_M31_a = [0]*nMultilocus
		shape_M31_b = [0]*nMultilocus

	if AD == '1':
		M14 = uniform(low = min(M_bound), high = max(M_bound), size = nMultilocus)
		M41 = uniform(low = min(M_bound), high = max(M_bound), size = nMultilocus)
		shape_M14_a = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
		shape_M14_b = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
		shape_M41_a = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
		shape_M41_b = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
	else:
		M14 = [0]*nMultilocus
		M41 = [0]*nMultilocus
		shape_M14_a = [0]*nMultilocus
		shape_M14_b = [0]*nMultilocus
		shape_M41_a = [0]*nMultilocus
		shape_M41_b = [0]*nMultilocus

	if BC == '1':
		M23 = uniform(low = min(M_bound), high = max(M_bound), size = nMultilocus)
		M32 = uniform(low = min(M_bound), high = max(M_bound), size = nMultilocus)
		shape_M32_a = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
		shape_M32_b = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
		shape_M23_a = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
		shape_M23_b = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
	else:
		M23 = [0]*nMultilocus
		M32 = [0]*nMultilocus
		shape_M32_a = [0]*nMultilocus
		shape_M32_b = [0]*nMultilocus
		shape_M23_a = [0]*nMultilocus
		shape_M23_b = [0]*nMultilocus

	if BD == '1':
		M24 = uniform(low = min(M_bound), high = max(M_bound), size = nMultilocus)
		M42 = uniform(low = min(M_bound), high = max(M_bound), size = nMultilocus)
		shape_M24_a = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
		shape_M24_b = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
		shape_M42_a = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
		shape_M42_b = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
	else:
		M24 = [0]*nMultilocus
		M42 = [0]*nMultilocus
		shape_M24_a = [0]*nMultilocus
		shape_M24_b = [0]*nMultilocus
		shape_M42_a = [0]*nMultilocus
		shape_M42_b = [0]*nMultilocus

	if CD == '1':
		M34 = uniform(low = min(M_bound), high = max(M_bound), size = nMultilocus)
		M43 = uniform(low = min(M_bound), high = max(M_bound), size = nMultilocus)
		shape_M34_a = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
		shape_M34_b = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
		shape_M43_a = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
		shape_M43_b = uniform(low = min(shape_bound), high=max(shape_bound), size = nMultilocus)
	else:
		M34 = [0]*nMultilocus
		M43 = [0]*nMultilocus
		shape_M34_a = [0]*nMultilocus
		shape_M34_b = [0]*nMultilocus
		shape_M43_a = [0]*nMultilocus
		shape_M43_b = [0]*nMultilocus


## times
epsilon = 0.00001

if target in ['opti1', 'opti2', 'opti3', 'opti4', 'opti5', 'gof']:
	if target == 'gof':
		Tsplit = [ posterior['Tsplit'] for i in range(nMultilocus) ]
	else:
		Tsplit = [ posterior['Tsplit'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]

else:
	Tsplit = uniform(low = T_bound[0], high = T_bound[1], size = nMultilocus)

if model == 'topo1':
	if target in ['opti1', 'opti2', 'opti3', 'opti4', 'opti5', 'gof']:
		if target == 'gof':
			Tsplit_24 = [ posterior['Tsplit_24']  for i in range(nMultilocus) ]
			Tsplit_24 = [ Tsplit_24[i] if Tsplit_24[i]<Tsplit[i] else Tsplit[i]-epsilon for i in range(nMultilocus) ]
			Tsplit_23 = [ posterior['Tsplit_23']  for i in range(nMultilocus) ]
			Tsplit_23 = [ Tsplit_23[i] if Tsplit_23[i]<Tsplit_24[i] else Tsplit_24[i]-epsilon for i in range(nMultilocus) ]
			if AB == '1':
				Tsc_12 = [ posterior['Tsc_12']  for i in range(nMultilocus) ]
				Tsc_12 = [ Tsc_12[i] if Tsc_12[i]<Tsplit_23[i] else Tsplit_23[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_12 = [0] * nMultilocus
			if AC == '1':
	#			Tsc_13 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_13 = [ posterior['Tsc_13']  for i in range(nMultilocus) ]
				Tsc_13 = [ Tsc_13[i] if Tsc_13[i]<Tsplit_23[i] else Tsplit_23[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_13 = [0] * nMultilocus
			if AD == '1':
	#			Tsc_14 = [ uniform(low = T_bound[0], high = time_space*Tsplit_24[i], size = None) for i in range(nMultilocus) ]
				Tsc_14 = [ posterior['Tsc_14']  for i in range(nMultilocus) ]
				Tsc_14 = [ Tsc_14[i] if Tsc_14[i]<Tsplit_23[i] else Tsplit_24[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_14 = [0] * nMultilocus
			if BC == '1':
	#			Tsc_23 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_23 = [ posterior['Tsc_23']  for i in range(nMultilocus) ]
				Tsc_23 = [ Tsc_23[i] if Tsc_23[i]<Tsplit_23[i] else Tsplit_23[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_23 = [0] * nMultilocus
			if BD == '1':
	#			Tsc_24 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ] # maintenant limite par Tsplit_34 et non plus Tsplit_24
				Tsc_24 = [ posterior['Tsc_24']  for i in range(nMultilocus) ]
				Tsc_24 = [ Tsc_24[i] if Tsc_24[i]<Tsplit_23[i] else Tsplit_23[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_24 = [0] * nMultilocus
			if CD == '1':
				Tsc_34 = [ posterior['Tsc_34']  for i in range(nMultilocus) ]
				Tsc_34 = [ Tsc_34[i] if Tsc_34[i]<Tsplit_23[i] else Tsplit_23[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_34 = [0] * nMultilocus
		else:
			Tsplit_24 = [ posterior['Tsplit_24'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ] # [ uniform(low = T_bound[0], high = time_space*Tsplit[i], size = None) for i in range(nMultilocus) ] # crx
			Tsplit_24 = [ Tsplit_24[i] if Tsplit_24[i]<Tsplit[i] else Tsplit[i]-epsilon for i in range(nMultilocus) ]
			Tsplit_23 = [ posterior['Tsplit_23'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ] # [ uniform(low = T_bound[0], high = time_space*Tsplit_24[i], size = None) for i in range(nMultilocus) ] # crx
			Tsplit_23 = [ Tsplit_23[i] if Tsplit_23[i]<Tsplit_24[i] else Tsplit_24[i]-epsilon for i in range(nMultilocus) ]
			if AB == '1':
				Tsc_12 = [ posterior['Tsc_12'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ] # Tsc_12 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_12 = [ Tsc_12[i] if Tsc_12[i]<Tsplit_23[i] else Tsplit_23[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_12 = [0] * nMultilocus
			if AC == '1':
	#			Tsc_13 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_13 = [ posterior['Tsc_13'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ] # Tsc_13 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_13 = [ Tsc_13[i] if Tsc_13[i]<Tsplit_23[i] else Tsplit_23[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_13 = [0] * nMultilocus
			if AD == '1':
	#			Tsc_14 = [ uniform(low = T_bound[0], high = time_space*Tsplit_24[i], size = None) for i in range(nMultilocus) ]
				Tsc_14 = [ posterior['Tsc_14'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ] # Tsc_14 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_14 = [ Tsc_14[i] if Tsc_14[i]<Tsplit_23[i] else Tsplit_24[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_14 = [0] * nMultilocus
			if BC == '1':
	#			Tsc_23 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_23 = [ posterior['Tsc_23'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ] # Tsc_23 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_23 = [ Tsc_23[i] if Tsc_23[i]<Tsplit_23[i] else Tsplit_23[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_23 = [0] * nMultilocus
			if BD == '1':
	#			Tsc_24 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ] # maintenant limite par Tsplit_34 et non plus Tsplit_24
				Tsc_24 = [ posterior['Tsc_24'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ] # Tsc_24 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_24 = [ Tsc_24[i] if Tsc_24[i]<Tsplit_23[i] else Tsplit_23[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_24 = [0] * nMultilocus
			if CD == '1':
				Tsc_34 = [ posterior['Tsc_34'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ] # Tsc_34 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_34 = [ Tsc_34[i] if Tsc_34[i]<Tsplit_23[i] else Tsplit_23[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_34 = [0] * nMultilocus
	else:
		Tsplit_24 = [ uniform(low = T_bound[0], high = time_space*Tsplit[i], size = None) for i in range(nMultilocus) ] # crx
		Tsplit_23 = [ uniform(low = T_bound[0], high = time_space*Tsplit_24[i], size = None) for i in range(nMultilocus) ] # crx
		if AB == '1':
			Tsc_12 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
		else:
			Tsc_12 = [0] * nMultilocus
		if AC == '1':
			Tsc_13 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
		else:
			Tsc_13 = [0] * nMultilocus
		if AD == '1':
			Tsc_14 = [ uniform(low = T_bound[0], high = time_space*Tsplit_24[i], size = None) for i in range(nMultilocus) ]
		else:
			Tsc_14 = [0] * nMultilocus
		if BC == '1':
			Tsc_23 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
		else:
			Tsc_23 = [0] * nMultilocus
		if BD == '1':
			Tsc_24 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ] # maintenant limite par Tsplit_34 et non plus Tsplit_24
		else:
			Tsc_24 = [0] * nMultilocus
		if CD == '1':
			Tsc_34 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
		else:
			Tsc_34 = [0] * nMultilocus

if model == 'topo2':
	if target in ['opti1', 'opti2', 'opti3', 'opti4', 'opti5', 'gof']:
		if target == 'gof':
			Tsplit_24 = [ posterior['Tsplit_24'] for i in range(nMultilocus) ]
			Tsplit_24 = [ Tsplit_24[i] if Tsplit_24[i]<Tsplit[i] else Tsplit[i]-epsilon for i in range(nMultilocus) ]
			Tsplit_34 = [ posterior['Tsplit_34'] for i in range(nMultilocus) ]
			Tsplit_34 = [ Tsplit_34[i] if Tsplit_34[i]<Tsplit_24[i] else Tsplit_24[i]-epsilon for i in range(nMultilocus) ]
			if AB == '1':
				Tsc_12 = [ posterior['Tsc_12'] for i in range(nMultilocus) ]
				Tsc_12 = [ Tsc_12[i] if Tsc_12[i]<Tsplit_24[i] else Tsplit_23[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_12 = [0] * nMultilocus
			if AC == '1':
	#			Tsc_13 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_13 = [ posterior['Tsc_13'] for i in range(nMultilocus) ]
				Tsc_13 = [ Tsc_13[i] if Tsc_13[i]<Tsplit_34[i] else Tsplit_23[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_13 = [0] * nMultilocus
			if AD == '1':
	#			Tsc_14 = [ uniform(low = T_bound[0], high = time_space*Tsplit_24[i], size = None) for i in range(nMultilocus) ]
				Tsc_14 = [ posterior['Tsc_14'] for i in range(nMultilocus) ]
				Tsc_14 = [ Tsc_14[i] if Tsc_14[i]<Tsplit_34[i] else Tsplit_24[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_14 = [0] * nMultilocus
			if BC == '1':
	#			Tsc_23 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_23 = [ posterior['Tsc_23'] for i in range(nMultilocus) ]
				Tsc_23 = [ Tsc_23[i] if Tsc_23[i]<Tsplit_34[i] else Tsplit_23[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_23 = [0] * nMultilocus
			if BD == '1':
	#			Tsc_24 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ] # maintenant limite par Tsplit_34 et non plus Tsplit_24
				Tsc_24 = [ posterior['Tsc_24'] for i in range(nMultilocus) ]
				Tsc_24 = [ Tsc_24[i] if Tsc_24[i]<Tsplit_34[i] else Tsplit_23[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_24 = [0] * nMultilocus
			if CD == '1':
				Tsc_34 = [ posterior['Tsc_34'] for i in range(nMultilocus) ]
				Tsc_34 = [ Tsc_34[i] if Tsc_34[i]<Tsplit_34[i] else Tsplit_23[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_34 = [0] * nMultilocus

		else:
			Tsplit_24 = [ posterior['Tsplit_24'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			Tsplit_24 = [ Tsplit_24[i] if Tsplit_24[i]<Tsplit[i] else Tsplit[i]-epsilon for i in range(nMultilocus) ]
			Tsplit_34 = [ posterior['Tsplit_34'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			Tsplit_34 = [ Tsplit_34[i] if Tsplit_34[i]<Tsplit_24[i] else Tsplit_24[i]-epsilon for i in range(nMultilocus) ]
			if AB == '1':
				Tsc_12 = [ posterior['Tsc_12'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ] # Tsc_12 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_12 = [ Tsc_12[i] if Tsc_12[i]<Tsplit_24[i] else Tsplit_23[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_12 = [0] * nMultilocus
			if AC == '1':
	#			Tsc_13 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_13 = [ posterior['Tsc_13'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ] # Tsc_13 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_13 = [ Tsc_13[i] if Tsc_13[i]<Tsplit_34[i] else Tsplit_23[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_13 = [0] * nMultilocus
			if AD == '1':
	#			Tsc_14 = [ uniform(low = T_bound[0], high = time_space*Tsplit_24[i], size = None) for i in range(nMultilocus) ]
				Tsc_14 = [ posterior['Tsc_14'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ] # Tsc_14 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_14 = [ Tsc_14[i] if Tsc_14[i]<Tsplit_34[i] else Tsplit_24[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_14 = [0] * nMultilocus
			if BC == '1':
	#			Tsc_23 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_23 = [ posterior['Tsc_23'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ] # Tsc_23 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_23 = [ Tsc_23[i] if Tsc_23[i]<Tsplit_34[i] else Tsplit_23[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_23 = [0] * nMultilocus
			if BD == '1':
	#			Tsc_24 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ] # maintenant limite par Tsplit_34 et non plus Tsplit_24
				Tsc_24 = [ posterior['Tsc_24'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ] # Tsc_24 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_24 = [ Tsc_24[i] if Tsc_24[i]<Tsplit_34[i] else Tsplit_23[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_24 = [0] * nMultilocus
			if CD == '1':
				Tsc_34 = [ posterior['Tsc_34'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ] # Tsc_34 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_34 = [ Tsc_34[i] if Tsc_34[i]<Tsplit_34[i] else Tsplit_23[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_34 = [0] * nMultilocus
	else:
		Tsplit_24 = [ uniform(low = T_bound[0], high = time_space*Tsplit[i], size = None) for i in range(nMultilocus) ]
		Tsplit_34 = [ uniform(low = T_bound[0], high = time_space*Tsplit_24[i], size = None) for i in range(nMultilocus) ]
		if AB == '1':
			Tsc_12 = [ uniform(low = T_bound[0], high = time_space*Tsplit_24[i], size = None) for i in range(nMultilocus) ]
		else:
			Tsc_12 = [0] * nMultilocus
		if AC == '1':
			Tsc_13 = [ uniform(low = T_bound[0], high = time_space*Tsplit_34[i], size = None) for i in range(nMultilocus) ]
		else:
			Tsc_13 = [0] * nMultilocus
		if AD == '1':
			Tsc_14 = [ uniform(low = T_bound[0], high = time_space*Tsplit_34[i], size = None) for i in range(nMultilocus) ]
		else:
			Tsc_14 = [0] * nMultilocus
		if BC == '1':
			Tsc_23 = [ uniform(low = T_bound[0], high = time_space*Tsplit_34[i], size = None) for i in range(nMultilocus) ]
		else:
			Tsc_23 = [0] * nMultilocus
		if BD == '1':
			Tsc_24 = [ uniform(low = T_bound[0], high = time_space*Tsplit_34[i], size = None) for i in range(nMultilocus) ] # maintenant limite par Tsplit_34 et non plus Tsplit_24
		else:
			Tsc_24 = [0] * nMultilocus
		if CD == '1':
			Tsc_34 = [ uniform(low = T_bound[0], high = time_space*Tsplit_34[i], size = None) for i in range(nMultilocus) ]
		else:
			Tsc_34 = [0] * nMultilocus

if model == 'topo3':
	if target in ['opti1', 'opti2', 'opti3', 'opti4', 'opti5', 'gof']:
		if target == 'gof':
			Tsplit_12 = [ posterior['Tsplit_12'] for i in range(nMultilocus) ]
			Tsplit_12 = [ Tsplit_12[i] if Tsplit_12[i]<Tsplit[i] else Tsplit[i]-epsilon for i in range(nMultilocus) ]
			Tsplit_34 = [ posterior['Tsplit_34'] for i in range(nMultilocus) ]
			Tsplit_34 = [ Tsplit_34[i] if Tsplit_34[i]<Tsplit[i] else Tsplit[i]-epsilon for i in range(nMultilocus) ]
			if AB == '1':
				Tsc_12 = [ posterior['Tsc_12'] for i in range(nMultilocus) ]
				Tsc_12 = [ Tsc_12[i] if Tsc_12[i]<Tsplit_24[i] else Tsplit_12[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_12 = [0] * nMultilocus
			if AC == '1':
	#			Tsc_13 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_13 = [ posterior['Tsc_13'] for i in range(nMultilocus) ]
				Tsc_13 = [ Tsc_13[i] if Tsc_13[i]<Tsplit_34[i] else Tsplit_34[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_13 = [0] * nMultilocus
			if AD == '1':
	#			Tsc_14 = [ uniform(low = T_bound[0], high = time_space*Tsplit_24[i], size = None) for i in range(nMultilocus) ]
				Tsc_14 = [ posterior['Tsc_14'] for i in range(nMultilocus) ]
				Tsc_14 = [ Tsc_14[i] if Tsc_14[i]<min([Tsplit_12[i], Tsplit_34[i]]) else min([Tsplit_12[i], Tsplit_34[i]])-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_14 = [0] * nMultilocus
			if BC == '1':
	#			Tsc_23 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_23 = [ posterior['Tsc_23'] for i in range(nMultilocus) ]
				Tsc_23 = [ Tsc_23[i] if Tsc_23[i]<min([Tsplit_12[i], Tsplit_34[i]]) else min([Tsplit_12[i], Tsplit_34[i]])-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_23 = [0] * nMultilocus
			if BD == '1':
	#			Tsc_24 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ] # maintenant limite par Tsplit_34 et non plus Tsplit_24
				Tsc_24 = [ posterior['Tsc_24'] for i in range(nMultilocus) ]
				Tsc_24 = [ Tsc_24[i] if Tsc_24[i]<min([Tsplit_12[i], Tsplit_34[i]]) else min([Tsplit_12[i], Tsplit_34[i]])-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_24 = [0] * nMultilocus
			if CD == '1':
				Tsc_34 = [ posterior['Tsc_34'] for i in range(nMultilocus) ]
				Tsc_34 = [ Tsc_34[i] if Tsc_34[i]<min([Tsplit_12[i], Tsplit_34[i]]) else min([Tsplit_12[i], Tsplit_34[i]])-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_34 = [0] * nMultilocus
		else:
			Tsplit_12 = [ posterior['Tsplit_12'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			Tsplit_12 = [ Tsplit_12[i] if Tsplit_12[i]<Tsplit[i] else Tsplit[i]-epsilon for i in range(nMultilocus) ]
			Tsplit_34 = [ posterior['Tsplit_34'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ]
			Tsplit_34 = [ Tsplit_34[i] if Tsplit_34[i]<Tsplit[i] else Tsplit[i]-epsilon for i in range(nMultilocus) ]
			if AB == '1':
				Tsc_12 = [ posterior['Tsc_12'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ] # Tsc_12 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_12 = [ Tsc_12[i] if Tsc_12[i]<Tsplit_24[i] else Tsplit_12[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_12 = [0] * nMultilocus
			if AC == '1':
	#			Tsc_13 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_13 = [ posterior['Tsc_13'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ] # Tsc_13 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_13 = [ Tsc_13[i] if Tsc_13[i]<Tsplit_34[i] else Tsplit_34[i]-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_13 = [0] * nMultilocus
			if AD == '1':
	#			Tsc_14 = [ uniform(low = T_bound[0], high = time_space*Tsplit_24[i], size = None) for i in range(nMultilocus) ]
				Tsc_14 = [ posterior['Tsc_14'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ] # Tsc_14 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_14 = [ Tsc_14[i] if Tsc_14[i]<min([Tsplit_12[i], Tsplit_34[i]]) else min([Tsplit_12[i], Tsplit_34[i]])-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_14 = [0] * nMultilocus
			if BC == '1':
	#			Tsc_23 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_23 = [ posterior['Tsc_23'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ] # Tsc_23 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_23 = [ Tsc_23[i] if Tsc_23[i]<min([Tsplit_12[i], Tsplit_34[i]]) else min([Tsplit_12[i], Tsplit_34[i]])-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_23 = [0] * nMultilocus
			if BD == '1':
	#			Tsc_24 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ] # maintenant limite par Tsplit_34 et non plus Tsplit_24
				Tsc_24 = [ posterior['Tsc_24'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ] # Tsc_24 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_24 = [ Tsc_24[i] if Tsc_24[i]<min([Tsplit_12[i], Tsplit_34[i]]) else min([Tsplit_12[i], Tsplit_34[i]])-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_24 = [0] * nMultilocus
			if CD == '1':
				Tsc_34 = [ posterior['Tsc_34'] * i/0.5 for i in beta(a=shapeOpti, b=shapeOpti, size=nMultilocus) ] # Tsc_34 = [ uniform(low = T_bound[0], high = time_space*Tsplit_23[i], size = None) for i in range(nMultilocus) ]
				Tsc_34 = [ Tsc_34[i] if Tsc_34[i]<min([Tsplit_12[i], Tsplit_34[i]]) else min([Tsplit_12[i], Tsplit_34[i]])-epsilon for i in range(nMultilocus) ]
			else:
				Tsc_34 = [0] * nMultilocus
	else:
		Tsplit_12 = [ uniform(low = T_bound[0], high = time_space*Tsplit[i], size = None) for i in range(nMultilocus) ]
		Tsplit_34 = [ uniform(low = T_bound[0], high = time_space*Tsplit[i], size = None) for i in range(nMultilocus) ]
		if AB == '1':
			Tsc_12 = [ uniform(low = T_bound[0], high = time_space*Tsplit_12[i], size = None) for i in range(nMultilocus) ]
		else:
			Tsc_12 = [0] * nMultilocus
		if CD == '1':
			Tsc_34 = [ uniform(low = T_bound[0], high = time_space*Tsplit_34[i], size = None) for i in range(nMultilocus) ]
		else:
			Tsc_34 = [0] * nMultilocus
		if AC == '1':
			Tsc_13 = [ uniform(low = T_bound[0], high = time_space*min([Tsplit_12[i] , Tsplit_34[i]]), size = None) for i in range(nMultilocus) ]
		else:
			Tsc_13 = [0] * nMultilocus
		if AD == '1':
			Tsc_14 = [ uniform(low = T_bound[0], high = time_space*min([Tsplit_12[i] , Tsplit_34[i]]), size = None) for i in range(nMultilocus) ]
		else:
			Tsc_14 = [0] * nMultilocus
		if BC == '1':
			Tsc_23 = [ uniform(low = T_bound[0], high = time_space*min([Tsplit_12[i] , Tsplit_34[i]]), size = None) for i in range(nMultilocus) ]
		else:
			Tsc_23 = [0] * nMultilocus
		if BD == '1':
			Tsc_24 = [ uniform(low = T_bound[0], high = time_space*min([Tsplit_12[i] , Tsplit_34[i]]), size = None) for i in range(nMultilocus) ]
		else:
			Tsc_24 = [0] * nMultilocus

if model == "topo1":
    #could be "topo1_SC_2M_2N" to define different models: msnsam tbs 10000 -t tbs -r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs  -m 1 3 tbs -m 3 1 tbs -m 2 4 tbs -m 4 2 tbs -m 1 2 tbs -m 2 1 tbs -m 2 3 tbs -m 3 2 tbs -m 1 4 tbs -m 4 1 tbs -m 3 4 tbs -m 4 3 tbs  -em tbs 1 3 0 -em tbs 3 1 0 -em tbs 2 4 0 -em tbs 4 2 0 -ej tbs 2 1 -en tbs 1 tbs -ej tbs 4 3 -en tbs 3 tbs -ej tbs 3 1 -eN tbs tbs\n" # no migration
    # param monolocus: values that will be read by ms
    priorfile = "N1\tN2\tN3\tN4\tNa_23\tNa_24\tNa\tshape_N_a\tshape_N_b\tTsplit_23\tTsplit_24\tTsplit\tTsc_23\tTsc_24\tTsc_14\tTsc_12\tTsc_13\tTsc_34\tM13\tshape_M13_a\tshape_M13_b\tM31\tshape_M31_a\tshape_M31_b\tM24\tshape_M24_a\tshape_M24_b\tM42\tshape_M42_a\tshape_M42_b\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\tM23\tshape_M23_a\tshape_M23_b\tM32\tshape_M32_a\tshape_M32_b\tM14\tshape_M14_a\tshape_M14_b\tM41\tshape_M41_a\tshape_M41_b\tM34\tshape_M34_a\tshape_M34_b\tM43\tshape_M43_a\tshape_M43_b\n"
    for sim in range(nMultilocus):
        priorfile += "{N1:.5f}\t{N2:.5f}\t{N3:.5f}\t{N4:.5f}\t{Na_23:.5f}\t{Na_24:.5f}\t{Na:.5f}\t{shape_N_a:.5f}\t{shape_N_b:.5f}\t{Tsplit_23:.5f}\t{Tsplit_24:.5f}\t{Tsplit:.5f}\t{Tsc_23:.5f}\t{Tsc_24:.5f}\t{Tsc_14:.5f}\t{Tsc_12:.5f}\t{Tsc_13:.5f}\t{Tsc_34:.5f}\t{M13:.5f}\t{shape_M13_a:.5f}\t{shape_M13_b:.5f}\t{M31:.5f}\t{shape_M31_a:.5f}\t{shape_M31_b:.5f}\t{M24:.5f}\t{shape_M24_a:.5f}\t{shape_M24_b:.5f}\t{M42:.5f}\t{shape_M42_a:.5f}\t{shape_M42_b:.5f}\t{M12:.5f}\t{shape_M12_a:.5f}\t{shape_M12_b:.5f}\t{M21:.5f}\t{shape_M21_a:.5f}\t{shape_M21_b:.5f}\t{M23:.5f}\t{shape_M23_a:.5f}\t{shape_M23_b:.5f}\t{M32:.5f}\t{shape_M32_a:.5f}\t{shape_M32_b:.5f}\t{M14:.5f}\t{shape_M14_a:.5f}\t{shape_M14_b:.5f}\t{M41:.5f}\t{shape_M41_a:.5f}\t{shape_M41_b:.5f}\t{M34:.5f}\t{shape_M34_a:.5f}\t{shape_M34_b:.5f}\t{M43:.5f}\t{shape_M43_a:.5f}\t{shape_M43_b}\n".format(N1=N1[sim], N2=N2[sim], N3=N3[sim], N4=N4[sim], Na_23=Na_23[sim], Na_24=Na_24[sim], Na=Na[sim], shape_N_a=shape_N_a[sim], shape_N_b=shape_N_b[sim], Tsplit_23=Tsplit_23[sim], Tsplit_24=Tsplit_24[sim], Tsplit=Tsplit[sim], Tsc_23=Tsc_23[sim], Tsc_24=Tsc_24[sim], Tsc_14=Tsc_14[sim], Tsc_12=Tsc_12[sim], Tsc_13=Tsc_13[sim], Tsc_34=Tsc_34[sim], M13=M13[sim], shape_M13_a=shape_M13_a[sim], shape_M13_b=shape_M13_b[sim], M31=M31[sim], shape_M31_a=shape_M31_a[sim], shape_M31_b=shape_M31_b[sim], M24=M24[sim], shape_M24_a=shape_M24_a[sim], shape_M24_b=shape_M24_b[sim], M42=M42[sim], shape_M42_a=shape_M42_a[sim], shape_M42_b=shape_M42_b[sim], M12=M12[sim], shape_M12_a=shape_M12_a[sim], shape_M12_b=shape_M12_b[sim], M21=M21[sim], shape_M21_a=shape_M21_a[sim], shape_M21_b=shape_M21_b[sim], M23=M23[sim], shape_M23_a=shape_M23_a[sim], shape_M23_b=shape_M23_b[sim], M32=M32[sim], shape_M32_a=shape_M32_a[sim], shape_M32_b=shape_M32_b[sim], M14=M14[sim], shape_M14_a=shape_M14_a[sim], shape_M14_b=shape_M14_b[sim], M41=M41[sim], shape_M41_a=shape_M41_a[sim], shape_M41_b=shape_M41_b[sim], M34=M34[sim], shape_M34_a=shape_M34_a[sim], shape_M34_b=shape_M34_b[sim], M43=M43[sim], shape_M43_a=shape_M43_a[sim], shape_M43_b=shape_M43_b[sim])
        
        # vectors of size 'nLoci' containing parameters
# vectors of size 'nLoci' containing parameters
        scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
        median_beta =  shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim])
        N1_vec = [ N1[sim]*i/median_beta for i in scalar_N ]
        N2_vec = [ N2[sim]*i/median_beta for i in scalar_N ]
        N3_vec = [ N3[sim]*i/median_beta for i in scalar_N ]
        N4_vec = [ N4[sim]*i/median_beta for i in scalar_N ]
        Na_vec = [ Na[sim]*i/median_beta for i in scalar_N ]
        Na_23_vec = [ Na_23[sim]*i/median_beta for i in scalar_N ]
        Na_24_vec = [ Na_24[sim]*i/median_beta for i in scalar_N ]
    
        if AB == '1':
                scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
                median_beta =  shape_M12_a[sim] / (shape_M12_a[sim] + shape_M12_b[sim])
                M12_vec = [ M12[sim] * i/median_beta for i in scalar_M12 ]
                
                scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
                median_beta =  shape_M21_a[sim] / (shape_M21_a[sim] + shape_M21_b[sim])
                M21_vec = [ M21[sim] * i/median_beta for i in scalar_M21 ]
        else:
                M12_vec = [0]*nLoci
                M21_vec = [0]*nLoci

        if AC == '1':
                scalar_M13 = beta(shape_M13_a[sim], shape_M13_b[sim], size = nLoci)
                median_beta =  shape_M13_a[sim] / (shape_M13_a[sim] + shape_M13_b[sim])
                M13_vec = [ M13[sim] * i/median_beta for i in scalar_M13 ]
		
                scalar_M31 = beta(shape_M31_a[sim], shape_M31_b[sim], size = nLoci)
                median_beta =  shape_M31_a[sim] / (shape_M31_a[sim] + shape_M31_b[sim])
                M31_vec = [ M31[sim] * i/median_beta for i in scalar_M31 ]
        else:
                M13_vec = [0]*nLoci
                M31_vec = [0]*nLoci

        if AD == '1':
                scalar_M14 = beta(shape_M14_a[sim], shape_M14_b[sim], size = nLoci)
                median_beta =  shape_M14_a[sim] / (shape_M14_a[sim] + shape_M14_b[sim])
                M14_vec = [ M14[sim] * i/median_beta for i in scalar_M14 ]
		
                scalar_M41 = beta(shape_M41_a[sim], shape_M41_b[sim], size = nLoci)
                median_beta =  shape_M41_a[sim] / (shape_M41_a[sim] + shape_M41_b[sim])
                M41_vec = [ M41[sim] * i/median_beta for i in scalar_M41 ]
        else:
                M14_vec = [0]*nLoci
                M41_vec = [0]*nLoci

        if BC == '1':
                scalar_M23 = beta(shape_M23_a[sim], shape_M23_b[sim], size = nLoci)
                median_beta =  shape_M23_a[sim] / (shape_M23_a[sim] + shape_M23_b[sim])
                M23_vec = [ M23[sim] * i/median_beta for i in scalar_M23 ]
		
                scalar_M32 = beta(shape_M32_a[sim], shape_M32_b[sim], size = nLoci)
                median_beta =  shape_M32_a[sim] / (shape_M32_a[sim] + shape_M32_b[sim])
                M32_vec = [ M32[sim] * i/median_beta for i in scalar_M32 ]
        else:
                M23_vec = [0]*nLoci
                M32_vec = [0]*nLoci

        if BD == '1':
                scalar_M24 = beta(shape_M24_a[sim], shape_M24_b[sim], size = nLoci)
                median_beta =  shape_M24_a[sim] / (shape_M24_a[sim] + shape_M24_b[sim])
                M24_vec = [ M24[sim] * i/median_beta for i in scalar_M24 ]
		
                scalar_M42 = beta(shape_M42_a[sim], shape_M42_b[sim], size = nLoci)
                median_beta =  shape_M42_a[sim] / (shape_M42_a[sim] + shape_M42_b[sim])
                M42_vec = [ M42[sim] * i/median_beta for i in scalar_M42 ]
        else:
                M24_vec = [0]*nLoci
                M42_vec = [0]*nLoci

        if CD == '1':
                scalar_M34 = beta(shape_M34_a[sim], shape_M34_b[sim], size = nLoci)
                median_beta =  shape_M34_a[sim] / (shape_M34_a[sim] + shape_M34_b[sim])
                M34_vec = [ M34[sim] * i/median_beta for i in scalar_M34 ]
		
                scalar_M43 = beta(shape_M43_a[sim], shape_M43_b[sim], size = nLoci)
                median_beta =  shape_M43_a[sim] / (shape_M43_a[sim] + shape_M43_b[sim])
                M43_vec = [ M43[sim] * i/median_beta for i in scalar_M43 ]
        else:
                M34_vec = [0]*nLoci
                M43_vec = [0]*nLoci
        
        for locus in range(nLoci):
            if mutation=='mu': # if mutations are conditionned on mu
                    print("{nsam}\t{theta:.5f}\t{rho:.5f}\t{L}\t{nsamA}\t{nsamB}\t{nsamC}\t{nsamD}\t{N1:.5f}\t{N2:.5f}\t{N3:.5f}\t{N4:.5f}\t{M13:.5f}\t{M31:.5f}\t{M24:.5f}\t{M42:.5f}\t{M12:.5f}\t{M21:.5f}\t{M23:.5f}\t{M32:.5f}\t{M14:.5f}\t{M41:.5f}\t{M34:.5f}\t{M43:.5f}\t{Tsc_23:.5f}\t{Tsc_23:.5f}\t{Tsc_13:.5f}\t{Tsc_13:.5f}\t{Tsc_34:.5f}\t{Tsc_34:.5f}\t{Tsc_24:.5f}\t{Tsc_24:.5f}\t{Tsc_12:.5f}\t{Tsc_12:.5f}\t{Tsc_14:.5f}\t{Tsc_14:.5f}\t{Tsplit_23:.5f}\t{Tsplit_23:.5f}\t{Na_23:.5f}\t{Tsplit_24:.5f}\t{Tsplit_24:.5f}\t{Na_24:.5f}\t{Tsplit:.5f}\t{Tsplit:.5f}\t{Na:.5f}".format(nsam=nsam_tot[locus], theta=theta[locus], rho=rho[locus], L=L[locus], nsamA=nsamA[locus], nsamB=nsamB[locus], nsamC=nsamC[locus], nsamD=nsamD[locus], N1=N1_vec[locus], N2=N2_vec[locus], N3=N3_vec[locus], N4=N4_vec[locus], M13=M13_vec[locus], M31=M31_vec[locus], M24=M24_vec[locus], M42=M42_vec[locus], M12=M12_vec[locus], M21=M21_vec[locus], M23=M23_vec[locus], M32=M32_vec[locus], M14=M14_vec[locus], M41=M41_vec[locus], M34=M34_vec[locus], M43=M43_vec[locus], Tsc_23=Tsc_23[sim], Tsc_13=Tsc_13[sim], Tsc_34=Tsc_34[sim], Tsc_24=Tsc_24[sim], Tsc_12=Tsc_12[sim], Tsc_14=Tsc_14[sim], Tsplit_23=Tsplit_23[sim], Na_23=Na_23_vec[locus], Tsplit_24=Tsplit_24[sim], Na_24=Na_24_vec[locus], Tsplit=Tsplit[sim], Na=Na_vec[locus]))
            else: # if mutation are conditionned on the number of SNPs
                    print("{nsam}\t{nSNPs}\t{rho:.5f}\t{L}\t{nsamA}\t{nsamB}\t{nsamC}\t{nsamD}\t{N1:.5f}\t{N2:.5f}\t{N3:.5f}\t{N4:.5f}\t{M13:.5f}\t{M31:.5f}\t{M24:.5f}\t{M42:.5f}\t{M12:.5f}\t{M21:.5f}\t{M23:.5f}\t{M32:.5f}\t{M14:.5f}\t{M41:.5f}\t{M34:.5f}\t{M43:.5f}\t{Tsc_23:.5f}\t{Tsc_23:.5f}\t{Tsc_13:.5f}\t{Tsc_13:.5f}\t{Tsc_34:.5f}\t{Tsc_34:.5f}\t{Tsc_24:.5f}\t{Tsc_24:.5f}\t{Tsc_12:.5f}\t{Tsc_12:.5f}\t{Tsc_14:.5f}\t{Tsc_14:.5f}\t{Tsplit_23:.5f}\t{Tsplit_23:.5f}\t{Na_23:.5f}\t{Tsplit_24:.5f}\t{Tsplit_24:.5f}\t{Na_24:.5f}\t{Tsplit:.5f}\t{Tsplit:.5f}\t{Na:.5f}".format(nsam=nsam_tot[locus], nSNPs=int(nSNPs[locus]), rho=rho[locus], L=L[locus], nsamA=nsamA[locus], nsamB=nsamB[locus], nsamC=nsamC[locus], nsamD=nsamD[locus], N1=N1_vec[locus], N2=N2_vec[locus], N3=N3_vec[locus], N4=N4_vec[locus], M13=M13_vec[locus], M31=M31_vec[locus], M24=M24_vec[locus], M42=M42_vec[locus], M12=M12_vec[locus], M21=M21_vec[locus], M23=M23_vec[locus], M32=M32_vec[locus], M14=M14_vec[locus], M41=M41_vec[locus], M34=M34_vec[locus], M43=M43_vec[locus], Tsc_23=Tsc_23[sim], Tsc_13=Tsc_13[sim], Tsc_34=Tsc_34[sim], Tsc_24=Tsc_24[sim], Tsc_12=Tsc_12[sim], Tsc_14=Tsc_14[sim], Tsplit_23=Tsplit_23[sim], Na_23=Na_23_vec[locus], Tsplit_24=Tsplit_24[sim], Na_24=Na_24_vec[locus], Tsplit=Tsplit[sim], Na=Na_vec[locus]))
    outfile = open("{simulationpath}/priorfile.txt".format(simulationpath=simulationpath), "w")
    outfile.write(priorfile)
    outfile.close()


if model == "topo2":
    #could be "topo1_SC_2M_2N" to define different models: msnsam tbs 10000 -t tbs -r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs  -m 1 3 tbs -m 3 1 tbs -m 2 4 tbs -m 4 2 tbs -m 1 2 tbs -m 2 1 tbs -m 2 3 tbs -m 3 2 tbs -m 1 4 tbs -m 4 1 tbs -m 3 4 tbs -m 4 3 tbs  -em tbs 1 3 0 -em tbs 3 1 0 -em tbs 2 4 0 -em tbs 4 2 0 -ej tbs 2 1 -en tbs 1 tbs -ej tbs 4 3 -en tbs 3 tbs -ej tbs 3 1 -eN tbs tbs\n" # no migration
    # param monolocus: values that will be read by ms
    priorfile = "N1\tN2\tN3\tN4\tNa_34\tNa_24\tNa\tshape_N_a\tshape_N_b\tTsplit_34\tTsplit_24\tTsplit\tTsc_23\tTsc_24\tTsc_14\tTsc_12\tTsc_13\tTsc_34\tM13\tshape_M13_a\tshape_M13_b\tM31\tshape_M31_a\tshape_M31_b\tM24\tshape_M24_a\tshape_M24_b\tM42\tshape_M42_a\tshape_M42_b\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\tM23\tshape_M23_a\tshape_M23_b\tM32\tshape_M32_a\tshape_M32_b\tM14\tshape_M14_a\tshape_M14_b\tM41\tshape_M41_a\tshape_M41_b\tM34\tshape_M34_a\tshape_M34_b\tM43\tshape_M43_a\tshape_M43_b\n"
    for sim in range(nMultilocus):
        priorfile += "{N1:.5f}\t{N2:.5f}\t{N3:.5f}\t{N4:.5f}\t{Na_34:.5f}\t{Na_24:.5f}\t{Na:.5f}\t{shape_N_a:.5f}\t{shape_N_b:.5f}\t{Tsplit_34:.5f}\t{Tsplit_24:.5f}\t{Tsplit:.5f}\t{Tsc_23:.5f}\t{Tsc_24:.5f}\t{Tsc_14:.5f}\t{Tsc_12:.5f}\t{Tsc_13:.5f}\t{Tsc_34:.5f}\t{M13:.5f}\t{shape_M13_a:.5f}\t{shape_M13_b:.5f}\t{M31:.5f}\t{shape_M31_a:.5f}\t{shape_M31_b:.5f}\t{M24:.5f}\t{shape_M24_a:.5f}\t{shape_M24_b:.5f}\t{M42:.5f}\t{shape_M42_a:.5f}\t{shape_M42_b:.5f}\t{M12:.5f}\t{shape_M12_a:.5f}\t{shape_M12_b:.5f}\t{M21:.5f}\t{shape_M21_a:.5f}\t{shape_M21_b:.5f}\t{M23:.5f}\t{shape_M23_a:.5f}\t{shape_M23_b:.5f}\t{M32:.5f}\t{shape_M32_a:.5f}\t{shape_M32_b:.5f}\t{M14:.5f}\t{shape_M14_a:.5f}\t{shape_M14_b:.5f}\t{M41:.5f}\t{shape_M41_a:.5f}\t{shape_M41_b:.5f}\t{M34:.5f}\t{shape_M34_a:.5f}\t{shape_M34_b:.5f}\t{M43:.5f}\t{shape_M43_a:.5f}\t{shape_M43_b}\n".format(N1=N1[sim], N2=N2[sim], N3=N3[sim], N4=N4[sim], Na_34=Na_34[sim], Na_24=Na_24[sim], Na=Na[sim], shape_N_a=shape_N_a[sim], shape_N_b=shape_N_b[sim], Tsplit_34=Tsplit_34[sim], Tsplit_24=Tsplit_24[sim], Tsplit=Tsplit[sim], Tsc_23=Tsc_23[sim], Tsc_24=Tsc_24[sim], Tsc_14=Tsc_14[sim], Tsc_12=Tsc_12[sim], Tsc_13=Tsc_13[sim], Tsc_34=Tsc_34[sim], M13=M13[sim], shape_M13_a=shape_M13_a[sim], shape_M13_b=shape_M13_b[sim], M31=M31[sim], shape_M31_a=shape_M31_a[sim], shape_M31_b=shape_M31_b[sim], M24=M24[sim], shape_M24_a=shape_M24_a[sim], shape_M24_b=shape_M24_b[sim], M42=M42[sim], shape_M42_a=shape_M42_a[sim], shape_M42_b=shape_M42_b[sim], M12=M12[sim], shape_M12_a=shape_M12_a[sim], shape_M12_b=shape_M12_b[sim], M21=M21[sim], shape_M21_a=shape_M21_a[sim], shape_M21_b=shape_M21_b[sim], M23=M23[sim], shape_M23_a=shape_M23_a[sim], shape_M23_b=shape_M23_b[sim], M32=M32[sim], shape_M32_a=shape_M32_a[sim], shape_M32_b=shape_M32_b[sim], M14=M14[sim], shape_M14_a=shape_M14_a[sim], shape_M14_b=shape_M14_b[sim], M41=M41[sim], shape_M41_a=shape_M41_a[sim], shape_M41_b=shape_M41_b[sim], M34=M34[sim], shape_M34_a=shape_M34_a[sim], shape_M34_b=shape_M34_b[sim], M43=M43[sim], shape_M43_a=shape_M43_a[sim], shape_M43_b=shape_M43_b[sim])
            
        # vectors of size 'nLoci' containing parameters
        scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
        median_beta =  shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim])
        N1_vec = [ N1[sim]*i/median_beta for i in scalar_N ]
        N2_vec = [ N2[sim]*i/median_beta for i in scalar_N ]
        N3_vec = [ N3[sim]*i/median_beta for i in scalar_N ]
        N4_vec = [ N4[sim]*i/median_beta for i in scalar_N ]
        Na_vec = [ Na[sim]*i/median_beta for i in scalar_N ]
        Na_34_vec = [ Na_34[sim]*i/median_beta for i in scalar_N ]
        Na_24_vec = [ Na_24[sim]*i/median_beta for i in scalar_N ]
	
        if AB == '1':
                scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
                median_beta =  shape_M12_a[sim] / (shape_M12_a[sim] + shape_M12_b[sim])
                M12_vec = [ M12[sim] * i/median_beta for i in scalar_M12 ]
                
                scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
                median_beta =  shape_M21_a[sim] / (shape_M21_a[sim] + shape_M21_b[sim])
                M21_vec = [ M21[sim] * i/median_beta for i in scalar_M21 ]
        else:
                M12_vec = [0]*nLoci
                M21_vec = [0]*nLoci

        if AC == '1':
                scalar_M13 = beta(shape_M13_a[sim], shape_M13_b[sim], size = nLoci)
                median_beta =  shape_M13_a[sim] / (shape_M13_a[sim] + shape_M13_b[sim])
                M13_vec = [ M13[sim] * i/median_beta for i in scalar_M13 ]
		
                scalar_M31 = beta(shape_M31_a[sim], shape_M31_b[sim], size = nLoci)
                median_beta =  shape_M31_a[sim] / (shape_M31_a[sim] + shape_M31_b[sim])
                M31_vec = [ M31[sim] * i/median_beta for i in scalar_M31 ]
        else:
                M13_vec = [0]*nLoci
                M31_vec = [0]*nLoci

        if AD == '1':
                scalar_M14 = beta(shape_M14_a[sim], shape_M14_b[sim], size = nLoci)
                median_beta =  shape_M14_a[sim] / (shape_M14_a[sim] + shape_M14_b[sim])
                M14_vec = [ M14[sim] * i/median_beta for i in scalar_M14 ]
		
                scalar_M41 = beta(shape_M41_a[sim], shape_M41_b[sim], size = nLoci)
                median_beta =  shape_M41_a[sim] / (shape_M41_a[sim] + shape_M41_b[sim])
                M41_vec = [ M41[sim] * i/median_beta for i in scalar_M41 ]
        else:
                M14_vec = [0]*nLoci
                M41_vec = [0]*nLoci

        if BC == '1':
                scalar_M23 = beta(shape_M23_a[sim], shape_M23_b[sim], size = nLoci)
                median_beta =  shape_M23_a[sim] / (shape_M23_a[sim] + shape_M23_b[sim])
                M23_vec = [ M23[sim] * i/median_beta for i in scalar_M23 ]
		
                scalar_M32 = beta(shape_M32_a[sim], shape_M32_b[sim], size = nLoci)
                median_beta =  shape_M32_a[sim] / (shape_M32_a[sim] + shape_M32_b[sim])
                M32_vec = [ M32[sim] * i/median_beta for i in scalar_M32 ]
        else:
                M23_vec = [0]*nLoci
                M32_vec = [0]*nLoci

        if BD == '1':
                scalar_M24 = beta(shape_M24_a[sim], shape_M24_b[sim], size = nLoci)
                median_beta =  shape_M24_a[sim] / (shape_M24_a[sim] + shape_M24_b[sim])
                M24_vec = [ M24[sim] * i/median_beta for i in scalar_M24 ]
		
                scalar_M42 = beta(shape_M42_a[sim], shape_M42_b[sim], size = nLoci)
                median_beta =  shape_M42_a[sim] / (shape_M42_a[sim] + shape_M42_b[sim])
                M42_vec = [ M42[sim] * i/median_beta for i in scalar_M42 ]
        else:
                M24_vec = [0]*nLoci
                M42_vec = [0]*nLoci

        if CD == '1':
                scalar_M34 = beta(shape_M34_a[sim], shape_M34_b[sim], size = nLoci)
                median_beta =  shape_M34_a[sim] / (shape_M34_a[sim] + shape_M34_b[sim])
                M34_vec = [ M34[sim] * i/median_beta for i in scalar_M34 ]
		
                scalar_M43 = beta(shape_M43_a[sim], shape_M43_b[sim], size = nLoci)
                median_beta =  shape_M43_a[sim] / (shape_M43_a[sim] + shape_M43_b[sim])
                M43_vec = [ M43[sim] * i/median_beta for i in scalar_M43 ]
        else:
                M34_vec = [0]*nLoci
                M43_vec = [0]*nLoci

        for locus in range(nLoci):
            # msnsam tbs nSIMULATIONS -t tbs -r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 3 tbs -m 3 1 tbs -m 2 4 tbs -m 4 2 tbs -m 1 2 tbs -m 2 1 tbs -m 2 3 tbs -m 3 2 tbs -m 1 4 tbs -m 4 1 tbs -m 3 4 tbs -m 4 3 tbs -em tbs 3 4 0 -em tbs 4 3 0 -em tbs 3 2 0 -em tbs 2 3 0 -em tbs 3 1 0 -em tbs 1 3 0 -em tbs 2 4 0 -em tbs 4 2 0 -em tbs 2 1 0 -em tbs 1 2 0 -em tbs 4 1 0 -em tbs 1 4 0 -ej tbs 3 4 -en tbs 4 tbs -ej tbs 4 2 -en tbs 2 tbs -ej tbs 2 1 -eN tbs tbs
            if mutation=='mu': # if mutations are conditionned on mu
                  print("{nsam}\t{theta:.5f}\t{rho:.5f}\t{L}\t{nsamA}\t{nsamB}\t{nsamC}\t{nsamD}\t{N1:.5f}\t{N2:.5f}\t{N3:.5f}\t{N4:.5f}\t{M13:.5f}\t{M31:.5f}\t{M24:.5f}\t{M42:.5f}\t{M12:.5f}\t{M21:.5f}\t{M23:.5f}\t{M32:.5f}\t{M14:.5f}\t{M41:.5f}\t{M34:.5f}\t{M43:.5f}\t{Tsc_34:.5f}\t{Tsc_34:.5f}\t{Tsc_23:.5f}\t{Tsc_23:.5f}\t{Tsc_13:.5f}\t{Tsc_13:.5f}\t{Tsc_24:.5f}\t{Tsc_24:.5f}\t{Tsc_12:.5f}\t{Tsc_12:.5f}\t{Tsc_14:.5f}\t{Tsc_14:.5f}\t{Tsplit_34:.5f}\t{Tsplit_34:.5f}\t{Na_34:.5f}\t{Tsplit_24:.5f}\t{Tsplit_24:.5f}\t{Na_24:.5f}\t{Tsplit:.5f}\t{Tsplit:.5f}\t{Na:.5f}".format(nsam=nsam_tot[locus], theta=theta[locus], rho=rho[locus], L=L[locus], nsamA=nsamA[locus], nsamB=nsamB[locus], nsamC=nsamC[locus], nsamD=nsamD[locus], N1=N1_vec[locus], N2=N2_vec[locus], N3=N3_vec[locus], N4=N4_vec[locus], M13=M13_vec[locus], M31=M31_vec[locus], M24=M24_vec[locus], M42=M42_vec[locus], M12=M12_vec[locus], M21=M21_vec[locus], M23=M23_vec[locus], M32=M32_vec[locus], M14=M14_vec[locus], M41=M41_vec[locus], M34=M34_vec[locus], M43=M43_vec[locus], Tsc_34=Tsc_34[sim], Tsc_23=Tsc_23[sim], Tsc_13=Tsc_13[sim], Tsc_24=Tsc_24[sim], Tsc_12=Tsc_12[sim], Tsc_14=Tsc_14[sim], Tsplit_34=Tsplit_34[sim], Na_34=Na_34_vec[locus], Tsplit_24=Tsplit_24[sim], Na_24=Na_24_vec[locus], Tsplit=Tsplit[sim], Na=Na_vec[locus]))
            else: 
                  print("{nsam}\t{nSNPs}\t{rho:.5f}\t{L}\t{nsamA}\t{nsamB}\t{nsamC}\t{nsamD}\t{N1:.5f}\t{N2:.5f}\t{N3:.5f}\t{N4:.5f}\t{M13:.5f}\t{M31:.5f}\t{M24:.5f}\t{M42:.5f}\t{M12:.5f}\t{M21:.5f}\t{M23:.5f}\t{M32:.5f}\t{M14:.5f}\t{M41:.5f}\t{M34:.5f}\t{M43:.5f}\t{Tsc_34:.5f}\t{Tsc_34:.5f}\t{Tsc_23:.5f}\t{Tsc_23:.5f}\t{Tsc_13:.5f}\t{Tsc_13:.5f}\t{Tsc_24:.5f}\t{Tsc_24:.5f}\t{Tsc_12:.5f}\t{Tsc_12:.5f}\t{Tsc_14:.5f}\t{Tsc_14:.5f}\t{Tsplit_34:.5f}\t{Tsplit_34:.5f}\t{Na_34:.5f}\t{Tsplit_24:.5f}\t{Tsplit_24:.5f}\t{Na_24:.5f}\t{Tsplit:.5f}\t{Tsplit:.5f}\t{Na:.5f}".format(nsam=nsam_tot[locus], nSNPs=int(nSNPs[locus]), rho=rho[locus], L=L[locus], nsamA=nsamA[locus], nsamB=nsamB[locus], nsamC=nsamC[locus], nsamD=nsamD[locus], N1=N1_vec[locus], N2=N2_vec[locus], N3=N3_vec[locus], N4=N4_vec[locus], M13=M13_vec[locus], M31=M31_vec[locus], M24=M24_vec[locus], M42=M42_vec[locus], M12=M12_vec[locus], M21=M21_vec[locus], M23=M23_vec[locus], M32=M32_vec[locus], M14=M14_vec[locus], M41=M41_vec[locus], M34=M34_vec[locus], M43=M43_vec[locus], Tsc_34=Tsc_34[sim], Tsc_23=Tsc_23[sim], Tsc_13=Tsc_13[sim], Tsc_24=Tsc_24[sim], Tsc_12=Tsc_12[sim], Tsc_14=Tsc_14[sim], Tsplit_34=Tsplit_34[sim], Na_34=Na_34_vec[locus], Tsplit_24=Tsplit_24[sim], Na_24=Na_24_vec[locus], Tsplit=Tsplit[sim], Na=Na_vec[locus]))
    outfile = open("{simulationpath}/priorfile.txt".format(simulationpath=simulationpath), "w")
    outfile.write(priorfile)
    outfile.close()


if model == "topo3":
    # msnsam tbs nSIMULATIONS -t tbs -r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 1 4 tbs -m 4 1 tbs -m 2 3 tbs -m 3 2 tbs -m 2 4 tbs -m 4 2 tbs -m 3 4 tbs -m 4 3 tbs -em tbs 1 2 0 -em tbs 2 1 0 -em tbs 1 3 0 -em tbs 3 1 0 -em tbs 1 4 0 -em tbs 4 1 0 -em tbs 2 3 0 -em tbs 3 2 0 -em tbs 2 4 0 -em tbs 4 2 0 -em tbs 3 4 0 -em tbs 4 3 0 -ej tbs 2 1 -en tbs 1 tbs -ej tbs 4 3 -en tbs 3 tbs -ej tbs 3 1 -eN tbs tbs
    # param monolocus: values that will be read by ms
    priorfile = "N1\tN2\tN3\tN4\tNa_12\tNa_34\tNa\tshape_N_a\tshape_N_b\tTsplit_12\tTsplit_34\tTsplit\tTsc_12\tTsc_13\tTsc_14\tTsc_23\tTsc_24\tTsc_34\tM13\tshape_M13_a\tshape_M13_b\tM31\tshape_M31_a\tshape_M31_b\tM24\tshape_M24_a\tshape_M24_b\tM42\tshape_M42_a\tshape_M42_b\tM12\tshape_M12_a\tshape_M12_b\tM21\tshape_M21_a\tshape_M21_b\tM23\tshape_M23_a\tshape_M23_b\tM32\tshape_M32_a\tshape_M32_b\tM14\tshape_M14_a\tshape_M14_b\tM41\tshape_M41_a\tshape_M41_b\tM34\tshape_M34_a\tshape_M34_b\tM43\tshape_M43_a\tshape_M43_b\n"
    for sim in range(nMultilocus):
        priorfile += "{N1:.5f}\t{N2:.5f}\t{N3:.5f}\t{N4:.5f}\t{Na_12:.5f}\t{Na_34:.5f}\t{Na:.5f}\t{shape_N_a:.5f}\t{shape_N_b:.5f}\t{Tsplit_12:.5f}\t{Tsplit_34:.5f}\t{Tsplit:.5f}\t{Tsc_12:.5f}\t{Tsc_13:.5f}\t{Tsc_14:.5f}\t{Tsc_23:.5f}\t{Tsc_24:.5f}\t{Tsc_34:.5f}\t{M13:.5f}\t{shape_M13_a:.5f}\t{shape_M13_b:.5f}\t{M31:.5f}\t{shape_M31_a:.5f}\t{shape_M31_b:.5f}\t{M24:.5f}\t{shape_M24_a:.5f}\t{shape_M24_b:.5f}\t{M42:.5f}\t{shape_M42_a:.5f}\t{shape_M42_b:.5f}\t{M12:.5f}\t{shape_M12_a:.5f}\t{shape_M12_b:.5f}\t{M21:.5f}\t{shape_M21_a:.5f}\t{shape_M21_b:.5f}\t{M23:.5f}\t{shape_M23_a:.5f}\t{shape_M23_b:.5f}\t{M32:.5f}\t{shape_M32_a:.5f}\t{shape_M32_b:.5f}\t{M14:.5f}\t{shape_M14_a:.5f}\t{shape_M14_b:.5f}\t{M41:.5f}\t{shape_M41_a:.5f}\t{shape_M41_b:.5f}\t{M34:.5f}\t{shape_M34_a:.5f}\t{shape_M34_b:.5f}\t{M43:.5f}\t{shape_M43_a:.5f}\t{shape_M43_b}\n".format(N1=N1[sim], N2=N2[sim], N3=N3[sim], N4=N4[sim], Na_12=Na_12[sim], Na_34=Na_34[sim], Na=Na[sim], shape_N_a=shape_N_a[sim], shape_N_b=shape_N_b[sim], Tsplit_12=Tsplit_12[sim], Tsplit_34=Tsplit_34[sim], Tsplit=Tsplit[sim], Tsc_12=Tsc_12[sim], Tsc_13=Tsc_13[sim], Tsc_14=Tsc_14[sim], Tsc_23=Tsc_23[sim], Tsc_24=Tsc_24[sim], Tsc_34=Tsc_34[sim], M13=M13[sim], shape_M13_a=shape_M13_a[sim], shape_M13_b=shape_M13_b[sim], M31=M31[sim], shape_M31_a=shape_M31_a[sim], shape_M31_b=shape_M31_b[sim], M24=M24[sim], shape_M24_a=shape_M24_a[sim], shape_M24_b=shape_M24_b[sim], M42=M42[sim], shape_M42_a=shape_M42_a[sim], shape_M42_b=shape_M42_b[sim], M12=M12[sim], shape_M12_a=shape_M12_a[sim], shape_M12_b=shape_M12_b[sim], M21=M21[sim], shape_M21_a=shape_M21_a[sim], shape_M21_b=shape_M21_b[sim], M23=M23[sim], shape_M23_a=shape_M23_a[sim], shape_M23_b=shape_M23_b[sim], M32=M32[sim], shape_M32_a=shape_M32_a[sim], shape_M32_b=shape_M32_b[sim], M14=M14[sim], shape_M14_a=shape_M14_a[sim], shape_M14_b=shape_M14_b[sim], M41=M41[sim], shape_M41_a=shape_M41_a[sim], shape_M41_b=shape_M41_b[sim], M34=M34[sim], shape_M34_a=shape_M34_a[sim], shape_M34_b=shape_M34_b[sim], M43=M43[sim], shape_M43_a=shape_M43_a[sim], shape_M43_b=shape_M43_b[sim])
        
        # vectors of size 'nLoci' containing parameters
        scalar_N = beta(shape_N_a[sim], shape_N_b[sim], size=nLoci)
        median_beta =  shape_N_a[sim] / (shape_N_a[sim] + shape_N_b[sim])
        N1_vec = [ N1[sim]*i/median_beta for i in scalar_N ]
        N2_vec = [ N2[sim]*i/median_beta for i in scalar_N ]
        N3_vec = [ N3[sim]*i/median_beta for i in scalar_N ]
        N4_vec = [ N4[sim]*i/median_beta for i in scalar_N ]
        Na_vec = [ Na[sim]*i/median_beta for i in scalar_N ]
        Na_12_vec = [ Na_12[sim]*i/median_beta for i in scalar_N ]
        Na_34_vec = [ Na_34[sim]*i/median_beta for i in scalar_N ]
        
        if AB == '1':
                scalar_M12 = beta(shape_M12_a[sim], shape_M12_b[sim], size = nLoci)
                median_beta =  shape_M12_a[sim] / (shape_M12_a[sim] + shape_M12_b[sim])
                M12_vec = [ M12[sim] * i/median_beta for i in scalar_M12 ]
                
                scalar_M21 = beta(shape_M21_a[sim], shape_M21_b[sim], size = nLoci)
                median_beta =  shape_M21_a[sim] / (shape_M21_a[sim] + shape_M21_b[sim])
                M21_vec = [ M21[sim] * i/median_beta for i in scalar_M21 ]
        else:
                M12_vec = [0]*nLoci
                M21_vec = [0]*nLoci

        if AC == '1':
                scalar_M13 = beta(shape_M13_a[sim], shape_M13_b[sim], size = nLoci)
                median_beta =  shape_M13_a[sim] / (shape_M13_a[sim] + shape_M13_b[sim])
                M13_vec = [ M13[sim] * i/median_beta for i in scalar_M13 ]
		
                scalar_M31 = beta(shape_M31_a[sim], shape_M31_b[sim], size = nLoci)
                median_beta =  shape_M31_a[sim] / (shape_M31_a[sim] + shape_M31_b[sim])
                M31_vec = [ M31[sim] * i/median_beta for i in scalar_M31 ]
        else:
                M13_vec = [0]*nLoci
                M31_vec = [0]*nLoci

        if AD == '1':
                scalar_M14 = beta(shape_M14_a[sim], shape_M14_b[sim], size = nLoci)
                median_beta =  shape_M14_a[sim] / (shape_M14_a[sim] + shape_M14_b[sim])
                M14_vec = [ M14[sim] * i/median_beta for i in scalar_M14 ]
		
                scalar_M41 = beta(shape_M41_a[sim], shape_M41_b[sim], size = nLoci)
                median_beta =  shape_M41_a[sim] / (shape_M41_a[sim] + shape_M41_b[sim])
                M41_vec = [ M41[sim] * i/median_beta for i in scalar_M41 ]
        else:
                M14_vec = [0]*nLoci
                M41_vec = [0]*nLoci

        if BC == '1':
                scalar_M23 = beta(shape_M23_a[sim], shape_M23_b[sim], size = nLoci)
                median_beta =  shape_M23_a[sim] / (shape_M23_a[sim] + shape_M23_b[sim])
                M23_vec = [ M23[sim] * i/median_beta for i in scalar_M23 ]
		
                scalar_M32 = beta(shape_M32_a[sim], shape_M32_b[sim], size = nLoci)
                median_beta =  shape_M32_a[sim] / (shape_M32_a[sim] + shape_M32_b[sim])
                M32_vec = [ M32[sim] * i/median_beta for i in scalar_M32 ]
        else:
                M23_vec = [0]*nLoci
                M32_vec = [0]*nLoci

        if BD == '1':
                scalar_M24 = beta(shape_M24_a[sim], shape_M24_b[sim], size = nLoci)
                median_beta =  shape_M24_a[sim] / (shape_M24_a[sim] + shape_M24_b[sim])
                M24_vec = [ M24[sim] * i/median_beta for i in scalar_M24 ]
		
                scalar_M42 = beta(shape_M42_a[sim], shape_M42_b[sim], size = nLoci)
                median_beta =  shape_M42_a[sim] / (shape_M42_a[sim] + shape_M42_b[sim])
                M42_vec = [ M42[sim] * i/median_beta for i in scalar_M42 ]
        else:
                M24_vec = [0]*nLoci
                M42_vec = [0]*nLoci

        if CD == '1':
                scalar_M34 = beta(shape_M34_a[sim], shape_M34_b[sim], size = nLoci)
                median_beta =  shape_M34_a[sim] / (shape_M34_a[sim] + shape_M34_b[sim])
                M34_vec = [ M34[sim] * i/median_beta for i in scalar_M34 ]
		
                scalar_M43 = beta(shape_M43_a[sim], shape_M43_b[sim], size = nLoci)
                median_beta =  shape_M43_a[sim] / (shape_M43_a[sim] + shape_M43_b[sim])
                M43_vec = [ M43[sim] * i/median_beta for i in scalar_M43 ]
        else:
                M34_vec = [0]*nLoci
                M43_vec = [0]*nLoci

            # msnsam tbs nSIMULATIONS -t tbs -r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 1 4 tbs -m 4 1 tbs -m 2 3 tbs -m 3 2 tbs -m 2 4 tbs -m 4 2 tbs -m 3 4 tbs -m 4 3 tbs -em tbs 1 2 0 -em tbs 2 1 0 -em tbs 1 3 0 -em tbs 3 1 0 -em tbs 1 4 0 -em tbs 4 1 0 -em tbs 2 3 0 -em tbs 3 2 0 -em tbs 2 4 0 -em tbs 4 2 0 -em tbs 3 4 0 -em tbs 4 3 0 -ej tbs 2 1 -en tbs 1 tbs -ej tbs 4 3 -en tbs 3 tbs -ej tbs 3 1 -eN tbs tbs        for locus in range(nLoci):
        for locus in range(nLoci):
            if mutation=='mu':
                  print("{nsam}\t{theta:.5f}\t{rho:.5f}\t{L}\t{nsamA}\t{nsamB}\t{nsamC}\t{nsamD}\t{N1:.5f}\t{N2:.5f}\t{N3:.5f}\t{N4:.5f}\t{M12:.5f}\t{M21:.5f}\t{M13:.5f}\t{M31:.5f}\t{M14:.5f}\t{M41:.5f}\t{M23:.5f}\t{M32:.5f}\t{M24:.5f}\t{M42:.5f}\t{M34:.5f}\t{M43:.5f}\t{Tsc_12:.5f}\t{Tsc_12:.5f}\t{Tsc_13:.5f}\t{Tsc_13:.5f}\t{Tsc_14:.5f}\t{Tsc_14:.5f}\t{Tsc_23:.5f}\t{Tsc_23:.5f}\t{Tsc_24:.5f}\t{Tsc_24:.5f}\t{Tsc_34:.5f}\t{Tsc_34:.5f}\t{Tsplit_12:.5f}\t{Tsplit_12:.5f}\t{Na_12:.5f}\t{Tsplit_34:.5f}\t{Tsplit_34:.5f}\t{Na_34:.5f}\t{Tsplit:.5f}\t{Tsplit:.5f}\t{Na:.5f}".format(nsam=nsam_tot[locus], theta=theta[locus], rho=rho[locus], L=L[locus], nsamA=nsamA[locus], nsamB=nsamB[locus], nsamC=nsamC[locus], nsamD=nsamD[locus], N1=N1_vec[locus], N2=N2_vec[locus], N3=N3_vec[locus], N4=N4_vec[locus], M12=M12_vec[locus], M21=M21_vec[locus], M13=M13_vec[locus], M31=M31_vec[locus], M14=M14_vec[locus], M41=M41_vec[locus], M23=M23_vec[locus], M32=M32_vec[locus], M24=M24_vec[locus], M42=M42_vec[locus], M34=M34_vec[locus], M43=M43_vec[locus], Tsc_12=Tsc_12[sim], Tsc_13=Tsc_13[sim], Tsc_14=Tsc_14[sim], Tsc_23=Tsc_23[sim], Tsc_24=Tsc_24[sim], Tsc_34=Tsc_34[sim], Tsplit_12=Tsplit_12[sim], Na_12=Na_12_vec[locus], Tsplit_34=Tsplit_34[sim], Na_34=Na_34_vec[locus], Tsplit=Tsplit[sim], Na=Na_vec[locus]))
            else:
                  print("{nsam}\t{nSNPs}\t{rho:.5f}\t{L}\t{nsamA}\t{nsamB}\t{nsamC}\t{nsamD}\t{N1:.5f}\t{N2:.5f}\t{N3:.5f}\t{N4:.5f}\t{M12:.5f}\t{M21:.5f}\t{M13:.5f}\t{M31:.5f}\t{M14:.5f}\t{M41:.5f}\t{M23:.5f}\t{M32:.5f}\t{M24:.5f}\t{M42:.5f}\t{M34:.5f}\t{M43:.5f}\t{Tsc_12:.5f}\t{Tsc_12:.5f}\t{Tsc_13:.5f}\t{Tsc_13:.5f}\t{Tsc_14:.5f}\t{Tsc_14:.5f}\t{Tsc_23:.5f}\t{Tsc_23:.5f}\t{Tsc_24:.5f}\t{Tsc_24:.5f}\t{Tsc_34:.5f}\t{Tsc_34:.5f}\t{Tsplit_12:.5f}\t{Tsplit_12:.5f}\t{Na_12:.5f}\t{Tsplit_34:.5f}\t{Tsplit_34:.5f}\t{Na_34:.5f}\t{Tsplit:.5f}\t{Tsplit:.5f}\t{Na:.5f}".format(nsam=nsam_tot[locus], nSNPs=int(nSNPs[locus]), rho=rho[locus], L=L[locus], nsamA=nsamA[locus], nsamB=nsamB[locus], nsamC=nsamC[locus], nsamD=nsamD[locus], N1=N1_vec[locus], N2=N2_vec[locus], N3=N3_vec[locus], N4=N4_vec[locus], M12=M12_vec[locus], M21=M21_vec[locus], M13=M13_vec[locus], M31=M31_vec[locus], M14=M14_vec[locus], M41=M41_vec[locus], M23=M23_vec[locus], M32=M32_vec[locus], M24=M24_vec[locus], M42=M42_vec[locus], M34=M34_vec[locus], M43=M43_vec[locus], Tsc_12=Tsc_12[sim], Tsc_13=Tsc_13[sim], Tsc_14=Tsc_14[sim], Tsc_23=Tsc_23[sim], Tsc_24=Tsc_24[sim], Tsc_34=Tsc_34[sim], Tsplit_12=Tsplit_12[sim], Na_12=Na_12_vec[locus], Tsplit_34=Tsplit_34[sim], Na_34=Na_34_vec[locus], Tsplit=Tsplit[sim], Na=Na_vec[locus]))
    outfile = open("{simulationpath}/priorfile.txt".format(simulationpath=simulationpath), "w")
    outfile.write(priorfile)
    outfile.close()

