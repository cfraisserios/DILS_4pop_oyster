#!/usr/bin/python
import sys
import os

# binpath: path where priorgen_4pop.py, msnsam and mscalc_4pop.py are located
#binpath='/home/croux/Documents/quentin_rougemont/4pop/reabc/example/bin/dev'

# datapath: path where the bpfile is present
#datapath='/home/croux/Documents/quentin_rougemont/4pop/reabc/example/bin/dev'

model=None
nMultilocus=None
mutation=None
AB=None
AC=None
AD=None
BC=None
BD=None
CD=None

N_bound = []
T_bound = []
M_bound = []
shape_bound = []

for arg in sys.argv:
	tmp = arg.split('=')
	if tmp[0]=='binpath': # where the code is.
		binpath=tmp[1]
	if tmp[0]=='datapath': # where the fasta and bpfile are.
		datapath=tmp[1]
	if tmp[0]=='simulationpath': # where the simulation has to be performed.
		simulationpath=tmp[1]
	if tmp[0]=='model':
		model=tmp[1]
	if tmp[0]=='nMultilocus':
		nMultilocus=int(tmp[1])
	if tmp[0]=='mutation':
		mutation=tmp[1]
	if tmp[0]=='AB':
		AB=tmp[1]
	if tmp[0]=='AC':
		AC=tmp[1]
	if tmp[0]=='AD':
		AD=tmp[1]
	if tmp[0]=='BC':
		BC=tmp[1]
	if tmp[0]=='BD':
		BD=tmp[1]
	if tmp[0]=='CD':
		CD=tmp[1]
	if tmp[0]=='N':
		N_bound.append(float(tmp[1]))
	if tmp[0]=='T':
		T_bound.append(float(tmp[1]))
	if tmp[0]=='M':
		M_bound.append(float(tmp[1]))
	if tmp[0]=='shape':
		shape_bound.append(float(tmp[1]))
	if tmp[0]=='target': # modelComp / posterior / opti1 / opti2 / opti3 / opti4 / opti5
		target=tmp[1]
	if tmp[0]=='posteriorFile': # name of the file containing estimated values (posterior or lower levels of optimization)
		posteriorFile=tmp[1]

if model==None:
	print("ERROR: please specify a correct model (i.e model=topo1, model=topo2 or model=topo3)\n")
	sys.exit()
if nMultilocus==None:
	print("ERROR: please specify a number of simulations to perform (for instance: nMultilocus=100)\n")
	sys.exit()
if mutation==None:
	print("ERROR: please specify a method to implement mutations: '=mu' (mutations conditionned on 4.N.mu) or '=nSNPs' (mutations conditionned on the number of nSNPs)\n")
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

if mutation=='mu':
	simulator='{binpath}/msnsam'.format(binpath=binpath)
	commandeMS = '-t tbs '
else:
	simulator='java -jar {binpath}/msms3.2rc-b163.jar'.format(binpath=binpath)
	commandeMS = '-s tbs '

if model in ['topo1', 'topo2', 'topo3']:
	if model=='topo1': # (A; (B;C); D)
		commandeMS += '-r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 3 tbs -m 3 1 tbs -m 2 4 tbs -m 4 2 tbs -m 1 2 tbs -m 2 1 tbs -m 2 3 tbs -m 3 2 tbs -m 1 4 tbs -m 4 1 tbs -m 3 4 tbs -m 4 3 tbs -em tbs 2 3 0 -em tbs 3 2 0 -em tbs 3 1 0 -em tbs 1 3 0 -em tbs 3 4 0 -em tbs 4 3 0 -em tbs 2 4 0 -em tbs 4 2 0 -em tbs 1 2 0 -em tbs 2 1 0 -em tbs 4 1 0 -em tbs 1 4 0 -ej tbs 3 2 -en tbs 2 tbs -ej tbs 2 4 -en tbs 4 tbs -ej tbs 4 1 -eN tbs tbs'
	if model=='topo2': # (A; B; (C; D))
		commandeMS += '-r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 3 tbs -m 3 1 tbs -m 2 4 tbs -m 4 2 tbs -m 1 2 tbs -m 2 1 tbs -m 2 3 tbs -m 3 2 tbs -m 1 4 tbs -m 4 1 tbs -m 3 4 tbs -m 4 3 tbs -em tbs 3 4 0 -em tbs 4 3 0 -em tbs 3 2 0 -em tbs 2 3 0 -em tbs 3 1 0 -em tbs 1 3 0 -em tbs 2 4 0 -em tbs 4 2 0 -em tbs 2 1 0 -em tbs 1 2 0 -em tbs 4 1 0 -em tbs 1 4 0 -ej tbs 3 4 -en tbs 4 tbs -ej tbs 4 2 -en tbs 2 tbs -ej tbs 2 1 -eN tbs tbs'
	if model=='topo3': # (A; B); (C; D))
		commandeMS += '-r tbs tbs -I 4 tbs tbs tbs tbs 0 -n 1 tbs -n 2 tbs -n 3 tbs -n 4 tbs -m 1 2 tbs -m 2 1 tbs -m 1 3 tbs -m 3 1 tbs -m 1 4 tbs -m 4 1 tbs -m 2 3 tbs -m 3 2 tbs -m 2 4 tbs -m 4 2 tbs -m 3 4 tbs -m 4 3 tbs -em tbs 1 2 0 -em tbs 2 1 0 -em tbs 1 3 0 -em tbs 3 1 0 -em tbs 1 4 0 -em tbs 4 1 0 -em tbs 2 3 0 -em tbs 3 2 0 -em tbs 2 4 0 -em tbs 4 2 0 -em tbs 3 4 0 -em tbs 4 3 0 -ej tbs 2 1 -en tbs 1 tbs -ej tbs 4 3 -en tbs 3 tbs -ej tbs 3 1 -eN tbs tbs'
else:
	print("ERROR: the specified model ({model}) is not correct (i.e; topo1, topo2 or topo3)".format(model=model))
	sys.exit()

# get number of loci
if os.path.isfile("{datapath}/bpfile".format(datapath=datapath)) == False:
	sys.exit("\n\t\033[1;31;40mERROR: bpfile was not found\n\033[0m\n")

infile = open('{datapath}/bpfile'.format(datapath=datapath), 'r')
tmp = infile.readline()
tmp = infile.readline().strip().split('\t')
nLoci = len(tmp)
infile.close()

# commande for the simulation
commande='python3 {binpath}/priorgen_4pop.py datapath={datapath} simulationpath={simulationpath} mutation={mutation} model={model} nMultilocus={nMultilocus} AB={AB} AC={AC} AD={AD} BC={BC} BD={BD} CD={CD} N={Nmin} N={Nmax} T={Tmin} T={Tmax} M={Mmin} M={Mmax} shape={shape_min} shape={shape_max} target={target} posteriorFile={posteriorFile} | {simulator} tbs {nSimulationsTot} {commandeMS} | pypy {binpath}/mscalc_4pop.py datapath={datapath} simulationpath={simulationpath}'.format(mutation=mutation, model=model, nMultilocus=nMultilocus, AB=AB, AC=AC, AD=AD, BC=BC, BD=BD, CD=CD, binpath=binpath, nSimulationsTot=nMultilocus*nLoci, commandeMS=commandeMS, datapath=datapath, Nmin=min(N_bound), Nmax=max(N_bound), Mmin=min(M_bound), Mmax=max(M_bound), Tmin=min(T_bound), Tmax=max(T_bound), shape_min=min(shape_bound), shape_max=max(shape_bound), simulationpath=simulationpath, simulator=simulator, target=target, posteriorFile=posteriorFile)

print(commande)
os.system(commande)
