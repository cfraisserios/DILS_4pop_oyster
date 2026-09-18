- [get the pipeline](#get-the-pipeline)
- [dependencies](#dependencies)
    + [R libraries](#r-libraries)
    + [Python libraries](#python-libraries)
    + [others](#others)
- [Snakefile](#snakefile)
- [config file](#config-file)
- [example](#example)
- [oyster updates](#oyster updates)

Pipeline using snakemake to perform demographic inferences in 4-population models. Three topologies are possible as well as different migration relationships making a maximum of 768 comparable models, depending on the user's specifications. The code is the same as Camille Roux's code available at https://github.com/popgenomics/DILS_4pop, except for the modifications reported in the "oyster updates" section.


# get the pipeline  
```
git clone https://github.com/cfraisserios/DILS_4pop_oyster
```  

# dependencies  
Like many modern bioinformatics tools, DILS relies heavily on external libraries and programs developed by real computer scientists.  Please install them before running DILS:   
 
### R libraries
- abcrf (_works with version: 1.9_)  
- data.table (_works with version: 1.14.2_)  
- FactoMineR (_works with version: 2.4_)  
- tidyverse (_works with version: 1.3.2_)  
- ggpubr (_works with version: 0.4.0_)  
- viridis (_works with version: 0.6.2_)  
  
In a R terminal:  
```
for(library in c('abcrf', 'data.table', 'FactoMineR', 'ggpubr', 'viridis')){
	install.packages(library)
}
```
  
### Python libraries  
- numpy (_works with version: 1.21.5_)  
- biopython (_works with version: 1.79_)  
  
### others  
- pypy ([fast implementation of Python](https://www.pypy.org/), _works with version: Python 3.6.9  [PyPy 7.3.1 with GCC 7.3.1 20180303 (Red Hat 7.3.1-5)]_)  
- snakemake ([installation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html), _works with version: 6.15.1_)  
- java (_works with version: openjdk 11.0.16 2022-07-19_)  
- msnsam  
- msms  
  
msnsam has to be compiled as follows:  
```
cd bin/msnsam_src
./clms
ln -s $PWD/msnsam ../msnsam
```  
  
  
The msms binary is located in:
> bin/msms3.2rc-b163.jar
  
# Snakefile  
Please, adapt the **line 4** of bin/Snakefile to your own system. This is simply the path to the scripts deposited in DILS_4pop/bin:  
- binpath: path to binaries for ABC inferences  
  
Currently, the line you should see and modify is:
```
binpath="/home/croux/Programmes/DILS_4pop/bin" # binpath: path to binaries for ABC inferences
```
   
# config file  
DILS requires a second input file (in addition to the fasta) to be executed. The second file, in yaml format, contains various information crucial to the analysis.  
Here is an example with **config.yaml**    
```
inputFile: /home/croux/Documents/zoe/data/RNAseq_nucl.fasta # full pathway to the input file (fasta format)
datapath: /home/croux/Documents/zoe/DILS/4pop_v1 # full pathway to the directory where the analysis will be performed, i.e, where the results will be stored
region: coding
nameA: E1
nameB: W3
nameC: W1
nameD: W2
maxN: 0.5
nMin: 8
Lmin: 200
mu: 0.00000000731
rec: 0.00000000731
N_bound_min: 0  # effective population sizes
N_bound_max: 10
T_bound_min: 1 # time of demographic events
T_bound_max: 10
M_bound_min: 0.4 # migration rates in 4.N.m
M_bound_max: 40
shape_bound_min: 1 # shape parameters of the Beta distributions
shape_bound_max: 20
topologies: topo2 #topo1 or topo2 or topo3 or topo1,topo2 or topo1,topo3 or topo2,topo3 or topo1,topo2,topo3
migAB: 0,1 #0 or 1 or 0,1
migAC: 0
migAD: 0
migBC: 0,1
migBD: 0
migCD: 0
```

generalities:  
- inputFile: path to the fasta input file  
- datapath: path to directory where inferences will be performed  
- region="coding" or "noncoding", depending on the nature of the data  
- nameA, nameB, nameC and nameD: names of the four populations already present in the identifiers of the sequences of the fasta file.  
- maxN=0.5 # loci with a proportion of missing data (i.e, non A, T, C or G) above maxN are excluded.  
- nMin=8 # loci with less than nMin gametes are excluded.  
- Lmin=200 # loci with less than Lmin valid positions (without missing data, and in at least nMin copies) are excluded.  
- mu=0.00000000731 # mutation rate per position and per generation. Used only if mutation='mu'.  
- rec=0.00000000731 # recombination rate per position and per generation.  

prior distributions in coalescent units:  
- N_bound_min=0  # effective population sizes  
- N_bound_max=10  
- T_bound_min=1 # time of demographic events  
- T_bound_max=10  
- M_bound_min=0.4 # migration rates in 4.N.m  
- M_bound_max=40  
- shape_bound_min=1 # shape parameters of the Beta distributions  
- shape_bound_max=20  
  
models to explore:  
 - topologies: **topo1** (only test (A, ((B, C), D)) ), **topo2** (only test (A, (B, (C, D))) ), **topo3** (only test ((A, B), (C, D))), **topo1,topo2** (topo1 _versus_ topo2), **topo1,topo2,topo3** (topo1 _versus_ topo2 _versus_ topo3), etc ...)   
 - migAB: **0** (assumes no migration between A and B), **1** (assumes migration between A and B) or **0,1** (will test migration _versus_ isolation between A and B)  
 - similarly for: migAC, migAD, migBC, migBD, migCD
  
# example  
Short example (dry run):  
```
cd DILS_4pop/example/
mkdir analysis
cd analysis
snakemake -n -p -s ../../bin/Snakefile
```
   
Full run:  
```
snakemake -p --snakefile ~/Programmes/DILS_4pop/bin/Snakefile --configfile ~/Programmes/DILS_4pop/example/analysis/config.yaml --cores 6
```
  
Of course, the file **config.yaml** has to be adapted by the user (i.e, by you), and can be located everywhere (not only in the DILS_4pop directory).  


# oyster updates
The original scripts have a "@true" extension. In the modified scripts, see tag "CFraisseRios" for modifications.
1. Python script **bin/priorgen_4pop.py** has been modified to correct a typo in the "topo3" with the times.
2. Snakemake script **bin/Snakefile** has been modified to increase the number of multilocus simulations. If inference is taking too long, decrease this number.
