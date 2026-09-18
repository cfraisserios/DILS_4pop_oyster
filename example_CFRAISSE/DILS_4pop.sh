#!/bin/bash
#SBATCH --partition=long
#SBATCH --job-name=dils4D
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=2G
#SBATCH --time=7-00:00:00

# last edit: 2026-05-27
# authors: Christelle Fraïsse

#-# usage: sbatch DILS_4pop.sh config.yaml
#-# prog required: DILS_4pop ; msnsam (it needs to be recompiled on the cluster)
#-# files required: 4pop.fasta
#-# description: Run 4D demographic inferences with DILS_4pop.
#	https://github.com/popgenomics/DILS_4pop

### START ###

#-# packages
export TCL_LIBRARY=/shared/software/miniconda/envs/pypy-2.7-5.10.0/lib/tcl
export R_LIBS_USER="/shared/home/cfraisserios/R/x86_64-conda-linux-gnu-library/4.5"
module load pypy/2.7-5.10.0
module load snakemake/7.7.0
module load r/4.5.1
module load python/3.12

#-# environment
ncpu=$SLURM_CPUS_PER_TASK
mem=$SLURM_MEM_PER_CPU 

#-# paths
binpath="/shared/home/cfraisserios/bin/dils4D/bin"
wdpath="/shared/home/cfraisserios/bin/dils4D/example_CFRAISSE"


#-# job start --------------------------------------------------------------------------------------------------------------------------------------------------------

TimeStart=$(echo "#START `date` ")

snakemake --snakefile ${binpath}/Snakefile -p -j 150 --cores 150 \
        --configfile ${1} \
        --cluster-config ${wdpath}/cluster.json \
        --latency-wait 300 \
        --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.n} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu}"

TimeEnd=$(echo "#STOP `date` ")

#-# job end ----------------------------------------------------------------------------------------------------------------------------------------------------------
sstat -j $SLURM_JOB_ID.batch --format=JobID,NodeList,MaxRSS
echo "${TimeStart} ${TimeEnd}"
