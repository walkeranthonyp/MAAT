#!/bin/bash

# bash script to excute MAAT on a cluster using slurm and sbatch script
# this example/default script excutes an MCMC MAAT run on 32 cores, but needs more arguments like met data and evaluation data 
# this script is copied to $PDIR and should only be edited there to customise a run 

# Optional single argument - the MAAT runid
#  - no argument will read 'init_MAAT.R' for initialisation
#  - with argument <ARG> will read 'init_MAAT_<ARG>.R' for initialisation from $PDIR 

# 1st argument - the MAAT runid
RUNID=$1

# directories - defined when 'setup_MAAT_project.bs' is run
SDIR="##SDIR##"
PDIR="##PDIR##"
MOBJ="##MOBJ##"



##########################################
### User defined variables 
#   - edit in this script once copied to $PDIR
#   - do not edit in repo

# met data directory
MDIR=""

# number of processors to use
NP=32

# command line arguments to pass to run_MAAT.R - argument names and options can be found in run_MAAT.R 
ARGS="srcdir<-'${SDIR}' pdir<-'${PDIR}' mod_obj<-'${MOBJ}' runid<-'${RUNID}' procs<-${NP} mdir<-'${MDIR}' multic<-T uq<-T procSA<-F salt<-F factorial<-F mcmc<-T"



##########################################
### DO NOT MODIFY ANYTHING BELOW THIS LINE

FORMATTED_ARGS=$(printf "%q" "$ARGS")

YMD=`date +%Y-%m-%d`

cd $PDIR
sbatch --job-name="${YMD}_MAAT_${RUNID}"  --export=ALL,SDIR=$SDIR,PDIR=$PDIR,MDIR=$MDIR,MOBJ=$MOBJ,NP=$NP,YMD=$YMD,FORMATTED_ARGS="'$FORMATTED_ARGS'" slurm_MAAT_arg_ORNLCADES.sbatch



### END ###
