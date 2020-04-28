#!/bin/bash

# bash script to excute MAAT on a cluster using slurm and sbatch script

# 1st argument - the MAAT runid
RUNID=$1

SDIR="blank"
PDIR="blank"
MDIR="blank"
MOBJ="blank"
NP=32

YMD=`date +%Y-%m-%d`

ARGS="srcdir<-'${SDIR}' pdir<-'${PDIR}' mod_obj<-'${MOBJ}' runid<-'${RUNID}' procs<-${NP} mdir<-'${MDIR}' multic<-T uq<-T procSA<-F salt<-F factorial<-F mcmc<-T"

FORMATTED_ARGS=$(printf "%q" "$ARGS")

cd $PDIR
sbatch --job-name="${YMD}_MAAT_${RUNID}_MCMC"  --export=ALL,SDIR=$SDIR,PDIR=$PDIR,MDIR=$MDIR,MOBJ=$MOBJ,NP=$NP,YMD=$YMD,FORMATTED_ARGS="'$FORMATTED_ARGS'" slurm_MAAT_arg_ORNLCADES.sbatch

### END ###
