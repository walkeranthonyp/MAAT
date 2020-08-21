#!/bin/bash
#===============================================================================

#SBATCH -A ccsi       ### account to run job under
#SBATCH -p batch      ### partition to run job in (sinfo)
#SBATCH --mem=0       ### <number>[unit], when 0, all available memory requested
#SBATCH -N 1          ### nodes requested for job
#SBATCH -n 32         ### tasks job will run per node
#SBATCH -c 1          ### cpu cores each task requires
#SBATCH -t 01:00:00   ### requested walltime
#SBATCH --no-requeue
#SBATCH --mail-type=ALL
#SBATCH --mail-user=${USER}@ornl.gov
#SBATCH -o ##PDIR##/logs/${SLURM_JOB_NAME}.o%j

#limit coredumpsize 1000000
#limit stacksize unlimited


#===============================================================================
# Set necessary environment variables - this can be made a hash variable or should already be an environment variable
export R_LIBS_USER="/home/alp/bin/Rlibs"


#===============================================================================
# Run script
cd "##SDIR##"

echo
echo $FORMATTED_ARGS
echo

Rscript run_MAAT.R $FORMATTED_ARGS

echo
echo "`date` -- MAAT EXECUTION HAS FINISHED"·



### END ###
