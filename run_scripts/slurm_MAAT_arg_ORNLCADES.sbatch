#!/bin/bash

#SBATCH --account=cades-ccsi
#SBATCH -A ccsi     ### account to run job under
#SBATCH -p testing  ### partition to run job in (sinfo)
#SBATCH --mem=0     ### <number>[unit], when 0, all available memory requested
#SBATCH -N 1        ### nodes requested for job
#SBATCH -n 32       ### tasks job will run per node
#SBATCH -c 1        ### cpu cores each task requires
#SBATCH -t 00:10:00 ### requested walltime
#SBATCH --no-requeue
#SBATCH --mail-type=ALL
#SBATCH --mail-user=1aj@ornl.gov
#SBATCH -o <insert_directory_here>/logs/${SLURM_JOB_NAME}.o%j


#limit coredumpsize 1000000
#limit stacksize unlimited

# Set necessary environment variables

export R_LIBS_USER="/home/alp/bin/Rlibs"

# Run script

cd "<insert_directory_here>/MAAT/src"

echo
echo $ARGS
echo

Rscript run_MAAT.R $ARGS

echo
echo "`date` -- MAAT EXECUTION HAS FINISHED"·


### END ###
