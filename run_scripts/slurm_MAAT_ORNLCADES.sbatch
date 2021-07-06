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
#SBATCH --mail-user=%u@ornl.gov
#SBATCH -o ##PDIR##/logs/%x.o%j

#limit coredumpsize 1000000
#limit stacksize unlimited


#===============================================================================
# Set necessary environment variables - this can be made a hash variable or should already be an environment variable
export R_LIBS_USER="/home/alp/bin/Rlibs"
#export R_LIBS_USER=/software/dev_tools/swtree/cs400_centos7.2_pe2016-08/R/3.5.0/centos7.2_gnu5.3.0/lib64/R/library/:/home/alp/bin/Rlibs

#===============================================================================
# Run script
cd "##SDIR##"

echo; echo; echo
echo "ARGS passed to batch script:"
echo $ARGS
echo

Rscript run_MAAT.R $ARGS

echo
echo "`date` -- MAAT EXECUTION HAS FINISHED"·



### END ###
