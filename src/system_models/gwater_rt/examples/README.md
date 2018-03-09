# Simple ground water model sensitivity analysis (SA) #
# Used in Walker et al. GMD submission #


This is an example of how to us the multi-assumption architecture and testbed (MAAT) to run a parocess and parameter SA using the simple groundwater model of Dai et al. (2017).


### set up the example ###

These instructions assume that MAAT code has been installed on your machine (as detailed on the main repo README).
 


* From within the highest level maat source code directory set up a MAAT project with the simple groundwater model:
```bash 
./run_scripts/setup_MAAT_project.bs gwater_rt <projectpath>
```
where `<projectpath>` is the full path of where the project is to be set up.
The lowest level directory in the path will be created if it does not already exist.


* Copy the simple groundwater example:
```bash 
./run_scripts/setup_MAAT_example.bs gwater_rt <projectpath>
```
Expects the path to exist.


* Run the process and parameter SA example:
Change directory to the project directory: 
```bash
cd <projectpath>
```  
and execute the job by either, 1) submitting to a que:  
```bash
./call_qsub_MAAT_SA.bs
```  
or 2) running locally:  
```bash
./call_run_MAAT_procSA.bs
./call_run_MAAT_saltSA.bs
```  
Either of these commands will run both the process and parameter SA ensemble for the the groundwater model in MAAT. 
Options that have been hard coded into the the above `call_*.bs` scripts are to run with a sample n of 1,000 for the process SA and a sample n of 1,000,000 for the parameter SA.
Also hard coded is the request to run these ensembles over multiple cores, 32 cores, for each ensemble. 
On a queued system this will run each SA ensemble on a separate node (if available) and each node will be required to have 32 cores.


* Post-process MAAT output and calculate sensitvity indices.
From the project directory open `analysis_MAAT_gwGMDms.R` using your favorite text editor (e.g. vim) or IDE (e.g. RStudio).
There are a number of configurable options but the only one needed is to set the date on which the ensemble was run replace `#OUTPUTDATE#` in `odate <- #OUTPUTDATE#` on line 34 with the ensemble run date in YYY-MM-DD format.
The analysis can then be run from the command line with:
```bash
Rscript analysis_MAAT_gwGMDms.R
```  
This will calculate the sensitivity indices, resample and bootstrap the sensitvity indices to test for convergence, write a pdf table of the process sensitivity indices, and plot the parameter sensitivity indices.
These will be placed in the `tables` and `plots` subdirectories of the results directory: `<projectpath>/results/#OUTPUTDATE#`.

