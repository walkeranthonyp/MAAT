# Leaf model photosynthesis solutions and temperature response curves #
# Used in Walker et al. GMD submission #


This is an example of how to us the multi-assumption architecture and testbed (MAAT) to run factorial ensembles. 
Two ensmebles are specified, 1) to run various numerical and analytical solutions to the photosynthesis equations with multiple representations of stomatal conductance, and 2) to run multiple combinations of temperature response functions for $V_cmax$.


### set up the example ###

These instructions assume that MAAT code has been installed on your machine (as detailed on the main repo README).
 


* From within the highest level maat source code directory set up a MAAT project with the leaf model:
```bash 
./run_scripts/setup_MAAT_project.bs leaf <projectpath>
```
where `<projectpath>` is the full path of where the project is to be set up.
The lowest level directory in the path will be created if it does not already exist.


* Copy the leaf example:
```bash 
./run_scripts/setup_MAAT_example.bs leaf <projectpath>
```
Expects the path to exist.


* Run the process and parameter SA example:
Change directory to the project directory: 
```bash
cd <projectpath>
```  
and execute the job by running locally:  
```bash
./call_run_MAAT_leafeg.bs
```  
This will run the ensemble for the the leaf model in MAAT. 
Options that have been hard coded into the the above `call_*.bs` script is the request to run these ensembles over multiple cores, 4 cores (the MAAT default), for each ensemble. 


* Generate figures from MAAT output.
From the project directory open `plot_figures.R` using your favorite text editor (e.g. vim) or IDE (e.g. RStudio).
There are a number of configurable options but the only one needed is to set the date on which the ensemble was run replace `##OUTPUTDATE##` in `odate <- '##OUTPUTDATE##'` on line 19 with the ensemble run date in YYYY-MM-DD format.
The figures can then be run from the command line with:
```bash
Rscript plot_figures.R
```  
This will generate the photosynthesis figures in Walker et al. manuscript submitted to GMD.
These will be placed in the `tables` and `plots` subdirectories of the results directory: `<projectpath>/results/#OUTPUTDATE#`.
Note that this plotting script depends on the R packages `latticeExtra`, `stringr`, `grid`, and `viridis` to function.



<!-- END -->
