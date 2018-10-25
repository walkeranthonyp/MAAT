# The Multi-Assumption Architecture and Testbed (MAAT) modelling system. #



### MAAT Version 1.0 ###

The multi-assumption architecture and testbed (MAAT) is a modelling framework designed to facilitate simple and rapid comparison of multiple modelling hypotheses and assumptions, i.e. ways in which to represent processes, in a systems context.
MAAT is designed to easily specify and run large ensembles of simulations.
Ensembles can vary in process representation, parameter values, and environmental conditions.
Built in sensitivity analysis and uncertainty quantification (SA/UQ) can be used to assess variability in model caused output by multiple ways in which to represent processes. 

The MAAT framework is a model wrapper and a model coding syntax that can be used to define a model of a system.
The wrapper and the model system are separated and the syntax allows the wrapper to pass system-model specific information to the model without any system-model specific information contained in the wrapper. 
All system-model specific information is contained within the system model and initialisation files. 
This separation of wrapper and system model allows rapid development of new system models and system model components (e.g. process representations) without the need to edit the wrapper.

Current systems models that come packaged with the MAAT are a leaf-scale model of photosynthesis, a simple ground-water model detailed in Dai et al (2017 WRR), and a simple 'Hello World' type example and template for building new system models. 

MAAT is written in R and uses the object oriented programming package 'proto' to develop the wrapper and model objects. 
MAAT has been developed on Linux Ubuntu 16.04 and is designed to be run from the command line (and in some instances from within the R console or IDE like Rstudio).


### Use guidelines ### 

If you want to use this code then please do!
If you plan to publish using this code the please be respectful of other people's work and their willingness to open source their code. 
If you are using recent modifications and code updates, please be aware that the author of these developments may not have yet published what they are intending to, please contact them to discuss.
Older code is fine to use without co-authorship, but again there is no harm in contacting code authors to clear this.
This is a balance, and likely no single rule really fits all situations, discussion and communication is the best policy.  


### MAAT set up ###

* Fork this repo and then clone your fork to your local machine. 


*All bash commands below are run from within the highest level maat source code directory* 

* Install R package dependencies ('proto', 'XML', 'xtable', 'randtoolbox'), from the command line (or manually run the script within the R console):
```bash 
Rscript install_dependencies.R
```


* Run unit tests (not necessary but will help confirm code is running on your system). 
Change directory to the system model source directory and open unit\_testing.R in RStudio or similar.
Make sure RStudio is closed before opening the unit testing script as this will allow the script to be opened in the correct working directory. 
```bash 
cd ./src/system_models/<modelobject>/
rstudio unit_testing.R
```
or in a Mac terminal
```bash
open -na Rstudio unit_testing.R
```
where `<modelobject>` is the name of the system model to be tested.
The unit testing script can be run line by line to run a number of tests of the model objects. 
An example from the leaf model object `unit_testing.R` script to run an ACi curve: 
```R 
source('leaf_object.R')
leaf_object$.test_aci(leaf.ca_conc=seq(0.1,2000,50))
```
The MAAT wrapper can also be tested to ensure the wrapper is working correctly and to test the various SA/UQ ensembles provided with MAAT. 
```bash 
cd ./src/
rstudio unit_testing.R
```
or in a Mac terminal
```bash
open -na Rstudio unit_testing.R
```


* Set up a MAAT project:
```bash 
./run_scripts/setup_MAAT_project.bs <modelobject> <projectpath>
```
where `<modelobject>` is the name of the system model to be used in the project and `<projectpath>` is the full path of where the project is to be set up.
The lowest level directory in the path will be created if it does not already exist.
Run the above command with `leaf` as `<modelobject>` and your prefered path to set up the MAAT project. 
Change directory to the project directory and a simple instance of MAAT can be run:  
```bash
cd <projectpath>
./call_run_MAAT.bs
```  
This should provide a simple simulation of Aci curves with three different electron transport functions. 
The variables of the run are defined in `init_MAAT.R` in the project directory.


* Initialisation files. 
Once the above steps have been completed and MAAT is working without error, the next step is to customise the run. 
A MAAT ensemble is defined by the process representations, the parameter values, and the environmental variables that a user defines. 
These values can be defined as static, i.e. values that are invariant across the whole ensemble, or dynamic values i.e. values that are varied across the ensemble. 
The values of the static variables and dynamic variables are defined by the user as either lists in an R script `init_MAAT.R` or as separate XML files `init_user_static.xml` and `init_user_dynamic.xml`. 
These are expected to be found in the highest level project directory. So that multiple simulations can be run from within the same project, these initialisation file names and be appended with `_<runid>` where `<runid>` is a character string that identifies the particular ensemble. 
For example, you can create a new init file, edit it, and then rerun MAAT with the runid as the first argument to the call script.
Assuming you are still in the project directory:
```bash
cp init_MAAT.R init_MAAT_yourrun.R
vim init_MAAT_yourrun.R   # change some of the values of the dynamic variables
./call_run_MAAT.bs yourrun
```  
  


* Options.
MAAT can be configured in many different ways.
The names of the variables and their various options can be found in `<model_object>_options.xml` in the source directory.  
Alternative command line options, their names, and how to specify them on the command line can be found on lines 33 - 110 of `run_MAAT.R`.

 
* Meteorological data files.
The path and filename for meteorological data can be passed as an option to run_MAAT.R through any of the `call_*.bs` scripts.
Currently all the meteorological data must reside in a single file in csv format with the first row (and only the first row) as a column header.
A file named `<modelobject>_user_met.xml` must also exist in the directory `<projectpath>`.
An example xml is copied across to that directory when the project is set up.
The xml is structured to represent the input data of the model object, i.e. the `env` list in the data struture.
To allow MAAT to read your meteorological data file, the values in the `<modelobject>_user_met.xml` list should be the name in the column header of the meteorological file of the variable that corresponds to the variable name in the MAAT `env` list.
Note that not all variables in the `env` list must be specified in the xml or the met data file, in this case the default or the value set in the static or dynamic input is used. 
In the case where a variable appears in the static or dynamic list AND the met data file, the values in the met data file will be used.        


### Further information ###
Further information can be found in the README files in sub-directories:
* How to set up and run examples can be found in `src/system_models/leaf/examples` and `src/system_models/gwater_rt/examples`.
* The MAAT formalism for defining model objects and how to start a template for coding a new model object can be found in `src/system_models`.



### Contribution guidelines ###

* For development please [fork this repo and sync to the main repo](https://help.github.com/articles/fork-a-repo/), create your own dev branch, and then [use the pull request functionality on GitHub](https://help.github.com/articles/creating-a-pull-request/).
Before starting new feature branches please make sure your fork and clone is up to date with the latest master branch on this repo.
 
* Code review - by walkeranthonyp 

* If you using this code and you are developing additional system models, please use pull requests to integrate your models back up with this central repo.
This will help to maintain a central code-base and the principles of open science, facilitating other people's work and benefitting everyone.  



### Sponsorship ###

We are grateful of the support from the U.S. Department of Energy (DOE) Office of Science, Biological and Environmental Research (BER). 

Elements of this work are supported through various DOE BER projects.
[NGEE Tropics](https://ngee-tropics.lbl.gov/) supports the development of the leaf and canopy scale photosynthesis models. 
Development of the MAAT framework and early versions of the photosynthesis models were supported by the [ORNL Terrestrial Ecosystem Science SFA](https://tes-sfa.ornl.gov/).    



### Please direct questions to ###

Anthony Walker (walkerap@ornl.gov)



<!-- END -->
