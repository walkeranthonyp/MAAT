# MAAT System Models Directory #


This directory contains each of the system models coded within the MAAT framework.
Each system model is contained within its own directory.
A number of READMEs and examples are provided to allow a user to familiarise themselves with a few applications of MAAT and to help code a new model.
Examples are set up so that runninng them should be straight forward with just a rudimentary knowledge of bash and R. Examples can be found in the `leaf/examples` and `gwater_rt/examples` directories and both should reproduce results presented in Walker et al. (2018 GMDD).
To create new system models in MAAT, a good working knowledge of R is assumed.
This README and the README in the template directory will help a user get started setting up and coding their own system model within MAAT. 

Unit testing functions are provided with every model object. Running these unit tests will help confirm that code is working without error messages and will provide the user with a very initial idea of what each system model can do.
To run the unit tests, change directory to the system model source directory and open unit\_testing.R in RStudio or similar.
Make sure RStudio is closed before opening the unit testing script as this will allow the script to be opened in the correct working directory. 
```bash 
cd ./src/system_models/<modelobject>/
rstudio unit_testing.R
```
where `<modelobject>` is the name of the system model to be tested.
The unit testing script can be run line by line to run a number of tests of the model objects. 
An example from the leaf model object `unit_testing.R` script to run an ACi curve: 
```R 
source('leaf_object.R')
leaf_object$.test_aci(leaf.ca_conc=seq(0.1,2000,50))
```


## MAAT system model structure ##

The MAAT framework is a model wrapper and a model coding syntax that can be used to define a model of a system.
The wrapper and the model system are separated and the syntax allows the wrapper to pass system-model specific information to the model without any system-model specific information contained in the wrapper. 
All system-model specific information is contained within the system model and initialisation files. 
This separation of wrapper and system model allows rapid development of new system models and system model components (e.g. process representations) without the need to edit the wrapper.
The remainder of this README is intended as a guide to the model coding struture and syntax to facilitate development of new system models. 

System models in MAAT are composed of three files of source code: `<modelobject>_object.R`, `<modelobject>_system_functions.R`, and `<modelobject>_functions.R`.
Respectively, these three files contain: the formal MAAT system model object structure including data structures (coded as a `proto` object), system representation functions (coded as R functions), and process representation functions for each process (also coded as R functions).       
Ancilliary files are also provided in the system model directory such as a unit testing R script and initialisation files.
The three main system model source code files are now described in a little more detail.  

### System model object (SMOs) ###
The system model object (SMO, an R proto object) contains data structures and a number of generic, semi-generic, and system model specific test functions.
In the proto package in R, the period refers to the object itself and is loosely equivalent to `self` in python objects. 
The data structure is composed of five named lists: 
* `fnames`, a named list of character strings 
* `pars`, a named list of numeric values 
* `env`, a named list of numeric values
* `state_pars`, a named list of numeric values
* `state`, a named list of numeric values

These lists are described in Walker et al. (2018 GMDD) and contain the names and values of, respectively: processes, parameters, environment or model input variables (e.g. meteorological data), dynamic parameters or secondary state variables, and primary state variables.
Making these data structures specific to a system model is a key task when setting up a new system model. 

The functions contained within the system model object are: 
* a generic `build` function 
* generic `configure` and `configure_sublist` functions
* a generic `run_met` function 
* a semi-generic `run` function 
* a semi-generic `output` function  
* system model specific test functions 

The `build` and `configure` functions are generic across system models, and vary only if a system model has a child model object embedded within it, e.g. a leaf model is embedded as a child within a canopy model.
The `configure_sublist` function is entirely generic.
The `run_met` function is also entirely generic and runs the model over an entire timeseries of environment variables, often a dataset of meteorological variables.
The `run` function calls the system representation function and the `output` function.
The `run` function is semi-generic beacuse it can be modified to suite a particular system model.
For example, a function could be called that calculates a number of parameter values that are parameters common to a number of the system representation functions, avoiding the need to set these parameters in every system representation function.
The `output` function provides the output from the model after a single iteration (e.g. timestep) of the model.
Usually this is the `state` list converted into a numeric vector, but various output options can be customised.
The test functions are used for unit testing and their names are prefixed with a period, e.g. `.test`, which means they are not copied when the model object is cloned.
The test functions are model specific and can be used to verify the model is operating as expected under various conditions.  


### System representation functions (SRFs)  ###
The system representation functions (SRFs, written as R functions) are the functions that represent the whole model system and are called by the `run` function of the SMO.
Multiple SRFs are possible where the system can be conceptualised in different ways.
For example, a leaf photosynthesis system can be conceptualised using principles of enzyme kinetics or more empirical concepts of light use efficiency.
These SRFs are entirely non-generic and are required to be coded for a specific system.
These SRFs access the named variables in the SMO data structure by combining the `.` and `$` syntax of R, e.g.:
```R 
.$fnames$<processname>
``` 
will return the value of the representation of process `<processname>` in the `fnames` list of the SMO data structure.
```R
.$pars$<parametername>
``` 
will return the value of the parameter `<parametername>` in the `pars` list of the SMO data structure.  

The first argument to the SRFs is the SMO and is named `.`.
Using the `.` as the only argument to this function call passes the SMO to the function allowing it to access the SMO data structure.
Passing the period as the first argument is necessary as the process representation functions are external to the SMO for reasons of cleaner code.

Local variables in these SRFs are discouraged.
Any variables used should reside in the the SMO data structure.
Explicit calculations in these SRFs are discouraged.
Where there is any chance that a calculation can be done with multiple alternative methods, the calculation should be considered a process and the calculation and any alternative should be given its own function, described below. 
To designate a calculation a process, the name of the process should be determined and added as an element of the `fnames` list in the SMO.
The value of this element in the `fnames` list is a character string that is the name of the function that is a specific representation of the process.
The character string allows the SRF to call the process representation function using the `get` function and the syntax:
```R
get(.$fnames$processname)(.)
``` 

### Process representation functions (PRFs) ###
The process representation functions (PRFs, written as R functions) are the functions that represent individual processes and are called by the SRF as described above.
As with the SRFs, the first argument to these functions is the SMO named `.` and the PRFs access the named variables in the SMO data structure using the `.` and `$` syntax described above.
As with the SRFs, local variables are discouraged, keeping all variables in the SMO data structure allows them to be easily varied in any SA/UQ type activities.
Different to the SRFs, explicit calculations in the PRFs are required to represent a particular hypothesis or assumption that represents a particular process.
These functions can be named anything a developer likes, remembering that the names of these PRFs are the values of the `fnames` list in the SMO data structure.
However, the use of a function naming convention is encouraged.
The following convention is used:
```
f_<processname>_<representationname>
```
where `f` identifies the function as a PRF; `<processname>` is the process name used in the `fnames` list in the SMO data structure; and `<representationname>` is a name that identifies the specific representation of the process, often a citation using first author name and year.
Using this naming convention will allow the `create_default_XMLs.R` script (described below) to generate an XML of all the representations for each process.


## Guide to integrating a new PRF in an existing system model ##

To integrate new process representations in an existing syetm model first familiarise yourself with the SMO and the SRF for the system model.
Find the name of the process of interest in the SMO `fnames` data struture or SRF.
Use this and the above naming convention to define a new PRF in the `<modelobject>_functions.R` file.
Code the calculations for the new PRF making sure to check for variables required by the new PRF that already exist in the SMO data structure.
These variables might exist in the `pars`, `env`, `state_pars`, or `state` lists.
This is all that is necessary.
Once the PRF has been coded, check the unit testing functions to see if there is one that can already test the process of interest and run the test function with the new PRF as an argument.
If a unit testing function does not exist, define one and test the new PRF. 
  

## Guide to setting up a template for a new system model ##

* Make sure you have forked this repo and cloned your fork as decribed in the README on the highest level MAAT directory. 


*All bash commands below are run from within the highest level maat source code directory* 

* Clone the template to create new `<modelobject>`:
```bash 
./run_scripts/create_new_model_object.bs <modelobject> 
```
where `<modelobject>` is the name of the system model you want to create.
This will create a system model directory named `<modelobject>` in this `system_models` directory.
All instances of the word template in the clones files and filenames will be replaced with `<modelobject>`, allowing the beginning of speficity for the new system model.  

* The new system model developer can then begin to develop their model. 
This involves:
    * developing the SRF(s), ensuring that each process is called by accessing the `fnames` list in the SMO data struture and is not hard coded in the SRF. `state` and `state_pars` variables that are set in the SRF should also be part of the SMO data structure. 
    * developing the PRFs to be called by the SRF, ensuring that parameters are accessed from the SMO data structure and are not local.
    * name the elements of the lists in the SMO data structure in the same way they have been labelled in the SRF(s) and PRFs
    * develop other non-generic elements of the SMO: potentially the `run` and `output` functions, and any unit testing functions
    * the `unit_testing.R` script can be modified to call various unit testing functions with various options to allow verification that the model is working correctly 

* Once the code has been verified to work correctly, initialisation files can be created.
These are needed for using the model object with the MAAT wrapper that generates SA / UQ type ensembles.
Initialisation files can be customised to the model object by running the XML generation script:  
```bash 
Rscript ./run_scripts/create_default_XMLs.R "mod_obj<-'<modelobject>'" 
```
This script will generate five XMLs that mimic various aspects of the SMO data structure in XML format.
The default XML is the most important as this is used everytime a model object is run within the MAAT wrapper (a model object can be run externally from the MAAT wrapper using the unit testing functions).
The default XML will mimic exactly the three list in the SMO that define a model simulation: the `fnames`, `pars`, and `env` lists.
Within the XML, the values of each named element will take the value of each named element in the SMO data structure in the SMO source code. 


