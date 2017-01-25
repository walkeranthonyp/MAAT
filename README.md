# README #

This README describes how to setup and use the Multi-Assumption Architecture and Testbed (MAAT) model.



### MAAT Version 0.1 ###

The multi-assumption architecture and testbed (MAAT) is a modelling framework designed to allow sensitivity analysis and uncertainty quantification (SA/UQ) in a model caused by multiple ways in which to represent a single process (e.g. Collatz et al 1991 vs Farquhar et al 1980 representation of the photosynthetic light repsonse), as well as variability in parameters and driving data. The MAAT is effectively a model wrapper that once initialised automates the SA/UQ. The model that comes pre-packaged with the MAAT is a multi-scale model of photosynthesis (currently leaf and canopy scale). MAAT is written in R and uses the object oriented programming package 'proto' to develop the MAAT and model objects.     



### MAAT set up ###

Fork this repo and then clone your fork to your local machine. From the directory 'maat/run_scripts' copy the 'example_*' run scripts to your <project> directory. Modify the directory names on lines 40-51 in 'example_run_MAAT.R'. An installation script will be developed shortly.  

* Configuration: The basic configuration is the MAAT wrapper object and a model object. The MAAT wrapper object consists of lists of variable names which are combined into dataframes and the model is run over these dataframes using a nested hierarchy of run functions each of which calls a configuration function embedded within the model object to allow automated configuration of the model for the SA/UQ. The final run function in the hierarchy then also runs the model. The model object consists of lists of variables and a primary run function. An associated set of process representation functions accompany the model object in a separate file.     
 
* Dependencies: R packages 'parallel' (now part of the 'base' package), 'proto', 'plyr', 'stringr'. 

* Database configuration: TBD

* How to run tests: Within the source code directory a 'unit_testing' script is provided which can be run line by line to run a number of tests of the model objects and the MAAT wrapper. 

* Deployment instructions: This code is under development and has not yet been used in publications. This code is entirely proprietary. Please do not share the code. Subsequent to publication the plan is to make this code available on the GNU licence (or similar).



### Contribution guidelines ###

* For development please fork this repo, create your own dev branch, and then [use the pull request functionality on BitBucket](https://confluence.atlassian.com/bitbucket/fork-a-teammate-s-repository-774243391.html). Before making a pull request please ensure that your development branch code is up-to-date by [syncing your forked repo with the 'master' branch in the original repo](https://confluence.atlassian.com/bitbucket/create-a-pull-request-774243413.html).
 
* Code review

* Other guidelines



### Please direct questions to ###

Anthony Walker (walkerap@ornl.gov)