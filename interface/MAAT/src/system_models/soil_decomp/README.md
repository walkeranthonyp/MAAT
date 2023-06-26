# MAAT Soil decomposition model object directory #


This directory contains the files needed for the soil decomposition model object: ```soil_decomp``` 


## MAAT use of SoilR functions ##

The soil decomposition model object uses modified versions of the [SoilR package](https://www.bgc-jena.mpg.de/TEE/software/soilr/) functions that allow generalisation of a pool and flux model.
A matrix model approach is used by SoilR and is solved using a function in the deSolve R package.
These functions represent a substantial amount of effort on the part of the SoilR designers - Carlos Sierra and Markus Müller - at the Max Planck Institute for Biogeochemistry.
We thank them for their open-access approach to science that has allowed us to use these functions and to modify them to work within MAAT.
If using output from this model object in any publication, please cite the following papers: 


Sierra, C.A., Müller, M., 2015. A general mathematical framework for representing soil organic matter dynamics. Ecological Monographs 85, 505–524. [doi.org/10.1890/15-0361.1](https://esajournals.onlinelibrary.wiley.com/doi/full/10.1890/15-0361.1)


Sierra, C.A., Müller, M., Trumbore, S.E., 2012. Models of soil organic matter decomposition: the SoilR package, version 1.0. Geosci. Model Dev. 5, 1045–1060. [doi.org/10.5194/gmd-5-1045-2012](https://www.geosci-model-dev.net/5/1045/2012/)

     








