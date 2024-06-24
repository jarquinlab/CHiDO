# CHiDO
__Characterization and Integration of Driven Omics__

## Description

CHiDO is a no-code application developed by the Jarquin Lab at the University of Florida's Department of Agronomy. The application
aims to remove technical and financial barriers associated with modeling Genotype-by-Environment (GxE) interactions. 

There is a demo version of CHiDO hosted as a free web application at: https://jarquinlab.shinyapps.io/chido/

## Authors

- Francisco González - University of Florida (UF)
- Julián García - Universidad Politécnica de Madrid (UPM)
- Dr. Diego Jarquin - University of Florida (UF)

## Installation

### Prerequisites

The following packages, listed in `requirements.txt`, are necessary to run CHiDO locally:
- shiny
- shinyjs
- shinyjqui
- shinydashboard
- dplyr
- ggplot2
- viridis
- scales
- gridExtra
- DT
- zip
- markdown
- BGLR

Once you clone/downloade the repository, ensure you are in the CHiDO directory and run the following command in the terminal:
```R
Rscript -e 'packages <- readLines("requirements.txt"); for(package in packages) if (!require(package, character.only = TRUE)) install.packages(package, dependencies = TRUE)'
```

This will ensure all prerequisite librarier and any of their dependencies are installed in your local environment.

## Contact

For questions about using CHiDO, please visit our [documentation](https://github.com/jarquinlab/CHiDO/docs) page for guidance.

If you encounter any issues with the application or have an idea to extend its capabilities,
please submit it to our GitHub issues page, found [here](https://github.com/jarquinlab/CHiDO/issues).

For any additional concerns or interest in collaborating with the project, please email Dr. Diego Jarquin (jhernandezjarqui@ufl.edu).
