## Introduction
This repository produces all outputs necessary for the CAID assessment starting from the first CAID edition, see
[Nature Methods article](https://www.nature.com/articles/s41592-021-01117-3).

In particular the **AlphaFold branch** (this branch) integrates AlphaFold predictions and update the assessment.

AlphaFold disorder predictions are generated starting from the predicted structures and using
[AlphaFold-disorder (GitHub)](https://github.com/BioComputingUP/AlphaFold-disorder)


## Index
* Requirements
* Installation
* Data
* Usage
* License


### Dependencies
CAID software packages wraps the vectorized_cls_metrics repository (https://github.com/marnec/vectorized_cls_metrics), 
which performs the calculations of the classification metrics used throughout CAID.

* `Python 3.6+`
* `vectorized_metrics` https://github.com/marnec/vectorized_cls_metrics
* `numpy`
* `matplotlib`
* `seaborn`
* `scipy`
* `pandas`

## Installation
Installation is only possible on Unix systems. In order to install the package follow these steps:
Typical install time is around 1 minute.

1. Clone or download the package from the GitHub repository. 
**N.B.** specify the branch (`-b alphafold --single-branch`) 


    git clone -b alphafold --single-branch https://github.com/BioComputingUP/CAID.git

2. CAID relies on `vectorized_cls_metrics` library. Clone or download the package from the GitHub repository


    git clone https://github.com/marnec/vectorized_cls_metrics

3. Add the `vectorized_cls_metrics` library to the PYTHONPATH environmental variable:


    export PYTHONPATH="${PYTHONPATH}:/path/where/the/library/was/cloned"


The library is successfully installed. In order to be able to copy-paste commands without the need of customize paths
the CAID package should be placed in this folder structure:

```
CAID-root
├── code                --> (THIS repository)
└── data
    ├── annotations
    ├── predictions
    │   ├── binding
    │   └── disorder
    ├── references
    │   ├── binding
    │   └── disorder
    └── results
```
  
## Data
Content of the `data` folder is available [here](https://idpcentral.org/caid/data/1_alphafold/)

## Run the assessment
The following command populates the `results` folder 

    bash launch_all.sh

The following Python Notebook generates data (figures and tables) in the `manuscript` folder

    jupyter notebook caid_paper.ipynb
        


## License
[CC BY 3.0](https://creativecommons.org/licenses/by/3.0/)



