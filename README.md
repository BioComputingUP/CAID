## CAID Assessment
This repository contains the code for [CAID challenge](https://caid.idpcentral.org/challenge) assessment. Upon having predictions and reference sets, you can use this repository to generate the evalutations and metrics. 
CAID software packages wraps the [vectorized_cls_metrics](https://github.com/marnec/vectorized_cls_metrics) repository (with small modifications), which performs the calculations of the classification metrics used throughout CAID. For the details of evaluations, please see the Papers section. 

If you use this code in your research, please **cite the following papers**:

- [CAID2](https://onlinelibrary.wiley.com/doi/full/10.1002/prot.26582) - Conte AD, Mehdiabadi M, Bouhraoua A, Miguel Monzon A, Tosatto SCE, Piovesan D. Critical assessment of protein intrinsic disorder prediction (CAID) - Results of round 2. Proteins. 2023; 91(12): 1925-1934 (2023)

- [CAID1](https://www.nature.com/articles/s41592-021-01117-3) - Necci, M., Piovesan, D., CAID Predictors. et al. Critical assessment of protein intrinsic disorder prediction. Nat Methods 18, 472–481 (2021)


## Installation
To run this package, you need to have `Python 3.8+` installed. 

```
git clone https://github.com/BioComputingUP/CAID.git        # clone the repository

pip install -r requirements.txt                             # install the requirements
```

The repository is structures as below (the demo-data just contains sample data from CAID3 and the results you get from the assessment).
```
CAID                     --> (CAID repository)
 ├── caid.py             --> the script to run the evaluaions
 ├── vectorized_metircs/ --> the assessment library 
 └── demo-data/          --> demo data directory, with sample data from CAID3 challenge  
   ├── predictions/      --> directory containing prediction of each method
   ├── references/       --> directory containing reference fasta file 
   └── results/          --> directory for saving results
```

## Input

### Predictions
In order to run the assessment, you have to have your predictions in CAID ouptut format (see https://caid.idpcentral.org/challenge), where columns correspond to position, residue type, disorder/binding score, and a binary state. If the state is not provided, it will be automatically calculated using a threshold by maximizing f1-score.  

```
>DP01234
1    M    0.892    1
2    E    0.813    1
...
```
Each file must be stored with .caid suffix. You can access and download all CAID challenge results from https://caid.idpcentral.org/challenge/results. 


### References
References must be provided as a single fasta file, includeing the sequence and the labels corresponding to each residue. In the labels, 0 indicates order, 1 indicates disorder/binding/linker, and - denotes that this residue is not included in the assessment. All the CAID challenge references can be downloaded from https://caid.idpcentral.org/challenge/results.

```
>DP01234
MNASDFRRRGKEMVDYMADYLE
000011111000----------
```

## Output

After running the assessment (see usage), the following files are generated.

```bash

# Score distribution for a given method. `rawscore` are all scores, `thresholds` is the unique list of thresholds
<method>.{rawscore,thresholds}.distribution.txt

# `dataset` the metrics for every considered threshold for a given `reference` and `method`
# `bootstrap` same as `dataset` but for every boostrap sample
# `target` same as `dataset` but for every predicted target
<reference>.analysis.<method>.{bootstrap,dataset,target}.metrics.csv

# Optimal thresholds for every calculated metric for a given `reference` and `method`
<reference>.analysis.<method>.thr.csv

# `ci` confidence intervals for all methods for a given `reference` and `optimization`
# `bootstrap` metrics for each method and each boostrap sample for every method for a given `reference` and `optimization`  
<reference>.all.{ci,bootstrap}.<optimization>.metrics.csv

# `cmat` confusion matrix for every method for a given `reference` and `optimization`
# `metrics` metrics for each method
<reference>.all.dataset.<optimization>.{cmat,metrics}.csv

# `cmat` confusion matrices for all methods and all thresholds for a given `reference`
# `pr` precision-recall data for all methods
# `roc` ROC data for all methods
# `predictions` scores and binary predictions for all methods at the residue level
<reference>.all.dataset._.{cmat,pr,predictions,roc}.csv

# metrics for all methods a the target level for a given `reference` and `optimization`
<reference>.all.target.<optimization>.metrics.csv

```
 
## Usage
To run the assessment, you can run the `caid.py` script with arguments explained as below:
```
python3 caid.py <path-to-reference-fasta> <directory-containing-predictions> -o <output-directory>
```
For example, the `demo-data/predictions` folder contains the predictions of 3 predictors from CAID3, and `demo-data/references/disorder_pdb.fasta` is the Disorder-PDB from [CAID3](https://caid.idpcentral.org/challenge/results). The script could be run by: 

```
python3 caid.py demo-data/references/disorder_pdb.fasta demo-data/predictions -o demo-data/results
```



## License
[CC BY 3.0](https://creativecommons.org/licenses/by/3.0/)