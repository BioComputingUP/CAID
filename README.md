## Introduction
The CAID software produces all outputs necessary for a CAID edition, including baselines, references, metrics and plots, 
starting from predictions and a reference (see Data Availability section to know how to obtain this data). 
CAID software packages wraps the vectorized_cls_metrics repository (https://github.com/marnec/vectorized_cls_metrics), 
which performs the calculations of the classification metrics used throughout CAID. More information at 
http://disprotcentral.org/caid


## Index
* Requirements
* Installation
* Demo
* Data
* Usage
* License

## Requirements
### Interpreter
Python 3.6+

### Dependencies
* `vectorized_metrics` https://github.com/marnec/vectorized_cls_metrics
* `numpy`
* `matplotlib`
* `seaborn`
* `scipy`
* `pandas`

## Installation
Installation is only possible on Unix systems. In order to install the package follow these steps:
Typical install time is around 1 minute.

1. Clone or download the package from the GitHub repository 

```
git clone https://github.com/BioComputingUP/CAID.git
```

2. CAID relies on `vectorized_cls_metrics` library. Clone or download the package from the GitHub repository

```
git clone https://github.com/marnec/vectorized_cls_metrics
```

3. Add the `vectorized_cls_metrics` library to the PYTHONPATH environmental variable:

```
export PYTHONPATH="${PYTHONPATH}:/path/where/the/library/was/cloned"
```

The library is successfully installed. In order to be able to copy-paste commands without the need of customize paths
the CAID package should be placed in this folder structure:

```
CAID-root
├── baseline
├── caid                --> (CAID repository)
├── data
│   ├── annotations
│   ├── predictions
│   │   ├── binding
│   │   └── disorder
│   ├── references
│   │   ├── binding
│   │   └── disorder
├── plots
└── results
```
  

## Demo
Once installed you can test if everything works fine by launching the `./demo.sh`. 
In around 1:30 minutes it should produce the following list of files in the `caid/demo/demo_output` folder:

```
D018_ESpritz-D.rawscores.distribution.txt
D018_ESpritz-D.thresholds.distribution.txt
D019_ESpritz-N.rawscores.distribution.txt
D019_ESpritz-N.thresholds.distribution.txt
D020_ESpritz-X.rawscores.distribution.txt
D020_ESpritz-X.thresholds.distribution.txt
demo-reference.analysis.all.bootstrap.bac.metrics.csv
demo-reference.analysis.all.bootstrap.csi.metrics.csv
demo-reference.analysis.all.bootstrap.default.metrics.csv
demo-reference.analysis.all.bootstrap.f1s.metrics.csv
demo-reference.analysis.all.bootstrap.f2s.metrics.csv
demo-reference.analysis.all.bootstrap.f05.metrics.csv
demo-reference.analysis.all.bootstrap.fnr.metrics.csv
demo-reference.analysis.all.bootstrap.fom.metrics.csv
demo-reference.analysis.all.bootstrap.fpr.metrics.csv
demo-reference.analysis.all.bootstrap.inf.metrics.csv
demo-reference.analysis.all.bootstrap.mcc.metrics.csv
demo-reference.analysis.all.bootstrap.mk.metrics.csv
demo-reference.analysis.all.bootstrap.npv.metrics.csv
demo-reference.analysis.all.bootstrap.ppv.metrics.csv
demo-reference.analysis.all.bootstrap.tnr.metrics.csv
demo-reference.analysis.all.bootstrap.tpr.metrics.csv
demo-reference.analysis.all.ci.bac.metrics.csv
demo-reference.analysis.all.ci.csi.metrics.csv
demo-reference.analysis.all.ci.default.metrics.csv
demo-reference.analysis.all.ci.f1s.metrics.csv
demo-reference.analysis.all.ci.f2s.metrics.csv
demo-reference.analysis.all.ci.f05.metrics.csv
demo-reference.analysis.all.ci.fnr.metrics.csv
demo-reference.analysis.all.ci.fom.metrics.csv
demo-reference.analysis.all.ci.fpr.metrics.csv
demo-reference.analysis.all.ci.inf.metrics.csv
demo-reference.analysis.all.ci.mcc.metrics.csv
demo-reference.analysis.all.ci.mk.metrics.csv
demo-reference.analysis.all.ci.npv.metrics.csv
demo-reference.analysis.all.ci.ppv.metrics.csv
demo-reference.analysis.all.ci.tnr.metrics.csv
demo-reference.analysis.all.ci.tpr.metrics.csv
demo-reference.analysis.all.dataset._.cmat.csv
demo-reference.analysis.all.dataset._.pr.csv
demo-reference.analysis.all.dataset._.predictions.csv
demo-reference.analysis.all.dataset._.roc.csv
demo-reference.analysis.all.dataset.bac.cmat.csv
demo-reference.analysis.all.dataset.bac.metrics.csv
demo-reference.analysis.all.dataset.csi.cmat.csv
demo-reference.analysis.all.dataset.csi.metrics.csv
demo-reference.analysis.all.dataset.default.cmat.csv
demo-reference.analysis.all.dataset.default.metrics.csv
demo-reference.analysis.all.dataset.f1s.cmat.csv
demo-reference.analysis.all.dataset.f1s.metrics.csv
demo-reference.analysis.all.dataset.f2s.cmat.csv
demo-reference.analysis.all.dataset.f2s.metrics.csv
demo-reference.analysis.all.dataset.f05.cmat.csv
demo-reference.analysis.all.dataset.f05.metrics.csv
demo-reference.analysis.all.dataset.fnr.cmat.csv
demo-reference.analysis.all.dataset.fnr.metrics.csv
demo-reference.analysis.all.dataset.fom.cmat.csv
demo-reference.analysis.all.dataset.fom.metrics.csv
demo-reference.analysis.all.dataset.fpr.cmat.csv
demo-reference.analysis.all.dataset.fpr.metrics.csv
demo-reference.analysis.all.dataset.inf.cmat.csv
demo-reference.analysis.all.dataset.inf.metrics.csv
demo-reference.analysis.all.dataset.mcc.cmat.csv
demo-reference.analysis.all.dataset.mcc.metrics.csv
demo-reference.analysis.all.dataset.mk.cmat.csv
demo-reference.analysis.all.dataset.mk.metrics.csv
demo-reference.analysis.all.dataset.npv.cmat.csv
demo-reference.analysis.all.dataset.npv.metrics.csv
demo-reference.analysis.all.dataset.ppv.cmat.csv
demo-reference.analysis.all.dataset.ppv.metrics.csv
demo-reference.analysis.all.dataset.tnr.cmat.csv
demo-reference.analysis.all.dataset.tnr.metrics.csv
demo-reference.analysis.all.dataset.tpr.cmat.csv
demo-reference.analysis.all.dataset.tpr.metrics.csv
demo-reference.analysis.all.target.bac.metrics.csv
demo-reference.analysis.all.target.csi.metrics.csv
demo-reference.analysis.all.target.default.metrics.csv
demo-reference.analysis.all.target.f1s.metrics.csv
demo-reference.analysis.all.target.f2s.metrics.csv
demo-reference.analysis.all.target.f05.metrics.csv
demo-reference.analysis.all.target.fnr.metrics.csv
demo-reference.analysis.all.target.fom.metrics.csv
demo-reference.analysis.all.target.fpr.metrics.csv
demo-reference.analysis.all.target.inf.metrics.csv
demo-reference.analysis.all.target.mcc.metrics.csv
demo-reference.analysis.all.target.mk.metrics.csv
demo-reference.analysis.all.target.npv.metrics.csv
demo-reference.analysis.all.target.ppv.metrics.csv
demo-reference.analysis.all.target.tnr.metrics.csv
demo-reference.analysis.all.target.tpr.metrics.csv
demo-reference.analysis.D018_ESpritz-D.bootstrap.metrics.csv
demo-reference.analysis.D018_ESpritz-D.dataset.metrics.csv
demo-reference.analysis.D018_ESpritz-D.target.metrics.csv
demo-reference.analysis.D019_ESpritz-N.bootstrap.metrics.csv
demo-reference.analysis.D019_ESpritz-N.dataset.metrics.csv
demo-reference.analysis.D019_ESpritz-N.target.metrics.csv
demo-reference.analysis.D020_ESpritz-X.bootstrap.metrics.csv
demo-reference.analysis.D020_ESpritz-X.dataset.metrics.csv
demo-reference.analysis.D020_ESpritz-X.target.metrics.csv
```

The content of `demo-reference.analysis.all.dataset.bac.metrics.csv` should look like this:

```
,bac,csi,f05,f1s,f2s,fnr,fom,fpr,inf,mcc,mk,npv,ppv,tnr,tpr,aucroc,aucpr,aps,thr
D020_ESpritz-X,0.692,0.261,0.329,0.414,0.56,0.268,0.074,0.349,0.383,0.287,0.215,0.926,0.289,0.651,0.732,0.739,0.304,0.303,0.048
D019_ESpritz-N,0.67,0.249,0.322,0.398,0.521,0.344,0.089,0.317,0.339,0.258,0.197,0.911,0.286,0.683,0.656,0.714,0.296,0.296,0.345
D018_ESpritz-D,0.704,0.27,0.338,0.426,0.576,0.248,0.068,0.345,0.407,0.305,0.229,0.932,0.297,0.655,0.752,0.774,0.41,0.409,0.248

```

## Data
CAID revolves around data obtained from different sources to build references and baselines

### Raw data
In order to use copy-pasted commands as they are (without customizing paths), the following files should be placed in
`CAID-root/data/annotations` 

disprot-2018-11-disorder.fasta obtained from: 
https://disprot.org/api/search?release=2018_11&show_ambiguous=false&show_obsolete=false&format=fasta&namespace=structural_state&get_consensus=true

disprot-2016-10-disorder.fasta obtained from: 
https://disprot.org/api/search?release=2016_10&show_ambiguous=false&show_obsolete=false&format=fasta&namespace=structural_state&get_consensus=true

disprot-2018-11-interaction.fasta obtained from: 
https://disprot.org/api/search?release=2018_11&show_ambiguous=false&show_obsolete=false&format=fasta&namespace=interaction_partner&get_consensus=true

disprot-2018-11.json obtained from: 
https://disprot.org/api/search?release=current&show_ambiguous=false&show_obsolete=false&format=json

InterProScan (5.38-76.0) output generated with the following command:

```
interproscan.sh -f tsv -dp -iprlookup -T /tmp/ -i disprot-2018-11_seq.fasta -o caid_interproscan -T /tmp/
python parse_interproscan.py ../data/annotations/disprot-2018-11-disorder.fasta ../data/annotations/disprot-2018-11.json caid_interproscan  > ../data/annotations/data/gene3d.fasta
```

PDB annotations obtained from MobiDB; downloaded on date 06/03/2020 from: 
https://mobidb.bio.unipd.it/mobidb3_datasets/latest/derived_disorder.mjson.gz

PDB-ateast definition

```
python parse_mobidb.py ../data/annotations/disprot-2018-11.json ../data/annotations/derived_disorder.mjson.gz > ../data/annotations/pdb-atleast.fasta
```

### Predictions
In order to use copy-pasted commands as they are (without customizing paths), the following files should be placed in
`CAID-root/data/predictions/{disorder|binding}`
 
predictions obtainable from: 
https://mobidb.org/caid/1/predictions

Disorder prediction filenames start with `D` character. Binding predictions filenames start with `B` character.


## Usage
Create reference

```
python datasets/make_references.py ../data/annotations/disprot-2018-11.json -d ../data/annotations/disprot-2018-11-disorder.fasta -e ../data/annotations/disprot-2016-10-disorder.fasta -s ../data/annotations/pdb-atleast.fasta ../data/annotations/gene3d.fasta -i ../data/annotations/disprot-2018-11-interaction.fasta
```

Calculate Reference statistics

```
python reference_stats.py ../data/annotations/disprot-2018-11.json ../data/references/disorder/disprot-disorder* ../data/references/binding/disprot-binding* -o ../data/dataset_stats/
```

Calculate Evaluation metrics

```
bash launch_all.sh
```

Draw plots

```
python plots.py ../results/ ../baseline/ ../data/references/disorder/ ../data/dataset_stats/ -o ../plots/ -n data/caid_names.json -ll DEBUG -g 'disprot-disorder*'
python plots.py ../results/ ../baseline/ ../data/references/binding/ ../data/dataset_stats/ -o ../plots/ -n data/caid_names.json -ll DEBUG -g 'disprot-binding*'
```

## License
[CC BY 3.0](https://creativecommons.org/licenses/by/3.0/)