## Introduction
The CAID software produces all outputs necessary for a CAID edition, including baselines, references, metrics and plots, 
starting from predictions and a reference (see Data Availability section to know how to obtain this data). 
CAID software packages wraps the vectorized_cls_metrics repository ([https://github.com/marnec/vectorized_cls_metrics]), 
which performs the calculations of the classification metrics used throughout CAID. More information at 
[http://disprotcentral.org/caid]


## Index
* Requirements
* Installation
* Data
* Usage
* License

## Requirements
### Interpreter
Python 3.6+

### Dependencies
* `vectorized_metrics` ([https://github.com/marnec/vectorized_cls_metrics])
* `numpy`
* `matplotlib`
* `seaborn`
* `scipy`
* `pandas`

## Installation
Installation is only possible on Unix systems. In order to install the package follow these steps:

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
  

## Data
CAID revolves around data obtained from different sources to build references and baselines

### Raw data
In order to use copy-pasted commands as they are (without customizing paths), the following files should be placed in
`CAID-root/data/annotations` 

disprot-2018-11-disorder.fasta obtained from: 
[https://disprot.org/api/search?release=2018_11&show_ambiguous=false&show_obsolete=false&format=fasta&namespace=structural_state&get_consensus=true]

disprot-2016-10-disorder.fasta obtained from: 
[https://disprot.org/api/search?release=2016_10&show_ambiguous=false&show_obsolete=false&format=fasta&namespace=structural_state&get_consensus=true]

disprot-2018-11-interaction.fasta obtained from: 
[https://disprot.org/api/search?release=2018_11&show_ambiguous=false&show_obsolete=false&format=fasta&namespace=interaction_partner&get_consensus=true]

disprot-2018-11.json obtained from: 
[https://disprot.org/api/search?release=current&show_ambiguous=false&show_obsolete=false&format=json]

InterProScan (5.38-76.0) output generated with the following command:

```
interproscan.sh -f tsv -dp -iprlookup -T /tmp/ -i disprot-2018-11_seq.fasta -o caid_interproscan -T /tmp/
python parse_interproscan.py ../data/annotations/disprot-2018-11-disorder.fasta ../data/annotations/disprot-2018-11.json caid_interproscan  > ../data/annotations/data/gene3d.fasta
```

PDB annotations obtained from MobiDB; downloaded on date 06/03/2020 from: 
[https://mobidb.bio.unipd.it/mobidb3_datasets/latest/derived_disorder.mjson.gz]

PDB-ateast definition

```
python parse_mobidb.py ../data/annotations/disprot-2018-11.json ../data/annotations/derived_disorder.mjson.gz > ../data/annotations/pdb-atleast.fasta
```

### Predictions
In order to use copy-pasted commands as they are (without customizing paths), the following files should be placed in
`CAID-root/data/predictions/{disorder|binding}`
 
predictions obtainable from: 
[https://mobidb.org/caid/1/predictions]

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