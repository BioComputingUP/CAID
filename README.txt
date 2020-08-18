CAID relies on vectorized_cls_metrics library.
At the moment this library needs to be cloned and manually added to the PYTHONPATH environmental variable. I plan to
publish it in pip in the future.

disprot-2018-11-disorder.fasta obtained from:
https://disprot.org/api/search?release=2018_11&show_ambiguous=false&show_obsolete=false&format=fasta&namespace=structural_state&get_consensus=true

disprot-2016-10-disorder.fasta obtained from:
https://disprot.org/api/search?release=2016_10&show_ambiguous=false&show_obsolete=false&format=fasta&namespace=structural_state&get_consensus=true

disprot-2018-11-interaction.fasta obtained from:
https://disprot.org/api/search?release=2018_11&show_ambiguous=false&show_obsolete=false&format=fasta&namespace=interaction_partner&get_consensus=true

disprot-2018-11.json obtained from:
https://disprot.org/api/search?release=current&show_ambiguous=false&show_obsolete=false&format=json

InterProScan (5.38-76.0) output generated with the following command:
interproscan.sh -f tsv -dp -iprlookup -T /tmp/ -i disprot-2018-11_seq.fasta -o caid_interproscan -T /tmp/

PDB annotations obtained from MobiDB; downloaded on date 06/03/2020 from:
https://mobidb.bio.unipd.it/mobidb3_datasets/latest/derived_disorder.mjson.gz

references were created with command
python datasets/make_references.py ../data/annotations/disprot-2018-11.json -d ../data/annotations/disprot-2018-11-disorder.fasta -e ../data/annotations/disprot-2016-10-disorder.fasta -s ../data/annotations/pdb-atleast.fasta ../data/annotations/gene3d.fasta -i ../data/annotations/disprot-2018-11-interaction.fasta

reference statistics calculated with command
python reference_stats.py ../data/annotations/disprot-2018-11.json ../data/references/disorder/disprot-disorder* ../data/references/binding/disprot-binding* -o ../data/dataset_stats/

results created with:
bash launch_all.sh

plots created with
python plots.py ../results/ ../baseline/ ../data/references/disorder/ ../data/dataset_stats/ -o ../plots/ -n data/caid_names.json -ll DEBUG -g 'disprot-disorder*'
python plots.py ../results/ ../baseline/ ../data/references/binding/ ../data/dataset_stats/ -o ../plots/ -n data/caid_names.json -ll DEBUG -g 'disprot-binding*'

REFERENCES
In the text we refer to proteins as targets, to disordered residues as positive labels and structured/ordered residues as negative labels. In CAID different references were built. References differ in the subset of DisProt used to define positive labels and in the definition of negatives labels. To define negatives we adopted two strategies that we called Simple and PDB.

Disorder Positive definition
To define positive labels DisProt annotation is used. All residues annotated as disordered in DisProt are considered positive labels in the reference, regardless of the experiment annotated in DisProt. Overlapping regions of annotation are merged in a single region and all residues are considered equally
regardless of the number of evidences they are covered by.

Disorder Negative definition: Simple
Simple refers to how negative labels are considered when building a reference. When using the Simple strategy DisProt is used to define positive labels for each target. All labels that are not positive are considered negative. In other words, when using the Simple strategy DisProt annotation defines regions of disorder while those parts of the sequence that are not covered by DisProt annotation are considered ordered.

Disorder Negative definition: PDB
Simple refers to how negative labels are considered when building a reference. When using the PDB strategy DisProt is used to define positive labels for each target while PDB is used to define negative labels. In other words, when using the PDB strategy DisProt annotation defines regions of disorder while PDB structures mapping on the proteins sequence define regions of order. When PDB and DisProt annotation cover the same stretch of sequence, the region is defined as disordered.

Binding Positive definition
To define positive labels DisProt annotation is used. All residues annotated in DisProt as binding region are considered positive labels in the reference, regardless of the experiment annotated in DisProt. Overlapping regions of annotation are merged in a single region and all residues are considered equally
regardless of the number of evidence they are covered by.

Binding Negative definition: Simple
Simple refers to how negative labels are considered when building a reference. When using the Simple strategy DisProt is used to define positive labels for each target. All labels that are not positive are considered negative. In other words, when using the Simple strategy DisProt annotation defines binding regions while those parts of the sequence that are not covered by DisProt annotation are considered outside of a binding region.

BASELINES
A number of baseline predictors have been built in order to be compared with actual predictors. Two are based on the reshuffling of a reference, one is completely random, one classifies labels randomly but with a fixed proportion of positive to negative labels, one is based on an estimate of residue conservation through evolution, two are based on taking the reverse of the annotation of structure as positive labels.
Shuffled baselines
Shuffled baselines are built by reshuffling a reference. A reference is reshuffled at the dataset level (shuffle-dataset) or at the target level (shuffle-target). The former is a random classification that preserves the proportion of positive to negative labels in the dataset. The latter is a random classification that preserves the proportion of positive to negative labels in each target (protein).

Shuffled baselines don't have prediction scores since they are just a reshuffling of the reference. Their prediction states are however treated as scores for the sake of computation. Default threshold is therefore set to 1.0 since thresholding is always applied following a greater-equal strategy (see threshold paragraph).

Random baselines
Random baselines are random classifiers in which the prediction score of each label is assigned randomly. Random baselines are built by randomly drawing floating point numbers out of a uniform distribution [0, 1]. Two of these baselines have been built: one that is completely random and one that has a fixed proportion of positive labels. Both assign scores from the same distribution. The former is thresholded as 0.5, the latter is thresholded at 1-fp, where fp is the fixed proportion of positive labels specified when building the baseline. Our fixed-positive-fraction baseline has fp equal to 0.347, which is the fraction of disordered amino acid residues in DisProt 7.0.
Conservation-based baseline
Conservation-based baseline is built by calculating the distance between the PSSM of each target in DisProt and the residue distribution of the BLOSUM62 substitution matrix. Target PSSMs are obtained by launching a 3 cycle PsiBLAST against UniRef90. The distance is calculated as the Jensens-Shannon divergence [https://doi.org/10.1093/bioinformatics/btm270] of the two frequencies.
Naive baselines
Naive baselines are based on the naive assumption that whatever is not annotated as structure is disorder. Two naive baselines are following this principle, one has the structure annotation defined by PDB structures, the other has structure annotations defined by Gene3d predictions. The former uses PDB data mapped on UniProt sequences by Mobi 2.0 [10.1093/bioinformatics/btx592], the latter uses Gene3D predictions generated by InterProScan.