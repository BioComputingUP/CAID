#!/usr/bin/env bash

# Execute it from inside the "code" folder

# Update this path
export PYTHONPATH="${PYTHONPATH}:/home/damiano/Projects/vectorized_cls_metrics/"

for fname in ../data/references/disorder/disprot-disorder* ;
   do
     echo ${fname}
     python3 baseline_rand.py ${fname} -o ../data/results -ll INFO
     python3 baseline_cons.py ${fname} ../data/pssm -o ../data/results -ll INFO;
     python3 baseline_naive.py ${fname} -o ../data/results --pRef ../data/references/disorder/pdb-atleast-reverse.txt -ll INFO;
     python3 baseline_naive.py ${fname} -o ../data/results --pRef ../data/references/disorder/gene3d-reverse.txt -ll INFO;
     python3 baseline_naive_pdb_homology.py ${fname} ../data/annotations/caid_646_seqres.blast_out -o ../data/results -n 20 -x 30 -ll INFO
     python3 baseline_naive_pdb_homology.py ${fname} ../data/annotations/caid_646_seqres.blast_out -o ../data/results -n 30 -x 100 -t 2018-11-30 -ll INFO
     python3 caid.py ${fname} ../data/predictions/disorder -ll INFO -o ../data/results
   done

for fname in ../data/references/binding/disprot-binding-all.txt*;
    do
      echo ${fname}
      python3 baseline_rand.py ${fname} -o ../data/results -ll INFO;
      python3 baseline_cons.py ${fname} ../data/pssm -o ../data/results -ll INFO;
      python3 baseline_naive.py ${fname} -o ../data/results --pRef ../data/references/disorder/pdb-atleast-reverse.txt -ll INFO;
      python3 baseline_naive.py ${fname} -o ../data/results --pRef ../data/references/disorder/gene3d-reverse.txt -ll INFO;
      python3 baseline_naive_pdb_homology.py ${fname} ../data/annotations/caid_646_seqres.blast_out -o ../data/results -n 20 -x 30 -ll INFO
      python3 baseline_naive_pdb_homology.py ${fname} ../data/annotations/caid_646_seqres.blast_out -o ../data/results -n 30 -x 100 -t 2018-11-30 -ll INFO
      python3 caid.py ${fname} ../data/predictions/binding -ll INFO -o ../data/results
    done
