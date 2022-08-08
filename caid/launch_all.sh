#!/usr/bin/env bash

for fname in ../data/references/*; do
    echo ${fname}
     python3 caid.py  ${fname} ../data/pred2 -o ../results2
#    python3 caid.py ${fname} ../data/predictions -o ../results -ll "WARNING"
#    python3 baseline_naive.py ${fname} -o ../baseline --pRef ../data/references/disorder/pdb-atleast-reverse.txt -ll INFO;
#    python3 baseline_naive.py ${fname} -o ../baseline --pRef ../data/references/disorder/gene3d-reverse.txt -ll INFO;
#    python3 baseline_rand.py ${fname} -o ../baseline -ll INFO
#    python3 baseline_cons.py ${fname} ../data/pssm -o ../baseline -ll INFO;
done

# for fname in ../data/references/binding/disprot-binding*;
# do
# echo ${fname}
#python3 caid.py ${fname} ../data/predictions/binding -ll INFO -o ../results -ll INFO;
#python3 baseline_rand.py ${fname} -o ../baseline -ll INFO;
#python3 baseline_cons.py ${fname} ../data/pssm -o ../baseline -ll INFO;
#python3 baseline_naive.py ${fname} -o ../baseline --pRef ../data/references/disorder/pdb-atleast-reverse.txt -ll INFO;
#python3 baseline_naive.py ${fname} -o ../baseline --pRef ../data/references/disorder/gene3d-reverse.txt -ll INFO;
# done
