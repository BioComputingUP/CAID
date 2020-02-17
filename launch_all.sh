#!/usr/bin/env bash
for fname in $(ls ../data/disorder/);
    do
#        pipenv run python3 caid.py ../data/disorder/${fname} ../data/results/disorder -ll INFO -o ../results
        pipenv run python3 baseline_rand.py ../data/disorder/${fname} -o ../baseline -ll DEBUG;
#        pipenv run python3 baseline_cons.py ../data/disorder/${fname} ../data/pssm -o ../baseline -ll DEBUG;
#        pipenv run python3 baseline_naive.py ../data/disorder/${fname} -o ../baseline --pRef ../data/disorder/new-pdb-r_simple.txt -ll DEBUG;
#        pipenv run python3 baseline_naive.py ../data/disorder/${fname} -o ../baseline --pRef ../data/disorder/new-gene3d-r_simple.txt -ll DEBUG;
    done

#for fname in $(ls ../data/binding/);
#    do
#      echo ../data/binding/${fname}
# 		  pipenv run python3 caid.py ../data/binding/${fname} ../data/results/binding -ll INFO -o ../results -ll DEBUG;
#      pipenv run python3 baseline_rand.py ../data/binding/${fname} -o ../baseline -ll DEBUG;
#      pipenv run python3 baseline_cons.py ../data/binding/${fname} -o ../data/pssm ../baseline -ll DEBUG;
#      pipenv run python3 baseline_naive.py ../data/binding/${fname} -o ../baseline --pRef ../data/disorder/new-disprot-all_gene3d.txt -ll DEBUG;
#      pipenv run python3 baseline_naive.py ../data/binding/${fname} -o ../baseline --pRef ../data/disorder/new-disprot-all_pdb.txt -ll DEBUG;
#    done
