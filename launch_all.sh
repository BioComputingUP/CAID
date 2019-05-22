#!/usr/bin/env bash
for fname in $(ls ../data/disorder/);
    do
#        pipenv run python3 evaluation.py disorder ../data/disorder/${fname} -p 'D[0-9]{1,3}' -ll INFO -o ../results;
#        pipenv run python3 baseline_rand.py ../data/disorder/${fname} -o ../baseline;
#        pipenv run python3 baseline_cons.py ../data/disorder/${fname} -o ../baseline;
#        pipenv run python3 baseline_naive.py ../data/disorder/${fname} -o ../baseline -p ../data/disorder/new-disprot-all_gene3d.txt
#        pipenv run python3 baseline_naive.py ../data/disorder/${fname} -o ../baseline -p ../data/disorder/new-disprot-all_pdb.txt
    done

for fname in $(ls ../data/binding/);
    do
#        pipenv run python3 evaluation.py binding ../data/binding/${fname} -p 'B[0-9]{1,3}' -ll INFO -o ../results;
        pipenv run python3 baseline_rand.py ../data/binding/${fname} -o ../baseline;
#        pipenv run python3 baseline_cons.py ../data/disorder/${fname} -o ../baseline;
#        pipenv run python3 baseline_naive.py ../data/disorder/${fname} -o ../baseline -p ../data/disorder/new-disprot-all_gene3d.txt
#        pipenv run python3 baseline_naive.py ../data/disorder/${fname} -o ../baseline -p ../data/disorder/new-disprot-all_pdb.txt
    done

