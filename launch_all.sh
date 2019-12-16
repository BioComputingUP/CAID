#!/usr/bin/env bash
for fname in $(ls ../data/references/disorder/);
    do
        pipenv run python3 caid.py disorder ../data/references/disorder/${fname} -p 'D[0-9]{1,3}' -ll INFO -o ../results;
        pipenv run python3 baseline_rand.py ../data/references/disorder/${fname} -o ../baseline;
        pipenv run python3 baseline_cons.py ../data/references/disorder/${fname} -o ../baseline;
        pipenv run python3 baseline_naive.py ../data/references/disorder/${fname} -o ../baseline -p ../data/references/disorder/new-pdb-r_simple.txt
        pipenv run python3 baseline_naive.py ../data/references/disorder/${fname} -o ../baseline -p ../data/references/disorder/new-gene3d-r_simple.txt
    done

#for fname in $(ls ../data/binding/);
#    do
#		pipenv run python3 caid.py binding ../data/binding/${fname} -p 'B[0-9]{1,3}' -ll INFO -o ../results;
#        pipenv run python3 baseline_rand.py ../data/binding/${fname} -o ../baseline;
#        pipenv run python3 baseline_cons.py ../data/disorder/${fname} -o ../baseline;
#        pipenv run python3 baseline_naive.py ../data/disorder/${fname} -o ../baseline -p ../data/disorder/new-disprot-all_gene3d.txt
#        pipenv run python3 baseline_naive.py ../data/disorder/${fname} -o ../baseline -p ../data/disorder/new-disprot-all_pdb.txt
#    done

