#!/usr/bin/env bash
for fname in ../data/references/disorder/disprot*;
    do
      echo ${fname}
      pipenv run python3 caid.py ${fname} ../data/predictions/disorder -ll INFO -o ../results
      pipenv run python3 baseline_naive.py ${fname} -o ../baseline --pRef ../data/references/disorder/pdb-atleast-reverse.txt -ll INFO;
      pipenv run python3 baseline_naive.py ${fname} -o ../baseline --pRef ../data/references/disorder/gene3d-reverse.txt -ll INFO;
      pipenv run python3 baseline_rand.py ${fname} -o ../baseline -ll INFO
      pipenv run python3 baseline_cons.py ${fname} ../data/pssm -o ../baseline -ll INFO;
    done

for fname in ../data/references/binding/disprot-binding*;
    do
      echo ${fname}
 		  pipenv run python3 caid.py ${fname} ../data/predictions/binding -ll INFO -o ../results -ll INFO;
      pipenv run python3 baseline_rand.py ${fname} -o ../baseline -ll INFO;
      pipenv run python3 baseline_cons.py ${fname} ../data/pssm -o ../baseline -ll INFO;
      pipenv run python3 baseline_naive.py ${fname} -o ../baseline --pRef ../data/references/disorder/pdb-atleast-reverse.txt -ll INFO;
      pipenv run python3 baseline_naive.py ${fname} -o ../baseline --pRef ../data/references/disorder/gene3d-reverse.txt -ll INFO;
    done
