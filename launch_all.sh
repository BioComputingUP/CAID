#!/usr/bin/env bash
module load python
for fname in $(ls ../data/disorder/);
    do
    python3 evaluation.py disorder ../data/disorder/${fname} -p 'D[0-9]{1,3}' -ll INFO;
    python3 baseline_rand.py ../data/disorder/${fname} ;
    python3 baseline_cons.py ../data/disorder/${fname};
    done

