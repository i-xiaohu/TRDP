# TRDP
This project explores the model of sequence alignment with tandem duplications, or DSI (Duplications, Substitutions, Indels) model. 
Tandem repeats are enriched in complex genomes, highly variable, and confuse downstream analysis.
Incorporating tandem copy events into sequence model will be instrumental for genomic studies.

`TRDP` can now identify tandem repeats within a sequence through self-alignment dynamic programming. 

Run the command to compile it.
```
mkdir build; cd build
cmake ..; make -j8
```

Its usage is described below.
```
Usage: TRDP [options] seq.fa
  Self-Alignment Options:
    -A [INT]  match score [1]
    -B [INT]  mismatch penalty [-4]
    -O [INT]  open gap(indel) penalty [-6]
    -E [INT]  extend gap penalty [-1]
  Duplication Alignment Options:
    -u [INT]  minimum repeat unit size [5]
    -d [INT]  open duplication penalty [-2]
    -p [INT]  close duplication penalty [-6]
    -a [INT]  match score [2]
    -b [INT]  mismatch penalty [-3]
    -o [INT]  open gap(indel) penalty [-3]
    -e [INT]  extend gap penalty [-1]
Note: duplication scoring matrix must reward more and/or
  penalize less than self-alignment matrix to drive out duplication events.
```