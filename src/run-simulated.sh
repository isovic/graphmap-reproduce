#! /bin/sh

currentpath=$PWD/../
### Create an output folder for the collected results.
mkdir -p $currentpath/results

### Install the Aligneval repository.
mkdir -p $currentpath/tools
cd $currentpath/tools
git clone https://github.com/isovic/aligneval.git
cd aligneval
./setup.py all

### Run all alignments.
./run-alignment.py
./run-alignment-graphmap_stages.py
./run-alignment-gridsearch.py

### Run all evaluations.
./run-evaluation.py v1
./run-evaluation.py v2
./run-evaluation.py v3

### Copy the relevant results to the collective results folder.
### Classic alignment results for all mappers:
cp $currentpath/tools/aligneval/results/precision_recall-50bp_distance-v1.csv $currentpath/results/
cp $currentpath/tools/aligneval/results/precision_recall-strict-correct_bases-v1.csv $currentpath/results/
cp $currentpath/tools/aligneval/results/memory-v1.csv $currentpath/results/
cp $currentpath/tools/aligneval/results/cpu_time-v1.csv $currentpath/results/

### Grid search results (running GraphMap on a range of read lengths and error rates):
cp $currentpath/tools/aligneval/results/precision_recall-50bp_distance-v2.csv $currentpath/results/
cp $currentpath/tools/aligneval/results/precision_recall-strict-correct_bases-v2.csv $currentpath/results/
cp $currentpath/tools/aligneval/results/memory-v2.csv $currentpath/results/
cp $currentpath/tools/aligneval/results/cpu_time-v2.csv $currentpath/results/

### GraphMap stages - evaluation of the influence of each algorithmic step of GraphMap.
cp $currentpath/tools/aligneval/results/precision_recall-50bp_distance-v3.csv $currentpath/results/
cp $currentpath/tools/aligneval/results/precision_recall-strict-correct_bases-v3.csv $currentpath/results/
cp $currentpath/tools/aligneval/results/memory-v3.csv $currentpath/results/
cp $currentpath/tools/aligneval/results/cpu_time-v3.csv $currentpath/results/
