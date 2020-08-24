#!/bin/bash
EXAMPLES_DIR=examples

##Run example1
DATASET=example1
DATA_DIR=./examples/example1/data/
RESULTS_DIR=./examples/example1/results

##Run example2
#DATASET=example2
#DATA_DIR=./examples/example2/data/
#RESULTS_DIR=./examples/example2/results

echo
echo Generating kernel matrices or vector of counts for each kernel method ...

### Vertex classification tasks

#Running Cumulative Random Hyperwalk Kernel 
./run_hypergraphkernel -p ${EXAMPLES_DIR}/${DATASET}/${DATASET}.pos -n ${EXAMPLES_DIR}/${DATASET}/${DATASET}.neg -g ${DATA_DIR} -l ${DATA_DIR} -e ${DATA_DIR} -t 0 -I 1000 -R 0.4 -k ${RESULTS_DIR}/${DATASET}_chrw.dat -v 

#Running Random Hyperwalk Kernel
./run_hypergraphkernel -p ${EXAMPLES_DIR}/${DATASET}/${DATASET}.pos -n ${EXAMPLES_DIR}/${DATASET}/${DATASET}.neg -g ${DATA_DIR} -l ${DATA_DIR} -e ${DATA_DIR} -t 1 -I 1000 -R 0.2 -k ${RESULTS_DIR}/${DATASET}_hrw.dat -v

#Running Standard Hypergraphlet Kernel
./run_hypergraphkernel -p ${EXAMPLES_DIR}/${DATASET}/${DATASET}.pos -n ${EXAMPLES_DIR}/${DATASET}/${DATASET}.neg -g ${DATA_DIR} -l ${DATA_DIR} -e ${DATA_DIR} -t 2 -k ${RESULTS_DIR}/${DATASET}_shgk.dat #Kernel Matrix

./run_hypergraphkernel -p ${EXAMPLES_DIR}/${DATASET}/${DATASET}.pos -n ${EXAMPLES_DIR}/${DATASET}/${DATASET}.neg -g ${DATA_DIR} -l ${DATA_DIR} -e ${DATA_DIR} -t 2 -s ${RESULTS_DIR}/${DATASET}_shgk.svml -v #SVML

#Running Label Substitutions Hypergraphlet Kernel
./run_hypergraphkernel -p ${EXAMPLES_DIR}/${DATASET}/${DATASET}.pos -n ${EXAMPLES_DIR}/${DATASET}/${DATASET}.neg -g ${DATA_DIR} -l ${DATA_DIR} -e ${DATA_DIR} -t 3 -K 0.34 -M 1 -V AC -E AB -k ${RESULTS_DIR}/${DATASET}_lshk.dat

./run_hypergraphkernel -p ${EXAMPLES_DIR}/${DATASET}/${DATASET}.pos -n ${EXAMPLES_DIR}/${DATASET}/${DATASET}.neg -g ${DATA_DIR} -l ${DATA_DIR} -e ${DATA_DIR} -t 3 -K 0.34 -M 1 -V AC -E AB -s ${RESULTS_DIR}/${DATASET}_lshk.svml -v

#Running Edge Indels Hypergraphlet Kernel
./run_hypergraphkernel -p ${EXAMPLES_DIR}/${DATASET}/${DATASET}.pos -n ${EXAMPLES_DIR}/${DATASET}/${DATASET}.neg -g ${DATA_DIR} -l ${DATA_DIR} -e ${DATA_DIR} -t 4 -m 1 -E AB -k ${RESULTS_DIR}/${DATASET}_eihk.dat

./run_hypergraphkernel -p ${EXAMPLES_DIR}/${DATASET}/${DATASET}.pos -n ${EXAMPLES_DIR}/${DATASET}/${DATASET}.neg -g ${DATA_DIR} -l ${DATA_DIR} -e ${DATA_DIR} -t 4 -m 1 -E AB -s ${RESULTS_DIR}/${DATASET}_eihk.svml -v

#Running Edit Distance Hypergraphlet Kernel
./run_hypergraphkernel -p ${EXAMPLES_DIR}/${DATASET}/${DATASET}.pos -n ${EXAMPLES_DIR}/${DATASET}/${DATASET}.neg -g ${DATA_DIR} -l ${DATA_DIR} -e ${DATA_DIR} -t 5 -m 1 -M 1 -K 0.34 -V AC -E AB -k ${RESULTS_DIR}/${DATASET}_edhk.dat

./run_hypergraphkernel -p ${EXAMPLES_DIR}/${DATASET}/${DATASET}.pos -n ${EXAMPLES_DIR}/${DATASET}/${DATASET}.neg -g ${DATA_DIR} -l ${DATA_DIR} -e ${DATA_DIR} -t 5 -m 1 -M 1 -K 0.34 -V AC -E AB -s ${RESULTS_DIR}/${DATASET}_edhk.svml -v


### Edge classification tasks

#Running Standard Hypergraphlet Kernel
./run_hypergraphkernel -p ${EXAMPLES_DIR}/${DATASET}/${DATASET}.pos -n ${EXAMPLES_DIR}/${DATASET}/${DATASET}.neg -g ${DATA_DIR} -l ${DATA_DIR} -e ${DATA_DIR} -t 2 -z 1 -s ${RESULTS_DIR}/${DATASET}_dual_shgk.svml -v #SVML

echo

date

