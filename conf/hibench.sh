#!/bin/bash

export APPDIR="./HiBench"
#export BINARY="./bin/run_all.sh"
export BINARYS="bin/workloads/graph/nweight/spark/run.sh
bin/workloads/micro/repartition/spark/run.sh
bin/workloads/micro/sleep/spark/run.sh
bin/workloads/micro/sort/spark/run.sh
bin/workloads/micro/terasort/spark/run.sh
bin/workloads/micro/wordcount/spark/run.sh
bin/workloads/ml/als/spark/run.sh
bin/workloads/ml/bayes/spark/run.sh
bin/workloads/ml/gbt/spark/run.sh
bin/workloads/ml/kmeans/spark/run.sh
bin/workloads/ml/lda/spark/run.sh
bin/workloads/ml/linear/spark/run.sh
bin/workloads/ml/lr/spark/run.sh
bin/workloads/ml/pca/spark/run.sh
bin/workloads/ml/rf/spark/run.sh
bin/workloads/ml/svd/spark/run.sh
bin/workloads/ml/svm/spark/run.sh
bin/workloads/sql/aggregation/spark/run.sh
bin/workloads/sql/join/spark/run.sh
bin/workloads/sql/scan/spark/run.sh
bin/workloads/websearch/pagerank/spark/run.sh"
export INPUT=""
export NumRunsTEST=3
export NumRunsBEST=10
export MAXTIME="10m"
export RUNSDE="yes"
export RUNPCM="no"
export RUNVTUNE="no"

if [ -n "${XEONHOST}" ]; then
	# on "normal" Xeon
	export TESTCONF="1|1"
	export BESTCONF="1|1"
	export SCALCONF=""
elif [ -n "${IKNLHOST}" ]; then
	# on one of the Phi (knl)
	export TESTCONF="1|1"
	export BESTCONF="1|1"
elif [ -n "${IKNMHOST}" ]; then
	# on one of the Phi (knm)
	export TESTCONF="1|1"
	export BESTCONF="1|1"
elif [ -n "${FUJIHOST}" ] || [ -n "${RFX7HOST}" ]; then
	export TESTCONF="1|1"
	export BESTCONF="1|1"
	export SCALCONF=""
fi
