#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
ulimit -s unlimited
ulimit -n 4096

# ============================ HPCG ===========================================
source conf/hpcg.sh
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/`hostname -s`/gem5run/hpcg/conf${1}.log"
mkdir -p ${LOG}_stat
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI=1
	#XXX: npx|y|z always 1 for 1==NumMPI: test to identify hpcg's internal dimensions
	X=$MAXXYZ; Y=$MAXXYZ; Z=$MAXXYZ
	#rm -f hpcg_log_* n*.yaml
	#rm -f hpcg20*.txt HPCG-Benchmark_3*.txt
	#mpiexec $MPIEXECOPT -n $NumMPI $BINARY -n 1 > /dev/null 2>&1
	#if [ ! "x$?" = "x0" ]; then continue; fi
	#if [ -f n*.yaml ]; then
	#	X=$(($MAXXYZ / `grep 'npx:' n*.yaml | awk -F 'npx:' '{print $2}'`))
	#	Y=$(($MAXXYZ / `grep 'npy:' n*.yaml | awk -F 'npy:' '{print $2}'`))
	#	Z=$(($MAXXYZ / `grep 'npz:' n*.yaml | awk -F 'npz:' '{print $2}'`))
	#elif [ -f HPCG-Benchmark_3*.txt ]; then
	#	#non-intel version needs to be div8 https://github.com/hpcg-benchmark/hpcg/issues/47
	#	X=$((($MAXXYZ / `grep 'npx=' HPCG-Benchmark_3*.txt | awk -F 'npx=' '{print $2}'` / 8) * 8))
	#	Y=$((($MAXXYZ / `grep 'npy=' HPCG-Benchmark_3*.txt | awk -F 'npy=' '{print $2}'` / 8) * 8))
	#	Z=$((($MAXXYZ / `grep 'npz=' HPCG-Benchmark_3*.txt | awk -F 'npz=' '{print $2}'` / 8) * 8))
	#else continue; fi
	rm -f hpcg_log_* n*.yaml
	rm -f hpcg20*.txt HPCG-Benchmark_3*.txt
	INPUT="`echo $DEFINPUT | sed -e \"s/NX/$X/\" -e \"s/NY/$Y/\" -e \"s/NZ/$Z/\"`"
	run_gem5_cmd "$@"
	cat hpcg_log_* >> $LOG 2>&1
	cat n*.yaml >> $LOG 2>&1
	rm -f hpcg_log_* n*.yaml
	cat hpcg20*.txt >> $LOG 2>&1
	cat HPCG-Benchmark_3*.txt >> $LOG 2>&1
	rm -f hpcg20*.txt HPCG-Benchmark_3*.txt
done
cd $ROOTDIR
