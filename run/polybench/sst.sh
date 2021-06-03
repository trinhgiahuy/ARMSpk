#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
source $ROOTDIR/dep/spack/share/spack/setup-env.sh
spack load gcc@8.4.0

export PIN_HOME=$ROOTDIR/dep/sst/pin-3.17-98314-g0c048d619-gcc-linux
export INTEL_PIN_DIRECTORY=$PIN_HOME
export SST_CORE_HOME=$ROOTDIR/dep/sst
export PATH=$SST_CORE_HOME/bin:$PATH
export LD_LIBRARY_PATH=$SST_CORE_HOME/lib:$LD_LIBRARY_PATH
export PYTHONPATH=$ROOTDIR/dep/sst:$PYTHONPATH
SIMTHREADS=1	# Gwen said that 2-4 is ideal for CPU sims
SST="sst -n $SIMTHREADS --verbose $ROOTDIR/dep/sst/bdw.py"

# ============================ PolyBench ====================================
source conf/polybench.sh
DEFLOG="$ROOTDIR/log/`hostname -s`/sstrun/polybench"
mkdir -p $DEFLOG
cd $APPDIR
for BEST in $BESTCONF; do
	for BMconf in $BINARYS; do
		NumMPI=1
		NumOMP=1
		BINARY="`echo ${BMconf} | tr '|' '_'`"
		BName="`basename $(echo ${BMconf} | cut -d'|' -f1)`"
		LOG="${DEFLOG}/${BName}.log"; mkdir -p ${LOG}_stat
		cp $ROOTDIR/dep/sst/bdwariel.cfg.in ${LOG}_stat/ariel.cfg
		sed -i -e "s#^executable:.*#executable: $BINARY#" ${LOG}_stat/ariel.cfg
		ARGC=0
		for ARG in $INPUT; do
			sed -i -e "/^executable:/a apparg$ARGC: $ARG" ${LOG}_stat/ariel.cfg
			ARGC=$((ARGC+1))
		done
		if [ $ARGC -ge 1 ]; then sed -i -e "/^executable:/a appargcount: $ARGC" ${LOG}_stat/ariel.cfg; fi
		echo "$SST --model-options=\"--verbose --config ${LOG}_stat/ariel.cfg --statfile ${LOG}_stat/stats.csv\"" >> $LOG 2>&1
		$SST --model-options="--verbose --config ${LOG}_stat/ariel.cfg --statfile ${LOG}_stat/stats.csv" >> $LOG 2>&1 &
	done
	wait
done
cd $ROOTDIR
