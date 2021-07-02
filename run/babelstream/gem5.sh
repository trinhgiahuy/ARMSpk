#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
ulimit -s unlimited
ulimit -n 4096

# ============================ BabelStream ====================================
source conf/babelstream.sh
DEFLOG="$ROOTDIR/log/`hostname -s`/gem5run/babelstream"
mkdir -p $DEFLOG
cd $APPDIR
DEFINPUT=$INPUT
for BEST in $BESTCONF; do
	for BINARY in $BINARYS; do
		NumMPI=1
		S="`echo $BINARY | cut -d '_' -f2`"
		BINARY="`echo $BINARY | cut -d '_' -f1`"
		LOG="${DEFLOG}/conf${1}.${S}gb.log"; mkdir -p ${LOG}_stat
		S=$((S*1024*1024*1024/8))
		run_gem5_cmd "$@"
	done
done
cd $ROOTDIR
