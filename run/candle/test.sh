#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

ulimit -s unlimited
ulimit -n 4096

# ============================ AMG ============================================
source conf/candle.sh
NumRUNS=3
LOG="$ROOTDIR/log/testrun/candle.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for TEST in $TESTCONF; do
	for BINARY in $BINARYS; do
		pushd "`find . -name $BINARY -exec dirname {} \;`"
		for i in `seq 1 $NumRUNS`; do
			python $BINARY >> $LOG 2>&1
		done
		popd
	done
done
cd $ROOTDIR
