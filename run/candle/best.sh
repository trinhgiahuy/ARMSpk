#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

ulimit -s unlimited
ulimit -n 4096

# ============================ CANDLE =========================================
source conf/candle.sh
NumRUNS=10
LOG="$ROOTDIR/log/bestrun/candle.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	for BINARY in $BINARYS; do
		pushd "`find . -name $BINARY -exec dirname {} \;`"
		for i in `seq 1 $NumRUNS`; do
			echo "Start at " `date --iso-8601=s` >> $LOG 2>&1
			python $BINARY >> $LOG 2>&1
			echo "Ended at " `date --iso-8601=s` >> $LOG 2>&1
		done
		popd
	done
done
cd $ROOTDIR
