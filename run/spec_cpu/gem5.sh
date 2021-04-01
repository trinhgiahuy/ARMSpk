#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
#source $ROOTDIR/conf/intel.cfg
#source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
#source $ROOTDIR/dep/spack/share/spack/setup-env.sh
#spack load gcc@8.4.0
### all static, no need for compilers and envs

export HOSTNAME="kiev0" # we don't care where gem runs
GEM5="$ROOTDIR/dep/gem5_riken/build/ARM/gem5.opt"
GEM5SE="$ROOTDIR/dep/gem5_riken/configs/example/se.py"

if [ $1 -eq 1 ]; then
	ARCHCONF="--cpu-type=O3_ARM_PostK_3 --caches --l2_size=16MB --mem_bus_width=64 --mem_resp_width=128 --mem-size=32GB"    # traditional
elif [ $1 -eq 2 ]; then
	ARCHCONF="--cpu-type=O3_ARM_PostK_3 --caches --l2_size=256MB --mem_bus_width=96 --mem_resp_width=192 --mem-size=32GB"   # aggressive
else echo "err: missing which arch conf to run"; exit 1; fi

# special core pinning? XXX: no, we run all simul on a node
#if [ -n $2 ] && which numactl >/dev/null 2>&1; then PIN="numactl -C $2"; else PIN=""; fi

SPECCMD="runcpu --config=nedo.cfg --nobuild --action=run --noreportable --use_submit_for_speed"

# ============================ SPEC CPU =======================================
source conf/spec_cpu.sh
LOGDIR="$ROOTDIR/log/`hostname -s`/gem5run/spec_cpu"
mkdir -p $LOGDIR
cd $APPDIR
for BENCH in $BINARY; do
	BM="`echo $BENCH | cut -d '|' -f1`"
	COMP="fuji"     #"`echo $BENCH | cut -d '|' -f2`"
	SIZE="`echo $BENCH | cut -d '|' -f3`"
	LOG="$LOGDIR/${BM}/conf${1}.log"; mkdir -p ${LOG}_stat
	NumMPI=1
	if [ $1 -eq 1 ]; then NumOMP="20"; elif [ $1 -eq 2 ]; then NumOMP="32"; fi
	echo -e "OMP_NUM_THREADS=$NumOMP\nOMP_NUM_PARALELL=$NumOMP\nFLIB_FASTOMP=FALSE\nFLIB_CNTL_BARRIER_ERR=FALSE" > `dirname $LOG`/omp${NumOMP}.txt
	cat <<EOF > ${LOG}_stat/gem5wrap.sh
#!/bin/bash
COMMAND="\$@"; BINARY="\`echo \$COMMAND | cut -d' ' -f1\`"; INPUT="\`echo \$COMMAND | cut -d' ' -f2-\`";
$GEM5 -d ${LOG}_stat $GEM5SE -c \$BINARY -o "\$INPUT" -n $NumOMP -e `dirname $LOG`/omp${NumOMP}.txt $ARCHCONF
EOF
	echo "=== gem5 run for $BENCH ===" >> $LOG 2>&1
	START="`date +%s.%N`"
	bash -c "source ./shrc; $SPECCMD --iterations=1 --size=$SIZE --threads=$NumOMP --define COMP=$COMP --define RESDIR=${LOG}_stat >> $LOG 2>&1
	ENDED="`date +%s.%N`"
	echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
done
cd $ROOTDIR
