#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.."
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096

source $ROOTDIR/dep/spack/share/spack/setup-env.sh
spack load gcc; spack load llvm

IACADIR=$ROOTDIR/dep/iaca-lin64
export PATH=$IACADIR:$PATH; export IACAINCL=$IACADIR
if [ ! -x "`which iaca 2>/dev/null`" ]; then echo "ERROR: IACA missing, please download from Intel iaca-version-v3.0-lin64.zip and unzip in ./dep folder"; exit; fi;

ALL=("amg" "babelstream2gb" "babelstream14gb" "comd" "ffvc" "ffb" "hpcg" "hpl" "laghos" "macsio" "miniamr" "minife" "minitri" "modylas" "mvmc" "nekbone" "ngsa" "nicam" "ntchem" "qcd" "sw4lite" "swfft" "xsbench" "dlproxy" "polybench" "spec_cpu" "spec_omp")
SUBBM=("polybench" "spec_cpu" "spec_omp")
#SMALL=("babelstream2gb" "babelstream14gb" "hpcg" "macsio" "minitri" "qcd" "xsbench" "dlproxy")
#MEDIUM=("amg" "ffvc" "ffb" "hpl" "laghos" "modylas" "ngsa" "nicam")
#LARGE=("comd" "miniamr" "minife" "mvmc" "nekbone" "ntchem" "sw4lite" "swfft")

BM=$1
ACTION=$2	# -1=preprocess only;  0=preprocess (and 1+2);   1=use previous asm.b (and 2);   2=analyze logs only
SELRANK=$3	# only a specific rank

if [ -z "${BM}" ] || [ -z "${ACTION}" ]; then echo 'ERR: missing parameter'; exit 1; fi
if ! [[ " ${ALL[@]} " =~ " ${BM} " ]]; then echo 'ERR: unknown BM'; exit 1; fi

if ! [[ " ${SUBBM[@]} " =~ " ${BM} " ]]; then
	D="`find -L ${ROOTDIR}/log/*/profrun -type d -name \"${BM}.log_*_sde\"`"
else
	D="`find -L ${ROOTDIR}/log/*/profrun -type d -name ${BM}`"
fi

if [ ${ACTION} -le 0 ]; then
	started=0
	for BDATA in `find -L ${D} -name '*.dcfg.json.bz2'`; do
		if [ -f ${BDATA%'.dcfg.json.bz2'}.asm.b ]; then continue; fi
		$ROOTDIR/util/parse_basic_blocks.py \
			-j ${BDATA} \
			-b ${BDATA%'.dcfg.json.bz2'}.bb.txt.bz2 \
			-s ${BDATA%'.dcfg.json.bz2'}.asm.b \
			>> ${BDATA%'.dcfg.json.bz2'}.log 2>&1 &
		started=$(($started + 1)); if [ $started -eq 10 ]; then wait; started=0; fi
	done
	wait
	if [ ${ACTION} -le -1 ]; then exit; fi
fi

if [ ${ACTION} -le 1 ]; then
	started=0
	for BDATA in `find -L ${D} -name '*.dcfg.json.bz2'`; do
		$ROOTDIR/util/parse_basic_blocks.py \
			-j ${BDATA} \
			-b ${BDATA%'.dcfg.json.bz2'}.bb.txt.bz2 \
			-l ${BDATA%'.dcfg.json.bz2'}.asm.b \
			>> ${BDATA%'.dcfg.json.bz2'}.log 2>&1 &
		started=$(($started + 1)); if [ $started -eq 10 ]; then wait; started=0; fi
	done
	wait
fi

if [ ${ACTION} -le 2 ]; then
	for BDATA in `find -L ${D} -name '*.dcfg.json.bz2'`; do
		echo ${BDATA%'.dcfg.json.bz2'}.log
		/bin/grep 'Total \|Converted ' ${BDATA%'.dcfg.json.bz2'}.log
		echo 'MAX:' `/bin/grep 'Converted ' ${BDATA%'.dcfg.json.bz2'}.log | cut -d'/' -f4 | sort -r -g | head -1`
	done
fi
