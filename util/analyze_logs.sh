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

ALL=("amg" "babelstream2gb" "babelstream14gb" "comd" "ffvc" "ffb" "hpcg" "hpl" "laghos" "macsio" "miniamr" "minife" "minitri" "modylas" "mvmc" "nekbone" "ngsa" "nicam" "ntchem" "qcd" "sw4lite" "swfft" "xsbench")
SMALL=("babelstream2gb" "babelstream14gb" "hpcg" "macsio" "minitri" "qcd" "xsbench")
MEDIUM=("amg" "ffvc" "ffb" "hpl" "laghos" "modylas" "ngsa" "nicam")
LARGE=("comd" "miniamr" "minife" "mvmc" "nekbone" "ntchem" "sw4lite" "swfft")

BM=$1
ACTION=$2	# 0=all, 1=use previous asm.b 2=analyze logs
SELRANK=$3	# only a specific rank

if [ -z "${BM}" ] || [ -z "${ACTION}" ]; then
	echo 'ERR: missing parameter'
	exit 1
fi

if ! [[ " ${ALL[@]} " =~ " ${BM} " ]]; then
	echo 'ERR: unknown BM'
	exit 1
fi

D="`find -L ${ROOTDIR}/log -type d -name ${BM}.log_sde`"
if [[ ${BM} = "ngsa" ]]; then
	N=`find -L ${D} -name '*.rank-*' | wc -l`
else
	N=`ls ${D}/*.bb.txt.bz2 | wc -l`
fi
NH=`echo $((N/2))`

if [ ${ACTION} -le 0 ]; then
	if [[ " ${LARGE[@]} " =~ " ${BM} " ]]; then
		for R in `seq 0 $((NH-1))`; do
			for BDATA in `find -L ${D} -name '*.dcfg.json.bz2' | /bin/grep -e "rank-${R}[\./]"`; do
				$ROOTDIR/util/parse_basic_blocks.py \
					-j ${BDATA} \
					-b ${BDATA%'.dcfg.json.bz2'}.bb.txt.bz2 \
					-s ${BDATA%'.dcfg.json.bz2'}.asm.b \
					>> ${BDATA%'.dcfg.json.bz2'}.log 2>&1 &
			done
		done
		wait
		for R in `seq ${NH} $((N-1))`; do
			for BDATA in `find -L ${D} -name '*.dcfg.json.bz2' | /bin/grep -e "rank-${R}[\./]"`; do
				$ROOTDIR/util/parse_basic_blocks.py \
					-j ${BDATA} \
					-b ${BDATA%'.dcfg.json.bz2'}.bb.txt.bz2 \
					-s ${BDATA%'.dcfg.json.bz2'}.asm.b \
					>> ${BDATA%'.dcfg.json.bz2'}.log 2>&1 &
			done
		done
		wait
	else
		for R in `seq 0 $((N-1))`; do
			for BDATA in `find -L ${D} -name '*.dcfg.json.bz2' | /bin/grep -e "rank-${R}[\./]"`; do
				$ROOTDIR/util/parse_basic_blocks.py \
					-j ${BDATA} \
					-b ${BDATA%'.dcfg.json.bz2'}.bb.txt.bz2 \
					-s ${BDATA%'.dcfg.json.bz2'}.asm.b \
					>> ${BDATA%'.dcfg.json.bz2'}.log 2>&1 &
			done
		done
		wait
	fi
fi

if [ ${ACTION} -le 1 ]; then
	if [[ " ${LARGE[@]} " =~ " ${BM} " ]]; then
		for R in `seq 0 $((NH-1))`; do
			for BDATA in `find -L ${D} -name '*.dcfg.json.bz2' | /bin/grep -e "rank-${R}[\./]"`; do
				$ROOTDIR/util/parse_basic_blocks.py \
					-j ${BDATA} \
					-b ${BDATA%'.dcfg.json.bz2'}.bb.txt.bz2 \
					-l ${BDATA%'.dcfg.json.bz2'}.asm.b \
					>> ${BDATA%'.dcfg.json.bz2'}.log 2>&1 &
			done
		done
		wait
		for R in `seq ${NH} $((N-1))`; do
			for BDATA in `find -L ${D} -name '*.dcfg.json.bz2' | /bin/grep -e "rank-${R}[\./]"`; do
				$ROOTDIR/util/parse_basic_blocks.py \
					-j ${BDATA} \
					-b ${BDATA%'.dcfg.json.bz2'}.bb.txt.bz2 \
					-l ${BDATA%'.dcfg.json.bz2'}.asm.b \
					>> ${BDATA%'.dcfg.json.bz2'}.log 2>&1 &
			done
		done
		wait
	else
		for R in `seq 0 $((N-1))`; do
			for BDATA in `find -L ${D} -name '*.dcfg.json.bz2' | /bin/grep -e "rank-${R}[\./]"`; do
				$ROOTDIR/util/parse_basic_blocks.py \
					-j ${BDATA} \
					-b ${BDATA%'.dcfg.json.bz2'}.bb.txt.bz2 \
					-l ${BDATA%'.dcfg.json.bz2'}.asm.b \
					>> ${BDATA%'.dcfg.json.bz2'}.log 2>&1 &
			done
		done
		wait
	fi
fi

if [ ${ACTION} -le 2 ]; then
	for R in `seq 0 $((N-1))`; do
		for BDATA in `find -L ${D} -name '*.dcfg.json.bz2' | /bin/grep -e "rank-${R}[\./]"`; do
			echo ${BDATA%'.dcfg.json.bz2'}.log
			/bin/grep '^Total\|Converted' ${BDATA%'.dcfg.json.bz2'}.log
			echo 'MAX:' `/bin/grep 'Converted' ${BDATA%'.dcfg.json.bz2'}.log | cut -d'/' -f4 | sort -r -g | head -1`
		done
	done
fi
