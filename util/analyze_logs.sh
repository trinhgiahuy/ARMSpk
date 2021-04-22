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

ALL=("amg" "babelstream2gb" "babelstream14gb" "comd" "ffvc" "ffb" "hpcg" "hpl" "laghos" "macsio" "miniamr" "minife" "minitri" "modylas" "mvmc" "nekbone" "ngsa" "nicam" "ntchem" "qcd" "sw4lite" "swfft" "xsbench" "dlproxy" "polybench" "spec_cpu" "spec_omp" "fs2020")
SUBBM=("polybench" "spec_cpu" "spec_omp" "fs2020")
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
	if [[ "${BM}" = "spec_"* ]]; then
		for E in `find -L ${D} -name '*.dcfg.json.bz2' -exec dirname {} \; | sort -u`; do
			echo ${E}
			jobSumLLVM=0; jobSumIACA=0; jobSumOSACA=0
			for BDATA in `find -L ${E} -name '*.dcfg.json.bz2'`; do
				echo -n "${BDATA%'.dcfg.json.bz2'}.log  ---  "
				for TOOL in "LLVM" "IACA" "OSACA"; do
					subWorkLoadMax="`/bin/grep \"$TOOL: Total \" ${BDATA%'.dcfg.json.bz2'}.log | cut -d':' -f3 | tr -s ' ' | tr 'e' 'E' | sort -r -g | head -1`"
					echo -n "$TOOL MAX: $subWorkLoadMax  "
					if [[ "$TOOL" = "LLVM" ]] && [ ! -z $subWorkLoadMax ]; then jobSumLLVM=$(echo "$jobSumLLVM + $subWorkLoadMax" | bc -l); fi
					if [[ "$TOOL" = "IACA" ]] && [ ! -z $subWorkLoadMax ]; then jobSumIACA=$(echo "$jobSumIACA + $subWorkLoadMax" | bc -l); fi
					if [[ "$TOOL" = "OSACA" ]] && [ ! -z $subWorkLoadMax ]; then jobSumOSACA=$(echo "$jobSumOSACA + $subWorkLoadMax" | bc -l); fi
				done
				echo
			done
			echo -e "\nSum for the (multi-workload) job: [LLVM; IACA; OSACA] = [ $jobSumLLVM ; $jobSumIACA ; $jobSumOSACA ]\n"
		done
	elif [[ "${BM}" = "ngsa"* ]]; then
		for E in `find -L ${D} -name '*.dcfg.json.bz2' -exec dirname {} \; | xargs dirname | xargs dirname | sort -u`; do
			rankMaxLLVM=0; rankMaxIACA=0; rankMaxOSACA=0
			for F in `find -L ${E} -name '*.dcfg.json.bz2' -exec dirname {} \; | xargs dirname | sort -u`; do
				echo ${E}
				jobSumLLVM=0; jobSumIACA=0; jobSumOSACA=0
				for BDATA in `find -L ${F} -name '*.dcfg.json.bz2'`; do
					echo -n "${BDATA%'.dcfg.json.bz2'}.log  ---  "
					for TOOL in "LLVM" "IACA" "OSACA"; do
						subWorkLoadMax="`/bin/grep \"$TOOL: Total \" ${BDATA%'.dcfg.json.bz2'}.log | cut -d':' -f3 | tr -s ' ' | tr 'e' 'E' | sort -r -g | head -1`"
						echo -n "$TOOL MAX: $subWorkLoadMax  "
						if [[ "$TOOL" = "LLVM" ]] && [ ! -z $subWorkLoadMax ]; then jobSumLLVM=$(echo "$jobSumLLVM + $subWorkLoadMax" | bc -l); fi
						if [[ "$TOOL" = "IACA" ]] && [ ! -z $subWorkLoadMax ]; then jobSumIACA=$(echo "$jobSumIACA + $subWorkLoadMax" | bc -l); fi
						if [[ "$TOOL" = "OSACA" ]] && [ ! -z $subWorkLoadMax ]; then jobSumOSACA=$(echo "$jobSumOSACA + $subWorkLoadMax" | bc -l); fi
					done
					echo
				done
				if (( $(echo "$jobSumLLVM > $rankMaxLLVM" |bc -l) )); then rankMaxLLVM=$jobSumLLVM; fi
				if (( $(echo "$jobSumIACA > $rankMaxIACA" |bc -l) )); then rankMaxIACA=$jobSumIACA; fi
				if (( $(echo "$jobSumOSACA > $rankMaxOSACA" |bc -l) )); then rankMaxOSACA=$jobSumOSACA; fi
			done
			echo -e "\nMax for the (multi-rank/worload) job: [LLVM; IACA; OSACA] = [ $rankMaxLLVM ; $rankMaxIACA ; $rankMaxOSACA ]\n"
		done
	else
		for E in `find -L ${D} -name '*.dcfg.json.bz2' -exec dirname {} \; | sort -u`; do
			echo ${E}
			jobMaxLLVM=0; jobMaxIACA=0; jobMaxOSACA=0
			for BDATA in `find -L ${E} -name '*.dcfg.json.bz2'`; do
				echo -n "${BDATA%'.dcfg.json.bz2'}.log  ---  "
				for TOOL in "LLVM" "IACA" "OSACA"; do
					#/bin/grep "$TOOL: Total CPU cycles on" ${BDATA%'.dcfg.json.bz2'}.log
					rankMax="`/bin/grep \"$TOOL: Total \" ${BDATA%'.dcfg.json.bz2'}.log | cut -d':' -f3 | tr -s ' ' | tr 'e' 'E' | sort -r -g | head -1`"
					echo -n "$TOOL MAX: $rankMax  "
					if [[ "$TOOL" = "LLVM" ]] && [ ! -z $rankMax ] && (( $(echo "$rankMax > $jobMaxLLVM" |bc -l) )); then jobMaxLLVM=$rankMax; fi
					if [[ "$TOOL" = "IACA" ]] && [ ! -z $rankMax ] && (( $(echo "$rankMax > $jobMaxIACA" |bc -l) )); then jobMaxIACA=$rankMax; fi
					if [[ "$TOOL" = "OSACA" ]] && [ ! -z $rankMax ] && (( $(echo "$rankMax > $jobMaxOSACA" |bc -l) )); then jobMaxOSACA=$rankMax; fi
				done
				echo
			done
			echo -e "\nMax for the (multi-rank) job: [LLVM; IACA; OSACA] = [ $jobMaxLLVM ; $jobMaxIACA ; $jobMaxOSACA ]\n"
		done
	fi
fi
