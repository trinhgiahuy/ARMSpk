#!/bin/bash

SELF="$(readlink -f "${BASH_SOURCE[0]}")"
ROOTDIR="$(readlink -f $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../)"
BenchID="$(basename $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) )"
cd ${ROOTDIR}

source ${ROOTDIR}/conf/host.cfg
source ${ROOTDIR}/conf/env.cfg
get_comp_env_name "${1}"
maybe_submit_job "${COMP}" "${SELF}" "${ROOTDIR}/conf/${BenchID}.sh"
load_compiler_env "${COMP}" "-1"

SPECCMD="runspec --config=nedo.cfg --nobuild --action=run --noreportable"

source ${ROOTDIR}/conf/${BenchID}.sh
LOGDIR="${ROOTDIR}/log/$(hostname -s)/testrun/${BenchID}"
mkdir -p ${LOGDIR}
move_to_scratch_area "${ROOTDIR}" "${APPDIR}"

for BENCH in ${BINARY}; do
	BM="$(echo ${BENCH} | cut -d '|' -f1)"
	SIZE="$(echo ${BENCH} | cut -d '|' -f2)"
	LOG="${LOGDIR}/${BM}.log"
	for TEST in ${TESTCONF}; do
		NumMPI="$(echo ${TEST} | cut -d '|' -f1)"; if skip_conf "${NumMPI}"; then continue; fi
		NumOMP="$(echo ${TEST} | cut -d '|' -f2)"
		echo -e "=== runing ${BENCH} ===\nsource ./shrc; ${SPECCMD} --iterations=${NumRunsTEST} --size=${SIZE} --threads=${NumOMP} --define COMP=${COMP} --define RESDIR=0 ${BM}" >> ${LOG} 2>&1
		START="$(date +%s.%N)"
		bash -c "source ./shrc; ${SPECCMD} --iterations=${NumRunsTEST} --size=${SIZE} --threads=${NumOMP} --define COMP=${COMP} --define RESDIR=0 ${BM}" >> ${LOG} 2>&1
		ENDED="$(date +%s.%N)"
		echo "Total running time: $(echo "${ENDED} - ${START}" | bc -l)" >> ${LOG} 2>&1
		REPORT="$(find ${APPROOT}/result -type f -name '*.log' | sort -g | tail -1)"
		cat ${REPORT} >> ${LOG} 2>&1
		/bin/grep 'Error .*runtime=\|Success .*runtime=' ${REPORT} >> ${LOG} 2>&1
		echo "Best ${BENCH} run: " "$(/bin/grep 'Success .*runtime=' ${REPORT} | awk -F 'runtime=' '{print $2}' |cut -d',' -f1 | sort -g | head -1)"
	done
done
cd ${ROOTDIR}
