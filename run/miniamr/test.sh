#!/bin/bash

SELF="$(readlink -f "${BASH_SOURCE[0]}")"
ROOTDIR="$(readlink -f $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../)"
BenchID="$(basename $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) )"
cd ${ROOTDIR}

source ${ROOTDIR}/conf/host.cfg
source ${ROOTDIR}/conf/env.cfg
get_comp_env_name "${1}"
maybe_submit_job "${COMP}" "${SELF}" "${ROOTDIR}/conf/${BenchID}.sh"
load_compiler_env "${COMP}"

source ${ROOTDIR}/conf/${BenchID}.sh
DEFINPUT=${INPUT}
LOG="${ROOTDIR}/log/$(hostname -s)/testrun/${BenchID}.log"
mkdir -p $(dirname ${LOG})
move_to_scratch_area "${ROOTDIR}" "${APPDIR}"

for TEST in ${TESTCONF}; do
	NumMPI="$(echo ${TEST} | cut -d '|' -f1)"; if skip_conf "${NumMPI}"; then continue; fi
	NumOMP="$(echo ${TEST} | cut -d '|' -f2)"
	X="$(echo ${TEST} | cut -d '|' -f3)"
	Y="$(echo ${TEST} | cut -d '|' -f4)"
	Z="$(echo ${TEST} | cut -d '|' -f5)"
	INPUT="$(echo ${DEFINPUT} | sed -e "s/PX/${X}/" -e "s/PY/${Y}/" -e "s/PZ/${Z}/")"
	echo "$(get_mpi_cmd "${NumMPI}" "${NumOMP}" "${LOG}" "") ${BINARY} ${INPUT}" >> ${LOG} 2>&1
	for i in $(seq 1 ${NumRunsTEST}); do
		START="$(date +%s.%N)"
		timeout --kill-after=30s ${MAXTIME} $(get_mpi_cmd "${NumMPI}" "${NumOMP}" "${LOG}" "") ${BINARY} ${INPUT} >> ${LOG} 2>&1
		clenup_after_mpi_cmd
		if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding ${MAXTIME} timeout" >> ${LOG} 2>&1; fi
		ENDED="$(date +%s.%N)"
		echo "Total running time: $(echo "${ENDED} - ${START}" | bc -l)" >> ${LOG} 2>&1
	done
done
echo "Best ${BenchID} run:"
WALLT="$(/bin/grep 'total GFLOPS' ${LOG} | awk -F 'GFLOPS:' '{print $2}' | sort -r -g | head -1)"
/bin/grep "${WALLT}\|^${MPIRUNCMD}" ${LOG} | /bin/grep -B1 "${WALLT}"
echo ""
cd ${ROOTDIR}
