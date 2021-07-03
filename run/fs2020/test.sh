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
LOGDIR="${ROOTDIR}/log/$(hostname -s)/testrun/${BenchID}"
mkdir -p ${LOGDIR}
move_to_scratch_area "${ROOTDIR}" "${APPDIR}"

for TEST in ${TESTCONF}; do
	NumMPI="$(echo ${TEST} | cut -d '|' -f1)"; if skip_conf "${NumMPI}"; then continue; fi
	NumOMP="$(echo ${TEST} | cut -d '|' -f2)"
	moreMPI="-x OMP_NUM_PARALELL=${NumOMP} -x OMP_PROC_BIND=close -x FLIB_FASTOMP=FALSE -x FLIB_CNTL_BARRIER_ERR=FALSE"
	for BMconf in ${BINARYS}; do
		BINARY="$(echo ${BMconf} | cut -d '|' -f1)"
		MAXOMP="$(echo ${BMconf} | cut -s -d '|' -f2)"
		if [ ! -z ${MAXOMP} ] && [ ${NumOMP} -gt ${MAXOMP} ]; then continue; fi
		if [[ "${BINARY}" = "23."* ]] && [ ${NumOMP} -lt 6 ]; then continue; fi	# needs at least 6 threads to avoid floating-point exception
		LOG="${LOGDIR}/$(echo ${BINARY} | cut -d'/' -f1).log"
		echo "$(get_mpi_cmd "${NumMPI}" "${NumOMP}" "${LOG}" "${moreMPI}") ${BINARY} ${INPUT}" >> ${LOG} 2>&1
		for i in $(seq 1 ${NumRunsTEST}); do
			START="$(date +%s.%N)"
			timeout --kill-after=30s ${MAXTIME} $(get_mpi_cmd "${NumMPI}" "${NumOMP}" "${LOG}" "${moreMPI}") ${BINARY} ${INPUT} >> ${LOG} 2>&1
			if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then clenup_after_mpi_cmd; echo "Killed after exceeding ${MAXTIME} timeout" >> ${LOG} 2>&1; fi
			ENDED="$(date +%s.%N)"
			echo "Total running time: $(echo "${ENDED} - ${START}" | bc -l)" >> ${LOG} 2>&1
		done
		WALLT="$(/bin/grep '^Walltime' ${LOG} | awk -F 'kernel:' '{print $2}' | sed -e 's/D/e/g' | sort -g | head -1)"
		echo "Best ${BINARY} NumOMP=${NumOMP} run: ${WALLT}"
	done
done
cd ${ROOTDIR}
