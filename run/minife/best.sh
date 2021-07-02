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
LOG="${ROOTDIR}/log/$(hostname -s)/bestrun/${BenchID}.log"
mkdir -p $(dirname ${LOG})
move_to_scratch_area "${ROOTDIR}" "${APPDIR}"

for BEST in ${BESTCONF}; do
	for BINARY in ${BBINARY}; do
		NumMPI="$(echo ${BEST} | cut -d '|' -f1)"; if skip_conf "${NumMPI}"; then continue; fi
		NumOMP="$(echo ${BEST} | cut -d '|' -f2)"
		echo "$(get_mpi_cmd "${NumMPI}" "${NumOMP}" "${LOG}" "") ${BINARY} ${INPUT}" >> ${LOG} 2>&1
		for i in $(seq 1 ${NumRunsBEST}); do
			START="$(date +%s.%N)"
			timeout --kill-after=30s ${MAXTIME} $(get_mpi_cmd "${NumMPI}" "${NumOMP}" "${LOG}" "") ${BINARY} ${INPUT} >> ${LOG} 2>&1
			if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then clenup_after_mpi_cmd; echo "Killed after exceeding ${MAXTIME} timeout" >> ${LOG} 2>&1; fi
			ENDED="$(date +%s.%N)"
			cat ./miniFE.*.yaml >> ${LOG} 2>&1
			rm -f ./miniFE.*.yaml
			echo "Total running time: $(echo "${ENDED} - ${START}" | bc -l)" >> ${LOG} 2>&1
		done
	done
done
echo "Best ${BenchID} run:"
BEST="$(/bin/grep 'Total CG Mflops' ${LOG} | awk -F 'Mflops:' '{print $2}' | sort -r -g | head -1)"
/bin/grep "${BEST}\|mpiexec" ${LOG} | /bin/grep -B1 "${BEST}"
echo ""
cd ${ROOTDIR}
