#!/bin/bash
exit 1 #ignore in this study

function PreprocessInput {
	PROCS=${1}
	INDIR=${2}
	for P in $(seq 0 $((PROCS - 1))); do
		mkdir -p ${INDIR}/00-read-rank/${P}
	done
	for N in $(seq 1 2); do
		P=0; C=0
		for S in $(seq -w 0 11); do
			if [ "${C}" = "$((12/PROCS))" ]; then
				P=$((P + 1))
				C=0
			fi
			cp ${INDIR}/00-read/part_${N}.${S} ${INDIR}/00-read-rank/${P}/
			C=$((C + 1))
		done
	done
}

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
LOG="${ROOTDIR}/log/$(hostname -s)/testrun/${BenchID}.log"
mkdir -p $(dirname ${LOG})
move_to_scratch_area "${ROOTDIR}" "${APPDIR}"

for TEST in ${TESTCONF}; do
	NumMPI="$(echo ${TEST} | cut -d '|' -f1)"; if skip_conf "${NumMPI}"; then continue; fi
	NumOMP="$(echo ${TEST} | cut -d '|' -f2)"
	echo "$(get_mpi_cmd ${NumMPI} ${NumOMP} ${LOG} "") ${BINARY} ${INPUT}" >> ${LOG} 2>&1
	for i in $(seq 1 ${NumRunsTEST}); do
		# prep input (dep on numMPI; up to 12 supported)
		PreprocessInput ${NumMPI} ${INPUTDIR}
		START="$(date +%s.%N)"
		timeout --kill-after=30s ${MAXTIME} $(get_mpi_cmd ${NumMPI} ${NumOMP} ${LOG} "") ${BINARY} ${INPUT} >> ${LOG} 2>&1
		if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then clenup_after_mpi_cmd; echo "Killed after exceeding ${MAXTIME} timeout" >> ${LOG} 2>&1; fi
		ENDED="$(date +%s.%N)"
		echo "Total running time: $(echo "${ENDED} - ${START}" | bc -l)" >> ${LOG} 2>&1
		# clean up
		rm -rf workflow_*
		if [ -d ${INPUTDIR}/00-read-rank ]; then
			rm -rf ${INPUTDIR}/00-read-rank
		fi
	done
done
echo "Best ${BenchID} run:"
TEST="$(/bin/grep '^Walltime' ${LOG} | awk -F 'kernel:' '{print $2}' | sort -g | head -1)"
/bin/grep "${TEST}\|mpiexec" ${LOG} | /bin/grep -B1 "${TEST}"
echo ""
cd ${ROOTDIR}
