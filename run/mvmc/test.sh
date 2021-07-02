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

if [ -n "${FUJIHOST}" ]; then module load Python2-CN; export FLIB_CNTL_BARRIER_ERR=FALSE; fi

source ${ROOTDIR}/conf/${BenchID}.sh
LOG="${ROOTDIR}/log/$(hostname -s)/testrun/${BenchID}.log"
mkdir -p $(dirname ${LOG})
move_to_scratch_area "${ROOTDIR}" "${APPDIR}"

for TEST in ${TESTCONF}; do
	NumMPI="$(echo ${TEST} | cut -d '|' -f1)"; if skip_conf "${NumMPI}"; then continue; fi
	NumOMP="$(echo ${TEST} | cut -d '|' -f2)"
	# prepare input for strong scaling (scale down a bit from the default 64 node run)
	if [ -d ./job_mpi${NumMPI} ]; then rm -rf job_mpi${NumMPI}; fi
	sed -i -e 's/^Lx = Ly = 12/Lx = Ly = 4 #12/' -e 's/^NTotalSample = 4096/NTotalSample = 512 #4096/' -e 's/^NOuterMPI = 64/NOuterMPI = 2 #64/' ./makeDef/makeDef_large.py
	python2 ./makeDef/makeDef_large.py ${NumMPI}
	cd ./job_mpi${NumMPI}
	echo "$(get_mpi_cmd "${NumMPI}" "${NumOMP}" "${LOG}" "") ${BINARY} ${INPUT}" >> ${LOG} 2>&1
	for i in $(seq 1 ${NumRunsTEST}); do
		START="$(date +%s.%N)"
		timeout --kill-after=30s ${MAXTIME} $(get_mpi_cmd "${NumMPI}" "${NumOMP}" "${LOG}" "") ${BINARY} ${INPUT} >> ${LOG} 2>&1
		if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then clenup_after_mpi_cmd; echo "Killed after exceeding ${MAXTIME} timeout" >> ${LOG} 2>&1; fi
		ENDED="$(date +%s.%N)"
		cat Lx*Ly*/zvo_HitachiTimer.dat >> ${LOG} 2>&1
		echo "Total running time: $(echo "${ENDED} - ${START}" | bc -l)" >> ${LOG} 2>&1
		rm -f Lx*Ly*/zvo_*
	done
	cd ../
done
echo "Best ${BenchID} run:"
TEST="$(/bin/grep '^Walltime' ${LOG} | awk -F 'kernel:' '{print $2}' | sort -g | head -1)"
/bin/grep "${TEST}\|mpiexec" ${LOG} | /bin/grep -B1 "${TEST}"
echo ""
cd ${ROOTDIR}
