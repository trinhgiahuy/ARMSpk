#!/bin/bash

SELF="$(readlink -f "${BASH_SOURCE[0]}")"
ROOTDIR="$(readlink -f $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../)"
BenchID="$(basename $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) )"
cd ${ROOTDIR}

source ${ROOTDIR}/conf/host.cfg
source ${ROOTDIR}/conf/env.cfg
get_comp_env_name "${1}"
maybe_submit_job "${COMP}" "${SELF}" "${ROOTDIR}/conf/${BenchID}.sh"
load_compiler_env "${COMP}" "128"

source ${ROOTDIR}/conf/${BenchID}.sh
LOG="${ROOTDIR}/log/$(hostname -s)/testrun/${BenchID}.log"
mkdir -p $(dirname ${LOG})
move_to_scratch_area "${ROOTDIR}" "${APPDIR}"

for TEST in ${TESTCONF}; do
	NumMPI="$(echo ${TEST} | cut -d '|' -f1)"; if skip_conf "${NumMPI}"; then continue; fi
	NumOMP="$(echo ${TEST} | cut -d '|' -f2)"
	pushd "../nek_mgrid_nummpi${NumMPI}"
	# prep input for strong scaling test
	NEPP=$((${ielN} / ${NumMPI}))
	if [ ! -f ./data.rea.bak ]; then cp ./data.rea ./data.rea.bak; fi
	sed -e "s/1   50  1 = iel0/${NEPP}  ${NEPP}  1 = iel0/" -e 's/8   10  2 = nx0/8    8  2 = nx0/' ./data.rea.bak > ./data.rea
	echo "$(get_mpi_cmd "${NumMPI}" "${NumOMP}" "${LOG}" "") ${BINARY} ${INPUT}" >> ${LOG} 2>&1
	for i in $(seq 1 ${NumRunsTEST}); do
		START="$(date +%s.%N)"
		timeout --kill-after=30s ${MAXTIME} $(get_mpi_cmd "${NumMPI}" "${NumOMP}" "${LOG}" "") ${BINARY} ${INPUT} >> ${LOG} 2>&1
		clenup_after_mpi_cmd
		if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding ${MAXTIME} timeout" >> ${LOG} 2>&1; fi
		ENDED="$(date +%s.%N)"
		echo "Total running time: $(echo "${ENDED} - ${START}" | bc -l)" >> ${LOG} 2>&1
	done
	popd
done

echo "Best ${BenchID} run:"
WALLT="$(/bin/grep '^Avg MFlops' ${LOG} | awk -F 'MFlops =' '{print $2}' | sort -g -r | head -1)"
/bin/grep "^Avg.*${WALLT}\|^${MPIRUNCMD}" ${LOG} | /bin/grep -B1 "${WALLT}"
echo ""
cd ${ROOTDIR}
