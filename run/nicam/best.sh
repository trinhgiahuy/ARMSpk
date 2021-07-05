#!/bin/bash

SELF="$(readlink -f "${BASH_SOURCE[0]}")"
ROOTDIR="$(readlink -f $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../)"
BenchID="$(basename $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) )"
cd ${ROOTDIR}

source ${ROOTDIR}/conf/host.cfg
source ${ROOTDIR}/conf/env.cfg
get_comp_env_name "${1}"
maybe_submit_job "${COMP}" "${SELF}" "${ROOTDIR}/conf/${BenchID}.sh"
load_compiler_env "${COMP}" "8"         # no clue but on Fu it needs much bigger stack

moreMPI="-x FORT_FMT_RECL=400"
if [ -n "${FUJIHOST}" ] || [ -n "${RFX7HOST}" ]; then
	#XXX: my love for fujitsu needs to be endless
	moreMPI+="${moreMPI} -x FORT90L='-Wl,-T'"
fi

source ${ROOTDIR}/conf/${BenchID}.sh
LOG="${ROOTDIR}/log/$(hostname -s)/bestrun/${BenchID}.log"
mkdir -p $(dirname ${LOG})
move_to_scratch_area "${ROOTDIR}" "${APPDIR}"

subOMP="$(for C in ${BESTCONF}; do echo ${C} | cut -d'|' -f2; done | sort -g -u)"
for NumOMP in ${subOMP}; do
	pushd "${APPROOT}/$(echo "${APPDIR}" | sed -e "s#NICAM#omp${NumOMP}#g")"
	# scale down #steps from 11 days to 1 day, and create input data set
	sed -i -e 's/^LSMAX  = 0/LSMAX  = 72/'  ../../test.conf
	make jobshell > /dev/null 2>&1
	if [ -d "${INPUT}" ] && [ -n "${INPUT}" ]; then pushd ${INPUT}; else exit; fi
	ln -s ../../../../bin/nhm_driver .
	ln -s ../../../../data/mnginfo/rl00-prc10.info .
	ln -s ../../../../data/grid/vgrid/vgrid40_24000-600m.dat .
	ln -s ../../../../data/grid/boundary/gl05rl00pe10/boundary_GL05RL00.pe000000 .
	ln -s ../../../../data/grid/boundary/gl05rl00pe10/boundary_GL05RL00.pe000001 .
	ln -s ../../../../data/grid/boundary/gl05rl00pe10/boundary_GL05RL00.pe000002 .
	ln -s ../../../../data/grid/boundary/gl05rl00pe10/boundary_GL05RL00.pe000003 .
	ln -s ../../../../data/grid/boundary/gl05rl00pe10/boundary_GL05RL00.pe000004 .
	ln -s ../../../../data/grid/boundary/gl05rl00pe10/boundary_GL05RL00.pe000005 .
	ln -s ../../../../data/grid/boundary/gl05rl00pe10/boundary_GL05RL00.pe000006 .
	ln -s ../../../../data/grid/boundary/gl05rl00pe10/boundary_GL05RL00.pe000007 .
	ln -s ../../../../data/grid/boundary/gl05rl00pe10/boundary_GL05RL00.pe000008 .
	ln -s ../../../../data/grid/boundary/gl05rl00pe10/boundary_GL05RL00.pe000009 .
	popd; popd
done

for BEST in ${BESTCONF}; do
	NumMPI="$(echo ${BEST} | cut -d '|' -f1)"; if skip_conf "${NumMPI}"; then continue; fi
	NumOMP="$(echo ${BEST} | cut -d '|' -f2)"
	pushd "${APPROOT}/$(echo "${APPDIR}" | sed -e "s#NICAM#omp${NumOMP}#g")/${INPUT}"
	echo "$(get_mpi_cmd "${NumMPI}" "${NumOMP}" "${LOG}" "${moreMPI}") ${BINARY}" >> ${LOG} 2>&1
	for i in $(seq 1 ${NumRunsBEST}); do
		START="$(date +%s.%N)"
		timeout --kill-after=30s ${MAXTIME} $(get_mpi_cmd "${NumMPI}" "${NumOMP}" "${LOG}" "${moreMPI}") ${BINARY} >> ${LOG} 2>&1
		clenup_after_mpi_cmd
		if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding ${MAXTIME} timeout" >> ${LOG} 2>&1; fi
		ENDED="$(date +%s.%N)"
		cat ./msg.pe00000 >> ${LOG} 2>&1
		echo "Total running time: $(echo "${ENDED} - ${START}" | bc -l)" >> ${LOG} 2>&1
	done
	popd
done

# clean up
for NumOMP in ${subOMP}; do
	pushd "${APPROOT}/$(echo "${APPDIR}" | sed -e "s#NICAM#omp${NumOMP}#g")"
	if [ -d "${INPUT}" ] && [ -n "${INPUT}" ]; then rm -rf ${INPUT}; fi
	popd
done

echo "Best ${BenchID} run:"
WALLT="$(/bin/grep '^Walltime' ${LOG} | awk -F 'kernel:' '{print $2}' | sort -g | head -1)"
/bin/grep "${WALLT}\|^${MPIRUNCMD}" ${LOG} | /bin/grep -B1 "${WALLT}"
echo ""
cd ${ROOTDIR}
