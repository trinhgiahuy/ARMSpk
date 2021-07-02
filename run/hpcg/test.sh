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

if [ -n "${XEONHOST}" ]; then                           moreMPI="-x KMP_AFFINITY=granularity=fine,compact,1,0"
elif [ -n "${IKNLHOST}" ] || [ -n "${IKNMHOST}" ]; then moreMPI="-x KMP_AFFINITY=compact"
else                                                    moreMPI=""; fi

source ${ROOTDIR}/conf/${BenchID}.sh
DEFINPUT=${INPUT}
LOG="${ROOTDIR}/log/$(hostname -s)/testrun/${BenchID}.log"
mkdir -p $(dirname ${LOG})
move_to_scratch_area "${ROOTDIR}" "${APPDIR}"

for TEST in ${TESTCONF}; do
	NumMPI="$(echo ${TEST} | cut -d '|' -f1)"; if skip_conf "${NumMPI}"; then continue; fi
	NumOMP="$(echo ${TEST} | cut -d '|' -f2)"
	# test to identify hpcg's internal dimensions
	rm -f hpcg_log_* n*.yaml
	rm -f hpcg20*.txt HPCG-Benchmark_3*.txt
	$(get_mpi_cmd "${NumMPI}" "1" "/dev/null" "${moreMPI}") ${BINARY} -n 1 > /dev/null 2>&1
	if [ ! "x$?" = "x0" ]; then continue; fi
	if [ -f n*.yaml ]; then
		X=$((${MAXXYZ} / $(/bin/grep 'npx:' n*.yaml | awk -F 'npx:' '{print $2}')))
		Y=$((${MAXXYZ} / $(/bin/grep 'npy:' n*.yaml | awk -F 'npy:' '{print $2}')))
		Z=$((${MAXXYZ} / $(/bin/grep 'npz:' n*.yaml | awk -F 'npz:' '{print $2}')))
	elif [ -f HPCG-Benchmark_3*.txt ]; then
		#non-intel version needs to be div8 https://github.com/hpcg-benchmark/hpcg/issues/47
		X=$(((${MAXXYZ} / $(/bin/grep 'npx=' HPCG-Benchmark_3*.txt | awk -F 'npx=' '{print $2}') / 8) * 8))
		Y=$(((${MAXXYZ} / $(/bin/grep 'npy=' HPCG-Benchmark_3*.txt | awk -F 'npy=' '{print $2}') / 8) * 8))
		Z=$(((${MAXXYZ} / $(/bin/grep 'npz=' HPCG-Benchmark_3*.txt | awk -F 'npz=' '{print $2}') / 8) * 8))
	else continue; fi
	rm -f hpcg_log_* n*.yaml
	rm -f hpcg20*.txt HPCG-Benchmark_3*.txt
	INPUT="$(echo ${DEFINPUT} | sed -e "s/NX/${X}/" -e "s/NY/${Y}/" -e "s/NZ/${Z}/")"
	echo "$(get_mpi_cmd "${NumMPI}" "${NumOMP}" "${LOG}" "${moreMPI}") ${BINARY} ${INPUT}" >> ${LOG} 2>&1
	for i in $(seq 1 ${NumRunsTEST}); do
		START="$(date +%s.%N)"
		timeout --kill-after=30s ${MAXTIME} $(get_mpi_cmd "${NumMPI}" "${NumOMP}" "${LOG}" "${moreMPI}") ${BINARY} ${INPUT} >> ${LOG} 2>&1
		if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then clenup_after_mpi_cmd; echo "Killed after exceeding ${MAXTIME} timeout" >> ${LOG} 2>&1; fi
		ENDED="$(date +%s.%N)"
		cat hpcg_log_* >> ${LOG} 2>&1
		cat n*.yaml >> ${LOG} 2>&1
		rm -f hpcg_log_* n*.yaml
		cat hpcg20*.txt >> ${LOG} 2>&1
		cat HPCG-Benchmark_3*.txt >> ${LOG} 2>&1
		rm -f hpcg20*.txt HPCG-Benchmark_3*.txt
		echo "Total running time: $(echo "${ENDED} - ${START}" | bc -l)" >> ${LOG} 2>&1
	done
done
echo "Best ${BenchID} run:"
TEST="$(/bin/grep '^Walltime' ${LOG} | awk -F 'kernel:' '{print $2}' | sort -g | head -1)"
/bin/grep "${TEST}\|mpiexec" ${LOG} | /bin/grep -B1 "${TEST}"
echo ""
cd ${ROOTDIR}
