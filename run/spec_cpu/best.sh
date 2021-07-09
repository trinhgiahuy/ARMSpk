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

SPEC0CMD="runcpu --config=nedo.cfg --nobuild --action=run --noreportable --use_submit_for_speed"

source ${ROOTDIR}/conf/${BenchID}.sh
LOGDIR="${ROOTDIR}/log/$(hostname -s)/bestrun/${BenchID}"
mkdir -p ${LOGDIR}
move_to_scratch_area "${ROOTDIR}" "${APPDIR}"

for BENCH in ${BINARY}; do
	BM="$(echo ${BENCH} | cut -d '|' -f1)"
	SIZE="$(echo ${BENCH} | cut -d '|' -f2)"
	SnowflakeNumOMP="$(echo ${BENCH} | cut -d '|' -f3)"
	#XXX: my love for fujitsu needs to be endless
	if [[ "${BM}" = @(*".wrf_"*|*".pop2_"*) ]]; then SPECCMD="export FORT90L='-Wl,-T'; ${SPEC0CMD}"; else SPECCMD="${SPEC0CMD}"; fi
	LOG="${LOGDIR}/${BM}.log"
	for BEST in ${BESTCONF}; do
		NumMPI="$(echo ${BEST} | cut -d '|' -f1)"; if skip_conf "${NumMPI}"; then continue; fi
		NumOMP="$(echo ${BEST} | cut -d '|' -f2)"
		#Some run better with pow2 or other OMP configs, so let them...
		if [ -n "${SnowflakeNumOMP}" ]; then NumOMP="${SnowflakeNumOMP}"; fi
		#For SPECspeed Integer, only the compression benchmark 657.xz_s includes OpenMP, but also does not scale with threads.
		if [[ "${BM}" = @(*".perlbench_"*|*".gcc_"*|*".mcf_"*|*".omnetpp_"*|*".xalancbmk_"*|*".x264_"*|*".deepsjeng_"*|*".leela_"*|*".exchange2_"*|*".xz_"*) ]]; then NumOMP=1; fi
		echo -e "=== runing ${BENCH} ===\nsource ./shrc; ${SPECCMD} --iterations=${NumRunsBEST} --size=${SIZE} --threads=${NumOMP} --define COMP=${COMP} --define RESDIR=0 ${BM}" >> ${LOG} 2>&1
		START="$(date +%s.%N)"
		bash -c "source ./shrc; ${SPECCMD} --iterations=${NumRunsBEST} --size=${SIZE} --threads=${NumOMP} --define COMP=${COMP} --define RESDIR=0 ${BM}" >> ${LOG} 2>&1
		ENDED="$(date +%s.%N)"
		echo "Total running time: $(echo "${ENDED} - ${START}" | bc -l)" >> ${LOG} 2>&1
		REPORT="$(find ${APPROOT}/result -type f -name '*.log' | sort -g | tail -1)"
		cat ${REPORT} >> ${LOG} 2>&1
		/bin/grep 'Error .*runtime=\|Success .*runtime=' ${REPORT} >> ${LOG} 2>&1
		echo "Best ${BENCH} run: " "$(/bin/grep 'Success .*runtime=' ${REPORT} | awk -F 'runtime=' '{print $2}' |cut -d',' -f1 | sort -g | head -1)"
		if [[ "${BM}" = @(*".perlbench_"*|*".gcc_"*|*".mcf_"*|*".omnetpp_"*|*".xalancbmk_"*|*".x264_"*|*".deepsjeng_"*|*".leela_"*|*".exchange2_"*|*".xz_"*) ]]; then break; fi
	done
done
cd ${ROOTDIR}
