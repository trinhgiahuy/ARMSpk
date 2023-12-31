#!/bin/bash
exit 1 #ignore in this study

SELF="$(readlink -f "${BASH_SOURCE[0]}")"
ROOTDIR="$(readlink -f $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../)"
BenchID="$(basename $( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd ) )"
cd ${ROOTDIR}

source ${ROOTDIR}/conf/host.cfg
source ${ROOTDIR}/conf/env.cfg
get_comp_env_name "${1}"
maybe_submit_job "${COMP}" "${SELF}" "${ROOTDIR}/conf/${BenchID}.sh"
load_compiler_env "${COMP}"

source ${ROOTDIR}/dep/spack/share/spack/setup-env.sh
spack load openjdk; spack load maven; spack load scala; spack load hadoop; spack load spark
export HADOOP_HOME="$(spack find -p | /bin/grep hadoop | cut -d' ' -f2- | tr -d ' ')"
export HADOOP_CONF_DIR=${HADOOP_HOME}/etc/hadoop
export HADOOP_COMMON_LIB_NATIVE_DIR=${HADOOP_HOME}/lib/native
export LD_LIBRARY_PATH=${HADOOP_COMMON_LIB_NATIVE_DIR}:${LD_LIBRARY_PATH}
export HADOOP_LOG_DIR=/scr0/hadoop-logs; export YARN_LOG_DIR=${HADOOP_LOG_DIR}
export HADOOP_MAPRED_LOG_DIR=${HADOOP_LOG_DIR}; export SPARK_LOG_DIR=${HADOOP_LOG_DIR}
export SPARK_HOME="$(spack find -p | /bin/grep spark | cut -d' ' -f2- | tr -d ' ')"
export JAVA_HEAP_MAX="-Xmx8g"
export YARN_HEAPSIZE="4096"
export YARN_RESOURCEMANAGER_HEAPSIZE="4096"
export YARN_TIMELINESERVER_HEAPSIZE="4096"
export HADOOP_HEAPSIZE="8192"
export HADOOP_NAMENODE_INIT_HEAPSIZE="4096"
export MAHOUT_HEAPSIZE="8192"
export NUTCH_HEAPSIZE="8192"

source ${ROOTDIR}/conf/${BenchID}.sh
LOG="${ROOTDIR}/log/$(hostname -s)/bestrun/${BenchID}.log"
mkdir -p $(dirname ${LOG})
move_to_scratch_area "${ROOTDIR}" "${APPDIR}"

for BEST in ${BESTCONF}; do
	NumMPI="$(echo ${BEST} | cut -d '|' -f1)"; if skip_conf "${NumMPI}"; then continue; fi
	NumOMP="$(echo ${BEST} | cut -d '|' -f2)"
	sed -i '/localhost/d' ~/.ssh/known_hosts; sed -i '/0\.0\.0\.0/d' ~/.ssh/known_hosts
	ssh -O exit localhost; ssh -o StrictHostKeyChecking=no localhost echo 0
	ssh -O exit 0.0.0.0;   ssh -o StrictHostKeyChecking=no 0.0.0.0 echo 0
	ssh -O exit 127.0.0.1; ssh -o StrictHostKeyChecking=no 127.0.0.1 echo 0
	# start all this java trash
	killall -9 java.org; killall -9 java
	rm -rf /scr0/hadoop-`whoami` /scr0/hadoop-logs /tmp/Jetty_*
	cd ${HADOOP_HOME}; echo 'Y' | ./bin/hdfs namenode -format; ./sbin/start-dfs.sh; ./bin/hdfs dfs -mkdir /user; ./bin/hdfs dfs -mkdir /user/$(logname); ./sbin/start-yarn.sh; cd -; sleep 10
	for BINARY in ${BINARYS}; do
		cd ${SPARK_HOME}; rm -f ./conf/spark-defaults.conf ./conf/log4j.properties; ./sbin/start-master.sh ; ./sbin/start-slave.sh spark://$(hostname):7077; cd -; sleep 10
		echo "Prepare ${BINARY} ${INPUT}" >> ${LOG} 2>&1
		`dirname ${BINARY}`/../prepare/prepare.sh >> ${LOG} 2>&1
		sleep 10
		echo "${BINARY} ${INPUT}" >> ${LOG} 2>&1
		for i in $(seq 1 ${NumRunsBEST}); do
			START="$(date +%s.%N)"
			timeout --kill-after=30s ${MAXTIME} ${BINARY} ${INPUT} >> ${LOG} 2>&1
			clenup_after_mpi_cmd
			if [ "x$?" = "x124" ] || [ "x$?" = "x137" ]; then echo "Killed after exceeding ${MAXTIME} timeout" >> ${LOG} 2>&1; fi
			ENDED="$(date +%s.%N)"
			echo "Total running time: $(echo "${ENDED} - ${START}" | bc -l)" >> ${LOG} 2>&1
			sleep 10; cat report/hibench.report >> ${LOG}; rm -rf report/
		done
		# shutdown shit again
		cd ${SPARK_HOME}; ./sbin/stop-slaves.sh; ./sbin/stop-master.sh; cd -
	done
	cd ${HADOOP_HOME}; ./sbin/stop-yarn.sh; ./sbin/stop-dfs.sh; killall -9 java.org; killall -9 java; cd -
done
echo ""
cd ${ROOTDIR}
