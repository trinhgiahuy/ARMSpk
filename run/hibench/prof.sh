#!/bin/bash
exit 1 #ignore in this study

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
source $ROOTDIR/dep/spack/share/spack/setup-env.sh
spack load openjdk; spack load maven; spack load scala; spack load hadoop; spack load spark
export HADOOP_HOME=`spack find -p | /bin/grep hadoop | cut -d' ' -f2- | tr -d ' '`
export HADOOP_CONF_DIR=$HADOOP_HOME/etc/hadoop
export HADOOP_COMMON_LIB_NATIVE_DIR=$HADOOP_HOME/lib/native
export LD_LIBRARY_PATH=$HADOOP_COMMON_LIB_NATIVE_DIR:$LD_LIBRARY_PATH
export HADOOP_LOG_DIR=/scr0/hadoop-logs; export YARN_LOG_DIR=$HADOOP_LOG_DIR
export HADOOP_MAPRED_LOG_DIR=$HADOOP_LOG_DIR; export SPARK_LOG_DIR=$HADOOP_LOG_DIR
export SPARK_HOME=`spack find -p | /bin/grep spark | cut -d' ' -f2- | tr -d ' '`
export JAVA_HEAP_MAX="-Xmx8g"
export YARN_HEAPSIZE="4096"
export YARN_RESOURCEMANAGER_HEAPSIZE="4096"
export YARN_TIMELINESERVER_HEAPSIZE="4096"
export HADOOP_HEAPSIZE="8192"
export HADOOP_NAMENODE_INIT_HEAPSIZE="4096"
export MAHOUT_HEAPSIZE="8192"
export NUTCH_HEAPSIZE="8192"

export PATH=$ROOTDIR/dep/sde-external-8.35.0-2019-03-11-lin:$PATH
#export PATH=$ROOTDIR/dep/sde-external-8.63.0-2021-01-18-lin:$PATH	#working not better either
if [ ! -x "`which sde64 2>/dev/null`" ]; then echo "ERROR: SDE missing, please download from Intel sde-external-8.35.0-2019-03-11-lin.tar.bz2 and untar in ./dep folder"; exit; fi;

# ============================ HiBench ====================================
source conf/hibench.sh
LOG="$ROOTDIR/log/`hostname -s`/profrun/hibench.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI=1
	NumOMP=$BEST
	sed -i '/localhost/d' ~/.ssh/known_hosts; sed -i '/0\.0\.0\.0/d' ~/.ssh/known_hosts
	ssh -O exit localhost; ssh -o StrictHostKeyChecking=no localhost echo 0
	ssh -O exit 0.0.0.0;   ssh -o StrictHostKeyChecking=no 0.0.0.0 echo 0
	ssh -O exit 127.0.0.1; ssh -o StrictHostKeyChecking=no 127.0.0.1 echo 0
	# start all this java trash
	killall -9 java.org; killall -9 java
	rm -rf /scr0/hadoop-`logname` /scr0/hadoop-logs /tmp/Jetty_*
	cd $HADOOP_HOME; echo 'Y' | ./bin/hdfs namenode -format; ./sbin/start-dfs.sh; ./bin/hdfs dfs -mkdir /user; ./bin/hdfs dfs -mkdir /user/`logname`; ./sbin/start-yarn.sh; cd -; sleep 10
	if [ "x$RUNSDE" = "xyes" ]; then
		cd $JAVA_HOME/bin; rm -f ./java; ln -s ./java.sh java; cd -
		for BINARY in $BINARYS; do
			cd $SPARK_HOME; echo -e "spark.network.timeout 72000s\nspark.executor.heartbeatInterval 36000s" > ./conf/spark-defaults.conf; cp conf/log4j.properties.template conf/log4j.properties; ./sbin/start-master.sh; ./sbin/start-slave.sh spark://`hostname`:7077; cd -; sleep 10
			echo "Prepare $BINARY $INPUT" >> $LOG 2>&1
			`dirname $BINARY`/../prepare/prepare.sh >> $LOG 2>&1
			sleep 10
			BM=`echo $BINARY | sed -e 's#bin/workloads/##' -e 's#/spark/run.sh##' | tr '/' '.'`
			mkdir -p ${LOG}_${NumMPI}_${NumOMP}_sde/$BM/
			echo -e "`which sde64`\n${LOG}_${NumMPI}_${NumOMP}_sde/$BM/" > /dev/shm/RUNNINGWITHSDE
			echo "$BINARY $INPUT" >> $LOG 2>&1
			START="`date +%s.%N`"
			$BINARY $INPUT >> $LOG 2>&1
			ENDED="`date +%s.%N`"
			echo "Total running time: `echo \"$ENDED - $START\" | bc -l`" >> $LOG 2>&1
			sleep 10; cat report/hibench.report >> $LOG; rm -rf report/; rm -f /dev/shm/RUNNINGWITHSDE; sleep 10
			# shutdown shit again
			cd $SPARK_HOME; ./sbin/stop-slaves.sh; ./sbin/stop-master.sh
			#and wait 10min ...
			c=0; while [ `ps aux | grep org.apache.spark.deploy.worker.Worker | grep java | wc -l` -eq 1 ]; do sleep 1m; c=$((c+1)); if [ $c -eq 10 ]; then break; fi; done
			cd -
		done
		cd $JAVA_HOME/bin; rm -f ./java; ln -s ./java.org java; cd -
	fi
	cd $HADOOP_HOME; ./bin/hdfs dfs -rm -r /user; ./bin/hdfs dfs -rm -r /tmp; ./sbin/stop-yarn.sh; ./sbin/stop-dfs.sh; killall -9 java.org; killall -9 java; cd -
done
echo ""
cd $ROOTDIR
