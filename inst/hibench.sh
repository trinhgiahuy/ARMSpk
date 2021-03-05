#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort
alias ar=`which xiar`
alias ld=`which xild`
export ADVISOR_2018_DIR=${ADVISOR_2019_DIR}

source $ROOTDIR/dep/spack/share/spack/setup-env.sh
spack load openmpi@3.1.6%intel@19.0.1.144
export OMPI_CC=$I_MPI_CC
export OMPI_CXX=$I_MPI_CXX
export OMPI_F77=$I_MPI_F77
export OMPI_FC=$I_MPI_F90

spack load openjdk; spack load maven; spack load scala; spack load hadoop; spack load spark
export HADOOP_HOME=`spack find -p | /bin/grep hadoop | cut -d' ' -f2- | tr -d ' '`
export HADOOP_CONF_DIR=$HADOOP_HOME/etc/hadoop
export HADOOP_COMMON_LIB_NATIVE_DIR=$HADOOP_HOME/lib/native
export LD_LIBRARY_PATH=$HADOOP_COMMON_LIB_NATIVE_DIR:$LD_LIBRARY_PATH
export HADOOP_LOG_DIR=/scr0/hadoop-logs; export YARN_LOG_DIR=$HADOOP_LOG_DIR
export HADOOP_MAPRED_LOG_DIR=$HADOOP_LOG_DIR; export SPARK_LOG_DIR=$HADOOP_LOG_DIR
export SPARK_HOME=`spack find -p | /bin/grep spark | cut -d' ' -f2- | tr -d ' '`

BM="HiBench"
VERSION="v7.1.1"
if [ ! -f $ROOTDIR/$BM/sparkbench/micro/target/sparkbench-micro-7.1.1.jar ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	# need some hack derived from https://www.themetabytes.com/2017/11/25/ld_preload-hacks/
	icc -c -fPIC -I${ADVISOR_2018_DIR}/include -I. -I$JAVA_HOME/include -I$JAVA_HOME/include/linux ../util/sde_java_hack.c
	ld -shared ./sde_java_hack.o -o sde_java_hack.so ${ADVISOR_2018_DIR}/lib64/libittnotify.a
	#fix a stupid bug...https://github.com/Intel-bigdata/HiBench/issues/534
	sed -i -e 's/^ *$CMD/eval $CMD/' bin/functions/workload_functions.sh
	mvn -Phadoopbench -Psparkbench -Dspark=2.4 -Dscala=2.11 clean package
	#change config options to fit our need
	cp conf/hadoop.conf.template conf/hadoop.conf
	cp conf/spark.conf.template conf/spark.conf
	sed -i -e "s#/PATH/TO/YOUR/HADOOP/ROOT#$HADOOP_HOME#" conf/hadoop.conf
	sed -i -e "s#^hibench.hdfs.master.*#hibench.hdfs.master hdfs://localhost:9000/user/`logname`#" conf/hadoop.conf
	sed -i -e "s#/PATH/TO/YOUR/SPARK/HOME#$SPARK_HOME#" conf/spark.conf
	sed -i -e "s#^hibench.spark.master.*#hibench.spark.master yarn#" conf/spark.conf
	## derive from https://www.delltechnologies.com/nl-be/collaterals/unauth/competitive-reports/products/servers/poweredge_fx2_apache_spark_tco_comparison.pdf
	sed -i -e "s#^hibench.yarn.executor.num.*#hibench.yarn.executor.num $((`lscpu | /bin/grep ^Socket | cut -d ':' -f2` * `lscpu | /bin/grep ^Core | cut -d ':' -f2` / 2 - 4))#" conf/spark.conf
	sed -i -e "s#^hibench.yarn.executor.cores.*#hibench.yarn.executor.cores 1#" conf/spark.conf
	sed -i -e "s#^spark.executor.memory.*#spark.executor.memory 16g#" conf/spark.conf
	sed -i -e "s#^spark.driver.memory.*#spark.driver.memory 16g#" conf/spark.conf
	sed -i -e "s#^hibench.streambench.spark.receiverNumber.*#hibench.streambench.spark.receiverNumber 1#" conf/spark.conf
	sed -i -e "s#^hibench.streambench.spark.storageLevel.*#hibench.streambench.spark.storageLevel 0#" conf/spark.conf
	sed -i -e "s#^hibench.streambench.spark.checkpointPath.*#hibench.streambench.spark.checkpointPath /scr0/hadoop-logs#" conf/spark.conf
	sed -i -e "s#^hibench.masters.hostnames.*#hibench.masters.hostnames localhost#" conf/hibench.conf
	sed -i -e "s#^hibench.slaves.hostnames.*#hibench.slaves.hostnames localhost#" conf/hibench.conf
	sed -i -e "s#^hibench.scale.profile.*#hibench.scale.profile tiny#" conf/hibench.conf
	sed -i -e 's/.*hadoop/#hadoop/' conf/frameworks.lst
	cd $ROOTDIR
fi

