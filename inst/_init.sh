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

echo -e '\nATTENTION!!! This script will ask for sudo/root access\nNote: MSR changes are not persistent.\n      Installing in NFS will not work when mounted with nosuid flag.\n\n(10s time to abort and manually install)\n'
sleep 10

echo -e '\nInstalling known dependencies'
if [[ -f /etc/redhat-release ]];then
	sudo yum -y install cmake autoconf automake libtool cpupowerutils gcc libstdc++ gcc-c++ gcc-gfortran bzip2 patch zlib-devel ncurses-devel kernel-devel bc gawk coreutils grep perf glibc-static libstdc++-static ncurses-static
	# vim screen
	# for Intel Parallel XE: gtk2 gtk3 pango xorg-x11-server-Xorg
else
	echo -e '\nNote: untested linux distro; please install -- cmake autoconf automake libtool cpupowerutils gcc libstdc++ gcc-c++ gcc-gfortran bzip2 patch zlib-devel ncurses-devel kernel-devel bc gawk coreutils grep perf -- yourself;\n      everything hereafter may break, no guarantees ;-)'
	sleep 10
fi

echo -e '\nUpdating ldconfig'
echo $LD_LIBRARY_PATH | sed -e 's/:/\n/g' > /dev/shm/precision.conf
sudo mv /dev/shm/precision.conf /etc/ld.so.conf.d/
sudo ldconfig

echo -e '\nFix problem with Vtune getting stuck it seems'
ssh -o StrictHostKeyChecking=no `hostname` echo ''

echo -e '\nDisabling NMI watchdog'
sudo sh -c "echo 0 > /proc/sys/kernel/nmi_watchdog"
sudo sh -c "echo 'kernel.nmi_watchdog=0' >> /etc/sysctl.conf"

echo -e '\nRemoving broken Intel VTune Sampling Drivers and use perf'
cd ${VTUNE_AMPLIFIER_2018_DIR}/sepdk/src
sudo ./rmmod-sep -R
sudo ./boot-script --uninstall
sudo sh -c 'echo 0 > /proc/sys/kernel/perf_event_paranoid'
#sudo ./insmod-sep -r -g `whoami` -pu
#sudo ./boot-script --install -g `whoami` -pu

echo -e '\nInstalling Likwid'
BM="likwid"
VERSION="4e1a04eed1371f82d04eb9c1d1d739706a633b4a"
if [ ! -f $ROOTDIR/dep/$BM/likwid-setFrequencies ]; then
	cd $ROOTDIR/dep/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	sed -i -e "s# /usr/local\#NO# $ROOTDIR/dep/$BM\#/usr/local\#NO#" ./config.mk
	make
	sudo make install
	for x in `ls bin/likwid-*`; do if [ -x $x ]; then sudo setcap cap_sys_admin,cap_sys_rawio+ep $x; fi; done
	cat /proc/cmdline | grep 'intel_pstate=disable'  > /dev/null
	# vim /etc/default/grub
	# grub2-mkconfig -o /boot/grub2/grub.cfg
	# grub2-mkconfig -o /boot/efi/EFI/centos/grub.cfg
	# reboot
	if [ ! "x$?" = "x0" ]; then echo -e 'Note: for likwid to work, please add 'intel_pstate=disable' to kernel parameter and reboot; Afterwards, run this script again'; fi
	echo -e "Please execute:\n  export PATH=$ROOTDIR/dep/$BM/bin:\$PATH\n  export LD_LIBRARY_PATH=$ROOTDIR/dep/$BM/lib:\$LD_LIBRARY_PATH"
	cd $ROOTDIR
fi

echo -e '\nInstalling Intel PCM'
BM="intel-pcm"
VERSION="33ba100b7694c5f3aecbb61dbc82507daa6c5b74"
if [ ! -f $ROOTDIR/dep/$BM/pcm-memory.x ]; then
	cd $ROOTDIR/dep/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	# many counters just zero if pcm compiled with perf
	sed -i -e 's/-DPCM_USE_PERF/#-DPCM_USE_PERF/' ./Makefile
	# no KNM suport yet, so "fake" it and hope for the best
	if [[ $HOSTNAME = *"${IKNMHOST}"* ]]; then
		sed -i -e 's/KNL = 87/KNL = 133/' ./cpucounters.h
	fi
	make -j CXX=icpc
	for x in `ls *.x`; do if [ -x $x ]; then sudo setcap cap_sys_admin,cap_sys_rawio+ep $x; fi; done
	cd $ROOTDIR
fi

# enable support for MSR and MSR-safe counters from user space
echo -e '\nInstalling MSR/MSR-Safe'
BM="msr-safe"
VERSION="b5bdf8b200db5a0bfa7e9ba2aadb85159f72c697"
WL="`printf 'wl_%.2x%x\n' $(lscpu | awk '/CPU family:/{print $3}') $(lscpu | awk '/Model:/{print $2}')`"
if [ ! -f $ROOTDIR/dep/$BM/msr-safe.ko ] || [ ! -r /dev/cpu/0/msr ]; then
	cd $ROOTDIR/dep/$BM/
	git checkout -b precision ${VERSION}
	make
	# use KNL whitelist for KNM nodes
	if [ ! -e ./whitelists/wl_0685 ]; then cp ./whitelists/wl_0657 ./whitelists/wl_0685; fi
	sudo insmod msr-safe.ko
	sudo chmod go+rw /dev/cpu/msr_whitelist
	sudo cat ./whitelists/${WL} > /dev/cpu/msr_whitelist
	sudo touch /sys/firmware/acpi/tables/MCFG
	sudo chmod go+rw /sys/firmware/acpi/tables/MCFG
	sudo chmod go+rw /dev/cpu/*/msr
	sudo chmod go+rw /dev/cpu/*/msr_safe
	sudo chmod go+rw /proc/bus/pci/*/*.*
	sudo chmod go+rw /dev/mem
	cd $ROOTDIR
fi

echo -e '\nInit spack and new LLVM and MPI alternative'
BM="spack"
VERSION="96fa6f0c1be4ab55ec6ba7cd5af059e1ed95351f"
cd $ROOTDIR/dep/$BM/
if [[ "`git branch`" = *"develop"* ]]; then
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	source $ROOTDIR/dep/$BM/share/spack/setup-env.sh
	spack compiler find
	# check system compiler and install one we like
	if [[ "`gcc --version | /usr/bin/grep 'gcc (GCC)' | cut -d ' ' -f3`" = "8.4.0" ]]; then
		spack install gcc@8.4.0
		spack load gcc@8.4.0
		spack compiler find
	fi
	# llvm
	spack install libpfm4@4.10.1%gcc@8.4.0
	spack load libpfm4
	spack install llvm@10.0.0%gcc@8.4.0 gold=False libpfm=True
	# openmpi
	#...issue with static external hwloc and no static verbs
	sed -i -e "s/with-hwloc=.*/with-hwloc=internal')/" -e "s/depends_on('hwloc/#depends_on('hwloc/" $ROOTDIR/dep/$BM/var/spack/repos/builtin/packages/openmpi/package.py
	spack install openmpi@3.1.6%gcc@8.4.0 fabrics=none thread_multiple=True vt=False cxx=True cxx_exceptions=False ^numactl@2.0.12/`spack find -l numactl@2.0.12%gcc@8.4.0 | /bin/grep numactl | cut -d' ' -f1`
	spack install openmpi@3.1.6%intel@19.0.1.144 fabrics=none thread_multiple=True vt=False cxx=True cxx_exceptions=False
	# hadoop/spark (reinstall spark from scratch because prebuild has no hive integration)
	spack install openjdk@1.8.0_222-b10%gcc@8.4.0
	spack install maven@3.6.3%gcc@8.4.0 ^openjdk@1.8.0_222-b10
	spack install scala@2.11.11%gcc@8.4.0 ^openjdk@1.8.0_222-b10
	spack install hadoop@2.10.0%gcc@8.4.0 ^openjdk@1.8.0_222-b10
	spack load openjdk; spack load maven; spack load scala; spack load hadoop
	export HADOOP_HOME=`spack find -p | /bin/grep hadoop | cut -d' ' -f2- | tr -d ' '`
	export HADOOP_CONF_DIR=$HADOOP_HOME/etc/hadoop
	export HADOOP_COMMON_LIB_NATIVE_DIR=$HADOOP_HOME/lib/native
	export LD_LIBRARY_PATH=$HADOOP_COMMON_LIB_NATIVE_DIR:$LD_LIBRARY_PATH
	export HADOOP_LOG_DIR=/scr0/hadoop-logs; export YARN_LOG_DIR=$HADOOP_LOG_DIR
	export HADOOP_MAPRED_LOG_DIR=$HADOOP_LOG_DIR; export SPARK_LOG_DIR=$HADOOP_LOG_DIR
	export MAVEN_OPTS="-Xmx4g -XX:ReservedCodeCacheSize=1024m"
	sed -i "/a7e29e78bd43aa6d137f0bb0afd54a3017865d471456c6d436ae79475bbeb161/i \    version('2.4.0', sha256='b1d6d6cb49d8253b36df8372a722292bb323bd16315d83f0b0bafb66a4154ef2')" $ROOTDIR/dep/spack/var/spack/repos/builtin/packages/spark/package.py
	spack install spark@2.4.0%gcc@8.4.0 hadoop=True ^hadoop@2.10.0%gcc@8.4.0 ^openjdk@1.8.0_222-b10
	cd /tmp; wget https://archive.apache.org/dist/spark/spark-2.4.0/spark-2.4.0.tgz
	tar xzf spark-2.4.0.tgz; cd spark-2.4.0
	./dev/make-distribution.sh --name custom-spark \
		--pip --tgz -Phadoop-provided -Phive -Phive-thriftserver -Pyarn -DskipTests
	export SPARK_HOME=`spack find -p | /bin/grep spark | cut -d' ' -f2- | tr -d ' '`
	export BN_SPARK_HOME=`basename $SPARK_HOME`
	cd $SPARK_HOME/../; mv $BN_SPARK_HOME .$BN_SPARK_HOME; mkdir $BN_SPARK_HOME
	tar xzf /tmp/spark-2.4.0/spark-2.4.0-bin-custom-spark.tgz --strip-components=1 -C $BN_SPARK_HOME
	spack load spark
	# config hadoop/spark for pseudo-distributed cluster mode (=> hibench/hadoop in local mode will not work https://github.com/Intel-bigdata/HiBench/issues/120; need hadoop even in spark-only mode for data preparation)
	cd $HADOOP_HOME
	sed -i '/<\/configuration>/i <property>\n<name>fs.defaultFS<\/name>\n<value>hdfs:\/\/localhost:9000<\/value>\n<\/property>\n<property>\n<name>hadoop.tmp.dir<\/name>\n<value>\/scr0\/hadoop-${user.name}<\/value>\n<\/property>' etc/hadoop/core-site.xml
	sed -i '/<\/configuration>/i <property>\n<name>dfs.replication<\/name>\n<value>1<\/value>\n<\/property>' etc/hadoop/hdfs-site.xml
	sed '/<\/configuration>/i <property>\n<name>mapreduce.framework.name<\/name>\n<value>yarn<\/value>\n<\/property>\n<property>\n<name>mapreduce.map.memory.mb<\/name>\n<value>2048<\/value>\n<\/property>\n<property>\n<name>mapreduce.reduce.memory.mb<\/name>\n<value>4096<\/value>\n<\/property>\n<property>\n<name>mapreduce.map.java.opts<\/name>\n<value>-Xmx1638m<\/value>\n<\/property>\n<property>\n<name>mapreduce.reduce.java.opts<\/name>\n<value>-Xmx3278m<\/value>\n<\/property>' etc/hadoop/mapred-site.xml.template > etc/hadoop/mapred-site.xml
	sed -i '/<\/configuration>/i <property>\n<name>yarn.nodemanager.aux-services<\/name>\n<value>mapreduce_shuffle<\/value>\n<\/property>\n<property>\n<name>yarn.nodemanager.disk-health-checker.enable<\/name>\n<value>false<\/value>\n<\/property>\n<property>\n<name>yarn.nodemanager.local-dirs<\/name>\n<value>\/scr0\/hadoop-${user.name}\/nm-local-dir<\/value>\n<\/property>\n<property>\n<name>yarn.nodemanager.log-dirs<\/name>\n<value>\/scr0\/hadoop-${user.name}\/containers<\/value>\n<\/property>\n<property>\n<name>yarn.nodemanager.remote-app-log-dir<\/name>\n<value>\/scr0\/hadoop-${user.name}\/apps<\/value>\n<\/property>\n<property>\n<name>yarn.nodemanager.resource.detect-hardware-capabilities<\/name>\n<value>true<\/value>\n<\/property>\n<property>\n<name>yarn.scheduler.maximum-allocation-mb<\/name>\n<value>MAXALLOC<\/value>\n<\/property>' etc/hadoop/yarn-site.xml
	sed -i -e "s/MAXALLOC/$(echo "`free -m | /bin/grep 'Mem:' | awk -F'[^0-9]*' '$0=$2'` / 2" | bc)/" etc/hadoop/yarn-site.xml
	sed -i -e "s#export JAVA_HOME.*#export JAVA_HOME=$JAVA_HOME#" -e "s#.*export HADOOP_LOG_DIR.*#export HADOOP_LOG_DIR=$HADOOP_LOG_DIR#" -e "s#.*export HADOOP_SECURE_DN_LOG_DIR.*#export HADOOP_SECURE_DN_LOG_DIR=$HADOOP_LOG_DIR#" -e "s#.*export HADOOP_PID_DIR.*#export HADOOP_PID_DIR=$HADOOP_LOG_DIR#" -e '/^export HADOOP_OPTS/a export HADOOP_NAMENODE_OPTS=-Xmx8g\nexport HADOOP_DATANODE_OPTS=-Xmx8g\nexport HADOOP_SECONDARYNAMENODE_OPTS=-Xmx8g' etc/hadoop/hadoop-env.sh
	sed -i -e "s#.*export YARN_NODEMANAGER_OPTS=.*#export YARN_NODEMANAGER_OPTS=-Xmx4g#" -e "s#.*export YARN_RESOURCEMANAGER_OPTS=.*#export YARN_RESOURCEMANAGER_OPTS=-Xmx4g#" etc/hadoop/yarn-env.sh
	sed -i -e 's/logger=INFO/logger=FATAL/' etc/hadoop/log4j.properties
	# run test
	cd $HADOOP_HOME; echo 'N' | ./bin/hdfs namenode -format; ./sbin/start-dfs.sh
	./bin/hdfs dfs -mkdir /user; ./bin/hdfs dfs -mkdir /user/jens; ./sbin/start-yarn.sh
	cd $SPARK_HOME; ./sbin/start-master.sh ; ./sbin/start-slave.sh spark://`hostname`:7077
	# shutdown again
	cd $SPARK_HOME; ./sbin/stop-slaves.sh; ./sbin/stop-master.sh
	cd $HADOOP_HOME; ./sbin/stop-yarn.sh; ./sbin/stop-dfs.sh; killall -9 java
	spack unload spark; spack unload hadoop; spack unload scala; spack unload maven; spack unload openjdk; spack unload libpfm4
fi
cd $ROOTDIR

#BM=""
#VERSION=""
#echo -e "\nInit $BM"
#if []; then
#	git clone https://github.com/GoogleCloudPlatform/PerfKitBenchmarker.git
#	git checkout v1.15.1
#	sed -i -e 's/^PyYAML==.*/PyYAML==5.2/' ./requirements.txt
#	python3 -m pip install --user --upgrade -r requirements.txt
#	ssh-keygen -P '' -f ~/.ssh/id_rsa_perfkit
#	cat ~/.ssh/id_rsa_perfkit.pub >> ~/.ssh/authorized_keys
#	#on kiev1 and kiev2 as ROOT (XXX)
#	#echo "jens ALL=(ALL) NOPASSWD: /usr/sbin/fdisk -l, /usr/sbin/sysctl vm.drop_caches=3, /usr/bin/mkdir -p /opt/pkb, /usr/bin/chmod a+rwxt /opt/pkb, /usr/sbin/sync, /usr/bin/tee" >> /etc/sudoers.d/pccstaff
#	#yum -y install iperf
#	mkdir -p /scr0/jens/input; mkdir -p /scr0/jens/output
#	#sed -i -e 's/^DEFAULT =.*/DEFAULT = CENTOS7/' perfkitbenchmarker/os_types.py
#	cp ~/.ssh/id_rsa_perfkit.pub /dev/shm
#	docker build --rm -t perfkit .
#	docker run --privileged -ti -p 2222:22 -p 3003:3003 -p 12865:12865 -p 20000:20000 -v /dev/shm/id_rsa_perfkit.pub:/root/.ssh/authorized_keys -v /scr0/jens/input:/data0 -v /scr0/jens/output:/data1 -e SSH_ENABLE_ROOT=true -e MOTD='' perfkit
#	ssh root@localhost -p2222 -i ~/.ssh/id_rsa_perfkit -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null
#fi

echo -e '\nInit OSACA'
BM="OSACA"
VERSION="768a90de103755fa995c9f4e23e1f498e763aff2"
cd $ROOTDIR/dep/$BM/
if [[ "`git branch`" = *"master"* ]]; then
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	#get arch info in different way, since pmu_name is intel-only
	python3 -m pip install --user --upgrade archspec==0.1.2
	python3 -m pip install --user --upgrade pyparsing==2.4.2
	python3 -m pip install --user --upgrade ruamel.yaml==0.16.12
	python3 -m pip install --user --upgrade ruamel.yaml.clib==0.2.2
	python3 -m pip install --user --upgrade networkx==2.4
	python3 -m pip install --user --upgrade decorator==4.4.0
	python3 setup.py bdist_wheel 2>&1 | tee comp
	python3 -m pip install --user --upgrade dist/osaca-0.3.14-py3-none-any.whl
fi
cd $ROOTDIR

echo -e '\nInit IACA'
# IACA https://software.intel.com/content/www/us/en/develop/articles/intel-architecture-code-analyzer.html
if [ -f $ROOTDIR/dep/iaca-version-v3.0-lin64.zip ]; then
	cd $ROOTDIR/dep/; rm -rf iaca-lin64/; unzip ./iaca-version-v3.0-lin64.zip; cd -
else
	echo "ERR: missing ./iaca-version-v3.0-lin64.zip in dep/ folder"
fi

echo -e "\nIn case of reboot, run again:
sudo sh -c 'echo 0 > /proc/sys/kernel/perf_event_paranoid'
sudo sh -c 'echo 0 > /proc/sys/kernel/nmi_watchdog'
cd $ROOTDIR/dep/msr-safe
sudo insmod msr-safe.ko
sudo chmod go+rw /dev/cpu/msr_whitelist
sudo cat ./whitelists/${WL} > /dev/cpu/msr_whitelist
sudo touch /sys/firmware/acpi/tables/MCFG
sudo chmod go+rw /sys/firmware/acpi/tables/MCFG
sudo chmod go+rw /dev/cpu/*/msr
sudo chmod go+rw /dev/cpu/*/msr_safe
sudo chmod go+rw /proc/bus/pci/*/*.*
sudo chmod go+rw /dev/mem
sudo tuned-adm off
sleep 5
sudo tuned-adm profile latency-performance
tuned-adm verify
export PATH=$ROOTDIR/dep/likwid/bin:\$PATH
export LD_LIBRARY_PATH=$ROOTDIR/dep/likwid/lib:\$LD_LIBRARY_PATH
#Xeon: likwid-setFrequencies -g performance --freq 2.2 --turbo 1 --umin 2.7 --umax 2.7
#KNL:  likwid-setFrequencies -g performance --freq 1.301 --turbo 1
#KNM:  likwid-setFrequencies -g performance --freq 1.501 --turbo 1
likwid-setFrequencies -p"

