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
	sudo yum -y install cmake autoconf automake libtool cpupowerutils gcc libstdc++ gcc-c++ gcc-gfortran bzip2 patch zlib-devel ncurses-devel kernel-devel bc gawk coreutils grep perf glibc-static libstdc++-static ncurses-static libtool-ltdl-devel
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
	cd $JAVA_HOME/bin
	mv ./java ./java.org; ln -s ./java.org java
	cat <<EOF > ./java.sh
#!/bin/bash
trap cleanup2 2
trap cleanup9 9
trap cleanup11 11
trap cleanup15 15
cleanup2() { kill -2 "\$SUBC"; exit 0; }
cleanup9() { kill -9 "\$SUBC"; exit 0; }
cleanup11() { kill -11 "\$SUBC"; exit 0; }
cleanup15() { kill -15 "\$SUBC"; exit 0; }

CDIR="\$( cd "\$( dirname "\${BASH_SOURCE[0]}" )" && pwd )"
EXID="\$(echo \\"\$@\\" | awk 'match(\$0, "executor-id[[:blank:]]+([0-9]+)", m) {print m[1]}')"
SUBC=""

if [ "x\$EXID" != "x" ] && [ -f /dev/shm/RUNNINGWITHSDE ]; then
	SDEBIN=\$(cat /dev/shm/RUNNINGWITHSDE | head -1);
	SDEOUT=\$(cat /dev/shm/RUNNINGWITHSDE | tail -1);
	echo "\$SDEBIN -log -log:mt -log:basename \$SDEOUT/dcfg-out.executor-id-\$EXID" \\
		" -sse-sde -disasm_att 1 -align_checker_prefetch 0 -align_correct 0 -emu_fast 1 -bdw" \\
		" -- \$CDIR/java.org \$@" >> /dev/shm/java.call.log
	"\$CDIR"/java.org "\$@" & SUBC=\$!
	sleep 4
	\$SDEBIN -log -log:mt -log:basename "\$SDEOUT"/dcfg-out.executor-id-"\$EXID" \\
		-sse-sde -disasm_att 1 -align_checker_prefetch 0 -align_correct 0 -emu_fast 1 -bdw \\
		-attach-pid \$SUBC
	wait \$SUBC
else
	echo "\$CDIR/java.org \$@" >> /dev/shm/java.call.log
	"\$CDIR"/java.org "\$@" & SUBC=\$!
	wait \$SUBC
fi
EOF
	cd -
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
#	cd $SPARK_HOME
#	sed -i -e 's/.*exec "${CMD\[@\]}"/#exec "${CMD[@]}"/' bin/spark-class
#	cat <<EOF >> bin/spark-class
#if [ -z \$RUNNINGWITHSDE ] || [ -z \$ROOTDIR ]; then
#  exec "\${CMD[@]}"
#else
#  export PATH=\$ROOTDIR/dep/sde-external-8.35.0-2019-03-11-lin:\$PATH
#  LD_PRELOAD=\$ROOTDIR/HiBench/sde_java_hack.so \\
#  exec \`which sde64\` -follow_subprocess 1 -sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 \\
#                       -align_checker_prefetch 0 -align_correct 0 -emu_fast 1 \\
#                       -start_ssc_mark 111:repeat -stop_ssc_mark 222:repeat -bdw -- "\${CMD[@]}"
#fi
#EOF
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
	./bin/hdfs dfs -mkdir /user; ./bin/hdfs dfs -mkdir /user/`logname`; ./sbin/start-yarn.sh
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

BM="mpistub"
VERSION="e3b0d504e6daba0116907589c944a3be39547057"
cd $ROOTDIR/dep/$BM/
if [ ! -f $ROOTDIR/dep/$BM/lib/mpistub/libmpi.a ]; then
	if [[ "`hostname -s`" = *"peach"* ]]; then
		git checkout -b precision ${VERSION}
		git apply --check $ROOTDIR/patches/*1-${BM}*.patch
		if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
		module load FujitsuCompiler/202007
		rm -rf $ROOTDIR/dep/$BM/build; mkdir -p $ROOTDIR/dep/$BM/build; cd $ROOTDIR/dep/$BM/build
		if ! CC=fccpx CXX=FCCpx FC=frtpx cmake .. -DCMAKE_INSTALL_PREFIX=$ROOTDIR/dep/$BM/ ; then
			# peach has too old cmake
			source $ROOTDIR/dep/spack/share/spack/setup-env.sh
			spack install cmake@3.4.3; spack load cmake@3.4.3
			CC=fccpx CXX=FCCpx FC=frtpx cmake .. -DCMAKE_INSTALL_PREFIX=$ROOTDIR/dep/$BM/
		fi
		make
		make install
	fi
fi
cd $ROOTDIR

BM="gem5_riken"
VERSION="19103648cfd8c720128f26b63923551bd043d287"
cd $ROOTDIR/dep/$BM/
if [ ! -f $ROOTDIR/dep/$BM/build/ARM/gem5.opt ]; then
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	source $ROOTDIR/dep/spack/share/spack/setup-env.sh
	spack load gcc@8.4.0; spack load python
	if ! [[ "`hostname -s`" = *"peach"* ]]; then
		if ! which python3 >/dev/null 2>&1; then
			if ! [ `python -V 2>&1 | cut -d'.' -f2` -ge 7 ]; then
				wget https://www.python.org/ftp/python/2.7/Python-2.7.tar.bz2
				tar xjf Python-2.7.tar.bz2
				cd Python-2.7
				./configure --prefix=$ROOTDIR/dep/$BM/py27
				make -j; make install
				export PATH=$ROOTDIR/dep/$BM/py27/bin:$PATH
				export LD_LIBRARY_PATH=$ROOTDIR/dep/$BM/py27/lib:$LD_LIBRARY_PATH
				export SPACK_PYTHON=`which python2.7`
				cd -
			fi
			wget https://downloads.sourceforge.net/project/scons/scons/1.3.1/scons-1.3.1.tar.gz
			tar xzf scons-1.3.1.tar.gz
			cd scons-1.3.1
			python2 setup.py install
			cd -
		else
			python3 -m pip install --user --upgrade scons==3.1.2
		fi
	fi
	sed -i -e "s#PREFIX=/opt/riken_simulator#PREFIX=$ROOTDIR/dep/$BM#g" ./util/gem5-o3
	sed -i "369,372s:^:#:" ./SConstruct
	sed -i -e 's/ exit(/ sys.exit(/g' ./util/cpt_upgrader.py
	sed -i -e 's/if NO_FALLOCATE.*/if NO_FALLOCATE==0/' ./src/sim/syscall_emul.cc
	sed -i 's/typedef uint64_t SnoopMask;/typedef unsigned __int128 SnoopMask;/' ./src/mem/snoop_filter.hh
	sed -i '38 i std::ostream& operator<<(std::ostream& d, const unsigned __int128 v);' ./src/base/cprintf_formats.hh
	sed -i '41 i ostream& operator<<(ostream& d, const unsigned __int128 v) { d << "128int Hi:" << (void*)v << ";Lo:" << (void*)(v >> 64); return d;}' ./src/base/cprintf.cc
	if ! which python3 >/dev/null 2>&1; then
		echo "" | scons build/ARM/gem5.opt -j $(nproc)
	else
		echo "" | SCONS_LIB_DIR=`find $HOME/.local/lib -type d -name scons | head -1` scons build/ARM/gem5.opt -j $(nproc)
	fi
fi
cd $ROOTDIR

BM="sst"
if [ ! -f $ROOTDIR/dep/$BM/bin/sst-info ]; then
	cd $ROOTDIR/dep/$BM
	source $ROOTDIR/dep/spack/share/spack/setup-env.sh
	spack load gcc@8.4.0 ; spack load cmake@3.17.3
	# http://sst-simulator.org/SSTPages/SSTBuildAndInstall_11dot0dot0_SeriesAdditionalExternalComponents/#intel-pin-tool-317-98314
	wget https://software.intel.com/sites/landingpage/pintool/downloads/pin-3.17-98314-g0c048d619-gcc-linux.tar.gz
	tar xzf pin-3.17-98314-g0c048d619-gcc-linux.tar.gz
	export PIN_HOME=$ROOTDIR/dep/$BM/pin-3.17-98314-g0c048d619-gcc-linux
	export INTEL_PIN_DIRECTORY=$PIN_HOME
	# http://sst-simulator.org/SSTPages/SSTBuildAndInstall10dot1dot0SeriesDetailedBuildInstructions/
	# 'master' is the latest STABLE versions of SST
	#git clone -b master https://github.com/sstsimulator/sst-core.git
	#git clone -b master https://github.com/sstsimulator/sst-elements.git
	#git clone -b master https://github.com/sstsimulator/sst-macro.git
	#git clone -b master https://github.com/sstsimulator/sst-tools.git
	#git clone -b v10.1.0_Final https://github.com/sstsimulator/sst-external-element.git
	#git clone -b master https://github.com/sstsimulator/sst-tutorials
	#git clone -b 1.0.0 https://github.com/umd-memsys/DRAMsim3.git
	for SUB in sst-core sst-elements sst-macro sst-tools sst-external-element sst-tutorials DRAMsim3; do
		cd $ROOTDIR/dep/$BM/$SUB
		git apply --check $ROOTDIR/patches/*1-${SUB}*.patch
		if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${SUB}*.patch; fi
	done
	#
	cd $ROOTDIR/dep/$BM/DRAMsim3
	export DRAMSIM3_HOME=$ROOTDIR/dep/$BM/DRAMsim3
	mkdir -p build; cd build; cmake ..; make -j
	#
	cd $ROOTDIR/dep/$BM/sst-core
	export SST_CORE_HOME=$ROOTDIR/dep/$BM
	export SST_CORE_ROOT=$ROOTDIR/dep/$BM/sst-core
	./autogen.sh
	./configure --prefix=$SST_CORE_HOME --disable-mpi
	make all -j V=1 && make install
	export PATH=$SST_CORE_HOME/bin:$PATH
	#
	cd $ROOTDIR/dep/$BM/sst-elements
	export SST_ELEMENTS_HOME=$ROOTDIR/dep/$BM
	export SST_ELEMENTS_ROOT=$ROOTDIR/dep/$BM/sst-elements
	./autogen.sh
	./configure --prefix=$SST_ELEMENTS_HOME --with-sst-core=$SST_CORE_HOME --with-pin=$PIN_HOME --with-dramsim3=$DRAMSIM3_HOME
	make all -j V=1 && make install
	cd $ROOTDIR/dep/$BM/sst-external-element/src
	make
	# for loading
	export SST_CORE_HOME=$ROOTDIR/dep/$BM
	export PATH=$SST_CORE_HOME/bin:$PATH
	export LD_LIBRARY_PATH=$SST_CORE_HOME/lib:$LD_LIBRARY_PATH
	if ! which sst-info >/dev/null 2>&1; then
		echo "ERR: something went wrong in SST install"
	else
		echo "INFO: check if we have all we need: http://sst-simulator.org/SSTPages/SSTBuildAndInstall10dot1dot0SeriesDetailedBuildInstructions/#configuration-selection-for-internal-elements-and-optionalrequired-external-components"
	fi
fi
cd $ROOTDIR

echo -e '\nInit OSACA'
BM="OSACA"
VERSION="v0.4.1"
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

