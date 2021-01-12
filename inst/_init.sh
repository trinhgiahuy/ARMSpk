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
	sudo yum -y install cmake autoconf automake libtool cpupowerutils gcc libstdc++ gcc-c++ gcc-gfortran bzip2 patch zlib-devel ncurses-devel kernel-devel bc gawk coreutils grep perf
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

echo -e '\nInit spack and new LLVM'
BM="spack"
VERSION="96fa6f0c1be4ab55ec6ba7cd5af059e1ed95351f"
cd $ROOTDIR/dep/$BM/
if [[ "`git branch`" = *"develop"* ]]; then
        git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
        source $ROOTDIR/dep/spack/share/spack/setup-env.sh
        spack compiler find
        # check system compiler and install one we like
        if [[ "`gcc --version | /usr/bin/grep 'gcc (GCC)' | cut -d ' ' -f3`" = "8.4.0" ]]; then
                spack install gcc@8.4.0
                spack load gcc@8.4.0
                spack compiler find
        fi
	spack install libpfm4@4.10.1%gcc@8.4.0
	spack load libpfm4
	spack install llvm@10.0.0%gcc@8.4.0 gold=False libpfm=True

	#spack install openjdk@1.8.0_222-b10%gcc@8.4.0
	#spack install maven@3.6.3%gcc@8.4.0 ^openjdk@1.8.0_222-b10
	#spack install scala@2.11.11%gcc@8.4.0 ^openjdk@1.8.0_222-b10
	#spack install hadoop@3.2.1%gcc@8.4.0 ^openjdk@1.8.0_222-b10
	#sed -i "/a7e29e78bd43aa6d137f0bb0afd54a3017865d471456c6d436ae79475bbeb161/i \    version('2.4.0', sha256='b1d6d6cb49d8253b36df8372a722292bb323bd16315d83f0b0bafb66a4154ef2')" $ROOTDIR/dep/spack/var/spack/repos/builtin/packages/spark/package.py
	#spack install spark@2.4.0%gcc@8.4.0 ^openjdk@1.8.0_222-b10
	#spack load maven; spack load scala; spack load hadoop; spack load spark
	#mvn -Dhadoop=3.2 -Dspark=2.4 -Dscala=2.12 clean package
	#cd `spack find -p | /bin/grep hadoop | cut -d' ' -f2-`; sed -i '/<\/configuration>/i \  <property>\n    <name>fs.default.name</name>\n    <value>hdfs://localhost:8020</value>\n  </property>' etc/hadoop/core-site.xml
	#./sbin/stop-all.sh
	#rm -rf /tmp/*
	#echo 'Y' | ./bin/hdfs namenode -format
	#./sbin/start-all.sh
	#cd `spack find -p | /bin/grep spark | cut -d' ' -f2-`
	#export SPARK_DIST_CLASSPATH=`hadoop classpath`
	#./sbin/start-master.sh ; sbin/start-slave.sh spark://`hostname`:7077 ; ###./sbin/start-all.sh not working
	#cd HiBench
	#cp conf/hadoop.conf.template conf/hadoop.conf
	#cp conf/spark.conf.template conf/spark.conf
	#sed -i -e "s#/PATH/TO/YOUR/HADOOP/ROOT#`spack find -p | /bin/grep hadoop | cut -d' ' -f2-`#" conf/hadoop.conf
	#sed -i -e "s#/PATH/TO/YOUR/SPARK/HOME#`spack find -p | /bin/grep spark | cut -d' ' -f2-`#" conf/spark.conf
	#sed -i -e "s#^hibench.spark.master.*#hibench.spark.master spark://`hostname`:7077#" conf/spark.conf
	#sed -i -e "s#^hibench.masters.hostnames.*#hibench.masters.hostnames localhost#" conf/hibench.conf
	#sed -i -e "s#^hibench.slaves.hostnames.*#hibench.slaves.hostnames localhost#" conf/hibench.conf
	#./bin/workloads/micro/wordcount/prepare/prepare.sh
	#./bin/workloads/micro/wordcount/hadoop/run.sh
	#./bin/workloads/micro/wordcount/spark/run.sh
fi
cd $ROOTDIR

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

