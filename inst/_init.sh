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
	make
	for x in `ls likwid-*`; do if [ -x $x ]; then sudo setcap cap_sys_admin,cap_sys_rawio+ep $x; fi; done
	cat /proc/cmdline | grep 'intel_pstate=disable'  > /dev/null
	# vim /etc/default/grub
	# grub2-mkconfig -o /boot/grub2/grub.cfg
	# grub2-mkconfig -o /boot/efi/EFI/fedora/grub.cfg
	# reboot
	if [ ! "x$?" = "x0" ]; then echo -e 'Note: for likwid to work, please add 'intel_pstate=disable' to kernel parameter and reboot; Afterwards, run this script again'; fi
	echo -e "Please execute:\n  export PATH=$ROOTDIR/dep/$BM:\$PATH\n  export LD_LIBRARY_PATH=$ROOTDIR/dep/$BM:\$LD_LIBRARY_PATH"
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
sudo chmod go+rw /dev/mem"

