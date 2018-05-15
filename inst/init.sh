#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source `cat $ROOTDIR/conf/intel.cfg` intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort
alias ar=`which xiar`
alias ld=`which xild`

echo '\nATTENTION!!! This script will ask for sudo/root access\nNote: MSR changes are not persistent.\n      Installing in NFS will not work when mounted with nosuid flag.\n\n(10s time to abort and manually install)\n'
sleep 10

echo '\nInstalling Likwid'
BM="dep/likwid"
VERSION="4e1a04eed1371f82d04eb9c1d1d739706a633b4a"
if [ ! -f $ROOTDIR/$BM/likwid-setFrequencies ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	sed -i -e 's/GCC#NO/ICC#NO/' -e 's/accessdaemon#NO/direct#NO/' ./config.mk
	make
	for x in `ls likwid-*`; do if [ -x $x ]; then sudo setcap cap_sys_admin,cap_sys_rawio+ep $x; fi; done
	cat /proc/cmdline | grep 'intel_pstate=disable'  > /dev/null
	# vim /etc/default/grub
	# grub2-mkconfig -o /boot/grub2/grub.cfg
	# grub2-mkconfig -o /boot/efi/EFI/fedora/grub.cfg
	# reboot
	if [ ! "x$?" = "x0" ]; then echo "Note: for likwid to work, please add 'intel_pstate=disable' to kernel parameter and reboot; Afterwards, run this script again"; fi
	echo "Please execute:\n  export PATH=$ROOTDIR/$BM:\$PATH\n  export LD_LIBRARY_PATH=$ROOTDIR/$BM:\$LD_LIBRARY_PATH"
	cd $ROOTDIR
fi

echo '\nInstalling Intel PCM'
BM="dep/intel-pcm"
VERSION="33ba100b7694c5f3aecbb61dbc82507daa6c5b74"
if [ ! -f $ROOTDIR/$BM/pcm-memory.x ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am < $ROOTDIR/patches/*1-${BM}*.patch; fi
	make CXX=icpc
	for x in `ls *.x`; do if [ -x $x ]; then sudo setcap cap_sys_admin,cap_sys_rawio+ep $x; fi; done
	cd $ROOTDIR
fi

# enable support for MSR and MSR-safe counters from user space
echo '\nInstalling MSR/MSR-Safe'
BM="dep/msr-safe"
VERSION="b5bdf8b200db5a0bfa7e9ba2aadb85159f72c697"
if [ ! -f $ROOTDIR/$BM/msr-safe.ko ] || [ ! -r /dev/cpu/0/msr ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	make
	WL=`printf 'wl_%.2x%x\n' $(lscpu | grep "CPU family:" | awk -F: '{print $3}') $(lscpu | grep "Model:" | awk -F: '{print $3}')`
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

