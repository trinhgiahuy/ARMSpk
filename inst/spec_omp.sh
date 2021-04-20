#!/bin/bash

function dump_omp_config {
cat <<EOF > $1
action         = validate
delay          = 1
iterations     = 1
output_format  = text
teeout         = no
teerunout      = yes
tune           = peak
verbose        = 99
env_vars       = 1
makeflags      = -j 32
ext            = %{COMP}
%if '%{COMP}' eq 'fuji'
check_md5      = 0
%endif

%if '%{COMP}' eq 'gnu'
submit         = ulimit -n 4096; ulimit -s unlimited; \$command
BOPTS          = -O3 -fopenmp -march=native -funroll-loops -ffast-math -ftree-vectorize
BLINK          = -static
%elif '%{COMP}' eq 'intel'
submit         = ulimit -n 4096; ulimit -s unlimited; \$command
BOPTS          = -O3 -qopenmp -xHOST -no-prec_div -fp-model fast=2 -fma
BLINK          = -static -static-intel -qopenmp-link=static
%elif '%{COMP}' eq 'sde'
submit         = ulimit -n 4096; ulimit -s unlimited; sde64 -sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 -dcfg:out_base_name %{RESDIR}/dcfg-out.wl-\${workload} -align_checker_prefetch 0 -align_correct 0 -emu_fast 1 -bdw -- \$command
BOPTS          = -O3 -qopenmp -xHOST -no-prec_div -fp-model fast=2 -fma
BLINK          = -static -static-intel -qopenmp-link=static
%elif '%{COMP}' eq 'fuji'
%if '%{HOST}' eq 'peach'
submit         = ulimit -n 4096; ulimit -s unlimited; bash %{RESDIR}/gem5wrap.sh \$command
BOPTS          = -Kfast,eval_concurrent,openmp -O3 -march=armv8.3-a+sve -funroll-loops
BLINK          = -Bstatic -lm
%else
submit         = ulimit -s unlimited; \$command
BOPTS          = -Nclang -Kfast,eval_concurrent,openmp -O3 -march=armv8.3-a+sve -funroll-loops
%endif
%else
%error wrong or unsupported COMP variable specified
%endif

preENV_OMP_STACKSIZE=256M
preENV_OMP_DYNAMIC=false
preENV_OMP_NESTED=false
preENV_KMP_AFFINITY=balanced
preENV_KMP_BLOCKTIME=infinite
preENV_KMP_LIBRARY=turnaround
preENV_KMP_SCHEDULE=static,balanced
preENV_KMP_STACKSIZE=256M
%ifndef %{RESDIR}
%error please define RESDIR for launch
%endif

default=default=default=default:
%if '%{COMP}' eq 'gnu'
sw_compiler    = GNU Compilers
FC             = gfortran
CC             = gcc
CXX            = g++
%elif '%{COMP}' eq 'intel' || '%{COMP}' eq 'sde'
sw_compiler    = Intel Compiler Suite
FC             = ifort
CC             = icc
CXX            = icpc
%elif '%{COMP}' eq 'fuji'
sw_compiler    = FUJITSU Software Technical Computing Suite
%if '%{HOST}' eq 'peach'
FC             = frtpx
CC             = fccpx
CXX            = FCCpx
%else
FC             = frt
CC             = fcc
CXX            = FCC
%endif
%endif
FOPTIMIZE      = \${BOPTS}
COPTIMIZE      = \${BOPTS}
CXXOPTIMIZE    = \${BOPTS}
OS_LIBS        = \${BLINK}

350.md=default=default=default:
%if '%{COMP}' eq 'gnu'
FPORTABILITY   = -ffree-form -fno-range-check
%elif '%{COMP}' eq 'fuji'
FPORTABILITY   = -Free
%else
FPORTABILITY   = -free
%endif

357.bt331=default=default=default:
%if '%{COMP}' eq 'intel' || '%{COMP}' eq 'sde'
FPORTABILITY   = -mcmodel=medium
LDOPT          = -shared-intel
PASS1_LDOPT    = -shared-intel
%elif '%{COMP}' eq 'fuji'
FPORTABILITY   = -mcmodel=large
%endif

359.botsspar=default=default=default:
%if '%{COMP}' eq 'intel' || '%{COMP}' eq 'sde'
strict_rundir_verify = 0
%endif

363.swim=default=default=default:
%if '%{COMP}' eq 'intel' || '%{COMP}' eq 'sde'
FPORTABILITY   = -mcmodel=medium
LDOPT          = -shared-intel
PASS1_LDOPT    = -shared-intel
%elif '%{COMP}' eq 'fuji'
FPORTABILITY   = -mcmodel=large
%endif

367.imagick=default=default=default:
CPORTABILITY   = -std=c99

372.smithwa=default=default=default:
%if '%{COMP}' eq 'fuji'
CPORTABILITY   = -fsigned-char
%endif
EOF
}

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
if [ -z $1 ]; then
	source $ROOTDIR/conf/intel.cfg
	source $INTEL_PACKAGE intel64 > /dev/null 2>&1
	export I_MPI_CC=icc
	export I_MPI_CXX=icpc
	export I_MPI_F77=ifort
	export I_MPI_F90=ifort
	alias ar=`which xiar`
	alias ld=`which xild`
	export ADVISOR_2018_DIR=${ADVISOR_2019_DIR}
elif [[ "$1" = *"gnu"* ]]; then
	source $ROOTDIR/dep/spack/share/spack/setup-env.sh
	spack load gcc@8.4.0
elif [[ "`hostname -s`" = *"fn01"* ]] && [[ "$1" = *"fuji"* ]]; then
	echo "does not compile on login node"; exit 1
elif [[ "$1" = *"fuji"* ]]; then
	if [[ "`hostname -s`" = *"peach"* ]]; then module load FujitsuCompiler/202007; fi
else
	echo 'wrong compiler'
	exit 1
fi

BM="SPEC_OMP"
if [ ! -f $ROOTDIR/$BM/bin/runcpu ]; then
	mkdir -p $ROOTDIR/$BM/
	if [ ! -f $ROOTDIR/dep/mnt_$BM/shrc ]; then
		if [ ! -f $ROOTDIR/dep/omp2012-1.1.iso ]; then
			echo -e "ERR: cannot find omp2012-1.1.iso under dep/; please fix it!"
			exit 1
		fi
		cd $ROOTDIR/dep/
		mkdir -p /tmp/mnt_$BM
		fuseiso ./omp2012-1.1.iso /tmp/mnt_$BM
		cd /tmp/mnt_$BM/
		if [[ "$1" = *"fuji"* ]] && ! [[ "`hostname -s`" = *"peach"* ]]; then
			./install.sh -f -d $ROOTDIR/$BM/ -u linux-apm-arm64
		else
			./install.sh -f -d $ROOTDIR/$BM/ -u linux-suse10-amd64
		fi
		cd -
		fusermount -u /tmp/mnt_$BM
	else
		cd $ROOTDIR/dep/mnt_$BM
		if [[ "$1" = *"fuji"* ]] && ! [[ "`hostname -s`" = *"peach"* ]]; then
			chmod -R +w $ROOTDIR/dep/mnt_$BM
			for FILE in `find . -name config.guess`; do wget 'http://savannah.gnu.org/cgi-bin/viewcvs/*checkout*/config/config/config.guess' -O $FILE; done
			FILE="tools/src/make-3.82/glob/glob.c"; sed -i -e 's@#if !defined __alloca && !defined __GNU_LIBRARY__@#if !defined __alloca@g' $FILE
			FILE="tools/src/perl-5.12.3/Configure"; sed -i -e 's@optimize="$ans"@optimize="-O1"@g' $FILE
			FILE="tools/src/specsum/gnulib/stdio.in.h"; sed -i -e '/use fgets instead/i #if defined(__GLIBC__) && !defined(__UCLIBC__) && !__GLIBC_PREREQ(2, 16)' -e '/use fgets instead/a #endif' $FILE
			FILE="tools/src/tar-1.25/gnu/stdio.in.h"; sed -i -e '/use fgets instead/i #if defined(__GLIBC__) && !defined(__UCLIBC__) && !__GLIBC_PREREQ(2, 16)' -e '/use fgets instead/a #endif' $FILE
			FILE="tools/src/TimeDate-1.20/t/getdate.t"; sed -i -e 's@timegm(0,0,0,1,0,70)@timegm(0,0,0,1,0,1970)@g' $FILE
			sed -i -e 's@MAKEFLAGS= $MYMAKE check; testordie "error testing tar"@#MAKEFLAGS= $MYMAKE check; testordie "error testing tar"@g' -e 's@MAKEFLAGS= $MYMAKE check; testordie "error testing specsum"@#MAKEFLAGS= $MYMAKE check; testordie "error testing specsum"@g' ./tools/src/buildtools
			yes | ./tools/src/buildtools
			bash -c "source ./shrc; runspec -V"
			bash -c "source ./shrc; ./bin/packagetools fugaku"
			bash -c "source ./shrc; specmd5sum ./fugaku-1.1.tar > ./fugaku-1.1.tar.md5"
			sed -i -e 's@runspec --test=dots@echo runspec --test=dots@g' ./install.sh
			SPEC_INSTALL_NOCHECK=1 ./install.sh -f -d $ROOTDIR/$BM/ -u fugaku
		else
			./install.sh -f -d $ROOTDIR/$BM/ -u linux-suse10-amd64
		fi
		cd -
	fi
        cd $ROOTDIR/$BM/
        #patch -p1 --forward < $ROOTDIR/patches/*1-${BM}*.patch
        dump_omp_config $ROOTDIR/$BM/config/nedo.cfg
	if [ -z $1 ]; then
		bash -c "source ./shrc; runspec --config=nedo.cfg --action=scrub --define COMP=intel --define RESDIR=0 gross"
		bash -c "source ./shrc; runspec --config=nedo.cfg --action=build --size=train --define COMP=intel --define RESDIR=0 gross"
		bash -c "source ./shrc; runspec --config=nedo.cfg --action=build --size=train --define COMP=sde --define RESDIR=0 gross"
	elif [[ "$1" = *"gnu"* ]]; then
		bash -c "source ./shrc; runspec --config=nedo.cfg --action=scrub --define COMP=gnu --define RESDIR=0 gross"
		bash -c "source ./shrc; runspec --config=nedo.cfg --action=build --size=train --define COMP=gnu --define RESDIR=0 gross"
	elif [[ "$1" = *"fuji"* ]]; then
		bash -c "source ./shrc; runspec --config=nedo.cfg --action=scrub --define COMP=fuji --define RESDIR=0 gross"
		if [[ "`hostname -s`" = *"peach"* ]]; then    #XXX: peach fccpx/static doesnt support mcmodel
			bash -c "source ./shrc; runspec --config=nedo.cfg --action=build --size=train --define COMP=fuji --define HOST=`hostname -s` --define RESDIR=0 gross ^bt331 ^swim"
		else
			bash -c "source ./shrc; runspec --config=nedo.cfg --action=build --size=train --define COMP=fuji --define HOST=`hostname -s` --define RESDIR=0 gross"
		fi
	fi
	# check that most/all are static
	find $ROOTDIR/$BM/benchspec/ -path '*/build_peak_*.0000/*' -executable -type f -exec echo {} \; -exec ldd {} \;
	if [[ "$1" = *"fuji"* ]]; then echo -e "\nWRN: if running gem5, then copy SPEC_OMP/benchspec/OMP2012/*/build/ and SPEC_OMP/benchspec/OMP2012/*/exe/ to server which runs gem"; fi
        cd $ROOTDIR
fi
