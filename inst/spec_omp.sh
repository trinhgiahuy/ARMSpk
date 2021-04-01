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
submit         = ulimit -n 4096; ulimit -s unlimited; bash %{RESDIR}/gem5wrap.sh \$command
BOPTS          = -Kfast,eval_concurrent,openmp -O3 -march=armv8.3-a+sve -funroll-loops
BLINK          = -Bstatic -lm
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
FC             = frtpx
CC             = fccpx
CXX            = FCCpx
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
elif [[ "$1" = *"fuji"* ]]; then
	if [[ "`hostname -s`" = *"peach"* ]]; then
		module load FujitsuCompiler/202007
		#alias fcc=fccpx
		#alias FCC=FCCpx
		#alias frt=frtpx	# XXX: not working
	else
		echo "check fugaku later"; exit 1
	fi
else
	echo 'wrong compiler'
	exit 1
fi

BM="SPEC_OMP"
dump_omp_config $ROOTDIR/$BM/config/nedo.cfg; exit
if [ ! -f $ROOTDIR/$BM/bin/runcpu ]; then
	if [ ! -f $ROOTDIR/dep/omp2012-1.1.iso ]; then
		echo -e "ERR: cannot find omp2012-1.1.iso under dep/; please fix it!"
		exit 1
	fi
	cd $ROOTDIR/dep/
	mkdir -p /tmp/mnt_$BM
	fuseiso ./omp2012-1.1.iso /tmp/mnt_$BM
	cd /tmp/mnt_$BM/
	mkdir -p $ROOTDIR/$BM/
	if [[ "$1" = *"fuji"* ]] && ! [[ "`hostname -s`" = *"peach"* ]]; then
		./install.sh -f -d $ROOTDIR/$BM/ -u linux-apm-arm64
	else
	        ./install.sh -f -d $ROOTDIR/$BM/ -u linux-suse10-amd64
	fi
        cd -
        fusermount -u /tmp/mnt_$BM
        cd $ROOTDIR/$BM/
        #patch -p1 < $ROOTDIR/patches/*1-${BM}*.patch
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
		if [[ "`hostname -s`" = *"peach"* ]]; then
			bash -c "source ./shrc; runspec --config=nedo.cfg --action=build --size=train --define COMP=fuji --define RESDIR=0 gross ^bt331 ^swim"	#XXX: peach fccpx doesnt support mcmodel
		else
			bash -c "source ./shrc; runspec --config=nedo.cfg --action=build --size=train --define COMP=fuji --define RESDIR=0 gross"
		fi
	fi
	# check that most/all are static
	find $ROOTDIR/$BM/benchspec/ -path '*/build_peak_*.0000/*' -executable -type f -exec echo {} \; -exec ldd {} \;
        cd $ROOTDIR
fi
