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
%if '%{COMP}' eq 'fujitrad' || '%{COMP}' eq 'fujiclang' || '%{COMP}' eq 'gem5' || '%{COMP}' eq 'llvm12'
check_md5      = 0
%endif

%if '%{COMP}' eq 'gnu'
submit         = ulimit -n 4096; ulimit -s unlimited; \$command
BOPTS          = -O3 -fopenmp -march=native -flto -funroll-loops -ffast-math -ftree-vectorize
BLINK          = -flto ${MAYBESTATIC}
%elif '%{COMP}' eq 'intel'
submit         = ulimit -n 4096; ulimit -s unlimited; \$command
BOPTS          = -O3 -qopenmp -xHOST -no-prec_div -fp-model fast=2 -fma
BLINK          = -static -static-intel -qopenmp-link=static
%elif '%{COMP}' eq 'sde'
submit         = ulimit -n 4096; ulimit -s unlimited; sde64 -sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 -dcfg:out_base_name %{RESDIR}/dcfg-out.wl-\${workload} -align_checker_prefetch 0 -align_correct 0 -emu_fast 1 -bdw -- \$command
BOPTS          = -O3 -qopenmp -xHOST -no-prec_div -fp-model fast=2 -fma
BLINK          = -static -static-intel -qopenmp-link=static
%elif '%{COMP}' eq 'fujitrad'
submit         = export XOS_MMM_L_PAGING_POLICY=demand:demand:demand; export XOS_MMM_L_ARENA_LOCK_TYPE=0; export OMP_PROC_BIND=close; ulimit -s unlimited; \$command
BOPTS          = -Kfast,openmp,ocl,eval_concurrent,largepage
BLINK          =
%elif '%{COMP}' eq 'fujiclang'
submit         = export XOS_MMM_L_PAGING_POLICY=demand:demand:demand; export XOS_MMM_L_ARENA_LOCK_TYPE=0; export OMP_PROC_BIND=close; ulimit -s unlimited; \$command
BOPTS          = -Nclang -mcpu=a64fx+sve -fopenmp
BLINK          =
%elif '%{COMP}' eq 'gem5'
submit         = ulimit -n 4096; ulimit -s unlimited; bash %{RESDIR}/gem5wrap.sh \$command
BOPTS          = -Nclang -mcpu=a64fx+sve -fopenmp
BLINK          =
%elif '%{COMP}' eq 'llvm12'
submit         = export XOS_MMM_L_PAGING_POLICY=demand:demand:demand; export XOS_MMM_L_ARENA_LOCK_TYPE=0; export OMP_PROC_BIND=close; ulimit -s unlimited; \$command
BOPTS          = -mcpu=a64fx+sve -mtune=a64fx+sve -fopenmp
BLINK          = -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)
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
FOPTIMIZE      = \${BOPTS} -fallow-argument-mismatch -fallow-invalid-boz
COPTIMIZE      = \${BOPTS}
CXXOPTIMIZE    = \${BOPTS}
OS_LIBS        = \${BLINK}
%elif '%{COMP}' eq 'intel' || '%{COMP}' eq 'sde'
sw_compiler    = Intel Compiler Suite
FC             = ifort
CC             = icc
CXX            = icpc
FOPTIMIZE      = \${BOPTS}
COPTIMIZE      = \${BOPTS}
CXXOPTIMIZE    = \${BOPTS}
OS_LIBS        = \${BLINK}
%elif '%{COMP}' eq 'fujitrad'
sw_compiler    = FUJITSU Software Technical Computing Suite
FC             = frt
CC             = fcc
CXX            = FCC
FOPTIMIZE      = \${BOPTS} -Klto
COPTIMIZE      = \${BOPTS}
CXXOPTIMIZE    = \${BOPTS}
OS_LIBS        = \${BLINK}
%elif '%{COMP}' eq 'fujiclang'
sw_compiler    = FUJITSU Software Technical Computing Suite
FC             = frt
CC             = fcc
CXX            = FCC
FOPTIMIZE      = \${BOPTS} -Kfast,ocl,eval_concurrent,largepage,lto
COPTIMIZE      = \${BOPTS} -Ofast -ffj-ocl -ffj-eval-concurrent -ffj-largepage -flto
CXXOPTIMIZE    = \${BOPTS} -Ofast -ffj-ocl -ffj-eval-concurrent -ffj-largepage -flto
OS_LIBS        = \${BLINK}
%elif '%{COMP}' eq 'gem5'
sw_compiler    = FUJITSU Software Technical Computing Suite
FC             = frt
CC             = fcc
CXX            = FCC
FOPTIMIZE      = \${BOPTS} -Kfast,ocl,eval_concurrent,nolargepage,nolto
COPTIMIZE      = \${BOPTS} -Ofast -ffj-ocl -ffj-eval-concurrent -ffj-no-largepage -fno-lto
CXXOPTIMIZE    = \${BOPTS} -Ofast -ffj-ocl -ffj-eval-concurrent -ffj-no-largepage -fno-lto
OS_LIBS        = \${BLINK}
%elif '%{COMP}' eq 'llvm12'
sw_compiler    = LLVM Compiler Infrastructure release 12.0.0
FC             = frt
CC             = clang
CXX            = clang++
FOPTIMIZE      = \${BOPTS} -Kfast,ocl,eval_concurrent,largepage,lto
COPTIMIZE      = \${BOPTS} -Ofast -ffast-math -mllvm -polly -mllvm -polly-vectorizer=polly -flto=thin
CXXOPTIMIZE    = \${BOPTS} -Ofast -ffast-math -mllvm -polly -mllvm -polly-vectorizer=polly -flto=thin
OS_LIBS        = \${BLINK}
%endif

350.md=default=default=default:
%if '%{COMP}' eq 'gnu'
FPORTABILITY   = -ffree-form -fno-range-check
%elif '%{COMP}' eq 'fujitrad' || '%{COMP}' eq 'fujiclang' || '%{COMP}' eq 'gem5' || '%{COMP}' eq 'llvm12'
FPORTABILITY   = -Free
%else
FPORTABILITY   = -free
%endif

357.bt331=default=default=default:
%if '%{COMP}' eq 'intel' || '%{COMP}' eq 'sde'
FPORTABILITY   = -mcmodel=medium
LDOPT          = -shared-intel
PASS1_LDOPT    = -shared-intel
%elif '%{COMP}' eq 'fujitrad' || '%{COMP}' eq 'fujiclang' || '%{COMP}' eq 'gem5' || '%{COMP}' eq 'llvm12'
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
%elif '%{COMP}' eq 'fujitrad' || '%{COMP}' eq 'fujiclang' || '%{COMP}' eq 'gem5' || '%{COMP}' eq 'llvm12'
FPORTABILITY   = -mcmodel=large
%endif

367.imagick=default=default=default:
CPORTABILITY   = -std=c99

372.smithwa=default=default=default:
%if '%{COMP}' eq 'fujitrad' || '%{COMP}' eq 'fujiclang' || '%{COMP}' eq 'gem5'
CPORTABILITY   = -fsigned-char
%endif
EOF
}

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/inst/_common.sh
load_compiler_env "$1"

BM="SPEC_OMP"
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
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
		if lscpu | grep 'aarch64' >/dev/null 2>&1 || [[ "$1" = *"gem5"* ]]; then
			./install.sh -f -d $ROOTDIR/$BM/ -u linux-apm-arm64
		else
			./install.sh -f -d $ROOTDIR/$BM/ -u linux-suse10-amd64
		fi
		cd -
		fusermount -u /tmp/mnt_$BM
	elif [ ! -f $ROOTDIR/$BM/shrc ]; then
		cd $ROOTDIR/dep/mnt_$BM
		if lscpu | grep 'aarch64' >/dev/null 2>&1 || [[ "$1" = *"gem5"* ]]; then
			chmod -R +w $ROOTDIR/dep/mnt_$BM
			if [ ! -f ./fugaku-1.1.tar ]; then
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
			fi
			SPEC_INSTALL_NOCHECK=1 ./install.sh -f -d $ROOTDIR/$BM/ -u fugaku
		else
			./install.sh -f -d $ROOTDIR/$BM/ -u linux-suse10-amd64
		fi
		cd -
	fi
        cd $ROOTDIR/$BM/
        #patch -p1 --forward < $ROOTDIR/patches/*1-${BM}*.patch
        dump_omp_config $ROOTDIR/$BM/config/nedo.cfg
	if [[ "$1" = *"intel"* ]]; then
		bash -c "source ./shrc; runspec --config=nedo.cfg --action=scrub --define COMP=intel --define RESDIR=0 gross"
		bash -c "source ./shrc; runspec --config=nedo.cfg --action=build --size=train --define COMP=intel --define RESDIR=0 gross --ignore_error"
		bash -c "source ./shrc; runspec --config=nedo.cfg --action=build --size=train --define COMP=sde --define RESDIR=0 gross --ignore_error"
	elif [[ "$1" = *"gnu"* ]]; then
		bash -c "source ./shrc; runspec --config=nedo.cfg --action=scrub --define COMP=gnu --define RESDIR=0 gross"
		bash -c "source ./shrc; runspec --config=nedo.cfg --action=build --size=train --define COMP=gnu --define RESDIR=0 gross --ignore_error"
	elif [[ "$1" = *"fujitrad"* ]]; then
		bash -c "source ./shrc; runspec --config=nedo.cfg --action=scrub --define COMP=fujitrad --define RESDIR=0 gross"
		bash -c "source ./shrc; runspec --config=nedo.cfg --action=build --size=train --define COMP=fujitrad --define RESDIR=0 gross --ignore_error"
	elif [[ "$1" = *"fujiclang"* ]]; then
		bash -c "source ./shrc; runspec --config=nedo.cfg --action=scrub --define COMP=fujiclang --define RESDIR=0 gross"
		bash -c "source ./shrc; runspec --config=nedo.cfg --action=build --size=train --define COMP=fujiclang --define RESDIR=0 gross --ignore_error"
	elif [[ "$1" = *"gem5"* ]]; then
		bash -c "source ./shrc; runspec --config=nedo.cfg --action=scrub --define COMP=gem5 --define RESDIR=0 gross"
		bash -c "source ./shrc; runspec --config=nedo.cfg --action=build --size=train --define COMP=gem5 --define RESDIR=0 gross --ignore_error"
		#XXX: peach fccpx/static doesnt support mcmodel
		#bash -c "source ./shrc; runspec --config=nedo.cfg --action=build --size=train --define COMP=gem5 --define HOST=`hostname -s` --define RESDIR=0 gross ^bt331 ^swim"
	elif [[ "$1" = *"llvm12"* ]]; then
		bash -c "source ./shrc; runspec --config=nedo.cfg --action=scrub --define COMP=llvm12 --define RESDIR=0 gross"
		bash -c "source ./shrc; runspec --config=nedo.cfg --action=build --size=train --define COMP=llvm12 --define RESDIR=0 gross --ignore_error"
	fi
	# check that most/all are static
	find $ROOTDIR/$BM/benchspec/ -path '*/build_peak_*.0000/*' -executable -type f -exec echo {} \; -exec ldd {} \;
	if [[ "$1" = *"gem5"* ]]; then echo -e "\nWRN: if running gem5, then copy SPEC_OMP/benchspec/OMP2012/*/build/ and SPEC_OMP/benchspec/OMP2012/*/exe/ to server which runs gem"; fi
        cd $ROOTDIR
fi
