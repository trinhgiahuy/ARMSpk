#!/bin/bash

function dump_cpu_config {
cat <<EOF > $1
default:
   action               = validate
   bench_post_setup     = sync
   command_add_redirect = 1
#   flagsurl1            = http://www.spec.org/cpu2017/flags/Intel-ic19.0u1-official-linux64.2019-07-09.xml
   iterations           = 1
   line_width           = 1020
   log_line_width       = 1020
   makeflags            = -j 32
   mean_anyway          = 1
   output_format        = txt
   preenv               = 1
   reportable           = 0
   tune                 = peak
   label                = %{COMP}
%if '%{COMP}' eq 'fujitrad' || '%{COMP}' eq 'fujiclang' || '%{COMP}' eq 'gem5' || '%{COMP}' eq 'llvm12'
   verify_binaries      = no
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

#speed: 1 copy of each benchmark in a suite is run. The tester may choose how many OpenMP threads to use.
#rate:  The tester chooses how many concurrent copies to run. OpenMP is disabled.
intspeed,fpspeed:
%if '%{COMP}' eq 'sde'
   submit               = ulimit -n 4096; ulimit -s unlimited; sde64 -sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 -dcfg:out_base_name %{RESDIR}/dcfg-out.wl-\${workload} -align_checker_prefetch 0 -align_correct 0 -emu_fast 1 -bdw -- \$command
%elif '%{COMP}' eq 'fujitrad' || '%{COMP}' eq 'fujiclang' || '%{COMP}' eq 'llvm12'
   submit               = export XOS_MMM_L_PAGING_POLICY=demand:demand:demand; export XOS_MMM_L_ARENA_LOCK_TYPE=0; export OMP_PROC_BIND=close; ulimit -s unlimited; \$command
%elif '%{COMP}' eq 'gem5'
   submit               = ulimit -n 4096; ulimit -s unlimited; bash %{RESDIR}/gem5wrap.sh \$command
%else
   submit               = ulimit -n 4096; ulimit -s unlimited; \$command
%endif

default:
%if '%{COMP}' eq 'gnu'
   CC                   = gcc -m64 -std=c11
   CXX                  = g++ -m64
   FC                   = gfortran -m64
   OPT_ROOT             = -O3 -march=native -flto -funroll-loops -ffast-math -ftree-vectorize
   LDCFLAGS             = -flto ${MAYBESTATIC}
   LDCXXFLAGS           = -flto ${MAYBESTATIC}
   LDFFLAGS             = -flto ${MAYBESTATIC}
%elif '%{COMP}' eq 'intel' || '%{COMP}' eq 'sde'
   CC                   = icc -m64 -std=c11
   CXX                  = icpc -m64
   FC                   = ifort -m64
   OPT_ROOT             = -O3 -xHOST -no-prec-div -qopt-mem-layout-trans=4 -qopt-prefetch
   EXTRA_FOPTIMIZE      = -nostandard-realloc-lhs
   LDCFLAGS             = -static -static-intel -qopenmp-link=static
   LDCXXFLAGS           = -static -static-intel -qopenmp-link=static
   LDFFLAGS             = -static -static-intel -qopenmp-link=static
%elif '%{COMP}' eq 'fujitrad'
   CC                   = fcc -m64 -std=c11
   CXX                  = FCC -m64
   FC                   = frt -m64
   OPT_ROOT             = -Kfast,ocl,eval_concurrent,largepage
   EXTRA_FOPTIMIZE      = -Klto
   LDCFLAGS             =
   LDCXXFLAGS           =
   LDFFLAGS             =
%elif '%{COMP}' eq 'fujiclang'
   CC                   = fcc -m64 -std=c11
   CXX                  = FCC -m64
   FC                   = frt -m64
   OPT_ROOT             = -Nclang -mcpu=a64fx+sve
   EXTRA_FOPTIMIZE      = -Kfast,ocl,eval_concurrent,largepage,lto
   EXTRA_COPTIMIZE      = -Ofast -ffj-ocl -ffj-eval-concurrent -ffj-largepage -flto
   EXTRA_CXXOPTIMIZE    = -Ofast -ffj-ocl -ffj-eval-concurrent -ffj-largepage -flto
   LDCFLAGS             =
   LDCXXFLAGS           =
   LDFFLAGS             =
%elif '%{COMP}' eq 'gem5'
   CC                   = fcc -m64 -std=c11
   CXX                  = FCC -m64
   FC                   = frt -m64
   OPT_ROOT             = -Nclang -mcpu=a64fx+sve
   EXTRA_FOPTIMIZE      = -Kfast,ocl,eval_concurrent,nolargepage,nolto
   EXTRA_COPTIMIZE      = -Ofast -ffj-ocl -ffj-eval-concurrent -ffj-no-largepage -fno-lto
   EXTRA_CXXOPTIMIZE    = -Ofast -ffj-ocl -ffj-eval-concurrent -ffj-no-largepage -fno-lto
   LDCFLAGS             =
   LDCXXFLAGS           =
   LDFFLAGS             =
%elif '%{COMP}' eq 'llvm12'
   CC                   = clang -m64 -std=c11
   CXX                  = clang++ -m64
   FC                   = frt -m64
   OPT_ROOT             = -ffast-math -mcpu=a64fx -mtune=a64fx
   EXTRA_FOPTIMIZE      = -Kfast,ocl,eval_concurrent,largepage,lto
   EXTRA_COPTIMIZE      = -Ofast -mllvm -polly -mllvm -polly-vectorizer=polly -flto=thin
   EXTRA_CXXOPTIMIZE    = -Ofast -mllvm -polly -mllvm -polly-vectorizer=polly -flto=thin
   LDCFLAGS             = -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)
   LDCXXFLAGS           = -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)
   LDFFLAGS             = -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)
%else
%error wrong or unsupported COMP variable specified
%endif
   CC_VERSION_OPTION    = --version
   CXX_VERSION_OPTION   = --version
   FC_VERSION_OPTION    = --version
   COPTIMIZE            = \$(OPT_ROOT)
   CXXOPTIMIZE          = \$(OPT_ROOT)
   FOPTIMIZE            = \$(OPT_ROOT)

%if '%{COMP}' eq 'intel' || '%{COMP}' eq 'sde'
intspeed:
   EXTRA_COPTIMIZE      = -qopenmp -DSPEC_OPENMP

fpspeed:
   EXTRA_OPTIMIZE       = -qopenmp -DSPEC_OPENMP
%elif '%{COMP}' eq 'gnu' || '%{COMP}' eq 'llvm12'
intspeed:
   EXTRA_COPTIMIZE      = -fopenmp -DSPEC_OPENMP -fno-strict-aliasing

fpspeed:
   EXTRA_OPTIMIZE       = -fopenmp -DSPEC_OPENMP
%else
intspeed:
   EXTRA_COPTIMIZE      = -fopenmp -DSPEC_OPENMP -fno-strict-aliasing -Knofp_relaxed

fpspeed:
   EXTRA_OPTIMIZE       = -fopenmp -DSPEC_OPENMP
%endif

intspeed,fpspeed:
   PORTABILITY          = -DSPEC_LP64

500.perlbench_r,600.perlbench_s:
   CPORTABILITY         = -DSPEC_LINUX_X64

502.gcc_r,602.gcc_s:
%if '%{COMP}' eq 'gnu' || '%{COMP}' eq 'fujitrad' || '%{COMP}' eq 'fujiclang' || '%{COMP}' eq 'gem5' || '%{COMP}' eq 'llvm12'
   CPORTABILITY         = -fgnu89-inline
%endif

521.wrf_r,621.wrf_s:
   CPORTABILITY         = -DSPEC_CASE_FLAG
%if '%{COMP}' eq 'gnu'
   FPORTABILITY         = -fconvert=big-endian -fallow-argument-mismatch
%elif '%{COMP}' eq 'llvm12'
   FPORTABILITY         = -fconvert=big-endian
%elif '%{COMP}' eq 'fujitrad' || '%{COMP}' eq 'fujiclang' || '%{COMP}' eq 'gem5'
   FPORTABILITY         = -fconvert=big-endian
   LDOUT_EXTRA_OPTIONS  = -lfj90f_sve
%else
   FPORTABILITY         = -convert big_endian
%endif

523.xalancbmk_r,623.xalancbmk_s:
   CXXPORTABILITY       = -DSPEC_LINUX

526.blender_r:
   CPORTABILITY         = -DSPEC_LINUX -funsigned-char

527.cam4_r,627.cam4_s:
   CPORTABILITY         = -DSPEC_CASE_FLAG
%if '%{COMP}' eq 'fujitrad' || '%{COMP}' eq 'fujiclang' || '%{COMP}' eq 'gem5'
   LDOUT_EXTRA_OPTIONS  = -lfj90f_sve
%elif '%{COMP}' eq 'gnu'
   FPORTABILITY         = -fallow-argument-mismatch
%endif

607.cactuBSSN_s:
%if '%{COMP}' eq 'llvm12'
   FPORTABILITY         = -Knohpctag
%endif

621.wrf_s,627.cam4_s,628.pop2_s:
%if '%{COMP}' eq 'fujiclang'
   EXTRA_COPTIMIZE      = -Ofast -ffj-ocl -ffj-eval-concurrent -ffj-largepage
   EXTRA_CXXOPTIMIZE    = -Ofast -ffj-ocl -ffj-eval-concurrent -ffj-largepage
%elif '%{COMP}' eq 'llvm12'
   EXTRA_COPTIMIZE      = -Ofast -mllvm -polly -mllvm -polly-vectorizer=polly
   EXTRA_CXXOPTIMIZE    = -Ofast -mllvm -polly -mllvm -polly-vectorizer=polly
%endif

625.x264_s:
%if '%{COMP}' eq 'gnu' || '%{COMP}' eq 'llvm12'
   #yes it's not see https://www.spec.org/cpu2017/Docs/benchmarks/625.x264_s.html but i need sleep
   PORTABILITY          = -fcommon
%endif

628.pop2_s:
   CPORTABILITY         = -DSPEC_CASE_FLAG
%if '%{COMP}' eq 'gnu'
   FPORTABILITY         = -fconvert=big-endian -fallow-argument-mismatch
%elif '%{COMP}' eq 'llvm12'
   FPORTABILITY         = -fconvert=big-endian
%elif '%{COMP}' eq 'fujitrad' || '%{COMP}' eq 'fujiclang' || '%{COMP}' eq 'gem5'
   FPORTABILITY         = -fconvert=big-endian
   LDOUT_EXTRA_OPTIONS  = -lfj90f_sve
%else
   FPORTABILITY         = -convert big_endian -assume byterecl
%endif

648.exchange2_s:
%if '%{COMP}' eq 'fujitrad' || '%{COMP}' eq 'fujiclang' || '%{COMP}' eq 'gem5'
   LDOUT_EXTRA_OPTIONS  = -lfj90i -lfj90fmt_sve -lfj90f -lfj90i -lfjsrcinfo
%endif
EOF
}

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/inst/_common.sh
load_compiler_env "$1"

BM="SPEC_CPU"
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if [ ! -f $ROOTDIR/$BM/bin/runcpu ]; then
	mkdir -p $ROOTDIR/$BM/
	if [ ! -f $ROOTDIR/dep/mnt_$BM/shrc ]; then
		if [ ! -f $ROOTDIR/dep/cpu2017-1.1.0.iso ]; then
			echo -e "ERR: cannot find cpu2017-1.1.0.iso under dep/; please fix it!"
			exit 1
		fi
		cd $ROOTDIR/dep/
		mkdir -p /tmp/mnt_$BM
		fuseiso ./cpu2017-1.1.0.iso /tmp/mnt_$BM
		cd /tmp/mnt_$BM/
		if lscpu | grep 'aarch64' >/dev/null 2>&1 || [[ "$1" = *"gem5"* ]]; then
			./install.sh -f -d $ROOTDIR/$BM/ -u linux-aarch64
		else
			./install.sh -f -d $ROOTDIR/$BM/ -u linux-x86_64
		fi
		cd -
		fusermount -u /tmp/mnt_$BM
	elif [ ! -f $ROOTDIR/$BM/shrc ]; then
		cd $ROOTDIR/dep/mnt_$BM
		if lscpu | grep 'aarch64' >/dev/null 2>&1 || [[ "$1" = *"gem5"* ]]; then
			./install.sh -f -d $ROOTDIR/$BM/ -u linux-aarch64
		else
			./install.sh -f -d $ROOTDIR/$BM/ -u linux-x86_64
		fi
		cd -
	fi
        cd $ROOTDIR/$BM/
        #patch -p1 --forward < $ROOTDIR/patches/*1-${BM}*.patch
        dump_cpu_config $ROOTDIR/$BM/config/nedo.cfg
	if [[ "$1" = *"intel"* ]]; then
		bash -c "source ./shrc; runcpu --config=nedo.cfg --action=scrub --define COMP=intel --define RESDIR=0 intspeed fpspeed intrate fprate"
		#XXX: scrub ignores compiler: bash -c "source ./shrc; runcpu --config=nedo.cfg --action=scrub --define COMP=sde --define RESDIR=0 intspeed fpspeed intrate fprate"
		bash -c "source ./shrc; runcpu --config=nedo.cfg --action=build --define COMP=intel --define RESDIR=0 intspeed fpspeed --ignore_error"
		bash -c "source ./shrc; runcpu --config=nedo.cfg --action=build --define COMP=sde --define RESDIR=0 intspeed fpspeed --ignore_error"
	elif [[ "$1" = *"gnu"* ]]; then
		if [ -n "$FJBLAS" ]; then sed -i -e 's/ -m64//' $ROOTDIR/$BM/config/nedo.cfg; fi
		bash -c "source ./shrc; runcpu --config=nedo.cfg --action=scrub --define COMP=gnu --define RESDIR=0 intspeed fpspeed intrate fprate"
		bash -c "source ./shrc; runcpu --config=nedo.cfg --action=build --define COMP=gnu --define RESDIR=0 intspeed fpspeed --ignore_error"
	elif [[ "$1" = *"fujitrad"* ]]; then
		bash -c "source ./shrc; runcpu --config=nedo.cfg --action=scrub --define COMP=fujitrad --define RESDIR=0 intspeed fpspeed intrate fprate"
		bash -c "source ./shrc; runcpu --config=nedo.cfg --action=build --define COMP=fujitrad --define RESDIR=0 intspeed fpspeed --ignore_error"
	elif [[ "$1" = *"fujiclang"* ]]; then
		bash -c "source ./shrc; runcpu --config=nedo.cfg --action=scrub --define COMP=fujiclang --define RESDIR=0 intspeed fpspeed intrate fprate"
		bash -c "source ./shrc; runcpu --config=nedo.cfg --action=build --define COMP=fujiclang --define RESDIR=0 intspeed fpspeed --ignore_error"
	elif [[ "$1" = *"gem5"* ]]; then
		bash -c "source ./shrc; runcpu --config=nedo.cfg --action=scrub --define COMP=gem5 --define RESDIR=0 intspeed fpspeed intrate fprate"
		bash -c "source ./shrc; runcpu --config=nedo.cfg --action=build --define COMP=gem5 --define RESDIR=0 intspeed fpspeed --ignore_error"
	elif [[ "$1" = *"llvm12"* ]]; then
		bash -c "source ./shrc; runcpu --config=nedo.cfg --action=scrub --define COMP=llvm12 --define RESDIR=0 intspeed fpspeed intrate fprate"
		bash -c "source ./shrc; runcpu --config=nedo.cfg --action=build --define COMP=llvm12 --define RESDIR=0 intspeed fpspeed --ignore_error"
		#bash -c "source ./shrc; runcpu --config=nedo.cfg --action=build --define COMP=llvm12 --define RESDIR=0 intspeed fpspeed ^603.bwaves_s ^621.wrf_s ^627.cam4_s ^628.pop2_s ^638.imagick_s ^654.roms_s" #XXX fix later
	fi
	# check that all are static
	find $ROOTDIR/$BM/benchspec/ -path '*/build_peak_*.0000/*' -executable -type f -exec echo {} \; -exec ldd {} \;
	if [[ "$1" = *"gem5"* ]]; then echo -e "\nWRN: if running gem5, then copy SPEC_CPU/benchspec/CPU/*/build/ and SPEC_CPU/benchspec/CPU/*/exe/ to server which runs gem"; fi
        cd $ROOTDIR
fi
