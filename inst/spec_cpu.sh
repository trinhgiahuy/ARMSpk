#!/bin/bash

function dump_cpu_config {
cat <<EOF > $1
default:
   action               = validate
   bench_post_setup     = sync
   command_add_redirect = 1
   flagsurl1            = http://www.spec.org/cpu2017/flags/Intel-ic19.0u1-official-linux64.2019-07-09.xml
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
%else
   submit               = ulimit -n 4096; ulimit -s unlimited; \$command
%endif

default:
%if '%{COMP}' eq 'gnu'
   CC                   = gcc -m64 -std=c11
   CXX                  = g++ -m64
   FC                   = gfortran -m64
   OPT_ROOT             = -O3 -march=native -funroll-loops -ffast-math -ftree-vectorize
   LDCFLAGS             = -static
   LDCXXFLAGS           = -static
   LDFFLAGS             = -static
%elif '%{COMP}' eq 'intel' || '%{COMP}' eq 'sde'
   CC                   = icc -m64 -std=c11
   CXX                  = icpc -m64
   FC                   = ifort -m64
   OPT_ROOT             = -O3 -xHOST -no-prec-div -qopt-mem-layout-trans=4 -qopt-prefetch
   EXTRA_FOPTIMIZE      = -nostandard-realloc-lhs
   LDCFLAGS             = -static -static-intel -qopenmp-link=static
   LDCXXFLAGS           = -static -static-intel -qopenmp-link=static
   LDFFLAGS             = -static -static-intel -qopenmp-link=static
%elif '%{COMP}' eq 'fuji'
   CC                   = fcc -m64 -std=c11
   CXX                  = FCC -m64
   FC                   = frt -m64
   OPT_ROOT             = -Kfast -O3 -march=native -funroll-loops -ffast-math -ftree-vectorize
   EXTRA_FOPTIMIZE      = -nostandard-realloc-lhs
   LDCFLAGS             = -Bstatic
   LDCXXFLAGS           = -Bstatic
   LDFFLAGS             = -Bstatic
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
%else
intspeed:
   EXTRA_COPTIMIZE      = -fopenmp -DSPEC_OPENMP

fpspeed:
   EXTRA_OPTIMIZE       = -fopenmp -DSPEC_OPENMP
%endif

intspeed,fpspeed:
   PORTABILITY          = -DSPEC_LP64

500.perlbench_r,600.perlbench_s:
   CPORTABILITY         = -DSPEC_LINUX_X64

502.gcc_r,602.gcc_s:
%if '%{COMP}' eq 'gnu' || '%{COMP}' eq 'fuji'
   CPORTABILITY         = -fgnu89-inline
%endif

521.wrf_r,621.wrf_s:
   CPORTABILITY         = -DSPEC_CASE_FLAG
%if '%{COMP}' eq 'gnu' || '%{COMP}' eq 'fuji'
   FPORTABILITY         = -fconvert=big-endian
%else
   FPORTABILITY         = -convert big_endian
%endif

523.xalancbmk_r,623.xalancbmk_s:
   CXXPORTABILITY       = -DSPEC_LINUX

526.blender_r:
   CPORTABILITY         = -DSPEC_LINUX -funsigned-char

527.cam4_r,627.cam4_s:
   CPORTABILITY         = -DSPEC_CASE_FLAG

628.pop2_s:
   CPORTABILITY         = -DSPEC_CASE_FLAG
%if '%{COMP}' eq 'gnu' || '%{COMP}' eq 'fuji'
   FPORTABILITY         = -fconvert=big-endian
%else
   FPORTABILITY         = -convert big_endian -assume byterecl
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
		alias fcc=fccpx
		alias FCC=FCCpx
		alias frt=frtpx
	else
		echo "check fugaku later"; exit 1
	fi
else
	echo 'wrong compiler'
	exit 1
fi

BM="SPEC_CPU"
if [ ! -f $ROOTDIR/$BM/bin/runcpu ]; then
        if [ ! -f $ROOTDIR/dep/cpu2017-1.1.0.iso ]; then
                echo -e "ERR: cannot find cpu2017-1.1.0.iso under dep/; please fix it!"
                exit 1
        fi
        cd $ROOTDIR/dep/
        mkdir -p /tmp/mnt_$BM
        fuseiso ./cpu2017-1.1.0.iso /tmp/mnt_$BM
        cd /tmp/mnt_$BM/
        mkdir -p $ROOTDIR/$BM/
	if [[ "$1" = *"fuji"* ]] && ! [[ "`hostname -s`" = *"peach"* ]]; then
		./install.sh -f -d $ROOTDIR/$BM/ -u linux-aarch64
	else
	        ./install.sh -f -d $ROOTDIR/$BM/ -u linux-x86_64
	fi
        cd -
        fusermount -u /tmp/mnt_$BM
        cd $ROOTDIR/$BM/
        #patch -p1 < $ROOTDIR/patches/*1-${BM}*.patch
        dump_cpu_config $ROOTDIR/$BM/config/nedo.cfg
	if [ -z $1 ]; then
		bash -c "source ./shrc; runcpu --config=nedo.cfg --action=scrub --define COMP=intel --define RESDIR=0 intspeed fpspeed intrate fprate"
		#XXX: scrub ignores compiler: bash -c "source ./shrc; runcpu --config=nedo.cfg --action=scrub --define COMP=sde --define RESDIR=0 intspeed fpspeed intrate fprate"
		bash -c "source ./shrc; runcpu --config=nedo.cfg --action=build --define COMP=intel --define RESDIR=0 intspeed fpspeed"
		bash -c "source ./shrc; runcpu --config=nedo.cfg --action=build --define COMP=sde --define RESDIR=0 intspeed fpspeed"
	elif [[ "$1" = *"gnu"* ]]; then
		bash -c "source ./shrc; runcpu --config=nedo.cfg --action=scrub --define COMP=gnu --define RESDIR=0 intspeed fpspeed intrate fprate"
		bash -c "source ./shrc; runcpu --config=nedo.cfg --action=build --define COMP=gnu --define RESDIR=0 intspeed fpspeed"
	elif [[ "$1" = *"fuji"* ]]; then
		bash -c "source ./shrc; runcpu --config=nedo.cfg --action=scrub --define COMP=fuji --define RESDIR=0 intspeed fpspeed intrate fprate"
		bash -c "source ./shrc; runcpu --config=nedo.cfg --action=build --define COMP=fuji --define RESDIR=0 intspeed fpspeed"
	fi
	# check that all are static
	find $ROOTDIR/$BM/benchspec/ -path '*/build_peak_*.0000/*' -executable -type f -exec echo {} \; -exec ldd {} \;
        cd $ROOTDIR
fi