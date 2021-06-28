#!/bin/bash

function load_compiler_env {
	if [ -z "${XEONHOST}" ] && [ -z "${IKNLHOST}" ] && [ -z "${IKNMHOST}" ] && [ -z "${FUJIHOST}" ] && [ -z "${RFX7HOST}" ] ;
	then
		echo "ERR: new env and/or host, no known compiler for it, please fix me"
	fi

	if [ -n "${XEONHOST}" ] || [ -n "${IKNLHOST}" ] || [ -n "${IKNMHOST}" ]; then
		if [[ "$1" = *"intel"* ]]; then
			source $ROOTDIR/conf/intel.cfg
			source $INTEL_PACKAGE intel64 > /dev/null 2>&1
			export I_MPI_CC=icc; export I_MPI_CXX=icpc
			export I_MPI_F77=ifort; export I_MPI_F90=ifort
			alias ar=`which xiar`
			alias ld=`which xild`
			export ADVISOR_2018_DIR=${ADVISOR_2019_DIR}
			source $ROOTDIR/dep/spack/share/spack/setup-env.sh
			spack load openmpi@3.1.6%intel@19.0.1.144
			export OMPI_CC=$I_MPI_CC; export OMPI_CXX=$I_MPI_CXX
			export OMPI_F77=$I_MPI_F77; export OMPI_FC=$I_MPI_F90
		elif [[ "$1" = *"gnu"* ]]; then
			source $ROOTDIR/conf/intel.cfg
			source $(echo $INTEL_PACKAGE | cut -d'/' -f-3)/mkl/bin/mklvars.sh intel64 > /dev/null 2>&1
			source $ROOTDIR/dep/spack/share/spack/setup-env.sh
			spack load gcc@8.4.0
			spack load openmpi@3.1.6%gcc@8.4.0
			export OMPI_CC=gcc; export OMPI_CXX=g++
			export OMPI_F77=gfortran; export OMPI_FC=gfortran
			export MAYBESTATIC="-static"
		fi
	elif [ -n "${FUJIHOST}" ]; then
		if ! lscpu | grep 'sve' >/dev/null 2>&1; then
			echo "ERR: does not compile on login node; please use compute node"; exit 1
		elif [[ "$1" = *"gnu"* ]]; then
			if ! which spack >/dev/null 2>&1; then . /vol0004/apps/oss/spack/share/spack/setup-env.sh; fi
			spack load gcc@10.2.0 arch=linux-rhel8-a64fx; export LD_LIBRARY_PATH=/lib64:$LD_LIBRARY_PATH
			spack load fujitsu-mpi%gcc arch=linux-rhel8-a64fx; export LD_LIBRARY_PATH=/lib64:$LD_LIBRARY_PATH
			spack load hwloc@2.2.0%gcc arch=linux-rhel8-a64fx; export LD_LIBRARY_PATH=/lib64:$LD_LIBRARY_PATH
			export OMPI_CC=gcc; export OMPI_CXX=g++
			export OMPI_F77=gfortran; export OMPI_FC=gfortran
			export FJMPI="Y"; export FJBLAS="Y"; export MAYBESTATIC=""
		elif [[ "$1" = *"fujitrad"* ]] || [[ "$1" = *"fujiclang"* ]]; then
			sleep 0
		elif [[ "$1" = *"gem5"* ]]; then
			#module load FujitsuCompiler/202007
			export LD_LIBRARY_PATH=$ROOTDIR/dep/mpistub/lib:$LD_LIBRARY_PATH
		elif [[ "$1" = *"arm"* ]]; then
			module load /opt/arm/modulefiles/A64FX/RHEL/8/arm-linux-compiler-20.3/armpl/20.3.0
			#module load /opt/arm/modulefiles/A64FX/RHEL/8/gcc-9.3.0/armpl/20.3.0
		elif [[ "$1" = *"llvm12"* ]]; then
			if ! which spack >/dev/null 2>&1; then . /vol0004/apps/oss/spack/share/spack/setup-env.sh; fi
			#spack load llvm@12%gcc; export LD_LIBRARY_PATH=/lib64:$LD_LIBRARY_PATH
			spack load fujitsu-mpi%fj; export LD_LIBRARY_PATH=/lib64:$LD_LIBRARY_PATH
			spack load hwloc@1.11.11%fj; export LD_LIBRARY_PATH=/lib64:$LD_LIBRARY_PATH
			#
			LLVMDIR=$HOME/llvm-v12.0.0
			export PATH=$LLVMDIR/bin:$PATH
			export LD_LIBRARY_PATH=$LLVMDIR/lib:$LD_LIBRARY_PATH
			#
			export OMPI_CC=clang; export OMPI_CXX=clang++
			#export OMPI_F77=flang; export OMPI_FC=flang
			#export F18_FC=frt; export FORT90C="-Kfast,openmp,largepage,ocl,lto"; export FORT90CPX=${FORT90C}
		fi
	elif [ -n "${RFX7HOST}" ]; then
		if ! lscpu | grep 'sve' >/dev/null 2>&1; then
			echo "ERR: does not compile on login node; please use compute node"; exit 1
		elif [[ "$1" = *"gnu"* ]]; then
			module load /opt/arm/modulefiles/A64FX/RHEL/8/gcc-9.3.0/armpl/20.3.0
		elif [[ "$1" = *"arm"* ]]; then
			module load /opt/arm/modulefiles/A64FX/RHEL/8/arm-linux-compiler-20.3/armpl/20.3.0
		fi
	else
		echo 'ERR: wrong compiler, only support [intel | gnu | fujitrad | fujiclang | gem5 | llvm12]'
		exit 1
	fi
}

function instrument_kernel {

	if [[ "$1" = *"intel"* ]]; then

		sleep 0

	elif [[ "$1" = *"gnu"* ]]; then

		for FILE in `/bin/grep 'include.*ittnotify' -r "$2" | cut -d':' -f1 | sort -u`; do
			sed -i  -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE;
		done

	elif [[ "$1" = *"fujitrad"* ]]; then

		for FILE in `/bin/grep 'include.*ittnotify' -r "$2" | cut -d':' -f1 | sort -u`; do
			sed -i  -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE;
			#sed -i  -e 's/.*include.*ittnotify\.h.*/#include "fj_tool\/fapp.h"\n#define __itt_resume() fapp_start("kernel",1,0);\n#define __itt_pause() fapp_stop("kernel",1,0);\n#define __SSC_MARK(hex)/' $FILE;
		done

	elif [[ "$1" = *"fujiclang"* ]]; then

		for FILE in `/bin/grep 'include.*ittnotify' -r "$2" | cut -d':' -f1 | sort -u`; do
			sed -i  -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE;
			#sed -i  -e 's/.*include.*ittnotify\.h.*/#include "fj_tool\/fapp.h"\n#define __itt_resume() fapp_start("kernel",1,0);\n#define __itt_pause() fapp_stop("kernel",1,0);\n#define __SSC_MARK(hex)/' $FILE;
		done

	elif [[ "$1" = *"gem5"* ]]; then

		for FILE in `/bin/grep 'include.*ittnotify' -r "$2" | cut -d':' -f1 | sort -u`; do
			sed -i  -e 's/.*include.*ittnotify\.h.*/#ifndef _POSIX_C_SOURCE\n#define _POSIX_C_SOURCE 199309L\n#endif\n#include <time.h>\n#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE;
			if ! [[ "$2" = *"DLproxy"* ]] && ! [[ "$2" = *"fs2020"* ]]; then
			sed -i  -e '/double mkrts, mkrte;/i struct timespec mkrtsclock;' \
				-e 's/mkrts = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrts = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' \
				-e 's/mkrte = MPI_Wtime();/clock_gettime(CLOCK_MONOTONIC, \&mkrtsclock); mkrte = (mkrtsclock.tv_sec + mkrtsclock.tv_nsec * .000000001);/' $FILE;
			fi
		done

	elif [[ "$1" = *"arm"* ]]; then

		for FILE in `/bin/grep 'include.*ittnotify' -r "$2" | cut -d':' -f1 | sort -u`; do
			sed -i  -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE;
		done

	elif [[ "$1" = *"llvm12"* ]]; then

		for FILE in `/bin/grep 'include.*ittnotify' -r "$2" | cut -d':' -f1 | sort -u`; do
			sed -i  -e 's/.*include.*ittnotify\.h.*/#define __itt_resume()\n#define __itt_pause()\n#define __SSC_MARK(hex)/' $FILE;
		done

	fi
}

