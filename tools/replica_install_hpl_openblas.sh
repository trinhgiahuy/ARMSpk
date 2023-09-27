#!/bin/bash
# =============================================================================
#
# Script Name: replica_install_hpl_openblas.sh
# Description: This script will first check/install OpenBLAS on $HOME directory
#              Then it will patch first the inst/hpl/sh script without installing HPL yet
#              Add exit 1 before 1 to exit
#              Modify compiler flags and link to OpenBLAS
#              Finally install HPL with OpenBLAS
# Author:      Huy Trinh
# Emai:        huy.trinh@a.riken.jp
# Date:        Sept 22, 2023
# Usage:
#
#       tools/replica_install_hpl_openblas.sh inst/hpl.sh gnu|llvm
# =============================================================================
#

ROOTDIR=$(cd "$(dirname "$0")/.." && pwd)
CPU_ARCH=NEOVERSEN1
ARCH=Linux_Intel64
BM="HPL"
BM_DIR=$ROOTDIR/$BM
MAKE_ARCH_FILE="Make.$ARCH"
INPUT_FILE=$1


if [[ -z "$1" || -z "$2" ]]; then
    echo "USAGE $1 <inst/hpl.sh> gnu|llvm"
    exit 1
fi

if [ ! -f $ROOTDIR/$BM/bin/$ARCH/xhpl ];then
    if [ -d ~/OpenBLAS ];then
      echo "FOUND OpenBLAS. SKIP INSTALLATION"
    else
      ##INSTALL OpenBLAS TARGETING ON SUPERCOMP01
      cd ~
      git clone https://github.com/OpenMathLib/OpenBLAS.git
      cd OpenBLAS

      ##CHECK THE KERNEL DIR FOR ARCH
      make TARGET=$CPU_ARCH
      make PREFIX=`pwd` install

      ##CHECK INSTALLATION COMPLETE BY EXISTANCE OF lib/libopenblas.a
      if [ ! -e ~/OpenBLAS/lib/libopenblas.a ];then
        echo "[ERROR] CANNOT FIND LIBRARY FILE FOR OpenBLAS. CHECK INSTALLATION AGAIN!"
        exit 1
      fi
    fi

    cd
    cd $ROOTDIR

    # Add exit 1 to before make line in Jens code
    last_make_line=$(grep -nE '\t*make\b' "$INPUT_FILE" | tail -n1 | cut -d: -f1)
    echo last_make_line $last_make_line
    TEMP_FILE=$(mktemp)
    org_permission=$(stat -c %a "${INPUT_FILE}")
    if ! rg 'Adding exit' $INPUT_FILE > /dev/null 2>&1; then
        echo "DO NOT FIND EXIT.ADD"
        echo last_make_line $last_make_line
        awk -v n="$last_make_line" 'NR == n {print "\techo \"Adding exit1\"; exit 1"} 1' "$INPUT_FILE" > "$TEMP_FILE" && mv "$TEMP_FILE" "$INPUT_FILE"
        # awk -v n="$last_make_line" 'NR == n && !/exit 1/ {print "\techo \"Adding exit1\"; exit 1"} 1' "$INPUT_FILE" > "$TEMP_FILE" && mv "$TEMP_FILE" "$INPUT_FILE"
        chmod $org_permission "$INPUT_FILE"
    else
        echo "FIND EXIT"
    fi

    # Call "FAKE" install to apply patches and get latest Jens' modification
    $1 $2


    # DIFFERENT FOR EACH BENCHMARK SPECS
    # REAL UPDATE COMPILE OPTIONS FOR SUPERCOMP07 NODE
    cd
    if [ ! -d ~/hpl ];then
        echo "CREATING HPL dir in home directory "
        mkdir ~/hpl
    fi

    # Soft symbolic link to ROOTDIR
    echo $MAKE_ARCH_FILE

    cd ~/hpl
    # ln -s $BM_DIR/$MAKE_ARCH_FILE ./$MAKE_ARCH_FILE
    [ ! -e "./$MAKE_ARCH_FILE" ] && ln -s "$BM_DIR/$MAKE_ARCH_FILE" "./$MAKE_ARCH_FILE" || echo "Symbolic link ./$MAKE_ARCH_FILE already exists."
    [ ! -d "./bin" ] && ln -s $BM_DIR/bin ./bin
    [ ! -d "./include" ] && ln -s $BM_DIR/include ./include
    [ ! -d "./lib" ] && ln -s $BM_DIR/lib ./lib


    # LINK TO OPENBLAS LIBRARY
    cd $BM_DIR
    pwd
    sed -i '/^LAlib\s*=/,/-Wl,--end-group -lpthread -ldl/c\LAlib        = -L$(LAdir)/lib -lopenblas' $MAKE_ARCH_FILE
    sed -i -e 's#\$(MKLROOT)#~/OpenBlas#g' -e 's#mkl/include#include#g' $MAKE_ARCH_FILE
    # REMOVE INTEL COMPILER FLAGS OPTION
    # -m64: generate code for a 64-bit environment
    # -ansi-alias: assume that the program adheres to the ANSI aliasability rules, no data aliasing, enable aggressive optimization
    # -i-static: statisticallyt link to Intel-provided libraries
    # -nocompchk: disable run-time check for compability between Intel Compiler version and Intel MKL
    # -mt_mpi: enable multi-thread libraries in MPI applicationwhen use Intel MPI Lib
    # -littnotify: associate with Intel's performance profiling tool (Intel VTune Amplifier
    sed -i 's/-m64\|-ansi-alias\|-i-static\|-nocompchk\|-mt_mpi\|-littnotify//g' $MAKE_ARCH_FILE
    source $HOME/spack/share/spack/setup-env.sh
    command -v mpicc >/dev/null || spack load openmpi@4.1.5%gcc@12.2.1 arch=linux-fedora37-neoverse_n1
    # NO NEED FOR CLEAN AS IT MAY RAISE INFINITE LOOP
    # make clean arch=$ARCH
    make arch=$ARCH
else
    echo "[LOG] XHPL binary exist! Nothing to do"
fi
cd $ROOTDIR

# COPY BINARY HPL TO OUR CUSTOM BIN DIR
if [ ! -e $ROOTDIR/bin/hpl/$2/xhpl ];then
    echo "COPYING XHPL TO bin/hpl/$2"
    echo "$BM_DIR/bin/$ARCH/xhpl"
    cp -p $BM_DIR/bin/$ARCH/xhpl $ROOTDIR/bin/hpl/$2/
fi
