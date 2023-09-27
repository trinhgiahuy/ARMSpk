#!/bin/bash
# =============================================================================
#
# Script Name: replica_install_mvmc.sh
# Description: This script will first xcheck/install OpenBLAS on $HOME directory.
#              Also check for installation of SCALAPACK (from Spack in my case)
#              Then it will patch first the inst/mvmc.sh script without installing MVMC yet
#              Add exit 1 before 1 to exit
#              Modify compiler flags and link to OpenBLAS, SCALAPACK
#              Finally install MVMC with OpenBLAS, SCALAPACK
# Author:      Huy Trinh
# Emai:        huy.trinh@a.riken.jp
# Date:        Sept 26, 2023
# Usage:
#
#       tools/replica_install_mvmc.sh inst/mvmc.sh gnu|llvm
# =============================================================================
#

ROOTDIR=$(cd "$(dirname "$0")/.." && pwd)
CPU_ARCH=NEOVERSEN1
ARCH=ARM64
BM="MVMC"
BM_DIR=$ROOTDIR/$BM
MAKE_ARCH_FILE="Makefile_$ARCH"
INPUT_FILE=$1

export BLAS_SCALA_CODE="\
LIB_PATH = ~/OpenBLAS/lib
INC_PATH = ~/OpenBLAS/include
BLAS = -lopenblas
SCALA_PATH := \$(shell spack location -i scalapack)/lib
"

if [[ -z "$1" || -z "$2" ]]; then
    echo "USAGE $1 <inst/hpl.sh> gnu|llvm"
    exit 1
fi
if [ ! -f $ROOTDIR/$BM/src/vmc.out ] || [ "x`ls -s $ROOTDIR/$BM/src/vmc.out | awk '{print $1}'`" = "x0" ]; then
    if [ -d ~/OpenBLAS ];then
        echo "FOUND OpenBLAS. SKIP INSTALLATION"
    else
        #INSTALL OpenBLAS TARGETING ON SUPERCOMP01
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
    last_make_line=$(grep -nE '\t*make\b' "${INPUT_FILE}" | tail -n1 | cut -d: -f1)
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

    cd $BM_DIR
    # DIFFERENT FOR EACH BENCHMARK SPECS
    # REAL UPDATE COMPILE OPTIONS FOR SUPERCOMP07 NODE
    cp -p src/Makefile_intel src/$MAKE_ARCH_FILE
    cp -p src/sfmt/Makefile_intel src/sfmt/$MAKE_ARCH_FILE
    cp -p src/pfapack/Makefile_intel src/pfapack/$MAKE_ARCH_FILE
    dos2unix src/$MAKE_ARCH_FILE

    sed -i '/^MKL\s*=/,/^-lpthread/{/^MKL\s*=/!d;s/-L.*/ /}' src/$MAKE_ARCH_FILE
    sed -i 's/^\(LIB\s*=\).*/\1 -L$(LIB_PATH) $(BLAS) -L$(SCALA_PATH) -lscalapack/' src/$MAKE_ARCH_FILE
    TEMP_FILE=$(mktemp)
    original_permission=$(stat -c %a src/"$MAKE_ARCH_FILE")
    # echo "TEST $(pwd)"
    awk -v blas_scala_code="$BLAS_SCALA_CODE" '
    {
        if ($0 ~ /^LIB/){
            print blas_scala_code
            print $0
        }else{
            print $0
        }
    }' src/${MAKE_ARCH_FILE} > "${TEMP_FILE}"
    mv "${TEMP_FILE}" src/"${MAKE_ARCH_FILE}"
    chmod $original_permission src/"${MAKE_ARCH_FILE}"

    TARGET_CODE="\
$ARCH :
\t\$(MAKE) -f Makefile_$ARCH
"
    if ! rg "ARM64" src/Makefile;then
        TEMP_FILE=$(mktemp)
        original_permission=$(stat -c %a src/Makefile)
        # Add our new target to Makefile
        awk -v target_code="$TARGET_CODE" '
        {
            if ($0 ~ /^clean/){
                print target_code
                print $0
            }else{
                print $0
            }
        }' src/Makefile > "$TEMP_FILE"
        mv "$TEMP_FILE" src/Makefile
        chmod $original_permission src/Makefile
    else
        echo "FOUND TARGET ARM64 CODE IN Makefile ALREADY"
    fi
    # sed -i '/^LAlib\s*=/,/-Wl,--end-group -lpthread -ldl/c\LAlib        = -L$(LAdir)/lib -lopenblas' $MAKE_ARCH_FILE
    # sed -i -e 's#\$(MKLROOT)#~/OpenBlas#g' -e 's#mkl/include#include#g' $MAKE_ARCH_FILE
    # sed -i 's/-m64\|-ansi-alias\|-i-static\|-nocompchk\|-mt_mpi\|-littnotify//g' $MAKE_ARCH_FILE
    source $HOME/spack/share/spack/setup-env.sh
    command -v mpicc >/dev/null || spack load openmpi@4.1.5%gcc@12.2.1 arch=linux-fedora37-neoverse_n1
    # NO NEED FOR CLEAN AS IT MAY RAISE INFINITE LOOP
    # make clean arch=$ARCH


    # HERE UPDATE ALL INCOMPATIBLE OPTIONS
    cd $ROOTDIR
    echo "CALL UPDATE"
    tools/update-options-gcc.sh "MVMC" "-xHost|-ipo|-opt-prefetch=3|-nofor-main|-vec-report|-DHAVE_SSE2|-xSSE2"

    cd $BM_DIR/src
    # -fopenmp is Intel specific flag, change to -fopenmp
    sed -i -e 's/Makefile_intel/Makefile_'${ARCH}'/g' -e 's/openmp/fopenmp/g' $MAKE_ARCH_FILE
    sed -i -e 's/-ipo -xHost/-march=native/g' -e 's/ ifort/ gfortran/' -e 's/-implicitnone/-fimplicit-none/' ./pfapack/Makefile_${ARCH}
    sed -i -e 's/-ipo -xHost/-march=native/g' -e 's/ icc/ gcc/' -e 's/-no-ansi-alias//' ./sfmt/Makefile_${ARCH}
    make $ARCH
else
    echo "[LOG] $BM binary exist! Nothing to do"
fi
cd $ROOTDIR
exit 0

#TODO: COPY BINARY VMC  TO OUR CUSTOM BIN DIR
# if [ ! -e $ROOTDIR/bin/hpl/$2/xhpl ];then
    # echo "COPYING XHPL TO bin/hpl/$2"
    # echo "$BM_DIR/bin/$ARCH/xhpl"
    # cp -p $BM_DIR/bin/$ARCH/xhpl $ROOTDIR/bin/hpl/$2/
# fi
