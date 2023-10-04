#!/bin/bash
# =============================================================================
#
# Script Name: replica_install_nicam.sh
# Description: This script will first
# Author:      Huy Trinh
# Emai:        huy.trinh@a.riken.jp
# Date:        Sept 27, 2023
# Usage:
#
#       tools/replica_install_nicam.sh inst/nicam.sh gnu|llvm
# =============================================================================
#

source $HOME/spack/share/spack/setup-env.sh
# ACTIVATE SPACK ENVIRONMENT FOR SUPERCOMP07 NODE
spack env activate supercomp07-env-ffb
command -v mpicc >/dev/null || spack load openmpi@4.1.5%gcc@12.2.1 arch=linux-fedora37-neoverse_n1



ROOTDIR=$(cd "$(dirname "$0")/.." && pwd)
CPU_ARCH=NEOVERSEN1
ARCH=ARM64
BM="NICAM"
BM_DIR=$ROOTDIR/$BM
# MAKE_ARCH_FILE="Makefile_$ARCH"
INPUT_FILE=$1
COMPILER_OPT=$2

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg

# IMPORTANT: FOR NICAM, APPEND CONFIG BEFORE INST
tools/append_config.sh conf/nicam.sh

if [[ -z "$1" || -z "$2" ]]; then
    echo "USAGE $1 <inst/hpl.sh> gnu|llvm"
    exit 1
fi
if [ ! -f $BM_DIR/bin/driver-dc ]; then
    cd $ROOTDIR

    # Add exit 1 to before make line in Jens code,
    # DEPENDS ON EACH BM, which is the SOONEST EXIT
    # FOR NICAM, we want the 2nd to last line
    last_make_line=$(grep -nE '\t*make\b' "${INPUT_FILE}" | tail -n2 | head -n1 | cut -d: -f1)
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
    # Copy from gnu Jens modified code instead of orignial Linux64-gnu-openmpi
    [ ! -e $BM_DIR/sysdep/Makedef.$ARCH ] && touch $BM_DIR/sysdep/Makedef.${ARCH}
    cp -p $BM_DIR/sysdep/Makedef.Linux64-intel-impi $BM_DIR/sysdep/Makedef.${ARCH}
    echo "[LOG] Finish call fake install"
    echo "BINARY ORG DOES NOT EXIST. MAKING!"
    cd $BM_DIR

    rg -l 'NICAM_SYS' --glob '!README.md' | while read -r file; do
        sed -i "s/\$(NICAM_SYS)/$ARCH/g" "$file"
        sed -i "s/\${NICAM_SYS}/$ARCH/g" "$file"
    done
    $ROOTDIR/tools/append_flag_opt.sh "-fallow-argument-mismatch" "FFLAGS_FAST" $BM_DIR/sysdep/Makedef.${ARCH}
    #
    # WARNING: Options: -heap-arrays if subtitute to , when compile when get segmentation error (kind of related to memory issue or something? Investigating on this"
    # $ROOTDIR/tools/update-options-gcc.sh "NICAM" 'mpiifort|mpiicc# -fpp3| -ip| -xHost| -convert big_endian| -fp-model precise| -heap-arrays| -fno-alias' "-assume|byterecl|-ftz|-pc\s80|-shared-intel"

    # KEEP THIS AS THIS COMMAND WORKS
    $ROOTDIR/tools/update-options-gcc.sh "NICAM" 'mpiifort|mpiicc# -fpp3| -ip| -xHost| -convert big_endian| -fp-model precise' "-assume|byterecl|-ftz|-pc\s80|-shared-intel|-heap-arrays|-fno-alias"
    cd src
    make ENABLE_OPENMP=1 -j
    cd '../test/case/jablonowski'
    make ENABLE_OPENMP=1 -j

    # ADDED THESE
    cp -p ../../../sysdep/Mkjobshell.Linux64-intel-impi.sh ../../../sysdep/Mkjobshell.${ARCH}.sh
    make jobshell
    if [ -e gl05rl00z40pe10/nhm_driver.cnf ]; then
        echo "Existing config in jablonowski test case"
        ln -s ./gl05rl00z40pe10/nhm_driver.cnf .
        # ln -s ./gl05rl00z40pe10/nh_driver .
        ln -s ../../../bin/nhm_driver .
     else
        echo "No exist conf file in jablonowski test case"
    fi
    cd $ROOTDIR

    echo "COMPIL:COMPILER_OPT $COMPILER_OPT"
    source $ROOTDIR/conf/nicam.sh $COMPILER_OPT
    echo "TESTCONF: $TESTCONF"
    subOMP="$(for C in ${TESTCONF}; do echo ${C} | cut -d'|' -f2;done | sort -g -u)"
    echo "subOMP: $subOMP"
    for NumOMP in $subOMP;do
        cd $BM_DIR
        if [ ! -f $BM_DIR/omp${NumOMP}/bin/driver-dc ]; then
            if [[ "$2" = *"gnu"* ]];then
                # pwd
                # no point of creating multiple for gnu
                if [ -n "$(find $ROOTDIR/$BM/ -type f -executable -path '*/bin/driver-dc')" ];then
                    tmpDIR=$ROOTDIR/$BM/omp${NumOMP}
                    echo $tmpDIR
                    [ ! -d $tmpDIR ] && mkdir -p $tmpDIR

                    # Copy the whole directory itself to omp
                    # head -1 here exclude binaries of omp dirs
                    source_dir="$(readlink -f $(dirname $(dirname $(find $ROOTDIR/$BM/ -type f -executable -path '*/bin/driver-dc' | head -1))))"
                    rsync -aq --exclude='omp[0-9]' --exclude='omp[0-9][0-9]' --exclude='omp[0-9][0-9][0-9]' "$source_dir/" "${tmpDIR}"
                    cd $tmpDIR
                    make jobshell > /dev/null 2>&1
                    cd -
                    # Copy Makedef.$ARCH to each sub directories
                    # cp -p $BM_DIR/sysdep/Makedef.${ARCH} $tmpDIR/sysdep/
                    continue
                fi
            fi
        fi
    done
fi


if [ ! -e $ROORDIR/bin/nicam/$2/nhm_driver ]; then
    echo "COPYING nhm_driver TO bin/nicam/$2"
    cp -p $BM_DIR/bin/nhm_driver $ROOTDIR/bin/nicam/$2/
fi
