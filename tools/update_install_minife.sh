#!/bin/bash
# =============================================================================
#
# Script Name: update_install_minife.sh
# Description: This script will update/modify necessary installing option for ARM chip
#							 on original inst file.
# Author:      Huy Trinh
# Emai:        huy.trinh@a.riken.jp
# Date:        Oct 5, 2023
# Usage:
#
#       tools/update_install_minife.sh inst/minife.sh gnu|llvm
# =============================================================================
#

source $HOME/spack/share/spack/setup-env.sh
# ACTIVATE SPACK ENVIRONMENT FOR SUPERCOMP07 NODE
spack env activate supercomp07-env-ffb
command -v mpicc >/dev/null || spack load openmpi@4.1.5%gcc@12.2.1 arch=linux-fedora37-neoverse_n1


ROOTDIR=$(cd "$(dirname "$0")/.." && pwd)
CPU_ARCH=NEOVERSEN1
ARCH=ARM64
BM="MiniFE"
BM_DIR=$ROOTDIR/$BM
INPUT_FILE=$1
COMPILER_OPT=$2
UPDATE_ARM_MSG="UPDATED INSTALL CODE FOR ARM"

# FOR GNU LAST MAKE, THIS DEPENDS ON BM
# CHANGE THIS FOR LLVM
MAKE_EXIT_REVERSE_ORDER=1
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg

# CODE TO ADD TO INSTALL FILE

# LINK AGAIST OPENBLAS .so FILES
# OpenBLAS provides both BLAS and LAPACK functionalities/implementations
# IN NTChem CASE: dcopy_ function (BLAS) and _ dpotrf_ and dtrtri_ functions (LAPACK)
LINK_BLAS_CODE="	sed -i -e 's|--lapack=|&'\${HOME}'/OpenBLAS/lib/libopenblas.so|' -e 's|--blas=|&'\$HOME'/OpenBLAS/lib/libopenblas.so|' config_mine"


PRIOR_MAKE_CODE="\
sed -i 's/-littnotify//g' src/mp2/GNUmakefile;
\$ROOTDIR/tools/update-options-gcc.sh \$BM '' \"-m64|-littnotify\""
# GET RID OF ERROR OF / \n ' ' in sed
PRIOR_MAKE_CODE=$(echo "${PRIOR_MAKE_CODE}" | sed 's|/|\\/|g' | tr '\n' ' ')

# The double backslashes \\\\ will be interpreted by the shell as \\, and then awk will interpret them as a single backslash \, which is what we want for the sed command.
GNU_UPDATE_CODE="\
            if [[ \"\$SUB\" = *\"mkl\"* ]] || [[ \"\$SUB\" = *\"knl\"* ]]; then continue; fi
                echo \"GO HERE \"
                sed -i -e 's/-ipo -x[a-zA-Z0-9\\\\-]*/-flto -march=native -fopenmp/g' -e 's# -I\${ADVISOR_2018_DIR}/include##g' -e 's#-L\${ADVISOR_2018_DIR}/lib64 -littnotify##g' -e 's# -mavx##g' ./Makefile"


MOVE_BINARY_CODE="\
        if [[ \"\$1\" = *\"gnu\"* ]]; then cp miniFE.x ../../mkl/src/; fi"
if [[ -z "$1" || -z "$2" ]]; then
    echo "USAGE $1 <inst/hpl.sh> gnu|llvm"
    exit 1
fi


cd $ROOTDIR
if ! rg 'UPDATED|ARM' $INPUT_FILE >> /dev/null;then
    gnu_elif_line=$(grep -nE '\t*elif' inst/minife.sh | grep "gnu" | cut -d: -f1)
    comment_start_line=$((gnu_elif_line + 1))

    if [ -n "$comment_start_line" ];then
        # 1 line before next elif line from elif [[ "$1" = *"gnu"* ]]
        comment_end_line=$(awk -v start="$comment_start_line" 'NR > start && /\t*elif \[\[ "\$1" = \*/ {print NR-1; exit}' $INPUT_FILE)
        echo "$comment_start_line, $comment_end_line"

        if [ -n "$comment_end_line" ];then
            # Comment out old elif code (which use MKL) and add update code ( should not rely on MKL)
            gawk -i inplace -v start="$comment_start_line" -v end="$comment_end_line" 'NR >= start && NR <= end {print "#"$0; next} 1' $INPUT_FILE


            TMP_FILE=$(mktemp)
            org_permission=$(stat -c %a $INPUT_FILE)
            # Append update code
            awk -v start="$comment_start_line" -v gnu_update_code="$GNU_UPDATE_CODE" 'BEGIN {FS="\n"}
                NR == start {
                    print gnu_update_code
                    print $0
                }
                NR != start {
                    print $0
                } ' $INPUT_FILE > $TMP_FILE
            mv "$TMP_FILE" "$INPUT_FILE"
            chmod $org_permission $INPUT_FILE
        fi

        last_make_line=$(grep -nE '\t*make\b' "${INPUT_FILE}" | tail -n${MAKE_EXIT_REVERSE_ORDER} | head -n1 | cut -d: -f1)
        insert_line=$((last_make_line + 1))

        TMP_FILE=$(mktemp)
        org_permission=$(stat -c %a $INPUT_FILE)
        awk -v insert_line="$insert_line" -v move_bin_code="$MOVE_BINARY_CODE" '
            NR == insert_line {
                print move_bin_code
                print $0
            }
            NR != insert_line {
                print $0
            } ' $INPUT_FILE > $TMP_FILE
        mv "$TMP_FILE" "$INPUT_FILE"
        chmod $org_permission "$INPUT_FILE"

        sed -i "${last_make_line}s/\(make\)/\1 ;echo \"UPDATED INSTALL CODE FOR ARM\"/" "$INPUT_FILE"
    fi
else
	echo "FIND UPDATE ARM MSG!"
fi

# DO NOT NEED FOR THIS UPDATE SCRIPT AS THE insert_copy_bin.sh WILL APPEND THE MOVING CODE
# BIN_NAME="rimp2.exe"
# if [ ! -e $ROORDIR/bin/ntchem/$2/$BIN_NAME ]; then
		# echo "COPYING $BIN_NAME TO bin/ntchem/$2"
		# cp -p $BM_DIR/bin/$BIN_NAME $ROOTDIR/bin/ntchem/$2/
# fi
