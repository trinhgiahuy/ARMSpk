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
#       tools/update_install_spec_cpu.sh inst/spec_cpu.sh gnu|llvm
# =============================================================================
#
# ADDITIONAL INFORMATION: It will coment out the fuseiso part as iso mounting requires root but only mount to /tmp directory
#                         Assuming the cluster has /tmp/mnt_SPEC_[CPU|OMP]

# [WORKING]


ROOTDIR=$(cd "$(dirname "$0")/.." && pwd)
CPU_ARCH=NEOVERSEN1
ARCH=ARM64
BM="SPEC_CPU"
BM_DIR=$ROOTDIR/$BM
INPUT_FILE=$1
COMPILER_OPT=$2
UPDATE_ARM_MSG="UPDATED INSTALL CODE FOR ARM"
UPDATE_ARM_CODE="\
        echo \"${UPDATE_ARM_MSG}\""
# FOR GNU LAST MAKE, THIS DEPENDS ON BM
# CHANGE THIS FOR LLVM
MAKE_EXIT_REVERSE_ORDER=1
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg

GNU_UPDATE_CODE="\
        # REMOVE -m64 option in GNU installation nedo config file
        grep -A 5 -nEx \"%if '%{COMP}' eq 'gnu'\" config/nedo.cfg | grep m64 | cut -d'-' -f1 | while read -r line;do
            sed -i \"\$line s/ -m64//\" config/nedo.cfg
        done"

cd $ROOTDIR
if ! rg 'UPDATED|ARM' $INPUT_FILE >> /dev/null;then
    gnu_elif_line=$(grep -nE '\t*elif.*\$1.*gnu' $INPUT_FILE | cut -d: -f1)
    echo "gnu_elif_line $gnu_elif_line"
    if [ -n "$gnu_elif_line" ];then
        gawk -i inplace -v gnu_elif_line="$gnu_elif_line" -v gnu_update_code="$GNU_UPDATE_CODE" 'NR==gnu_elif_line {print $0; print gnu_update_code ; next} 1' $INPUT_FILE
        # 1 line before next elif line from elif [[ "$1" = *"gnu"* ]]
        gnu_clause_end_line=$(awk -v start="$gnu_elif_line" 'NR > start && /\t*elif \[\[ "\$1" = \*/ {print NR-1; exit}' $INPUT_FILE)
        gawk -i inplace -v gnu_clause_end_line="$gnu_clause_end_line" -v update_arm_msg="$UPDATE_ARM_CODE" 'NR==gnu_clause_end_line {print $0; print update_arm_msg; next} 1' $INPUT_FILE
    fi
fi

exit 1

#TODO: CHECK THIS: SOMEHOW INSTALL THE SPEC_CPU IS ONE HIRECHY UP FROM, NEED COPY


