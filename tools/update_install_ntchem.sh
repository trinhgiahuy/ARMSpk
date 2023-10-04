#!/bin/bash
# =============================================================================
#
# Script Name: replica_install_nicam.sh
# Description: This script will update/modify necessary installing option for ARM chip
#							 on original inst file.
# Author:      Huy Trinh
# Emai:        huy.trinh@a.riken.jp
# Date:        Oct 4, 2023
# Usage:
#
#       tools/update_install_ntchem.sh inst/ntchem.sh gnu|llvm
# =============================================================================
#

source $HOME/spack/share/spack/setup-env.sh
# ACTIVATE SPACK ENVIRONMENT FOR SUPERCOMP07 NODE
spack env activate supercomp07-env-ffb
command -v mpicc >/dev/null || spack load openmpi@4.1.5%gcc@12.2.1 arch=linux-fedora37-neoverse_n1


ROOTDIR=$(cd "$(dirname "$0")/.." && pwd)
CPU_ARCH=NEOVERSEN1
ARCH=ARM64
BM="NTChem"
BM_DIR=$ROOTDIR/$BM
INPUT_FILE=$1
COMPILER_OPT=$2
UPDATE_ARM_MSG="UPDATED INSTALL CODE FOR ARM"

# FOR GNU LAST MAKE, THIS DEPENDS ON BM
# CHANGE THIS FOR LLVM
MAKE_EXIT_REVERSE_ORDER=4
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

if [[ -z "$1" || -z "$2" ]]; then
    echo "USAGE $1 <inst/hpl.sh> gnu|llvm"
    exit 1
fi

cd $ROOTDIR
if ! rg 'UPDATED|ARM' $INPUT_FILE >> /dev/null;then 
	last_make_line=$(grep -nE '\t*make\b' "${INPUT_FILE}" | tail -n${MAKE_EXIT_REVERSE_ORDER} | head -n1 | cut -d: -f1)
	echo last_make_line $last_make_line
	TEMP_FILE=$(mktemp)
	org_permission=$(stat -c %a "${INPUT_FILE}")
	echo "DO NOT FIND EXIT.ADD"
	echo last_make_line $last_make_line
	sed -i "${last_make_line}s/\(else\)/\1 echo \"UPDATED INSTALL CODE FOR ARM\"; ${PRIOR_MAKE_CODE}; /" "$INPUT_FILE"

	# ADD LINK BLAS CODE
	sed -i "/\.\/config_mine/ i ${LINK_BLAS_CODE}" $INPUT_FILE
else
	echo "FIND UPDATE ARM MSG!"
fi

# DO NOT NEED FOR THIS UPDATE SCRIPT AS THE insert_copy_bin.sh WILL APPEND THE MOVING CODE
# BIN_NAME="rimp2.exe"
# if [ ! -e $ROORDIR/bin/ntchem/$2/$BIN_NAME ]; then
		# echo "COPYING $BIN_NAME TO bin/ntchem/$2"
		# cp -p $BM_DIR/bin/$BIN_NAME $ROOTDIR/bin/ntchem/$2/
# fi
