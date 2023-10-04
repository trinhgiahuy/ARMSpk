#!/bin/bash
# =============================================================================
# Script name: insert_copy_bin_minitri.sh
# Description: This file will modify inst/$bm file by appending code of copy
#               binary from $BM dir to bin custom dir
#
# Usage:
#       tools/insert_copy_bin_minitri inst/minitri.shh

# Author:       Huy Trinh
# Emai:	        huy.trinh@a.riken.jp
# Date:         Sept 19, 2023
# =============================================================================

if [ -z "$1" ]; then
    echo "USAGE: $0 <inst/\$bm"
    exit 1
fi

INPUT_FILE=$1

LOG_CODE="\
if [ -z \"\$1\" ]; then
    echo \"USAGE: source \$0 gnu|llvm-arm\"
    # exit for standalone script while return for source script which executes within the current shell context
    exit 1
fi

export LOG_P=\"[LOG][\$0] \"
E_LOG(){
    echo \"\$LOG_P\"
}

"
EXPORT_CODE="\
BENCHID=\$(basename \"\${BASH_SOURCE}\" .sh)
#echo BENCHID \$BENCHID
BENCHFILENAME=\$(basename \"\${BASH_SOURCE}\")
#echo \$BENCHFILENAME
CONF_FILE=\$ROOTDIR/conf/\$BENCHFILENAME
#echo \$CONF_FILE
BINARY_DIR=\$ROOTDIR/bin
#echo BINARY DIR \$BINARY_DIR

DOUBLE_QUOTE=(\$(extract_double_quote_minitri \"\$CONF_FILE\"))
MPI_BINARY_VALUE=\${DOUBLE_QUOTE[0]}
OMP_BINARY_VALUE=\${DOUBLE_QUOTE[1]}

APP_DIR=\${DOUBLE_QUOTE[2]}
MPI_BIN_ARR=\$(parse_binary \"\$MPI_BINARY_VALUE\")
OMP_BIN_ARR=\$(parse_binary \"\$OMP_BINARY_VALUE\")


echo MPI_BIN_ARR \$MPI_BIN_ARR
echo OMP_BIN_ARR \$OMP_BIN_ARR

MPI_BINEXE=\"\$ROOTDIR/\$APP_DIR/\$MPI_BIN_ARR\"
OMP_BINEXE=\"\$ROOTDIR/\$APP_DIR/\$OMP_BIN_ARR\"


if [[ ! -d \"\$BINARY_DIR/\$BENCHID/\$1\" ]]; then
    echo \"\$(E_LOG) BIN DIR FOR \$1 NOT FOUND! CREATING..\"
    mkdir -p \"\$BINARY_DIR/\$BENCHID/\$1\"
fi

"
if ! grep -q "E_LOG" "$INPUT_FILE"; then
    WORKING_DIR="$(dirname "$INPUT_FILE")"

    #TODO: Check for backup first! If not backup auto creating
    orginal_permission=$(stat -c  %a "$INPUT_FILE")
    TEMP_FILE=$(mktemp)
    awk -v logcode="$LOG_CODE" -v exportcode="$EXPORT_CODE" '
    {
        if ($0 ~ /^ROOTDIR=/) {
            print logcode
            print $0
        } else if ($0 ~ /^BM="/) {
            print exportcode
            print $0
        } else {
            print $0
        }
    }' "${INPUT_FILE}" >"${TEMP_FILE}"


    cat <<EOF >> "${TEMP_FILE}"

# CHECK IF DIRECTORY TO BINARY IS VALID
# #
# # ANOTHER USE IS ~= //
if [[ "\$MPI_BINEXE" == *//* ]];then
    echo "\$(E_LOG) BINARY \$MPI_BINEXE MAY NOT BE VALID. CHECK AGAIN "
    exit 1
elif [[ "$OMP_BINEXE" == *//* ]];then
    echo "\$(E_LOG) BINARY \$OMP_BINEXE MAY NOT BE VALID. CHECK AGAIN "
    exit 1
else
    echo "\$(E_LOG) BOTH MPI & OMP BINARIES MAY BE VALID AT THIS TIME."
fi

# CHECK IF CMAKE SUCCESS AND EXECUTABLE BM EXIST
echo \$MPI_BINEXE
echo \$OMP_BINEXE

if [[ ! -x \$MPI_BINEXE || ! -x \$OMP_BINEXE ]];then
    echo "\$(E_LOG) EXECUTABLE IN BENCHMARK FAIL. EXITING.."
    exit 1
else
    echo "\$(E_LOG) FOUND EXECUTABLE OR BUILD SUCCESS."
fi

# Check MPI and OpenMP Dir
if [ ! -d \$BINARY_DIR/\$BENCHID/\$1/MPI ];then
    mkdir -p \$BINARY_DIR/\$BENCHID/\$1/MPI
fi

if [ ! -d \$BINARY_DIR/\$BENCHID/\$1/OMP ];then
    mkdir -p \$BINARY_DIR/\$BENCHID/\$1/OMP
fi

if  [ ! -e \$BINARY_DIR/\$BENCHID/\$1/MPI/\$(basename \$MPI_BINEXE) ];then
    echo "\$(E_LOG) REPLICA BIN FOR \$(basename \$MPI_BINEXE) NOT FOUND! COPYING TO bin/\$BENCHID/\$1"
    cp -p \$MPI_BINEXE \$BINARY_DIR/\$BENCHID/\$1/MPI
else
    echo "\$(E_LOG) NO NEED REPLICA BIN FOR MPI \$(basename \$MPI_BINEXE). DONE!!"
fi


if  [ ! -e \$BINARY_DIR/\$BENCHID/\$1/OMP/\$(basename \$OMP_BINEXE) ];then
    echo "\$(E_LOG) REPLICA BIN FOR \$(basename \$OMP_BINEXE) NOT FOUND! COPYING TO bin/\$BENCHID/\$1"
    cp -p \$OMP_BINEXE \$BINARY_DIR/\$BENCHID/\$1/OMP
else
    echo "\$(E_LOG) NO NEED REPLICA BIN FOR OMP \$(basename \$OMP_BINEXE). DONE!!"
fi

EOF
    # Replace the original file with the modified one
    mv "${TEMP_FILE}" "${INPUT_FILE}"
    chmod "$orginal_permission" "$INPUT_FILE"
else
    echo "CODE COPY BINARY TO BIN DIR ALREADY EXIST!"
fi


