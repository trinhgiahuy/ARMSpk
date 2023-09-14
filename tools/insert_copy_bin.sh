#!/bin/bash
# =============================================================================
# Script name: insert_copy_bin.sh
# Description: This file will modify inst/$bm file by appending code of copy
#               binary from $BM dir to bin custom dir
#
# Usage:
#       tools/insert_copy_bin inst/amg.sh

# Author:       Huy Trinh
# Emai:	        huy.trinh@a.riken.jp
# Date:         Sept 9, 2023
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


DOUBLE_QUOTE=(\$(extract_double_quote \"\$CONF_FILE\"))
BINARY_VALUE=\${DOUBLE_QUOTE[0]}
APP_DIR=\${DOUBLE_QUOTE[1]}
BIN_ARR=\$(parse_binary \"\$BINARY_VALUE\")
echo BIN_ARR \$BIN_ARR
BINEXE=\"\$ROOTDIR/\$APP_DIR/\$BIN_ARR\"

# ?? NECESSARY Check for an existing executable binary in BM dir
# if [[ ! -x \"\$BINEXE\" ]]; then
#     echo \"[ERROR][\$0] NOT FOUND/EXECUTABLE \$BINEXE. CALL BUILDING..\"
#         # BUILD_BM=true
#         fi
#
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

if  [ ! -e \$BINARY_DIR/\$BENCHID/\$1/\$(basename \$BINEXE) ];then
    echo "\$(E_LOG) REPLICA BIN FOR \$(basename \$BINEXE) NOT FOUND! COPYING TO bin/\$BENCHID/\$1"
    cp -p \$BINEXE \$BINARY_DIR/\$BENCHID/\$1
else
    echo "\$(E_LOG) NO NEED REPLICA BIN FOR \$(basename \$BINEXE). DONE!!"
fi

EOF
    # Replace the original file with the modified one
    mv "${TEMP_FILE}" "${INPUT_FILE}"
    chmod "$orginal_permission" "$INPUT_FILE"
else
    echo "CODE COPY BINARY TO BIN DIR ALREADY EXIST!"
fi


