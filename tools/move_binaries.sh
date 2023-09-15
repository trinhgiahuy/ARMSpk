#/bin/bash
# =============================================================================
#
# Description: This script will
#   * Look for file $ROOTDIR/run/$BENCH_ID/test.sh tthe line `source $ROOTID/conf/${BENCH_ID}.sh` originally
#   + If have \$1, continue
#   + If exist but missing \$1, append to it
#   + If not exist the source line, output error for manual checkk [WON'T HAPPEN]
#
#   * Look for file $ROOTDIR/conf/${BENCH_ID}.sh the line `export BINARY=` originally.
#   + Check the corresponding binary with compiler options if exist
#   + Append the move code to binary directory
#
# Usage:
#
#     ./move_binaries.sh amg gnu
#
# Look for line export BINARY =
# MOVE THE LINE TO $ROOTDIR/bin
# The bin is checked in execute_top.sh, then create $ROOTDIR/bin/$bm/gnu & $ROOTDIR/bin/$bm/llvm-arm/
# Author:       Huy Trinh
# Emai:         huy.trinh@a.riken.jp
# Datue:        Sept 9, 2023
# =============================================================================
#
#
# Then look for $2 is gnu|llvm-arm then move binaries to it

#
E_LOG() {
    echo "$(eval echo $LOG_P)"
}

# Check for required input
if [[ -z "$1" || -z "$2" ]]; then
  echo "Usage: $0 <input-benchmark> gnu|llvm-arm"
  exit 1
fi

ROOTDIR=$(cd "$(dirname ${BASH_SOURCE})/../" && pwd)
# Input file
# Should look like amg
BENCH_ID=$1
CONF_FILE="$ROOTDIR/conf/${BENCH_ID}.sh"

# EXTRACT BINARY NAM EINSIDE "" INSIDE export BINARY= in CONF FILE
#TODO: ADD CASE OF MULTIPLE BINARIES NAME
source $ROOTDIR/conf/env.cfg
echo conffile $CONF_FILE
# DOUBLE_QUOTE=($(extract_double_quote "$CONF_FILE"))
DOUBLE_QUOTE=($(extract_double_quote "$CONF_FILE"))
# echo DOUBLE_QUOTE $DOUBLE_QUOTE
# echo ${DOUBLE_QUOTE[0]}
# echo ${DOUBLE_QUOTE[1]}
# echo ${DOUBLE_QUOTE[2]}
BINARYNAME=$(basename ${DOUBLE_QUOTE[0]})

# echo BINNAME $BINARYNAEM
# echo $BENCH_ID
COMPILER=$2

# HERE BIN_DIR is bin/b$BENCH_ID/$1 IS CORRECT!
# THIS BIN_DIR is WRITE TO CONFIG FILE WHEN SOURCE
# WILL TAKE $1 AS COMPILER AND GET BINARY FROM OUR BIN DIR
#
BIN_DIR=$ROOTDIR/bin/$BENCH_ID/\$1
BIN_DIR_CHECK=$ROOTDIR/bin/$BENCH_ID/$COMPILER

echo ${BIN_DIR}/$BINARYNAME
echo ${BIN_DIR_CHECK}

export MOVE_CODE="\
export BINARY=\"${BIN_DIR}/$BINARYNAME\"
echo \"NEW BINARY IS \$BINARY\""

export CHECK_PARAM_CODE="\
if [ -z \"\$1\" ]; then
    echo \"SOURCE CFG FILE REQUIRES INPUT gnu|llvm-arm PARAM\"
    return 1
fi

"


RUN_TEST_FILE="$ROOTDIR/run/${BENCH_ID}/test.sh"
# echo "RUN_TEST_FILE: $RUN_TEST_FILE"

# ALMOST NOT HAPPEN BUT STILL NEED TO CHECK
if ! rg --no-ignore 'source \$\{ROOTDIR\}/conf/\$\{BenchID\}.sh' $RUN_TEST_FILE > /dev/null 2>&1;then
    echo [ERROR] MISSING THE SOURCE CONFIG FILE IN RUN FILE! MANUALLY CHECK.
    exit 1
else
    echo ""
fi

# Check exist $1 option in line source conf file in run file
# if rg --no-ignore 'source \$\{ROOTDIR\}/conf/\$\{BenchID\}.sh "\$1"' $RUN_TEST_FILE > /dev/null 2>&1; then
if rg --no-ignore 'source \$\{ROOTDIR\}/conf/\$\{BenchID\}.sh ..' $RUN_TEST_FILE > /dev/null 2>&1; then
    echo "$(E_LOG) Run file for ${BENCH_ID} already has option \$1. Skipping.."
elif rg --no-ignore 'source \$\{ROOTDIR\}/conf/\$\{BenchID\}\.sh' $RUN_TEST_FILE > /dev/null 2>&1; then
   # sed -i 's|source \$\{ROOTDIR\}/conf/\$\{BenchID\}\.sh|& "$1"|' $RUN_TEST_FILE
  sed -i -E 's|source \$\{ROOTDIR\}/conf/\$\{BenchID\}\.sh|& "$1"|' $RUN_TEST_FILE
  echo "$(E_LOG) [run/${BENCH_ID}/test.sh] ADDING \$1 TO SOURCE LINE"
else
  echo "[ERROR] Pattern 'source \$\{ROOTDIR\}/conf/\$\{BenchID\}\.sh' not found in $CONF_FILE. Exiting."
  exit 1
fi

# Check if the bin/$BenchID/$2 directory exists. If not, create it.
if [ ! -d "$BIN_DIR_CHECK" ]; then
  echo "[LOG] Directory ${BIN_DIR} does not exist. Creating..."
  mkdir -p "$BIN_DIR"
else
  echo "[LOG] Directory ${BIN_DIR_CHECK} exists."
fi


original_permission=$(stat -c %a "$CONF_FILE")
# Check if exit param check code at the begining
if ! rg --no-ignore 'SOURCE CFG FILE' $CONF_FILE > /dev/null 2>&1;then
    echo "$(E_LOG) CODE PARAM CHECK NOT ADDED! ADDING CODE CHECK PARAM TO CONF FILE"
    TEMP_FILE=$(mktemp)
    awk -v PARAM_CHECK="$CHECK_PARAM_CODE" '
    {
        if ($0 ~ /^export APPDIR=/){
            print PARAM_CHECK
            print $0
        }else{
            print $0
        }
    }' "${CONF_FILE}" > "${TEMP_FILE}"
    mv "${TEMP_FILE}" "${CONF_FILE}"
    chmod $original_permission "$CONF_FILE"
else
    echo "PARAM CODE EXIST! NO ADD"
fi
# Check if the file already exist the lines,
#   export BINARY= source $ROOTDIR/conf/${BenchID}.sh
# If yes, then ignore,
# If no, append

if ! rg --no-ignore 'NEW BINARY' $CONF_FILE > /dev/null 2>&1;then
    echo "$(E_LOG) CODE BINARY NOT UPDATED TO COMPILERS! ADDING CODE MOVE BINARY TO CONF FILE.."
    TEMP_FILE=$(mktemp)
    awk -v MOVE_CODE="$MOVE_CODE" '
    {
        if ($0 ~ /^export BINARY=/){
            print "#" $0
            print MOVE_CODE
        } else {
            print $0
        }
    }' "${CONF_FILE}" > "${TEMP_FILE}"
    mv "${TEMP_FILE}" "${CONF_FILE}"
    chmod $original_permission "$CONF_FILE"
else
    echo "$(E_LOG) CONF FILE IS UPDATED TO COMPILERS! NO ADD"
fi


# Verify that those file has expected output
if rg --no-ignore 'source \$\{ROOTDIR\}/conf/\$\{BenchID\}\.sh "\$1"' $RUN_TEST_FILE > /dev/null 2>&1; then
  echo "[SUCCESS] RUN FILE $RUN_TEST_FILE HAS \$1"
else
  echo "[ERROR] RUN FILE $RUN_TEST_FILE DOES NOT HAVE \$1!!"
  exit 1
fi


if rg --no-ignore 'NEW BINARY' $CONF_FILE > /dev/null 2>&1; then
  echo "[SUCCESS] CONF FILE ${CONF_FILE} BINARY UPDATED TO COMPILERS"
else
  echo "[ERROR] CONF FILE ${CONF_FILE} BINARY NOT UPDATED TO COMPILERS"
  exit 1
fi
