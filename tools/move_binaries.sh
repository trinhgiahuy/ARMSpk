#/bin/bash
# =============================================================================
#
# Description: This script will
#   * Look for file $ROOTDIR/run/$BenchID the line `source $ROOTID/conf/${BenchID}.sh` originally
#   + If have \$1, continue
#   + If exist but missing \$1, append to it
#   + If not exist the source line, output error for manual checkk [WON'T HAPPEN]
#
#   * Look for file $ROOTDIR/conf/${BenchID}.sh the line `export BINARY=` originally.
#   + Check the corresponding binary with $2 compiler options if exist
#   + Append the move code to binary directory
#   
#
#
# Usage:
#     
#     ./move_binaries.sh conf/amg.sh gnu
#
#     
# Look for line export BINARY =
# MOVE THE LINE TO $ROOTDIR/bin
# The bin is checked in execute_top.sh, then create $ROOTDIR/bin/$bm/gnu & $ROOTDIR/bin/$bm/llvm-arm/
# r:      Huy Trinh
# Emai:        huy.trinh@a.riken.jp 
# Datue:        Sept 7, 2023
# =============================================================================
#
#
# Then look for $2 is gnu|llvm-arm then move binaries to it
#

# Check for required input
if [[ -z "$1" || -z "$2" ]]; then
  echo "Usage: $0 <input-benchmark> gnu|llvm-arm"
  exit 1
fi

ROOTDIR=$(cd "$(dirname ${BASH_SOURCE})/../" && pwd)
# Input file
# Should look like conf/amg.sh
INPUT_FILE=$1
BenchID=$(basename $INPUT_FILE .sh)
echo $BenchID
COMPILER=$2
BIN_DIR=$ROOTDIR/bin/$BM/$COMPILER/
# Temporary file for t
TEMP_FILE=$(mktemp)


export MOVE_CODE="\
mv \$BINARY ${BIN_DIR}
export BINARY=${BIN_DIR}/\$(basename \$BINARY)
echo \"NEW BINARY IS \$BINARY\""

RUN_TEST_FILE="$ROOTDIR/run/${BenchID}/test.sh"
echo "RUN_TEST_FILE: $RUN_TEST_FILE"

# Check exist $1 option in line source conf file in run file
if rg --no-ignore 'source \$\{ROOTDIR\}/conf/\$\{BenchID\}\.sh "\$1"' $RUN_TEST_FILE > /dev/null 2>&1; then
  echo "[LOG] Run file for ${BenchID} already has option \$1. Skipping.."
elif rg --no-ignore 'source \$\{ROOTDIR\}/conf/\$\{BenchID\}\.sh' $RUN_TEST_FILE > /dev/null 2>&1; then
   # sed -i 's|source \$\{ROOTDIR\}/conf/\$\{BenchID\}\.sh|& "$1"|' $RUN_TEST_FILE
  sed -i -E 's|source \$\{ROOTDIR\}/conf/\$\{BenchID\}\.sh|& "$1"|' $RUN_TEST_FILE
  echo "[LOG] Added \$1 to run file for ${BenchID}"
else
  echo "[ERROR] Pattern 'source \$\{ROOTDIR\}/conf/\$\{BenchID\}\.sh' not found in $INPUT_FILE. Exiting."
  exit 1
fi

# Check if the bin/$BenchID/$2 directory exists. If not, create it.
if [ ! -d "$BIN_DIR" ]; then
  echo "[LOG] Directory ${BIN_DIR} does not exist. Creating..."
  mkdir -p "$BIN_DIR"
else
  echo "[LOG] Directory ${BIN_DIR} exists."
fi


# Check if the file already exist the lines, 
#   export BINARY= source $ROOTDIR/conf/${BenchID}.sh
# If yes, then ignore, 
# If no, append

if ! rg --no-ignore 'NEW BINARY' $INPUT_FILE > /dev/null 2>&1;then
  echo "[LOG] COLD BINARY NOT UPDATED TO COMPILERS! ADDING CODE TO CONF FILE.."
  awk -v MOVE_CODE="$MOVE_CODE" '
  {
    if ($0 ~ /^export BINARY=/){
      print $0
      print MOVE_CODE
    } else {
      print $0
    }
  }' "${INPUT_FILE}" > "${TEMP_FILE}"

  mv "${TEMP_FILE}" "${INPUT_FILE}"
fi


# Verify that those file has expected output
if rg --no-ignore 'source \$\{ROOTDIR\}/conf/\$\{BenchID\}\.sh "\$1"' $RUN_TEST_FILE > /dev/null 2>&1; then
  echo "[SUCCESS] RUN FILE $RUN_TEST_FILE HAS \$1"
else
  echo "[ERROR] RUN FILE $RUN_TEST_FILE DOES NOT HAVE \$1!!"
  exit 1
fi

if rg --no-ignore 'NEW BINARY' $INPUT_FILE > /dev/null 2>&1; then
  echo "[SUCCESS] CONF FILE $(basename $INPUT_FILE) BINARY UPDATED TO COMPILERS"
else
  echo "[ERROR] CONF FILE $(basename $INPUT_FILE) BINARY NOT UPDATED TO COMPILERS"
  exit 1
fi
