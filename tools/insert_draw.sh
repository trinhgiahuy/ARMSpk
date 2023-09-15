#!/bin/bash
# =============================================================================
# Script Name: insert_draw.sh
# Description: This script will add progress bar code for Benchmark running test
# Additional add:
#
#       The script will add checl log file code to run file. In some cases,
#       run file cannot automatically create new log file
#
# Usage:
#
# 	    insert_draw.sh $ROOTDIR/run/amg/test.sh
#
# Example:
#
#       tools/insert_draw run/ffvc/test.sh
#
# Author:       Huy Trinh
# Emai:			huy.trinh@a.riken.jp
# Date:         Sept 7, 2023
# =============================================================================


if [ -z "$1" ]; then
    echo "Usage: $0 <input-file>"
    exit 1
fi

INPUT_FILE=$1

CHECK_LOG_CODE="\
if [ ! -e \$LOG ];then
    echo \"[\$0] DOES NOT HAVE LOG FILE. CREATING...\"
    touch \$LOG
fi

"



if ! grep -q "COUNTER" "$INPUT_FILE"; then

  # WORKING_DIR="$(cd $(dirname $1)/../../ && pwd)"
  # base_name="${INPUT_FILE%.sh}"
  # backup_file="$WORKING_DIR/${base_name}_bku.sh"

  # cp -p "$INPUT_FILE" "$backup_file"
  original_permissions=$(stat -c %a "$INPUT_FILE")
  TEMP_FILE=$(mktemp)
  awk -v check_log_code="${CHECK_LOG_CODE}" '
  {
      if ($0 ~ /\$\{TESTCONF\}/) {
          print check_log_code
          print "TOTAL_TESTS=$(echo $TESTCONF | tr \" \" \"\\n\" | wc -l)\necho \"TOTAL_TESTS: $TOTAL_TESTS\"\nCOUNTER=0"
          print $0
          print "\tCOUNTER=$((COUNTER+1))\n\tdraw_progress_bar $COUNTER $TOTAL_TESTS"
      } else {
          print $0
      }
  }' "${INPUT_FILE}" > "${TEMP_FILE}"
  mv "${TEMP_FILE}" "${INPUT_FILE}"
  chmod "$original_permissions" "$INPUT_FILE"
else
  echo "The file $1 already contain progress bar!"
fi
