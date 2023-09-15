#!/bin/bash
# =============================================================================
# Script Name: append_make_err_handling.sh
# Description: This script will append the error handling for the last make in inst files
#
# Usage:
#
#       tools/append_make_err_handling.sh inst/$bm.sh
#
# Author:       Huy Trinh
# Emai:         huy.trinh@a.riken.jp
# Date:         Sept 14, 2023
# =============================================================================


if [ -z $1 ];then
    echo "USAGE: $0 inst/\$bm.sh"
    exit 1
fi

INPUT_FILE=$1

last_make_line=$(grep -nE '\t*make\b' "$INPUT_FILE" | tail -n1 | cut -d: -f1)

if rg 'MAKE FAIL' $INPUT_FILE > /dev/null 2>&1; then
    echo "[$0] MAKE HAS ERR HANDLING ALREADY!"
else
    echo "last_make_line: ${last_make_line}"
    if [[ -n $last_make_line && "$last_make_line" =~ ^[0-9]+$ ]]; then
        echo "[$0] ADD MAKE ERROR HANDLING TO FILE $1"
        sed -i "${last_make_line}s/$/ || { echo \"MAKE FAIL! EXITING..\"; exit 1; }/" "${INPUT_FILE}"
    else
        echo "[ERROR] [$0] NO MAKE FOUND IN LINE... CHECK $INPUT_FILE MANUALLY"
        exit 1
    fi
fi
