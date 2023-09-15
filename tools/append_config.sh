#!/bin/bash
# =============================================================================
#
# Script Name: append_confg.sh
# Description: This script will add CONF and BESTCONF setting to all conf/$bm.sh
#							 files. Look for last `fi` line and append before it.
# Author:      Huy Trinh
# Emai:        huy.trinh@a.riken.jp
# Date:        Sept 7, 2023
# Usage:
#
#       tools/append_config.sh conf/$bm.sh
#TODO: Seperate $INPUT FOR CODE_DEFINE, some $bm does not have domain decomposition
# =============================================================================
export CODE_DEFINE="\
elif [ -n \"\${ARMHOST}\" ]; then
	export TESTCONF=\"1|64|1|1|1 1|128|1|1|1
                  4|16|2|2|1 4|32|2|2|1
                  8|1|2|2|2 8|2|2|2|2 8|4|2|2|2 8|8|2|2|2 8|16|2|2|2
                  16|4|4|2|2 16|8|4|2|2
                  32|2|4|4|2 32|4|4|4|2
                  64|1|4|4|4 64|2|4|4|4
                  128|1|8|4|4\"
    export BESTCONF=\"\""

INPUT_SCRIPT="$1"
org_permissions=$(stat -c %a "${INPUT_SCRIPT}")
TEMP_FILE=$(mktemp)

# awk requires $ENVIRON["CODE_DEFINE"] to access the global variable defined,
# that global var need to be export
# However, we can use `awk -v code="$CODE_DEFINE" ` and `print code`: A more common to do this
if ! grep -q "ARMHOST" $INPUT_SCRIPT; then
    echo "[$0] FOR $(basename "$INPUT_SCRIPT"), TESTCONF/BESTCONF NOT FOUND FOR ARMHOST. ADDING.."
    awk -v code_define="$CODE_DEFINE" '
    BEGIN {last_fi_line=0}
    {
        if ($0 - /^fi$/) {
            last_fi_line=NR;
        }
        lines[NR]=$0;
    }
    END {
        for (i=1; i<=NR; i++) {
            if (i == last_fi_line){
                printf "%s\n",code_define;
            }
            print lines[i];
        }
    } ' "$INPUT_SCRIPT" > "${TEMP_FILE}"

    mv "${TEMP_FILE}" "${INPUT_SCRIPT}"
    chmod "${org_permissions}" "$INPUT_SCRIPT"
else
    echo "[$0] TESTCONF/BESTCONF FOUND IN $(basename "$INPUT_SCRIPT")"
fi
