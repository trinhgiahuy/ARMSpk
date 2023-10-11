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
TESTCONFVAL="1|64|1|1|1 1|128|1|1|1
            4|16|2|2|1 4|32|2|2|1
            8|1|2|2|2 8|2|2|2|2 8|4|2|2|2 8|8|2|2|2 8|16|2|2|2
            16|4|4|2|2 16|8|4|2|2
            32|2|4|4|2 32|4|4|4|2
            64|1|4|4|4 64|2|4|4|4
            128|1|8|4|4"


# DO NOT USE THIS TESTCONFVAL for HPL, XHPL cannot run with 1 MPI process
# TESTCONFVAL="4|16|2|2|1 4|32|2|2|1
            # 8|1|2|2|2 8|2|2|2|2 8|4|2|2|2 8|8|2|2|2 8|16|2|2|2
            # 16|4|4|2|2 16|8|4|2|2
            # 32|2|4|4|2 32|4|4|4|2
            # 64|1|4|4|4 64|2|4|4|4
            # 128|1|8|4|4"

BM_ARR_NO_DOMAIN=(
    'babelstream'
    'candle'
    'dlproxy'
    'fs2020'
    'hpcg'
    'laghos'
    'macsio'
    'minife'
    'minitri'
    'modylas'
    'mvmc'
    'nekbone'
    'ngsa'
    # 'nicam'
    'ntchem'
    'sw4lite'
    'swfft'
    'xsbench'
)

INPUT_SCRIPT="$1"
BASENAME_SCRIPT=$(basename $INPUT_SCRIPT)
BENCH_ID="${BASENAME_SCRIPT%.sh}"


# If benchmark exist in an array of no domain
if [[ "${BM_ARR_NO_DOMAIN[@]}" =~ "$BENCH_ID" ]];then
    echo "[$0] BENCHMARK $BENCH_ID DO NOT REQUIRE DOMAIN IN TESTCONF"
    # export TESTCONFVAL=$(echo "${TESTCONFVAL}" | sed 's/\([0-9]*|[0-9]*\)|[0-9]*|[0-9]*|[0-9]*//g')
    export TESTCONFVAL=$(echo "${TESTCONFVAL}" | sed 's/\([0-9]*|[0-9]*\)|[0-9]*|[0-9]*|[0-9]*/\1/g')
    # echo "NEW $TESTCONFVAL"
elif [[ "$BENCH_ID" =~ hpl ]];then
    echo "BENCHMARK HPL DIFFERENT CONFIG"
    # GET 2 DOMAINS ONLY
    export TESTCONFVAL=$(echo "${TESTCONFVAL}" | sed 's/\([0-9]*|[0-9]*|[0-9]*|[0-9]*\)|[0-9]*/\1/g')
elif [[ "$BENCH_ID" =~ nicam ]];then
    echo "BENCHMARK NICAM DIFFERENT CONFIG"
    export TESTCONFVAL=$(echo "${TESTCONFVAL}" | tr ' ' '\n' | sed 's/^[0-9]*\([|][0-9]*\).*$/10\1/g' | sort | uniq | paste -sd' ' -)
else
    echo "[$0] BENCHMARK $BENCH_ID REQUIRE DOMAIN IN TESTCONF"
fi

#TODO: Add code hibench,polybench is exception with only test 1|1
if [[ "$BENCH_ID" =~ hpl ]];then
    CODE_DEFINE="\
elif [ -n \"\${ARMHOST}\" ]; then
    export TESTCONF=\"$TESTCONFVAL\"
    export BESTCONF=\"\"
    export HPLNB=\"192\""
elif [[ "$BENCH_ID" =~ minife ]];then
    echo "GOT MINIFE"
    CODE_DEFINE="\
elif [ -n \"\${ARMHOST}\" ]; then
    export TESTCONF=\"$TESTCONFVAL\"
    export BESTCONF=\"\"
    export BINARY=\"./mkl/src/miniFE.x\""
elif [[ "$BENCH_ID" =~ hibench || "$BENCH_ID" =~ polybench ]];then
    #TODO: ADD HERE
    CODE_DEFINE=""
else
    CODE_DEFINE="\
elif [ -n \"\${ARMHOST}\" ]; then
    export TESTCONF=\"$TESTCONFVAL\"
    export BESTCONF=\"\""
fi

# echo $CODE_DEFINE

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
