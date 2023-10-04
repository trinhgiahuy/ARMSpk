#!/bin/bash
# =============================================================================
# Script Name: insert_draw.sh
# Description: This script will add progress bar code for Benchmark running test
# Additional add:
#
#       The script will add checl log file code to run file. In some cases,
#       run file cannot automatically create new log file
#
#       For NICAM benchmark, it will conditional check for existing directory of
#       data symbolic link
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
    last_test_conf_line=$(grep -nE 'TESTCONF' "$INPUT_FILE" | tail -n1 | cut -d: -f1)
    echo last_test_conf_line $last_test_conf_line

    original_permissions=$(stat -c %a "$INPUT_FILE")
    TEMP_FILE=$(mktemp)
    awk -v check_log_code="${CHECK_LOG_CODE}" -v line_num="${last_test_conf_line}" '
    {
        if ( NR == line_num){
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


LINK_SYSDEP_CODE="\
    cd ../../../sysdep
    [ ! -e Makedef.ARM64 ] && ln -s ../../sysdep/Makedef.ARM64 .
    cd -"

# Add link sysdep to files test of NICAM
bm_id=$(echo $INPUT_FILE | cut -d/ -f2 )
echo bm_id $bm_id

if [[ $bm_id =~ nicam ]];then
    echo "Got nicam"
    make_jobshell_line=$(grep -nE 'jobshell' "$INPUT_FILE" | cut -d: -f1)
    echo "jobshell_line $make_jobshell_line"

    original_permissions=$(stat -c %a "$INPUT_FILE")
    TEMP_FILE=$(mktemp)
    if ! grep -q "sysdep" "$INPUT_FILE";then
        echo "Adding link sysdep"
        awk -v link_sysdep_code="$LINK_SYSDEP_CODE" -v line_num="${make_jobshell_line}" '
        {
            if ( NR == line_num){
                print $0
            }else{
                print $0
            }
        }' "${INPUT_FILE}" > "${TEMP_FILE}"
        mv "${TEMP_FILE}" "${INPUT_FILE}"
        chmod "$original_permissions" "${INPUT_FILE}"
    else
        echo "Already has link code"
    fi
    # SHOULD BE
    # Makedef.ARM64 nhm_driver rl00-prc10.info vgrid40_24000-600m.dat boundary_GL05RL00.pe000000 boundary_GL05RL00.pe000001
    # boundary_GL05RL00.pe000002 boundary_GL05RL00.pe000003 boundary_GL05RL00.pe000004 boundary_GL05RL00.pe000005
    # boundary_GL05RL00.pe000006 boundary_GL05RL00.pe000007 boundary_GL05RL00.pe000008 boundary_GL05RL00.pe000009
    linking_files_name_arr=($(grep -n 'ln' run/nicam/test.sh | rev |  cut -d/ -f1| cut -d' ' -f2 | rev))
    echo "${linking_files_name_arr[@]}"

    # TEMP_FILE=$(mktemp)
    # cp -p $INPUT_FILE $TEMP_FILE
    for fn in "${linking_files_name_arr[@]}";do
        echo $fn
        # sed -i '/ln -s .* '"$fn"'$/s#^#[ ! -e '"$fn"' ] \&\& #' $TEMP_FILE
        sed -i '/ln -s .*'"$fn"'.*/s#^\t#\t[ ! -e '"./$fn"' ] \&\& #' $INPUT_FILE
    done

    # mv $TEMP_FILE $INPUT_FILE
fi



