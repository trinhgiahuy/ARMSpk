#!/bin/bash
#
# =============================================================================
# #
# Script Name: update-options-gcc.sh
# Description: [TEMPORARY SOLUTION] This script will check for incompatible options (which raises conflict
#              during compilation for ARMv8 (gcc) and remove it
#              TODO: Look for interchanging options while stil keep same performance
#
# Information: It will ignore directories org/, patches/, and files name including
#              intel, kei, kashiwa, [ADDMORE]
#
# Author:      Huy Trinh
# Emai:        huy.trinh@a.riken.jp
# Date:        Sept 26, 2023
# Usage:
#
#              tools/update-options-gcc.sh "MVMC" "-xHost|-nofor-main|-vec-report"
#TODO:
# =============================================================================
#
ROOTDIR=$(cd "$(dirname $0)/.." && pwd)
# echo $ROOTDIR

# BM name like "MVMC"
BM="$1"
BMDIR="$ROOTDIR/$BM"
declare -A script_name_dict=(
    ["MVMC"]="mvmc"
)

INPUT_FILE="inst/${script_name_dict["$1"]}.sh"

echo $INPUT_FILE

# Add exit 1 to before make line in Jens code
last_make_line=$(grep -nE '\t*make\b' "$INPUT_FILE" | tail -n1 | cut -d: -f1)
echo last_make_line $last_make_line
TEMP_FILE=$(mktemp)
org_permission=$(stat -c %a "${INPUT_FILE}")
if ! rg 'Adding exit' $INPUT_FILE > /dev/null 2>&1; then
    echo "DO NOT FIND EXIT.ADD"
    echo last_make_line $last_make_line
    awk -v n="$last_make_line" 'NR == n {print "\techo \"Adding exit1\"; exit 1"} 1' "$INPUT_FILE" > "$TEMP_FILE" && mv "$TEMP_FILE" "$INPUT_FILE"
    # awk -v n="$last_make_line" 'NR == n && !/exit 1/ {print "\techo \"Adding exit1\"; exit 1"} 1' "$INPUT_FILE" > "$TEMP_FILE" && mv "$TEMP_FILE" "$INPUT_FILE"
    chmod $org_permission "$INPUT_FILE"
else
    echo "FIND EXIT"
fi

# OPTS: -DHAVE_SSE2 IS NOT AVAILABLE ON SUPERCOMP07.
# PUT INCOMPATIBLE COMPILE OPTS TO GCC AS $2 ARG
INCOPT="$2"

# SHOULD BE SOMETHING LIKE "\s-ipo|\s-qopt-prefetch=3|\s-nofor-main|\s-vevc-report"
SEARCHOPT=$(echo $INCOPT | awk -F'|' '{ for(i=1;i<=NF;i++) printf (i==NF) ? ((substr($i,1,1) == "-") ? "\\s"$i:$i) : ((substr($i,1,1) == "-") ? "\\s"$i"|":$i"|")}')
echo $SEARCHOPT
# THIS WORKS BUT USE AWK
# SUBTITUTEOPT=$(echo $INCOPT | awk -F'|' '{for(i=1;i<=NF;i++) printf (i==NF) ? $i : "\\|"$i}')
SUBTITUTEOPT=$(echo $INCOPT | sed 's/|/\\|/g')
echo "$SUBTITUTEOPT"

# Exclude directory first, then files, otherwise unexpected behavior
search_cmd() {
    rg -i "$SEARCHOPT" --glob '!{*org/,*patches/,}' --glob '!{*intel*,*kei*,*kashiwa*,}' "$@"
}

if search_cmd $BMDIR; then
# if rg "$SEARCHOPT" --glob '!{*org,*patches}/' $BMDIR; then
# if rg '\s-ipo|\s-qopt-prefetch=4|\s-nofor-main|\s-vevc-report' --glob '!{*org,*patches}/' $BMDIR; then
    echo "[LOG] EXIST OPTIONS NOT COMPATIBLE WITH GCC. MODYFYING.."
    subfiles=($(search_cmd -l $BMDIR))
    for file in "${subfiles[@]}";do
        echo $file
        sed -i "s/${SUBTITUTEOPT}//g" $file
    done

    exit 0
else
    echo "NO FIND INCOMPATIBLE FOUNDS!"
fi

# Get the path to scalapack using spack
SCALAPACK_PATH=$(spack location -i scalapack)

# Check if SCALAPACK_PATH exists
if [ ! -d "$SCALAPACK_PATH" ]; then
    echo "Error: Scalapack path not found!"
    exit 1
fi

# Check if LD_LIBRARY_PATH already contains SCALAPACK_PATH
if [[ ":$LD_LIBRARY_PATH:" != *":$SCALAPACK_PATH:"* ]]; then
    export LD_LIBRARY_PATH=$SCALAPACK_PATH:$LD_LIBRARY_PATH
    echo "Added Scalapack path to LD_LIBRARY_PATH."
else
    echo "Scalapack path already in LD_LIBRARY_PATH."
fi

