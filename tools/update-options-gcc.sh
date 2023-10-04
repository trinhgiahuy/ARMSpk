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
#

# Exclude directory first, then files, otherwise unexpected behavior
search_cmd() {
    rg -i "$1" --glob '!{*org/,*patches/,*src*}' --glob '!{*Linux*,*linux*,*intel*,*kei*,*kashiwa*,*src*}' "${@:2}"
}
declare -A ARM_INTERCHANGE=(
    ["mpiifort"]="mpifort"
    ["mpiicc"]="mpicc"

    # OTHER OPTION INTERCHANGE
    [" -fpp3"]=" -cpp"
    [" -ip"]=" -flto"
    [" -convert big_endian"]=" -fconvert=big-endian"
    [" -fp-model precise"]=" -ffloat-store"
    [" -heap-arrays"]=" -fstack-arrays"
    [" -fno-alias"]=" -fno-strict-aliasing"
    [" -xHost"]=" -march=native"
)

EXCHANGE_OPT=$2
COMPILER_INCOMPS=$(echo $EXCHANGE_OPT| cut -d'#' -f1)
echo $COMPILER_INCOMPS
FLAG_INCOMPS=$(echo $EXCHANGE_OPT| cut -d'#' -f2)
echo $FLAG_INCOMPS
# FLAG_INCOMPS=$(echo $FLAG_INCOMPS | awk -F'|' '{ for(i=1;i<=NF;i++) printf (i==NF) ? ((substr($i,1,1) == "-") ? "\\s"$i:$i) : ((substr($i,1,1) == "-") ? "\\s"$i"|":$i"|")}')
# echo $FLAG_INCOMPS


# SHOULD BE mpiifort|mpiicc|\s-fpp3|\s-ip
EXCHANGE_OPT_MERGED="$COMPILER_INCOMPS|$FLAG_INCOMPS"

echo $EXCHANGE_OPT_MERGED


echo ${inc_opt_arr[@]}

ROOTDIR=$(cd "$(dirname $0)/.." && pwd)

# echo $ROOTDIR

# BM name like "MVMC"
BM="$1"
BMDIR="$ROOTDIR/$BM"
declare -A script_name_dict=(
    ["MVMC"]="mvmc"
    ["NICAM"]="nicam"
)

INPUT_FILE="inst/${script_name_dict["$1"]}.sh"


# echo $INPUT_FILE

# subfiles=$(search_cmd "$inc_opt" $BMDIR -l)
# if [ -n "$subfiles" ]; then
    # for file in $subfiles;do
        # IFS='|' read -ra inc_opt_arr <<< "$EXCHANGE_OPT_MERGED"
        # for inc_opt in $inc_opt_arr;do
            # echo $inc_opt
            # if [[ -v ARM_INTERCHANGE["$inc_opt"] ]];then
                # echo "FOUND INTERCHANGE OPTION FOR $inc_opt"
                # chg_opt=${ARM_INTERCHANGE["$inc_opt"]}
                # echo "+++++++++++++++++++++++++++++ $file"
                # exit 1
                # sed -i "s/${inc_opt}/${chg_opt}/g" $file
            # fi
        # done
    # done
# fi

subfiles=$(search_cmd "$EXCHANGE_OPT_MERGED" $BMDIR -l)
# echo "$(type search_cmd) "$EXCHANGE_OPT_MERGED" $BMDIR -l"
echo "GET:${subfiles[@]}"
if [ -e "$subfiles" ]; then
    IFS='|' read -ra inc_opt_arr <<< "$EXCHANGE_OPT_MERGED"
    # inc_opt_arr=(${EXCHANGE_OPT_MERGED//|/ })

    for inc_opt in "${inc_opt_arr[@]}";do
        echo "1:  $inc_opt"

        if [[ ${ARM_INTERCHANGE["$inc_opt"]+isset} ]];then
        # if [[ -v ARM_INTERCHANGE["$inc_opt"] ]];then
            echo "FOUND INTERCHANGE OPTION FOR $inc_opt"
            chg_opt=${ARM_INTERCHANGE["$inc_opt"]}
            # echo "$(type search_cmd) $inc_opt $BMDIR -l"
            echo "CHANGE:  $chg_opt"
            subfiles=$(search_cmd "$inc_opt" $BMDIR -l)
            # echo "=== $subfiles"
            for file in "${subfiles[@]}";do
                echo "+++++++++++++++++++++++++++++ $file"
                sed -i "s/${inc_opt}/${chg_opt}/g" $file
                # sed -i "s/'${inc_opt}'/'${chg_opt}'/g" $file
            done
        fi
    done
fi


# Add exit 1 to before make line in Jens code
# last_make_line=$(grep -nE '\t*make\b' "$INPUT_FILE" | tail -n1 | cut -d: -f1)
# echo last_make_line $last_make_line
# TEMP_FILE=$(mktemp)
# org_permission=$(stat -c %a "${INPUT_FILE}")
# if ! rg 'Adding exit' $INPUT_FILE > /dev/null 2>&1; then
    # echo "DO NOT FIND EXIT.ADD"
    # echo last_make_line $last_make_line
    # awk -v n="$last_make_line" 'NR == n {print "\techo \"Adding exit1\"; exit 1"} 1' "$INPUT_FILE" > "$TEMP_FILE" && mv "$TEMP_FILE" "$INPUT_FILE"
    # #####awk -v n="$last_make_line" 'NR == n && !/exit 1/ {print "\techo \"Adding exit1\"; exit 1"} 1' "$INPUT_FILE" > "$TEMP_FILE" && mv "$TEMP_FILE" "$INPUT_FILE"
    # chmod $org_permission "$INPUT_FILE"
# else
    # echo "FIND EXIT"
# fi

# OPTS: -DHAVE_SSE2 IS NOT AVAILABLE ON SUPERCOMP07.
# PUT INCOMPATIBLE COMPILE OPTS TO GCC AS $3 ARG
INCOPT="$3"
# INCOPT=$(echo "$INCOPT" | sed 's/ /\ /g')
echo "HERE: $INCOPT"
# SHOULD BE SOMETHING LIKE "\s-ipo|\s-qopt-prefetch=3|\s-nofor-main|\s-vevc-report"
SEARCHOPT=$(echo $INCOPT | awk -F'|' '{ for(i=1;i<=NF;i++) printf (i==NF) ? ((substr($i,1,1) == "-") ? "\\s"$i:$i) : ((substr($i,1,1) == "-") ? "\\s"$i"|":$i"|")}')
echo $SEARCHOPT
# THIS WORKS BUT USE AWK
# SUBTITUTEOPT=$(echo $INCOPT | awk -F'|' '{for(i=1;i<=NF;i++) printf (i==NF) ? $i : "\\|"$i}')
SUBTITUTEOPT=$(echo $INCOPT | sed 's/|/\\|/g')
echo "SUBTITUTEOPT: $SUBTITUTEOPT"

if search_cmd $SEARCHOPT $BMDIR; then
# if rg "$SEARCHOPT" --glob '!{*org,*patches}/' $BMDIR; then
# if rg '\s-ipo|\s-qopt-prefetch=4|\s-nofor-main|\s-vevc-report' --glob '!{*org,*patches}/' $BMDIR; then
    echo "[LOG] EXIST OPTIONS NOT COMPATIBLE WITH GCC. MODYFYING.."
    subfiles=($(search_cmd $SEARCHOPT $BMDIR -l))
    for file in "${subfiles[@]}";do
        echo $file
        sed -i "s/${SUBTITUTEOPT}//g" $file
    done
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

