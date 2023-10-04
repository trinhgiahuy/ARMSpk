#!/bin/bash


#TODO: Advaced FLAG_OPTS shoud include multiple flags seperated by |
#For e.g: -fallow-argument-mismatch|-openmp
FLAG_OPTS=$1
VAR_SET=$2
FILE_NAME=$3


if [[ -z "$1" || -z "$2" || -z "$3" ]];then\
    echo "EXAMPLE USAGE: tools/append_flag_opt.sh \"-fallow-argument-mismatch|-openmp|...\" \"FFLAG_FAST\" NICAM/sysdep/Makedef.ARM64"
    exit 1
fi

if [[ ! -f $FILE_NAME ]];then
    echo  "[ERROR] $FILE_NAME NOT FOUND"
    exit 1
fi

echo "${FLAG_OPTS}"
echo "${VAR_SET}"
awk -v flag_opt="${FLAG_OPTS}" \
-v var_set="${VAR_SET}" \
    '
    $0 ~ var_set " *=.*" {
    print
    getline
    while ($0 ~ /\\$/){
        print
        getline
    }
    print $0, flag_opt
    next
}
{ print }' $FILE_NAME > "${FILE_NAME}.tmp" && mv "${FILE_NAME}.tmp" "${FILE_NAME}"
