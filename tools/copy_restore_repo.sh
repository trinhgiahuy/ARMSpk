#!/bin/bash
#

BENCHID_INPUTS=($@)

echo "${BENCHID_INPUTS[@]}"

# NOT RECOMMEND TO CHANGE THIS UNLESS JENS's REPO CHANGES
RESTORE_DIR=$HOME/restore/a64fxCvC
MY_DIR=$HOME/riken/a64fxCvC

for benchid in "${BENCHID_INPUTS[@]}";do

    # COPY CONF FILE
    cp -p $RESTORE_DIR/conf/"${benchid}.sh" $MY_DIR/conf/"${benchid}.sh"

    # COPY INST FILE
    echo $RESTORE_DIR/inst/"${benchid}.sh"
    echo $MY_DIR/inst/"${benchid}.sh"
    cp -p $RESTORE_DIR/inst/"${benchid}.sh" $MY_DIR/inst/"${benchid}.sh"

    # COPY TEST FILE
    cp -p $RESTORE_DIR/run/${benchid}/test.sh $MY_DIR/run/${benchid}/test.sh

done

