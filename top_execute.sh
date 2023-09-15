#!/bin/bash

# =============================================================================
# Script Name: top_execute.sh
#
# Description: This script will rerun all checked benchmarks by following steps
#
# *  Run with GNU first then LLVM later
# *  For each compilers, binary would be maked in directories GNU/LLVM  under benchmarks directory
# *
# *  Before run, it would change the log file with date and move old log file and move it to ~/archive_log/
# whilke
# *
#
#
# Usage: Change the benchmarks array want to run
# Author:       Huy Trinh
# Emai:         huy.trinh@a.riken.jp
# Date:         Sept 7, 2023
# =============================================================================


#TODO: It check benchmark type to append_config properly (BM with domain decomposition)
#TODO: After append_
#
#
#Multiple binaries:
#'babelstream'
#'fs2020'
#'hibench'
#'polybench'
#'spec_cpu'
#'spec_omp'
#'sw4lite'
#'swfft'
#

benchmark_id_arr=(
    # 'amg'
    # 'candle' # REQUIRE CONDA ENV
    # 'comd'
    # 'dlproxy' # INSTALL ERROR
    'ffb'
    'ffvc'
    # 'hpcg'
    # 'hpl' Current error due to make use patches, with arch=Linux64 not suitable
    #  Study patches and modify it properly

    # 'laghos'
    'macsio'
    'miniamr'
    'minife'
    'minitri'
    'modylas'
    'mvmc'
    'nekbone'
    'ngsa'
    'nicam'
    'ntchem'
    'qcd'
    'xsbench'
)
compiler_opts=(
  'gnu'
  # 'llvm-arm'
)

# echo "${benchmark_id_arr[@]}"
# echo "${compiler_opts[@]}"
ROOTDIR=$(dirname "$0")
TERM_LOG_DIR=$ROOTDIR/term-log
# echo $ROOTDIR
# echo $TERM_LOG_DIR

if [ ! -d $TERM_LOG_DIR ];then
  echo "Terminal log directory not found. Creating.."
  mkdir $TERM_LOG_DIR
fi


# Call copy from restoration for clean files
tools/copy_restore_repo.sh ${benchmark_id_arr[@]}

# Call check backup for all array first
tools/check_create_backup.sh ${benchmark_id_arr[@]}

total_iterations=$(( ${#benchmark_id_arr[@]} * ${#compiler_opts[@]} ))
# echo total $total_iterations
current_iteration=0

echo ""

source $ROOTDIR/conf/env.cfg

for compiler_opt in "${compiler_opts[@]}";do

  for bm_id in "${benchmark_id_arr[@]}";do

    # echo "RUN"
    sh_fn="${bm_id}.sh"
    # echo $sh_fn
    # python3 tools/subtitute_llvmarm.py --inst_file inst/$sh_fn
    tools/insert_copy_bin.sh inst/$sh_fn
    tools/append_make_err_handling.sh inst/$sh_fn

    # Call install benchmark
    inst/$sh_fn $compiler_opt
    if [ $? -ne 0 ]; then
        echo "[ERROR] FILE inst/$sh_fn RETURN ERROR. EXITING.."
        exit 1
    fi

    tools/append_config.sh conf/$sh_fn
    tools/move_binaries.sh $bm_id $compiler_opt
    tools/insert_draw.sh run/$bm_id/test.sh

    term_log_file=$TERM_LOG_DIR/"${bm_id}_term.log"
    if [ ! -e "$term_log_file" ];then
      echo "Creating terminal log file for benchmark $bm_id"
      touch $term_log_file
    fi

    # Update and print progress bar
    echo "RUNNING ON BENCHMARK $bm_id"
    current_iteration=$(( $current_iteration + 1 ))

    # echo "$total_iterations $current_iteration"
    draw_red_progress_bar $current_iteratio $total_iterations
    # Run benchmark
    CURR_DATETIME=$(date '+%Y-%m-%d %H:%M:%S')
    echo "" >> $term_log_file
    echo "==[${CURR_DATETIME}]=============================================================" >> $term_log_file
    run/$bm_id/test.sh $compiler_opt >> $term_log_file

  done
  current_iteration=$(( $current_iteration + 1 ))
done

echo ""
