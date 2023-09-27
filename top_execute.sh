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
    # 'ffb'
    # 'ffvc'
    # 'hpcg'
    # 'hpl'
    # 'laghos'
    # 'macsio'
    # 'miniamr'
    # 'minife'
    # 'minitri'
    # 'modylas'
    'mvmc'
    # 'nekbone'
    # 'ngsa'
    # 'nicam'
    # 'ntchem'
    # 'qcd'
    # 'xsbench'
)
compiler_opts=(
  'gnu'
  # 'llvm-arm'
)

# Deactivate if previous run MVMC
# replica_instsall_mvmc itself will activate the py27_env
conda deactivate


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
    if [[ $bm_id =~ minitri ]];then
        echo "Got minitri"
        tools/insert_copy_bin_minitri.sh inst/$sh_fn
    else
        tools/insert_copy_bin.sh inst/$sh_fn
    fi
    tools/append_make_err_handling.sh inst/$sh_fn

    # Call install benchmark
    # If hpl script, call replica_install_openblas.sh
    if [[ $bm_id =~ hpl ]];then
        echo "Got hpl"
        # rm -rf HPL/
        tools/replica_install_hpl_openblas.sh inst/$sh_fn $compiler_opt
    elif [[ $bm_id =~ mvmc ]];then
        echo "Got MVMC"
        tools/replica_install_mvmc.sh inst/$sh_fn $compiler_opt
    else
        # echo "SKIP"
        inst/$sh_fn $compiler_opt
    fi

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

    # CHECK IF BINARY EXIST! IF NOT EXIST
    BIN_DIR="$ROOTDIR/bin/$bm_id/$compiler_opt"
    if ! find "$BIN_DIR" -maxdepth 1 -type f -executable -print | grep -q .;then
        echo "NO EXECUTABLE"
        exit 1
    else
        echo "FOUND EXECUTABLE"
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
