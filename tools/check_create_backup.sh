#!/bin/bash
#
# Description: This script will
#   * Check for org directories in conf/, inst/, run/ directories,if not create
#   * Check for backup files in conf/org/, inst/org/, run/$Benchmark_ID/org/ directories, if not create
# 
# FIRST EVER RUN: it will check for Jens' code version and replica those files into org directories with _old prefix.
# OTHER RUN     : it will check for the exist of 2 versions of benchmark in each directories
#                 _old in org refer to Jens' code, normal refer to my maybe modified version
#
# Usage:
# 
#   check_crewate_backup.sh $Benchmark_ID_Arr
#
ROOTDIR=$(cd "$(dirname $BASH_SOURCE[0])/.." && pwd)
# echo $ROOTDIR

# would be like amg, babelstream, comd
Benchmark_ID_Arr=("$@")
# echo ${Benchmark_ID_Arr[@]}


# Directories to check, dir run handle different inside
CHECK_DIRS=("conf" "inst")
#

# Check and create 'org' directories for conf and inst
for dir in "${CHECK_DIRS[@]}"; do
  org_dir="$ROOTDIR/$dir/org"
  if [ ! -d "$org_dir" ]; then
    echo "[LOG] [$dir/org] BACK UP DIR INSIDE NOT FOUND! CREATING.."
    mkdir -p "$org_dir"
  fi

  for bm_id in ${Benchmark_ID_Arr[@]};do
    bm_script="${bm_id}.sh"
    if [ ! -e "$ROOTDIR/$dir/${bm_script}" ]; then
    # if ! rg "${bm_script}" "$ROOTDIR/$dir" > /dev/null 2>&1; then
      echo "[ERROR] [$dir/$bm_id] NORMAL BENCHMARK NOT FOUND!  MANUALLY CHECK AGAIN!"
      exit 1
    else
      if [ ! -f "${org_dir}/${bm_id}_old.sh" ]; then
      # if ! find "${org_dir}" -name "${bm_id}_old.sh"; then
      # if ! rg "${bm_id}_old.sh" "${org_dir}";then
        echo "[BACKUP] [$dir/org/${bm_id}_old] NOT FOUND! REPLICA..."
        # Exist Jens' files in $dir, replica with _old prefix
        cp "$ROOTDIR/$dir/${bm_script}" "$org_dir/${bm_id}_old.sh"
      fi
    fi 
  done
done

# Check part for backup files
for bm_id in ${Benchmark_ID_Arr[@]};do
  # Special handling for the run directory
  bm_run_org_dir="$ROOTDIR/run/$bm_id/org"
  if [ ! -d "$bm_run_org_dir" ]; then
    echo "[LOG] [run/$bm_id/org] BACK UP DIR INSIDE NOT FOUND! CREATING.."
    mkdir -p "$bm_run_org_dir"
  fi
  
  bm_id_bku="${bm_id}_old.sh"
  # echo $bm_id_bku
  run_dir="$ROOTDIR/run/$bm_id"
  # echo $run_dir
  if [[ ! -e "${run_dir}/test.sh" || \
        ! -e "${run_dir}/best.sh" ]];then
    echo "[ERROR] [run/$bm_id] NORMAL TEST/BEST FILE  NOT FOUND! MANUALLY CHECK AGAIN!"
    exit 1
  else
    # echo $bm_run_org_dir
    if [[ ! -e "${bm_run_org_dir}/test.sh" || \
          ! -e "${bm_run_org_dir}/best.sh" ]];then
      # Copy the whole run directory to bku
      echo "[BACKUP] [run/${bm_id}/org] REPLICA WHOLE DIRECTORY run/$bm_id..."
      cp "${run_dir}/"*.sh "${bm_run_org_dir}" 
    fi
  fi
done

echo "[LOG] BACKUP CHECK AND CREATE FINISH"
