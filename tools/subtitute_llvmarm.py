# =============================================================================

# Script Name: subtitute_llvmarm.py
# Description: This script will automatically subtitute the elif part for [[ "$1" = *"llvm-arm"* ]]
# to ALMOST all benchmarks in inst/$bm.sh based on the elif part of [[ "$1" = *"llvm12"* ]] options
# The options can be replaced manually by adding in the REPLACEMENTS dictionary.
# These options is added in order of from very specific (long one) to very general (short).
# For example:
# * '$(which mpifcc))/../lib64': '$(which clang))/../lib' BEFORE 'mpifcc': 'clang'
#
# Usage:
#
#   python3 subtitute_llvmarm.py --inst_file $ROOTDIR/inst/amg.sh
#
# Safe Features:
#   * The script checks for if there is exist llvm-arm code block, if yes then skip & exit
#   * The script checks for if exist the backup files ( format $bm_bku.sh). If not then create before make changes
# Author:      Huy Trinh
# Emai:        huy.trinh@a.riken.jp
# Date:        Sept 8, 2023
# =============================================================================

import os
import sys
import re
import argparse
import shutil


# Global replacements dictionary for the 'llvm-arm' block
REPLACEMENTS = {

#===AMG
    '$(which mpifcc))/../lib64': '$(which clang))/../lib',
    ' -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)': '',


#===General
    '-mcpu=a64fx': '-mcpu=native',
    'mtune=a64fx': 'mtune=native',
    'mpifcc': 'clang',
}


def is_set(arg_name):
    if arg_name in sys.argv:
        return True
    return False

def ensure_backup(file_path):
    directory, filename = os.path.split(file_path)
    basename, ext = os.path.splitext(filename)
    backup_file = os.path.join(directory, f"{basename}_bku{ext}")

    #  print(f"backup_file is {backup_file}")
    if not os.path.exists(backup_file):
        print(f"[LOG] Creating backup file {backup_file}...")
        shutil.copy2(file_path,backup_file)
    else:
        print(f"[LOG] Backup file {backup_file} already exists!")
    return backup_file

def subtitute_llvmarm(file_path):
    with open(file_path, 'r') as file:
        content = file.read()


    # Ensure backup file before make changes
    #  ensure_backup(file_path)

    # Check if the "llvm-arm" block already exist
    llvm_arm_patterm = r'elif\s*\[\[\s*"\$1"\s*=\s*.*"llvm-arm".*\]\]'
    if re.search(llvm_arm_patterm,content):
        print("[LOG] Code block for 'llvm-arm' was already subtituted")
        sys.exit(0)

    # Find the specific elif block for 'llvm12'
    # pattern = r'(\s*elif \[\[ "\$1" = \*"llvm12"\* \]\]; then[\s\S]*?\n\s*fi)(?=\s*elif|\s*fi|$)'
    # pattern = r'(\s*elif \[\[ "\$1" = \*"llvm12"\* \]\]; then[\s\S]*?\n\tfi)'
    pattern = r'(\s*elif \[\[ "\$1" = \*"llvm12"\* \]\]; then[\s\S]*?)(?=\n\tfi)'
    match = re.search(pattern, content, re.DOTALL)

    if match:
        llvmarm_block = match.group(1).replace('llvm12', 'llvm-arm')
        replacement_made = False
        for old_str, new_str in REPLACEMENTS.items():
            if old_str in llvmarm_block:
                llvmarm_block = llvmarm_block.replace(old_str, new_str)
                replacement_made=True
        if replacement_made:
            # Insert the modified llvmarm block after the original 'llvm12' block
            new_content = content.replace(match.group(0), match.group(0).rstrip() + '\n' + llvmarm_block)
            print("[LOG] Finish subtituting the llvm-arm code part")

        with open(file_path, 'w') as file:
            file.write(new_content)
    else:
        print("[ERROR] Pattern for 'llvm12' not found!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Subtitute the elif part for the \"llvm-arm\" based on \"llvm12\" code part")
    parser.add_argument(
        "--inst_file",
        type=str,
        required=True,
        dest="inst_file",
        help="Benchmark Installation Files: inst/$bm.sh",
    )

    args=parser.parse_args()
    if is_set(args.inst_file):
        FILE_PATH = args.inst_file
        #  print(f"FILE PATH {FILE_PATH}")

    subtitute_llvmarm(FILE_PATH)

