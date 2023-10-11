# =============================================================================

# Script Name: subtitute_config.py
#
# Usage:
#
#   python3 subtitute_config.py --config_file $ROOTDIR/config/fs2020.sh
#
# Safe Features:
# Author:      Huy Trinh
# Emai:        huy.trinh@a.riken.jp
# Date:        Oct 11 , 2023
# =============================================================================

import os
import sys
import re
import argparse
import shutil


def is_set(arg_name):
    """
    @Desc: The function check if an argument is setted or not
    """
    if arg_name in sys.argv:
        return True
    return False


def build_arm_host_config(file_path):
    """
    @Desc: This function will return the ARMHOST block config code based on XEONHOST block config code
    """

    with open(file_path, 'r') as file:
        content = file.read()

    xeon_host_pattern = r'((el)?if\s+\[\s*-n\s*"\${XEONHOST}"\s*\];\s*then[\s\S]*?)(?=\n(el)?if)'
    match = re.search(xeon_host_pattern, content)
    #  print(f"match {match}")

    if match:
        arm_host_block = match.group(1).replace('XEONHOST','ARMHOST').replace('if','elif')
        arm_host_block = '\n'.join(line for line in arm_host_block.splitlines() if not line.strip().startswith('#'))
        arm_host_block = re.sub(r'(export BESTCONF=)"[^"]*"', r'\1""', arm_host_block)
        #  print(arm_host_block)

    else:
        print("[ERROR] NO XEON CONFIG CODE FOUND! CHECK AGAIN..")
        arm_host_block=None

    file.close()


    return arm_host_block

def append_before_last_fi(file_path, code_define):

    with open(file_path, 'r') as file:
        lines = file.readlines()
        content = ''.join(lines)

    arm_host_pattern = r'elif\s*\[\s*-n\s*"\${ARMHOST}"\s*.*'

    if re.search(arm_host_pattern, content):
        print("[LOG] CONFIG BLOCK FOR ARMHOST ALREADY EXIST")
        sys.exit(0)

    # Find the index of the last occurrence of 'fi\n'
    last_fi_index = None
    for idx, line in reversed(list(enumerate(lines))):
        if line == "fi\n":
            last_fi_index = idx
            break

    # Insert the code_define right before the last 'fi\n'
    if last_fi_index is not None:
        lines.insert(last_fi_index, code_define + '\n')

    # Write the modified lines back to the file
    with open(file_path, 'w') as file:
        file.writelines(lines)

    file.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Subtitute the elif part for the \"ARMHOST\" based on \"XEONHOST\" code part")
    parser.add_argument(
        "--conf_file",
        type=str,
        required=True,
        dest="conf_file",
        help="Benchmark Config Files: conf/$bm.sh",
    )

    args=parser.parse_args()

    if is_set(args.conf_file):
        FILE_PATH = args.conf_file
        #  print(f"FILE PATH {FILE_PATH}")

    arm_host_block=build_arm_host_config(FILE_PATH)
    if arm_host_block is not None:
        append_before_last_fi(FILE_PATH, arm_host_block)


