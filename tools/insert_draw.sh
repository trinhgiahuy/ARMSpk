#!/bin/bash

# Check for required input
if [ -z "$1" ]; then
    echo "Usage: $0 <input-file>"
    exit 1
fi

# Input file
INPUT_FILE=$1

# Temporary file for the modified content
TEMP_FILE=$(mktemp)

# Use awk to process the input file
awk '
{
    if ($0 ~ /\$\{TESTCONF\}/) {
        print "TOTAL_TESTS=$(echo $TESTCONF | tr \" \" \"\\n\" | wc -l)\necho \"TOTAL_TESTS: $TOTAL_TESTS\"\nCOUNTER=0"
        print $0
        print "\tCOUNTER=$((COUNTER+1))\n\tdraw_progress_bar $COUNTER $TOTAL_TESTS"
    } else {
        print $0
    }
}' "${INPUT_FILE}" > "${TEMP_FILE}"

# Replace the original file with the modified one
mv "${TEMP_FILE}" "${INPUT_FILE}"

