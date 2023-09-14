## WORK FLOW


### [inst/conf/run] Check backup benchmark files inside inst/conf/run.

`tools/check_create_backup.sh 'amg' ['comd' ['laghos']]`

### [inst] Append the LLVM-ARM subtitution for inst benchmark file.

` python3 tools/subtitute_llvmarm.py --inst_file inst/amg.sh`

### [inst] Append the move binary code (from $BM dir to dir `bin/$compilers/` when call `inst/$bm.sh gnu|llvm-arm`) forr inst benchmark file.

` tools/insert_copy_bin.py inst/amg.sh`

[CALL INST/$BM.sh file]

### [conf] Append the TESTCONF and empty RUNCONF for conf benchmark files.

`tools/append_config.sh conf/amg.sh`

### [conf/run] Check `$1` in `source conf/${BENCH_ID}.sh` in run benchmark file, check `export BINARY=""` updated to our binary directory.

`tools/move_binaries.sh amg gnu`

### [run] Append drawing progress bar code for run benchmark files

`tools/insert_draw.sh run/`

[CALL RUN/$BM/test.sh file]

## DESCRIPTION

## TODO & DEVELOPMENT TRENDS

+ Inside each tools, check for $1 param inside specific directory

+ `subtitute_llvmarm` check backup, change to check `org` instead of `_bku.sh` file

+ [SAFE FEATURE], check for accidentally  multiple calls

**DONE**:`check_create_backup`, `subtitute_llvmarm.py`, `insert_copy_bin.sh`, `append_config`, `insert_draw.sh`

## NOTE:

- [HIDDNE BUG] `python3 tools/subtitute_llvmarm.py --inst_file inst/amg.sh` sometime yields `[ERROR] Pattern for 'llvm12' not found!`. [UNDER INVESTIGATION],
    related to swap? 2 files modified at same time?

## CURRENT ERROR

