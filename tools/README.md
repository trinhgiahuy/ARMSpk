## WORK FLOW

## THIS STEP IS OPTIONAL/USE WHEN OBTAINING CLEAN ORIGINAL CODE
### copy_restore_repo.sh
#### [inst/conf/run] Copy the original version of inst/conf/run file to our current directory

`tools/copy_restore_repo.sh 'amg' ['comd' ['laghos']]`


### check_create_backup.sh
#### [inst/conf/run] Check backup benchmark files inside inst/conf/run.

`tools/check_create_backup.sh 'amg' ['comd' ['laghos']]`


### subtitute_llvmarm.py
#### [inst] Append the LLVM-ARM subtitution for inst benchmark file.

` python3 tools/subtitute_llvmarm.py --inst_file inst/amg.sh`


### insert_copy_bin.sh
#### [inst] Append the move binary code (from $BM dir to dir `bin/$compilers/` when call `inst/$bm.sh gnu|llvm-arm`) for inst benchmark file.

` tools/insert_copy_bin.sh inst/amg.sh`

### append_make_err_handling.sh
#### [inst] Append the error handling exit in last make liune for inst benchmark file. If `make` fail, it will exit for manually chekck

`tools/append_make_err_handling.sh inst/amg.sh`

[CALL INST/$BM.sh file]


### append_config.sh
#### [conf] Append the TESTCONF and empty RUNCONF for conf benchmark files.

`tools/append_config.sh conf/amg.sh`


### move_binaries.sh
#### [conf/run] Check `$1` in `source conf/${BENCH_ID}.sh` in run benchmark file, check `export BINARY=""` updated to our binary directory.

`tools/move_binaries.sh amg gnu`


### insert_draw.sh
#### [run] Append drawing progress bar code for run benchmark files

`tools/insert_draw.sh run/`

[CALL RUN/$BM/test.sh file]



## TODO & DEVELOPMENT TRENDS

+ Deal with multiple binaries output
+ Inside each tools, check for $1 param inside specific directory


+ [SAFE FEATURE], check for accidentally  multiple calls

**DONE**:`check_create_backup`, `subtitute_llvmarm.py`, `insert_copy_bin.sh`, `append_config`, `insert_draw.sh`

## NOTE:

- [HIDDNE BUG] `python3 tools/subtitute_llvmarm.py --inst_file inst/amg.sh` sometime yields `[ERROR] Pattern for 'llvm12' not found!`. [UNDER INVESTIGATION],
    related to swap? 2 files modified at same time?

## CURRENT ERROR

