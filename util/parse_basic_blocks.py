#!/usr/bin/env python3

# usage:
# <script> -j ./dcfg-out.dcfg.json.bz2 -b ./dcfg-out.bb.txt.bz2
#  ( create files: ./dcfg-out.dcfg.json.bz2 and ./dcfg-out.bb.txt.bz2 via
#      `sde64 -dcfg 1 -dcfg:write_bb 1 ...` )

from os import path, getpid, remove
from sys import exit
from bz2 import open as bz2open
from json import dumps as jsondumps     # TODO: take this out
from re import compile, sub, IGNORECASE
from collections import defaultdict
from subprocess import run, PIPE
from networkx import DiGraph, has_path, find_cycle, NetworkXNoCycle, kamada_kawai_layout
from copy import deepcopy


def _get_OBJDUMP_ASSEMBLY(sde_files=None):
    assert(isinstance(sde_files, dict))

    objdp_asm = {}

    # offset, fn name
    fn_hdr = compile(r'^(\w+)\s+<(.+)>\s+\(File\s+Offset:\s+(\w+)\):$')
    # offset, instruction (+comment)
    fn_asm_part = compile(r'^\s+(\w+):\s+([\.\(]?\w{2,}.*)$')
    fn_asm_igno = compile(
        r'^(.*)\s+(<.*>)\s+\(\s*(File\s+Offset:\s+)(\w+)\s*\)$')
    __trash__ = compile(
        r'^$|^.*:\s+file format elf.*$|^Disassembly of section.*$')

    for fid in sde_files:
        FILE_NAME = sde_files[fid]

        objdp_asm[fid] = {'File': FILE_NAME, 'FileID': fid, 'FNs': {}, 'LO': 0}

        # ignore kernel lib
        if FILE_NAME == '[vdso]' or not path.exists(FILE_NAME):
            print('WRN 01: skipping missing file %s' % FILE_NAME)
            continue

        # get objdump for each file (main binary and all shared libs)
        p = run(['objdump',
                 '--disassemble',
                 '--disassembler-options=att-mnemonic',
                 '--no-show-raw-insn',
                 '--wide',
                 '--disassemble-zeroes',
                 '--file-offsets', FILE_NAME], stdout=PIPE)

        curr_fn_os, load_os = None, None
        for line in p.stdout.decode().splitlines():
            if fn_hdr.match(line):

                fn = fn_hdr.match(line)
                fn_os, fn_name, fn_real_os = \
                    fn.group(1).strip(), fn.group(2).strip(), \
                    fn.group(3).strip()

                # XXX: executable ELF64 binaries linked by ld have canonical
                #      base address offset of 0x400000, but SDE subtracts this
                if load_os is None:
                    load_os = int(fn_os, 16) - int(fn_real_os, 16)
                    if load_os == int('0x400000', 16):
                        objdp_asm[fid]['LO'] = load_os
                    elif load_os > 0:
                        exit('ERR: never seen before base address offset')

                curr_fn_os = int(fn_os, 16) - load_os
                os_str = '0x' + format(curr_fn_os, 'x')

                objdp_asm[fid]['FNs'][curr_fn_os] = {'Offset': os_str,
                                                     'Func': fn_name,
                                                     'ASM': []}
            elif fn_asm_part.match(line):
                assert(curr_fn_os is not None)

                fn = fn_asm_part.match(line)
                fn_asm_os, fn_asm_in = fn.group(1).strip(), fn.group(2).strip()
                # FIXME: strip comments at end as well??? => yes, llvm confused
                if fn_asm_igno.match(fn_asm_in):
                    instr, func, file_os, addr = \
                        fn_asm_igno.match(fn_asm_in).group(1), \
                        fn_asm_igno.match(fn_asm_in).group(2), \
                        fn_asm_igno.match(fn_asm_in).group(3), \
                        fn_asm_igno.match(fn_asm_in).group(4)
                    # mca doesn't like jump addresses without 0x, and jmpq is
                    # worse requiring a *0x before the address ... snowflakes
                    if instr.startswith('jmpq'):
                        prefix = '*0x'
                    else:
                        prefix = '0x'
                    instr = sub(r'(\s+)(%x)' % (load_os + int(addr, 16)),
                                r'\g<1>%s\g<2>' % prefix, instr,
                                count=0, flags=IGNORECASE)
                    # and add '#' before meaning data
                    fn_asm_in = '%s #%s (%s %s)' % (instr, func, file_os, addr)
                # XXX handle exceptions:
                # llvm hates: `bnd jmpq *%r11`, check llvm-objdump to replace
                fn_asm_in = sub(r'bnd\s+jmpq', r'repne\njmpq', fn_asm_in,
                                count=0, flags=IGNORECASE)
                # llvm hates: `fs addr32 nop`, check llvm-objdump to replace
                fn_asm_in = sub(r'fs\s+addr32\s+nop', r'nop', fn_asm_in,
                                count=0, flags=IGNORECASE)
                # llvm hates: `ds jmpq ...`, check llvm-objdump to replace
                fn_asm_in = sub(r'ds\s+jmpq\s+\*', r'jmpq *', fn_asm_in,
                                count=0, flags=IGNORECASE)

                fn_asm_os = '0x' + format(int(fn_asm_os, 16) - load_os, 'x')

                objdp_asm[fid]['FNs'][curr_fn_os]['ASM'].append([fn_asm_os,
                                                                 fn_asm_in])
            elif __trash__.match(line):
                continue
            else:
                exit('ERR: unknown line (%s) in objdump (%s)' %
                     (line, FILE_NAME))

    return objdp_asm


def _get_single_block_from_OBJDUMP(sde_file=None, from_addr=None,
                                   num_byte=None, load_os=0):
    assert(isinstance(sde_file, str) and isinstance(load_os, int)
           and isinstance(from_addr, int) and isinstance(num_byte, int))
    assert(sde_file is not '[vdso]')

    asm = []
    start_addr = load_os + from_addr

    if not path.exists(sde_file):
        print('WRN 10: %s missing on "this" node, or was deleted' % sde_file)
        return [['0x' + format(x, 'x'), 'nop']
                for x in range(start_addr, start_addr + num_byte)]

    fn_asm_part = compile(r'^\s+(\w+):\s+([\w\(\.].*)$')
    fn_asm_igno = compile(
        r'^(.*)\s+(<.*>)\s+\(\s*(File\s+Offset:\s+)(\w+)\s*\)$')

    # open interval [.,.), so --stop-address=0x... is excluded from output
    p = run(['objdump',
             '--disassemble',
             '--disassembler-options=att-mnemonic',
             '--no-show-raw-insn',
             '--wide',
             '--disassemble-zeroes',
             '--start-address=0x' + format(start_addr, 'x'),
             '--stop-address=0x' + format(start_addr + num_byte, 'x'),
             '--file-offsets', sde_file], stdout=PIPE)

    for line in p.stdout.decode().splitlines():
        if fn_asm_part.match(line):
            fn = fn_asm_part.match(line)
            fn_asm_os, fn_asm_in = fn.group(1).strip(), fn.group(2).strip()
            # FIXME: strip comments at end as well??? => yes, llvm confused
            if fn_asm_igno.match(fn_asm_in):
                instr, func, file_os, addr = \
                    fn_asm_igno.match(fn_asm_in).group(1), \
                    fn_asm_igno.match(fn_asm_in).group(2), \
                    fn_asm_igno.match(fn_asm_in).group(3), \
                    fn_asm_igno.match(fn_asm_in).group(4)
                # mca doesn't like jump addresses without 0x, and jmpq is
                # worse requiring a *0x before the address ... snowflakes
                if instr.startswith('jmpq'):
                    prefix = '*0x'
                else:
                    prefix = '0x'
                instr = sub(r'(\s+)(%x)' % (load_os + int(addr, 16)),
                            r'\g<1>%s\g<2>' % prefix, instr,
                            count=0, flags=IGNORECASE)
                # and add '#' before meaning data
                fn_asm_in = '%s #%s (%s %s)' % (instr, func, file_os, addr)
            # XXX handle exceptions:
            # llvm hates: `bnd jmpq *%r11`, check llvm-objdump to replace
            fn_asm_in = sub(r'bnd\s+jmpq', r'repne\njmpq', fn_asm_in,
                            count=0, flags=IGNORECASE)

            fn_asm_os = '0x' + format(int(fn_asm_os, 16) - load_os, 'x')
            asm.append([fn_asm_os, fn_asm_in])

    return asm


def _read_backup_OBJDUMP_ASSEMBLY(backup=None):
    assert(isinstance(backup, str))

    from pickle import load

    return load(open(backup, 'rb'))


def _store_backup_OBJDUMP_ASSEMBLY(objdp_asm=None, backup=None):
    assert(isinstance(objdp_asm, dict) and isinstance(backup, str))

    from pickle import dump

    return dump(objdp_asm, open(backup, 'wb'))


def _parse_FILE_NAMES(sde_bb_json=None):
    assert(isinstance(sde_bb_json, dict))

    # "FILE_NAMES" :
    #   [ [ "FILE_NAME_ID", "FILE_NAME" ],
    #     [ 2, "\/usr\/joe\/src\/misc\/hello-world" ],
    #     [ 5, "\/lib64\/libgcc_s.so" ], ...
    #   ]
    return dict(sde_bb_json['FILE_NAMES'][1:])


def _parse_SPECIAL_NODES(sde_bb_json=None):
    assert(isinstance(sde_bb_json, dict))

    # "SPECIAL_NODES" :
    #   [ [ "NODE_ID", "NODE_NAME" ],
    #     [ 1, "START" ],
    #     [ 2, "END" ],
    #     [ 3, "UNKNOWN" ]
    #   ]
    return dict(sde_bb_json['SPECIAL_NODES'][1:])


def _parse_EDGE_TYPES(sde_bb_json=None):
    assert(isinstance(sde_bb_json, dict))

    # "EDGE_TYPES" :
    #   [ [ "EDGE_TYPE_ID", "EDGE_TYPE" ],
    #     [ 1, "ENTRY" ],
    #     [ 2, "EXIT" ], ...
    #   }
    return dict(sde_bb_json['EDGE_TYPES'][1:])


def _parse_PROCESSES(sde_bb_json=None, sde_files=None, objdp_asm=None):
    assert(isinstance(sde_bb_json, dict) and isinstance(sde_files, dict)
           and isinstance(objdp_asm, dict))

    sde_procs = {}

    # "PROCESSES" :
    #   [ [ "PROCESS_ID", "PROCESS_DATA" ],
    #     [ 22814, dict{...} ],
    #     [ 958, dict{...} ], ...
    #   ]
    for PROCESS_ID, PROCESS_DATA in sde_bb_json['PROCESSES'][1:]:
        assert(PROCESS_ID not in sde_procs)

        # "PROCESS_DATA" :
        #   { "INSTR_COUNT" : 2134576,
        #     "INSTR_COUNT_PER_THREAD" : [ 1997750, 51676, 19794,18381, ...],
        #     "IMAGES" : [...],
        #     "EDGES" : [...]
        #   }
        # ==> "IMAGES" :    // basically just: binaries and shared objects
        #       [ [ "IMAGE_ID", "LOAD_ADDR", "SIZE", "IMAGE_DATA" ],
        #         [ 1, "0x400000", 2102216, {...} ],
        #         [ 2, "0x2aaaaaaab000", 1166728, {...} ]
        #       ]
        sde_procs[PROCESS_ID] = {'FileIDs': {}, 'Edges': {}}
        for _, LOAD_ADDR, _, IMAGE_DATA in PROCESS_DATA['IMAGES'][1:]:
            fid = IMAGE_DATA['FILE_NAME_ID']
            if sde_files[fid] == '[vdso]':
                _add_fake_vdso_info_to_asm(IMAGE_DATA, objdp_asm[fid]['FNs'])

            _, symbols, _, bb, rout = _parse_IMAGE_DATA(IMAGE_DATA, sde_files,
                                                        objdp_asm[fid]['FNs'])

            laddr = '0x' + format(int(LOAD_ADDR, 16), 'x')
            sde_procs[PROCESS_ID]['FileIDs'][fid] = {'FileID': fid,
                                                     'Offset': laddr,
                                                     'FNs': symbols,
                                                     'Blocks': bb,
                                                     'Routines': rout}

        _parse_EDGES(PROCESS_DATA['EDGES'], sde_procs[PROCESS_ID]['Edges'])

    # XXX: shouldn't have more than 1 process in here
    assert(1 == len(sde_procs))

    return sde_procs


def _parse_IMAGE_DATA(sde_image_data=None, sde_files=None, objdp_asm=None):
    assert(isinstance(sde_image_data, dict) and isinstance(sde_files, dict)
           and isinstance(objdp_asm, dict))

    # "IMAGE_DATA" :
    #   { "FILE_NAME_ID" : 4,
    #     "SYMBOLS" : [...],        // eq. functions (names and offsets)
    #     "SOURCE_DATA" : [...],    // XXX: SOURCE_DATA seems optional
    #     "BASIC_BLOCKS" : [...],
    #     "ROUTINES" : [...]        // FIXME: not sure, but its not full func
    #   }
    symbols, src_data, bb, routines = {}, {}, {}, {}

    if 'FILE_NAME_ID' in sde_image_data:
        FILE_NAME_ID = sde_image_data['FILE_NAME_ID']
    else:
        exit('ERR: no FILE_NAME_ID found in IMAGE_DATA')

    if 'SYMBOLS' in sde_image_data:
        _parse_SYMBOLS(sde_image_data, objdp_asm, symbols)
    if 'SOURCE_DATA' in sde_image_data:
        _parse_SOURCE_DATA(sde_image_data, src_data)
    if 'BASIC_BLOCKS' in sde_image_data:
        _parse_BASIC_BLOCKS(sde_image_data, sde_files, objdp_asm, bb)
    if 'ROUTINES' in sde_image_data:
        _parse_ROUTINES(sde_image_data, routines)

    return (FILE_NAME_ID, symbols, src_data, bb, routines)


def _parse_SYMBOLS(sde_image_data=None, objdp_asm=None, symbols=None):
    assert(isinstance(sde_image_data, dict) and isinstance(objdp_asm, dict)
           and isinstance(symbols, dict))

    # "SYMBOLS" :   // XXX: SDE incomplete, some symbols missing, check objdump
    #   [ [ "NAME", "ADDR_OFFSET", "SIZE" ],
    #     [ "free", "0x127f0", 60],
    #     [ "malloc", "0x12940", 13], ...
    #   ]
    for NAME, ADDR_OFFSET, SIZE in sde_image_data['SYMBOLS'][1:]:
        if NAME is not '.text':
            OFFSET, first = int(ADDR_OFFSET, 16), True
            if OFFSET in objdp_asm:
                NAME = objdp_asm[OFFSET]['Func']
            symbols[OFFSET] = {'Func': NAME, 'Size': SIZE,
                               'Offset': '0x' + format(OFFSET, 'x')}
        else:
            # some bigger functions seem to be split in 200k chunks and
            # the initial size is incorrect, but a sum of the .text chunks
            if first:
                symbols[OFFSET]['Size'] = 0
            symbols[OFFSET]['Size'] += SIZE


def _add_fake_vdso_info_to_asm(sde_image_data=None, objdp_asm=None):
    assert(isinstance(sde_image_data, dict) and isinstance(objdp_asm, dict))

    # can't find assembly of [vdso], and hence assume all is 'nop' instruction
    for NAME, ADDR_OFFSET, SIZE in sde_image_data['SYMBOLS'][1:]:
        curr_fn_os = int(ADDR_OFFSET, 16)
        os_str = '0x' + format(curr_fn_os, 'x')
        objdp_asm[curr_fn_os] = {'Offset': os_str,
                                 'Func': NAME,
                                 'ASM': [['0x' + format(x, 'x'), 'nop']
                                         for x in range(curr_fn_os,
                                                        curr_fn_os + SIZE)]}


def _parse_SOURCE_DATA(sde_image_data=None, src_data=None):
    assert(isinstance(sde_image_data, dict) and isinstance(src_data, dict))

    # "SOURCE_DATA" :
    #   [ [ "FILE_NAME_ID", "LINE_NUM", "ADDR_OFFSET", "SIZE", "NUM_INSTRS" ],
    #     [ 8, 25, "0x7a8", 4, 1 ],
    #     [ 8, 26, "0x7ac", 5, 1 ], ...
    #   ] -> FIXME: do we need this? maybe if we insert start/stop SDE cmd?
    return


def _parse_BASIC_BLOCKS(sde_image_data=None, sde_files=None, objdp_asm=None,
                        bb=None):
    assert(isinstance(sde_image_data, dict) and isinstance(sde_files, dict)
           and isinstance(objdp_asm, dict) and isinstance(bb, dict))

    # "BASIC_BLOCKS" :
    #   [ [ "NODE_ID", "ADDR_OFFSET", "SIZE", "NUM_INSTRS",
    #       "LAST_INSTR_OFFSET", "COUNT" ],
    #     [ 4, "0x7a8", 9, 2, 4, 1 ],
    #     [ 5, "0x7b1", 5, 1, 0, 9 ], ...
    #   ]
    # NODE_ID : is unique across the entire process
    # COUNT : total nr of times this block was executed across all threads
    for NODE_ID, ADDR_OFFSET, SIZE, NUM_INSTRS, _, COUNT in sde_image_data[
            'BASIC_BLOCKS'][1:]:
        fn_offset = _get_fn_os_for_block(objdp_asm, int(ADDR_OFFSET, 16))
        fn_name = objdp_asm[fn_offset]['Func']
        bb[NODE_ID] = {
            'Offset': '0x' + format(int(ADDR_OFFSET, 16), 'x'), 'Bytes': SIZE,
            'Func': fn_name, 'FuncOffset': fn_offset,
            'NumInst': NUM_INSTRS, 'ExecCnt': COUNT}
        if bb[NODE_ID]['Func'] is None:
            print('WRN 02: Cannot find block w/ offset %s in %s' %
                  (ADDR_OFFSET, sde_files[sde_image_data['FILE_NAME_ID']]))
            # can happen sometimes, SDE's SYMBOLS list seems incomplete, so
            # we should record the offset just in case we get ID/name later


def _get_fn_os_for_block(objdp_asm=None, bb_offset=None):
    assert(isinstance(objdp_asm, dict) and isinstance(bb_offset, int))

    for fn_offset in sorted(objdp_asm.keys(), reverse=False):
        if bb_offset >= fn_offset and \
                bb_offset <= int(objdp_asm[fn_offset]['ASM'][-1][0], 16):
            return fn_offset  # objdp_asm[fn_offset]['Func']
    else:
        exit('ERR: symbol (0x%s) not found basic block' %
             format(bb_offset, 'x'))
        # return None


def _parse_ROUTINES(sde_image_data=None, routines=None):
    assert(isinstance(sde_image_data, dict) and isinstance(routines, dict))

    # "ROUTINES" :  // XXX: those are NOT full functions... FUCK THIS S!*$
    #   [ [ "ENTRY_NODE_ID", "EXIT_NODE_IDS", "NODES", "LOOPS" ],
    #     [ 4, [ 4, 5, 6, 7 ],
    #       [ [ "NODE_ID", "IDOM_NODE_ID" ],
    #         [ 4, 4 ], ...
    #       ]]
    #     [ 63, [ 65, 67 ], ... ]
    #   ]
    # ENTRY_NODE_ID : first block of a function
    # NODES : list of all blocks belonging to same function (and their
    #         immediate dominators)
    for fn in sde_image_data['ROUTINES'][1:]:
        # some have the "LOOPS" array, while others don't
        if len(fn) == 4:
            ENTRY_NODE_ID, EXIT_NODE_IDS, NODES, _ = fn
        else:
            ENTRY_NODE_ID, EXIT_NODE_IDS, NODES = fn
        routines[ENTRY_NODE_ID] = sorted([nid for nid, _ in NODES[1:]],
                                         reverse=False)


def _parse_EDGES(sde_edge_data=None, edges=None):
    assert(isinstance(sde_edge_data, list) and isinstance(edges, dict))

    # "EDGES" :     // XXX: self-edges are possible
    #   [ ["EDGE_ID", "SOURCE_NODE_ID", "TARGET_NODE_ID", "EDGE_TYPE_ID",
    #      "COUNT_PER_THREAD" ],
    #     [ 4810, 1683, 3002, 16, [ 0, 0, 1, 1, 1, 0 ] ],
    #     [ 3460, 1597, 1598, 18, [ 1, 0, 0, 0, 0, 0 ] ], ...
    #   ]
    for dataset in sde_edge_data[1:]:
        EDGE_ID, SRC_BLOCK_ID, TARGET_BLOCK_ID, EDGE_TYPE_ID = dataset[0:4]
        COUNT_PER_THREAD = dataset[4]

        if SRC_BLOCK_ID not in edges:
            edges[SRC_BLOCK_ID] = {}
        else:
            assert(TARGET_BLOCK_ID not in edges[SRC_BLOCK_ID])

        edges[SRC_BLOCK_ID][TARGET_BLOCK_ID] = {
            'EdgeTypeID': EDGE_TYPE_ID,
            'ThreadExecCnts': COUNT_PER_THREAD}


def parse_SDE_JSON(args=None, sde_data=None):
    assert(isinstance(args, dict) and isinstance(sde_data, dict))

    from json import load as jsonload

    # JSON information follows the DCFG specification from:
    # from software.intel.com/sites/default/files/managed/32/2d/DCFG-format.pdf
    sde_bb_json = jsonload(bz2open(path.realpath(args.get('__sde_json_f__')),
                                   'rt'))
    assert(sde_bb_json is not None and
           6 == len(sde_bb_json.keys()) and
           'MAJOR_VERSION' in sde_bb_json and
           'MINOR_VERSION' in sde_bb_json and
           'FILE_NAMES' in sde_bb_json and
           'EDGE_TYPES' in sde_bb_json and
           'SPECIAL_NODES' in sde_bb_json and
           'PROCESSES' in sde_bb_json)

    sde_data['Files'] = _parse_FILE_NAMES(sde_bb_json)
    # first need to get the assembly, from either real files or previous runs
    if args.get('__load_objdump__') is None:
        sde_data['ObjDumpAsm'] = _get_OBJDUMP_ASSEMBLY(sde_data['Files'])
    else:
        sde_data['ObjDumpAsm'] = _read_backup_OBJDUMP_ASSEMBLY(
            path.realpath(args.get('__load_objdump__')))
    if args.get('__store_objdump__') is not None:
        _store_backup_OBJDUMP_ASSEMBLY(
            sde_data['ObjDumpAsm'],
            path.realpath(args.get('__store_objdump__')))
        exit("INFO: intermediate store of objdump data complete;\n"
             "      rerun with --load_objd %s to complete the workflow"
             % args.get('__store_objdump__'))

    # and then continue with the rest
    sde_data['SpecialBlockIDs'] = _parse_SPECIAL_NODES(sde_bb_json)
    sde_data['EdgeTypes'] = _parse_EDGE_TYPES(sde_bb_json)
    sde_data['Processes'] = _parse_PROCESSES(sde_bb_json, sde_data['Files'],
                                             sde_data['ObjDumpAsm'])


def parse_SDE_BB(args=None, sde_bb=None):
    assert(isinstance(args, dict) and isinstance(sde_bb, dict))

    # EXAMPLE BLOCK:
    #  BLOCK: 2586 PC: 530 ICOUNT: 14 EXECUTIONS: 1 #BYTES: 81 FN: main IMG: X
    #  XDIS 530: BASE     55                   push rbp
    #  XDIS 531: BASE     4889e5               mov rbp, rsp
    #  XDIS 534: BASE     4883e480             and rsp, 0xffffffffffffff80
    #  XDIS 538: BASE     4881ec80030000       sub rsp, 0x380

    with bz2open(path.realpath(args.get('__sde_block_f__')), 'rt') as bb_txt:
        # block number, program counter, #instruction (inst per block * #exec),
        # #{exec} of this block, #bytes of block, fn the block belongs to,
        # contained in binary/lib
        block_hdr = compile(
            r'^BLOCK:\s+(\d+)\s+PC:\s+(\w+)\s+ICOUNT:\s+(\d+)\s+' +
            r'EXECUTIONS:\s+(\d+)\s+#BYTES:\s+(\d+)(\s+FN:)?((?(6).*|\s+))' +
            r'IMG:\s+(.+)$')
        # program counter, instr.type (like AVX,FMA,etc), hex of instr,
        # and decoded instr./assembly (as seen by SDE)
        block_part = compile(r'^XDIS\s+(\w+):\s+(\w+)\s+(\w+)\s+(.*)$')
        no_block = compile('^$')

        curr_bb = None
        for line in bb_txt:
            if block_hdr.match(line):
                bb = block_hdr.match(line)
                bb_id, bb_pc, _, bb_ex, _, bb_fi = \
                    bb.group(1).strip(), bb.group(2).strip(), \
                    bb.group(3).strip(), bb.group(4).strip(), \
                    bb.group(5).strip(), bb.group(8).strip()
                bb_fn = bb.group(6)
                if bb_fn is not None:
                    bb_fn = bb.group(7).strip()
                curr_bb = int(bb_id)

                sde_bb[curr_bb] = {'BlockID': curr_bb, 'ProgCnt': bb_pc,
                                   'NumExec': bb_ex, 'Func': bb_fn,
                                   'File': bb_fi, 'ASM': []}

            elif block_part.match(line):
                bb = block_part.match(line)
                _, _, _, bb_in = \
                    bb.group(1).strip(), bb.group(2).strip(), \
                    bb.group(3).strip(), bb.group(4).strip()

                sde_bb[curr_bb]['ASM'].append(bb_in)

            elif no_block.match(line):
                continue
            else:
                exit('ERR: unknown line: %s' % line)


def _get_helping_mappers(sde_data=None):
    assert(isinstance(sde_data, dict))

    # 'sde_data' :
    #   { 'Files': {},              // executable and libs belonging to it
    #     'SpecialBlockIDs': {},    // block IDs without real instr
    #     'EdgeTypes': {},          // types of control flow edges
    #     'Processes': {},          // data for each process ID
    #     'ObjDumpAsm': {},         // assembly of files collected from objdump
    #   }

    # create quick lookup tables:
    #       file ID -to- file name
    fid2fn, fn2fid = sde_data['Files'], {}
    for file_id, file_name in fid2fn.items():
        fn2fid[file_name] = file_id
    #       special block ID -to- special block name
    sbid2sbn, sbn2sbid = sde_data['SpecialBlockIDs'], {}
    for spec_blk_id, spec_blk_name in sbid2sbn.items():
        sbn2sbid[spec_blk_name] = spec_blk_id
    #       edge type ID -to- edge type name
    etid2etn, etn2etid = sde_data['EdgeTypes'], {}
    for et_id, et_name in etid2etn.items():
        etn2etid[et_name] = et_id

    # 'sde_data'/'Processes' :
    #   { 'FileIDs': {ID: { 'FileID': ID,
    #                       'Offset': 0xYZ,
    #                       'FNs': {},
    #                       'Blocks': {},
    #                       'Routines': {} } },
    #     'Edges': {SRC: {DST: {'EdgeTypeID': ID,
    #                           'ThreadExecCnts': list} } }
    #   }

    #       basic block ID -to- process ID
    #       basic block ID -to- file ID
    #       offset/file -to- basic block ID
    bbid2pid, bbid2fid, offs2bbid = {}, {}, {}
    for pid, pid_data in sde_data['Processes'].items():
        for fid in pid_data['FileIDs']:
            offs2bbid[fid] = {}
            for bid in pid_data['FileIDs'][fid]['Blocks']:
                assert(bid not in bbid2pid)
                bbid2pid[bid] = pid
                assert(bid not in bbid2fid)
                bbid2fid[bid] = fid
                offset = pid_data['FileIDs'][fid]['Blocks'][bid]['Offset']
                assert(offset not in offs2bbid[fid])
                offs2bbid[fid][offset] = bid
    #       basic block ID -to- "routine" primary basic block
    bbid2rpbb = {}
    for pid, pid_data in sde_data['Processes'].items():
        for fid in pid_data['FileIDs']:
            for bid in pid_data['FileIDs'][fid]['Blocks']:
                for rpbbid in pid_data['FileIDs'][fid]['Routines']:
                    if bid in pid_data['FileIDs'][fid]['Routines'][rpbbid]:
                        assert(bid not in bbid2rpbb)
                        bbid2rpbb[bid] = rpbbid
                        break
                else:
                    print('WRN 03: BB ID (%s) missing from "routines" of %s' %
                          (bid, fid2fn[fid]))
    #       "routine" primary basic block -to- function primary basic block
    rpbb2fpbb, fpbb2rpbb = {}, defaultdict(list)
    for pid, pid_data in sde_data['Processes'].items():
        for fid in pid_data['FileIDs']:
            for rpbbid in pid_data['FileIDs'][fid]['Routines']:
                rpbbid_os = \
                    pid_data['FileIDs'][fid]['Blocks'][rpbbid]['Offset']
                found = False
                for fn_os in sde_data['ObjDumpAsm'][fid]['FNs']:
                    fn_asm = sde_data['ObjDumpAsm'][fid]['FNs'][fn_os]
                    for offset, _ in fn_asm['ASM']:
                        if offset == rpbbid_os:
                            # XXX: SDE doesn't have some blocks, like .plt.got
                            if fn_asm['Offset'] in offs2bbid[fid]:
                                rpbb2fpbb[rpbbid] = \
                                    offs2bbid[fid][fn_asm['Offset']]
                            else:
                                rpbb2fpbb[rpbbid] = offs2bbid[fid][offset]
                                print('WRN 04:', pid, fid, rpbbid, offset)
                            found = True
                            break
                    if found:
                        break
                else:
                    print('WRN 05: BB ID (%s) missing from objdump of %s' %
                          (rpbbid, fid2fn[fid]))
    for rpbb, fpbb in rpbb2fpbb.items():
        fpbb2rpbb[fpbb].append(rpbb)
    #       function/file -to- func primary block  // XXX: fn names not unique
    func2fpbb, fpbb2func = {}, {}
    for pid, pid_data in sde_data['Processes'].items():
        for fid in pid_data['FileIDs']:
            func2fpbb[fid] = defaultdict(list)
            for rpbbid in pid_data['FileIDs'][fid]['Routines']:
                fpbb = rpbb2fpbb[rpbbid]
                # diff. routine BB IDs can point to same function primary BB ID
                if fpbb not in fpbb2func:
                    fpbb2func[fpbb] = \
                        pid_data['FileIDs'][fid]['Blocks'][fpbb]['Func']
                # stupid WRN 04 situation can lead to instance where two (or
                # more) rpbbIDs belong to same function but the function itself
                # has no primary block in our list (eg. __cpu_indicator_init)
                func2fpbb[fid][fpbb2func[fpbb]].append(fpbb)

    return (fid2fn, fn2fid, sbid2sbn, sbn2sbid, etid2etn, etn2etid, bbid2pid,
            bbid2fid, offs2bbid, bbid2rpbb, rpbb2fpbb, fpbb2rpbb, func2fpbb,
            fpbb2func)


def convert_sde_data_to_something_usable(sde_data=None):
    assert(isinstance(sde_data, dict))

    data = {}

    # ___Naming Scheme___
    # fid2fn/fn2fid        :  file names & file IDs
    # sbid2sbn/sbn2sbid    :  special block ID & special block name
    # etid2etn/etn2etid    :  edge type ID & edge type name
    # bbid2pid             :  basic block ID & process ID
    # bbid2fid             :  basic block ID & file ID
    # offs2bbid            :  offsets & basic blocks
    # bbid2rpbb            :  basic block ID & "routine" primary basic block
    # rpbb2fpbb/fpbb2rpbb  :  "routine" p. basic block & func. p. basic block
    # func2fpbb/fpbb2func  :  function name/file & func. p. basic block

    fid2fn, fn2fid, sbid2sbn, sbn2sbid, etid2etn, etn2etid, bbid2pid, \
        bbid2fid, offs2bbid, bbid2rpbb, rpbb2fpbb, fpbb2rpbb, func2fpbb, \
        fpbb2func = _get_helping_mappers(sde_data)
    mapper = {'fid2fn': fid2fn, 'fn2fid': fn2fid,
              'sbid2sbn': sbid2sbn, 'sbn2sbid': sbn2sbid,
              'etid2etn': etid2etn, 'etn2etid': etn2etid,
              'bbid2pid': bbid2pid, 'bbid2fid': bbid2fid,
              'offs2bbid': offs2bbid,
              'bbid2rpbb': bbid2rpbb,
              'rpbb2fpbb': rpbb2fpbb, 'fpbb2rpbb': fpbb2rpbb,
              'func2fpbb': func2fpbb, 'fpbb2func': fpbb2func}

    maybe_sinks = set()
    for pid, pid_data in sde_data['Processes'].items():
        for source_bb in pid_data['Edges']:

            # if no edge in and no edge out, then ignore this block
            in_edges, out_edges = {}, {}
            for block, edges in pid_data['Edges'].items():
                # XXX: self-edges (source_bb == block ID) are possible!!!!!!!
                if source_bb == block and not source_bb not in edges:
                    continue
                elif source_bb in edges and \
                        sum(edges[source_bb]['ThreadExecCnts']) > 0:
                    in_edges[block] = deepcopy(edges[source_bb])
            for target_bb, edge_data in pid_data['Edges'][source_bb].items():
                if sum(edge_data['ThreadExecCnts']) > 0:
                    out_edges[target_bb] = deepcopy(edge_data)

            if len(in_edges) == 0 and len(out_edges) == 0:
                continue

            # get some metadata for block, but careful with special snowflakes
            if source_bb in sbid2sbn:
                func, file_name = sbid2sbn[source_bb], 'N/A'
            else:
                func = fpbb2func[rpbb2fpbb[bbid2rpbb[source_bb]]]
                file_name = fid2fn[bbid2fid[source_bb]]

            # get AT&T assembly for the block
            if source_bb in sbid2sbn:
                asm = [['0xffffffffffffffff', 'nop']]
            else:
                fid = bbid2fid[source_bb]
                blk_data = pid_data['FileIDs'][fid]['Blocks'][source_bb]
                objdp_asm = \
                    sde_data['ObjDumpAsm'][fid]['FNs'][blk_data['FuncOffset']]
                # sometimes blocks start in middle of other instructions, e.g.
                #    XDIS 4835ff: BINARY    BASE       00C3    add %al, %bl
                # and the PC jumps to 483600 and executes the C3 (meaning a
                # "ret" instruction...) -> these blocks won't be found in
                # native objdump data, so we have to get it manually ;-(
                try:
                    idx = [offset
                           for offset, _
                           in objdp_asm['ASM']].index(blk_data['Offset'])
                    asm = objdp_asm['ASM'][idx:idx + blk_data['NumInst']]
                except BaseException:
                    print('WRN 06: odd block; could not find [bb=%s] %s in %s'
                          % (source_bb, blk_data['Offset'], file_name)
                          + ' (perform fuzzy _get_single_block_from_OBJDUMP)')
                    asm = _get_single_block_from_OBJDUMP(
                        file_name, int(blk_data['Offset'], 16),
                        blk_data['Bytes'], sde_data['ObjDumpAsm'][fid]['LO'])
                    pass
                if len(asm) != blk_data['NumInst']:
                    print('WRN 07: strange block ID %s in fn %s (%s) with a'
                          % (source_bb, func, file_name)
                          + ' mismatch in len(asm) and instructions (%s vs %s)'
                          % (len(asm), blk_data['NumInst']) + '; please check')

            data[source_bb] = {'ID': source_bb,
                               'File': file_name,
                               'Func': func,
                               'ASM': asm,
                               'NumASM': len(asm),
                               'in_edges': in_edges,
                               'out_edges': out_edges}
            maybe_sinks.update(out_edges.keys())

            #if len(in_edges)>0 and len(out_edges)>0:
            #    print('jens:', source_bb)
            #    print('IN:', in_edges)
            #    for b in in_edges:
            #        if b in data: print(data[b])
            #    print('SELF:', data[source_bb])
            #    print('out:', out_edges)
            #    for b in out_edges:
            #        if b in data: print(data[b])


    for sink_bb in maybe_sinks:
        # only process remaining real sinks
        if sink_bb not in data:
            in_edges = {}
            for block, edges in pid_data['Edges'].items():
                if sink_bb in edges and \
                        sum(edges[sink_bb]['ThreadExecCnts']) > 0:
                    in_edges[block] = deepcopy(edges[sink_bb])

            # get some metadata for block, but careful with special snowflakes
            if sink_bb in sbid2sbn:
                func, file_name = sbid2sbn[sink_bb], 'N/A'
            else:
                func = fpbb2func[rpbb2fpbb[bbid2rpbb[sink_bb]]]
                file_name = fid2fn[bbid2fid[sink_bb]]

            # get AT&T assembly for the block
            if sink_bb in sbid2sbn:
                asm = [['0xffffffffffffffff', 'nop']]
            else:
                fid = bbid2fid[sink_bb]
                blk_data = pid_data['FileIDs'][fid]['Blocks'][sink_bb]
                objdp_asm = \
                    sde_data['ObjDumpAsm'][fid]['FNs'][blk_data['FuncOffset']]
                # sometimes blocks start in middle ... (see complaint above)
                try:
                    idx = [offset
                           for offset, _
                           in objdp_asm['ASM']].index(blk_data['Offset'])
                    asm = objdp_asm['ASM'][idx:idx + blk_data['NumInst']]
                except BaseException:
                    print('WRN 08: odd block; could not find [bb=%s] %s in %s'
                          % (sink_bb, blk_data['Offset'], file_name)
                          + ' (perform fuzzy _get_single_block_from_OBJDUMP)')
                    asm = _get_single_block_from_OBJDUMP(
                        file_name, int(blk_data['Offset'], 16),
                        blk_data['Bytes'], sde_data['ObjDumpAsm'][fid]['LO'])
                    pass
                if len(asm) != blk_data['NumInst']:
                    print('WRN 09: strange block ID %s in fn %s (%s) with a'
                          % (sink_bb, func, file_name)
                          + ' mismatch in len(asm) and instructions (%s vs %s)'
                          % (len(asm), blk_data['NumInst']) + '; please check')

            data[sink_bb] = {'ID': sink_bb,
                             'File': file_name,
                             'Func': func,
                             'ASM': asm,
                             'NumASM': len(asm),
                             'in_edges': in_edges,
                             'out_edges': {}}

    return (data, mapper)


def simulate_cycles_with_LLVM_MCA(blockdata=None, mapper=None):
    assert(isinstance(blockdata, dict) and isinstance(mapper, dict))

    # make the llvm-mca simulation more realistic / first step from gut feeling
    # assuming we have 3 basic blocks (random index fn of ld-linux-x86-64.so.2)
    # 1.
    #   mov    %rdi,%rdx
    #   and    $0x7,%edx
    #   mov    %rdi,%rax
    #   je     0x19c2a
    # 2.
    #   neg    %edx
    #   add    $0x8,%edx
    # 3.
    #   mov    (%rax),%cl
    #   cmp    %cl,%sil
    #   je     0x19da0
    # getting the timeline for each of them individually looks like this:
    #   [0,0]     DeER .  movq  %rdi, %rdx
    #   [0,1]     D=eER.  andl  $7, %edx
    #   [0,2]     DeE-R.  movq  %rdi, %rax
    #   [0,3]     D==eER  je   105514
    #  and
    #   [0,0]     DeER.   negl  %edx
    #   [0,1]     D=eER   addl  $8, %edx
    #  and
    #   [0,0]     DeeeeeER .   movb     (%rax), %cl
    #   [0,1]     D=====eER.   cmpb     %cl, %sil
    #   [0,2]     D======eER   je       105888
    # so 1. needs 6 cycles; second needs 5 cycles; last needs 10 cycles
    # but this is awfully pessimistic, so what if we add the previous block as
    # "preheating" phase and check the instruction retirement instead of the
    # pure cycles we can get more reasonable results maybe
    # assuming we combine 2. and 3. block for example:
    #   [0,0]     DeER .   .   negl     %edx
    #   [0,1]     D=eER.   .   addl     $8, %edx
    #   [0,2]     DeeeeeER .   movb     (%rax), %cl
    #   [0,3]     D=====eER.   cmpb     %cl, %sil
    #   [0,4]     .D=====eER   je       105888
    # and check the difference from last R of 2. block to last R of 3. block,
    # which is 5, or 10 respectively... meaning the "runtime for the 3. block
    # is 10-5=5 cycles and not 10 as predicted for standalone kernel

    asm_cnt = compile(r'^Instructions:\s+(\d+)$')
    cyc_cnt = compile(r'^Total Cycles:\s+(\d+)$')
    tl_pipe = compile(r'^(\[[\d,]+)\]\s+([\.D]+[DReE=\-\.\s]+)\s+([\.\(]?\w{2,}.*)$')

    for bbid, bdata in blockdata.items():
        for sink_bbid, sink_data in bdata['out_edges'].items():
            # llvm-mca read assembly from file, so dump a copy to ramdisk
            mca_in_fn = '/dev/shm/asm_%s_%s_%s.s' % (getpid(), bbid, sink_bbid)

            selfedge, icnt = False, 1
            num_asm_bbid, num_asm_sink_bbid = \
                blockdata[bbid]['NumASM'], blockdata[sink_bbid]['NumASM']

            with open(mca_in_fn, 'w') as mca_in_file:
                mca_in_file.write('\n'.join([instr
                                             for offset, instr
                                             in blockdata[bbid]['ASM']]))
                # self-edges indicate long running loops and compute phases
                # hence we skip adding sink_block and increase iteration count
                if bbid == sink_bbid:
                    selfedge = True
                    icnt = 100
                    #FIXME: have to change to account for ThreadExecCnts array
                else:
                    mca_in_file.write('\n')
                    mca_in_file.write('\n'.join([instr
                                                 for offset, instr
                                                 in blockdata[sink_bbid]['ASM']]))

            # mtriple: see http://clang.llvm.org/docs/CrossCompilation.html
            # get 'Total Cycles' from llvm-mca for each basic block
            p = run(['llvm-mca',
                     '--mtriple=x86_64-unknown-linux-gnu',
                     '--mcpu=native',
                     '--iterations=%s' % icnt,
                     '-timeline',
                     '-timeline-max-iterations=1',
                     '-timeline-max-cycles=2147483648', # call instr long
                     '-all-stats',
                     '-print-imm-hex', mca_in_fn], stdout=PIPE, stderr=PIPE)

            stderr = p.stderr.decode()
            stderr = sub(r'warning: found a return instruction in the input' +
                         ' assembly sequence.\nnote: program counter updates' +
                         ' are ignored.\n',
                         r'', stderr, count=0, flags=IGNORECASE)
            stderr = sub(r'warning: found a call in the input assembly' +
                         ' sequence.\nnote: call instructions are not' +
                         ' correctly modeled. Assume a latency of 100cy.\n',
                         r'', stderr, count=0, flags=IGNORECASE)
            if len(stderr):
                print(stderr)

            num_instr, num_cycles, timeline_data = 0, 0, []
            for line in p.stdout.decode().splitlines():
                if asm_cnt.match(line):
                    num_instr = int(asm_cnt.match(line).group(1).strip())

                elif cyc_cnt.match(line):
                    num_cycles = int(cyc_cnt.match(line).group(1).strip())

                elif tl_pipe.match(line):
                    tl = tl_pipe.match(line)
                    index, timeline, asm = \
                        tl.group(1).strip(), tl.group(2).strip(), \
                        tl.group(3).strip()
                    timeline_data.append([index, timeline, asm])

            assert(num_instr == len(timeline_data)              # merge 2 blk
                   or num_instr == icnt * len(timeline_data))   # self-edge

            # getting avg. estimated cycles per block is easy for self-edges
            if selfedge:
                cycles_per_iter = num_cycles / icnt
            # but for others we have to compare retirement (R) time difference
            # between the last instruction of the two blocks
            else:
                cycles_per_iter = \
                    timeline_data[num_instr - 1][1].find('R') \
                    - timeline_data[num_asm_bbid - 1][1].find('R')
            assert(cycles_per_iter >= 0)

            blockdata[bbid]['out_edges'][sink_bbid]['CyclesPerIteration'] = \
                cycles_per_iter

            remove(mca_in_fn)


def bb_graph(blockdata=None, mapper=None, thread_id=0):
    assert(isinstance(blockdata, dict) and isinstance(mapper, dict)
           and isinstance(thread_id, int) and thread_id > -1)

    # do depth-first traversal from START to build graph while avoiding cycles
    start_bbid = mapper['sbn2sbid']['START']
    if thread_id >= len(next(iter(blockdata[start_bbid]['out_edges'].values()))['ThreadExecCnts']):
        return None

    G = DiGraph()

    #for bbid, bdata in blockdata.items():
    #    for sink_bbid, sink_data in bdata['out_edges'].items():
    #        for thread in range(len(sink_data['ThreadExecCnts'])):
    #            w = sink_data['ThreadExecCnts'][thread] \
    #                * sink_data['CyclesPerIteration']
    #            G.add_edge(bbid, sink_bbid, cpu_cycles=w)
    #    if bbid==1:
    #        print(bdata)

    branches = [[start_bbid, out_edge]
                for out_edge in blockdata[start_bbid]['out_edges'].keys()]
    # build spanning tree
    while len(branches) > 0:
        curr_bbid, next_bbid = branches.pop()
        curr2next_bbid_edge = blockdata[curr_bbid]['out_edges'][next_bbid]

        # if the edge (curr_bbid->next_bbid) isn't used by the thread?
        iter_cnt = curr2next_bbid_edge['ThreadExecCnts'][thread_id]
        if iter_cnt < 1:
            continue
        else:
            edge_weight = iter_cnt * curr2next_bbid_edge['CyclesPerIteration']

        # check if we would close a (self-)loop -> convert to leaf
        if curr_bbid == next_bbid or G.has_node(next_bbid):
            G.add_edge(curr_bbid, '%s->%s->|' % (curr_bbid, next_bbid),
                       cpu_cycles=edge_weight)
            continue

        # or simply add the edge and continue
        G.add_edge(curr_bbid, next_bbid, cpu_cycles=edge_weight)

        branches += [[next_bbid, out_edge]
                     for out_edge in blockdata[next_bbid]['out_edges'].keys()]

    return G


def _id_in_range(interval, testid):
    if interval[0] > 0 and testid < interval[0]:
        return False
    if interval[1] > 0 and testid > interval[1]:
        return False
    return True


def fn_graph(data=None, mapper=None, block_id_range=None, fn_name=None):
    assert(isinstance(data, dict) and isinstance(mapper, dict)
           and isinstance(block_id_range, list) and isinstance(fn_name, str))

    import plotly.graph_objs as go

    G = DiGraph()

    # show all functions and their call-dependency
    for fid in mapper['func2fpbb']:
        for func, fpbbs in mapper['func2fpbb'][fid].items():
            for fpbb in fpbbs:
                added = []
                for block_id, block_data in data.items():
                    if block_data['Func'] != func:
                        continue
                    if mapper['rpbb2fpbb'][mapper['bbid2rpbb'][block_id]] \
                            != fpbb:
                        continue

                    for adj_func in set([data[blk]['Func']
                                         for blk in block_data['out_edges']]):
                        if [func, adj_func] in added:
                            continue
                        added.append([func, adj_func])
                        G.add_edge(func, adj_func)

    pos = kamada_kawai_layout(G)
    # pos = spectral_layout(G)
    for node in G.nodes:
        G.nodes[node]['pos'] = list(pos[node])

    plt_data, anno = [], []
    for edge in G.edges:
        x0, y0 = G.nodes[edge[0]]['pos']
        x1, y1 = G.nodes[edge[1]]['pos']
        trace = go.Scattergl(x=tuple([x0, x1, None]), y=tuple([y0, y1, None]),
                             mode='lines', line={'color': 'black'}, opacity=1)
        plt_data.append(trace)
        anno.append(
            dict(
                ax=x0,
                ay=y0,
                axref='x',
                ayref='y',
                x=(x1 - x0) / 2 + x0,
                y=(y1 - y0) / 2 + y0,
                xref='x',
                yref='y',
                showarrow=True,
                arrowhead=3,
                arrowsize=1.5,
                arrowwidth=1.5,
                opacity=1))

    trace = go.Scattergl(
        x=[],
        y=[],
        hovertext=[],
        text=[],
        mode='markers+text',
        textposition='bottom center',
        hoverinfo='text',
        marker={
            'size': 10,
            'color': 'LightSkyBlue'})
    for node in G.nodes():
        x, y = G.nodes[node]['pos']
        hovertext = 'NodeID: %s' % node
        trace['x'] += tuple([x])
        trace['y'] += tuple([y])
        trace['hovertext'] += tuple([hovertext])
        trace['text'] += tuple([hovertext])
    plt_data.append(trace)

    figure = {
        'data': plt_data,
        'layout': go.Layout(title='', showlegend=False, hovermode='closest',
                            margin={'b': 0, 'l': 0, 'r': 0, 't': 0},
                            xaxis={'showgrid': False, 'zeroline': False,
                                   'showticklabels': False},
                            yaxis={'showgrid': False, 'zeroline': False,
                                   'showticklabels': False},
                            height=800, clickmode='event+select',
                            annotations=anno,
                            )}

    return figure


def init_dash_for_vis(data=None, mapper=None):
    assert(isinstance(data, dict) and isinstance(mapper, dict))

    import dash
    import dash_core_components as dcc
    import dash_html_components as html

    block_id_range = ['', '']
    search_fn_name = ''

    # import a css template, and pass it into dash
    css = ['https://codepen.io/chriddyp/pen/bWLwgP.css']
    app = dash.Dash(__name__, external_stylesheets=css)
    app.title = 'Basic Block Graph Visualizer'

    styles = {
        'data': {
            'border': 'thin lightgrey solid',
            'overflowX': 'scroll',
            'overflowY': 'scroll'}}

    app.layout = html.Div([
        # Title
        html.Div([html.H1('%s' % app.title)],
                 className='row',
                 style={'textAlign': 'center'}),
        # define the full page layout
        html.Div(
            className='row',
            children=[
                # left side : input components on top and info stuff at bottom
                html.Div(
                    className='two columns',
                    children=[
                        html.Div(
                            className='twelve columns',
                            children=[
                                dcc.Markdown('**RANGE: From Block ID**'),
                                dcc.Input(
                                    id='my-block-from', type='text',
                                    value=block_id_range[0], debounce=True,
                                    placeholder='Block ID'),
                                html.Br(),
                                html.Div(id='output-text-block-id-from')
                            ],
                            style={'height': '100px'}
                        ),
                        html.Div(
                            className='twelve columns',
                            children=[
                                dcc.Markdown('**RANGE: To Block ID**'),
                                dcc.Input(
                                    id='my-block-to', type='text',
                                    value=block_id_range[1], debounce=True,
                                    placeholder='Block ID'),
                                html.Br(),
                                html.Div(id='output-text-block-id-to')
                            ],
                            style={'height': '100px'}
                        ),
                        html.Div(
                            className='twelve columns',
                            children=[
                                dcc.Markdown('**Search for Function Name**'),
                                dcc.Input(
                                    id='my-fn-input', type='text',
                                    value=search_fn_name, debounce=True,
                                    placeholder='Function Name'),
                                html.Br(),
                                html.Div(id='output-text-func-name')
                            ],
                            style={'height': '100px'}
                        ),
                        html.Div(
                            className='twelve columns',
                            children=[
                                dcc.Markdown(
                                    '**Hover-over Data**\n\n' +
                                    'Mouse over values in the graph.'),
                                html.Pre(id='hover-data', style=styles['data'])
                            ],
                            style={'height': '500px'}
                        ),
                    ]
                ),
                # middle : basic block graph component
                html.Div(
                    className='eight columns',
                    children=[dcc.Graph(id='my-graph',
                                        figure=fn_graph(data, mapper,
                                                        block_id_range,
                                                        search_fn_name))],
                ),
                html.Div(
                    className='two columns',
                    children=[
                        html.Div(
                            className='twelve columns',
                            children=[
                                dcc.Markdown(
                                    '**Full Block Data**\n\n' +
                                    'Click on vertices in the graph.'),
                                html.Pre(id='click-data', style=styles['data'])
                            ],
                            style={'height': '800px'}
                        ),
                    ],
                ),
            ],
        ),
    ])
    ###########################################################################
    # callback subroutines to control the interactivity in dash: ##############
    #   callback for left side components

    @app.callback(
        dash.dependencies.Output('my-graph', 'figure'),
        [dash.dependencies.Input('my-block-from', 'value'),
         dash.dependencies.Input('my-block-to', 'value'),
         dash.dependencies.Input('my-fn-input', 'value')])
    def update_block_range(input1, input2, input3):
        nonlocal data, mapper, block_id_range, search_fn_name
        block_id_range = [int(input1), int(input2)]
        search_fn_name = input2
        return fn_graph(data, mapper, [int(input1), int(input2)], input3)

    @app.callback(
        dash.dependencies.Output('hover-data', 'children'),
        [dash.dependencies.Input('my-graph', 'hoverData')])
    def display_hover_data(hoverData):
        return jsondumps(hoverData, indent=2)       # TODO: need correct update

    @app.callback(
        dash.dependencies.Output('click-data', 'children'),
        [dash.dependencies.Input('my-graph', 'clickData')])
    def display_click_data(clickData):
        return jsondumps(clickData, indent=2)       # TODO: need correct update
    ###########################################################################

    return app


def main():
    from psutil import cpu_freq
    from argparse import ArgumentParser

    arg_parser = ArgumentParser()
    arg_parser.add_argument('-j', '--sde_json', dest='__sde_json_f__',
                            help='name of existing *.dcfg.json.bz2 file',
                            type=str, metavar='<input file>', default=None,
                            required=True)
    arg_parser.add_argument('-b', '--sde_block', dest='__sde_block_f__',
                            help='name of existing *.bb.txt.bz2 file',
                            type=str, metavar='<input file>', default=None,
                            required=True)
    arg_parser.add_argument('-v', '--vis', dest='__visualize__',
                            help='run dash server and show basic block graph',
                            action='store_true', default=False)
    arg_parser.add_argument('-l', '--load_objd', dest='__load_objdump__',
                            help='load previously stored objdump data from' +
                            ' file instead of accessing all objects locally',
                            type=str, metavar='<input file>', default=None)
    arg_parser.add_argument('-s', '--store_objd', dest='__store_objdump__',
                            help='backup objdump data to file to postprocess' +
                            ' on a different computer or get speedup locally',
                            type=str, metavar='<output file>', default=None,)
    args = vars(arg_parser.parse_args())

    assert(args.get('__sde_json_f__') is not None
           and args.get('__sde_block_f__') is not None)
    assert(path.isfile(path.realpath(args.get('__sde_json_f__')))
           and path.isfile(path.realpath(args.get('__sde_block_f__'))))
    if args.get('__load_objdump__') is not None:
        assert(path.isfile(path.realpath(args.get('__load_objdump__'))))

    sde_data = {'BasicBlocks': {}}
    parse_SDE_JSON(args, sde_data)
    parse_SDE_BB(args, sde_data['BasicBlocks'])

    data, mapper = convert_sde_data_to_something_usable(sde_data)
    del sde_data

    cpufreq = cpu_freq()
    simulate_cycles_with_LLVM_MCA(data, mapper)
    total_cycles_per_thread, thread_id = [], 0
    while True:
        G = bb_graph(data, mapper, thread_id)
        if G == None:
            break
        thread_id += 1

        #bbids= set(data.keys())
        #nodes= set(list(G))
        #print('jens G:', G)
        #print('jens bbids:', bbids)
        #print('jens nodes:', nodes)
        #print('jens diff1:', nodes.difference(bbids))
        #print('jens diff2:', bbids.difference(nodes))
        #print("Edges of graph: ")
        #print(G.edges.data())

        total_cycles = G.size(weight='cpu_cycles')
        print('Total CPU cycles on rank %s and thread ID %s : %s\n'
              % (0, thread_id, total_cycles) +
              '(Converted to time (with min/curr/max freq.): %ss / %ss / %ss)'
              % (total_cycles / (cpufreq.min * pow(10, 6)),
                 total_cycles / (cpufreq.current * pow(10, 6)),
                 total_cycles / (cpufreq.max * pow(10, 6))))

    #for x in nx.simple_cycles(G):
    #    if 3 not in x:
    #        print(x)

#    total_cycles = nx.dag_longest_path_length(G, weight='cycles')
#    print('Extrapolated runtime: %s', total_cycles)

    if args.get('__visualize__'):
        dash_app = init_dash_for_vis(data, mapper)
        dash_app.run_server(debug=False)


if __name__ == '__main__':
    main()
