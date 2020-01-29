#!/usr/bin/env python3

# usage:
# <script> ./dcfg-out.dcfg.json.bz2 ./dcfg-out.bb.txt.bz2
#  ( create files: ./dcfg-out.dcfg.json.bz2 and ./dcfg-out.bb.txt.bz2 via
#      `sde64 -dcfg 1 -dcfg:write_bb 1 ...` )

from os import path
from sys import argv, exit
from bz2 import open as bz2open
from json import load as jsonload
from re import compile
import subprocess as subp

if not (len(argv) == 3 and path.isfile(path.realpath(
        argv[1])) and path.isfile(path.realpath(argv[2]))):
    exit("ERR: usage: <script> ./*.dcfg.json.bz2 ./*.bb.txt.bz2")

# basic block information from SDE
basic_blocks = {}
# basic blocks from objdump
objdp_blocks = {}


def parse_FILE_NAMES(sde_bb_json=None, sde_files=None, objdp_blocks=None):
    assert(isinstance(sde_bb_json, dict)
           and isinstance(sde_files, dict)
           and isinstance(objdp_blocks, dict))
    ###########################################################################
    # offset, fn name
    block_hdr = compile(r'^(\w+)\s+<(.+)>:$')
    # offset, instruction (+comment)
    block_part = compile(r'^\s+\w+:\s+([\w\(\.].*)$')
    no_block = compile(
        r'^$|^.*:\s+file format elf.*$|^Disassembly of section.*$')
    # "FILE_NAMES" :
    #   [ [ "FILE_NAME_ID", "FILE_NAME" ],
    #     [ 2, “\/usr\/joe\/src\/misc\/hello-world" ],
    #     [ 5, "\/lib64\/libgcc_s.so" ], ...
    #   ]
    for FILE_NAME_ID, FILE_NAME in sde_bb_json['FILE_NAMES'][1:]:
        sde_files[FILE_NAME_ID] = FILE_NAME
        continue    # TODO: take out

        objdp_blocks[FILE_NAME_ID] = {
            'File': FILE_NAME,
            'FileID': FILE_NAME_ID,
            'Blocks': {}}

        # ignore kernel lib
        if FILE_NAME == '[vdso]':
            continue

        if not path.exists(FILE_NAME):
            print('WRN: skipping missing file %s' % FILE_NAME)
            continue

        # get objdump for each file (main binary and all shared libs)
        p = subp.run(['objdump',
                      '--disassemble',
                      '--disassembler-options=att-mnemonic',
                      '--no-show-raw-insn',
                      '--wide',
                      '--disassemble-zeroes',
                      FILE_NAME],
                     stdout=subp.PIPE)

        curr_bb = None
        for line in p.stdout.decode().splitlines():
            if block_hdr.match(line):
                bb = block_hdr.match(line)
                bb_os, bb_fn = bb.group(1).strip(), bb.group(2).strip()
                curr_bb = format(int(bb_os, 16), 'x')
                objdp_blocks[FILE_NAME_ID]['Blocks'][curr_bb] = {
                    'BlockID': curr_bb,
                    'Offset': bb_os,
                    'Func': bb_fn,
                    'Inst': []}
            elif block_part.match(line):
                bb = block_part.match(line)
                bb_in = bb.group(1).strip()  # FIXME: strip comment?
                objdp_blocks[FILE_NAME_ID]['Blocks'][curr_bb]['Inst'].append(
                    bb_in)
            elif no_block.match(line):
                continue
            else:
                exit(
                    'ERR: unknown line (%s) in objdump (%s)' %
                    (line, FILE_NAME))
        print(
            FILE_NAME_ID,
            FILE_NAME,
            objdp_blocks[FILE_NAME_ID]['Blocks'][curr_bb])  # print last


def parse_SPECIAL_NODES(sde_bb_json=None, sde_special_nodes=None):
    assert(
        isinstance(
            sde_bb_json,
            dict) and isinstance(
            sde_special_nodes,
            dict))
    ###########################################################################
    # "SPECIAL_NODES" :
    #   [ [ "NODE_ID", "NODE_NAME" ],
    #     [ 1, "START" ],
    #     [ 2, "END" ],
    #     [ 3, "UNKNOWN" ]
    #   ]
    for NODE_ID, NODE_NAME in sde_bb_json['SPECIAL_NODES'][1:]:
        sde_special_nodes[NODE_ID] = NODE_NAME


def parse_EDGE_TYPES(sde_bb_json=None, sde_edge_types=None):
    assert(isinstance(sde_bb_json, dict) and isinstance(sde_edge_types, dict))
    ###########################################################################
    # "EDGE_TYPES" :
    #   [ [ "EDGE_TYPE_ID", "EDGE_TYPE" ],
    #     [ 1, "ENTRY" ],
    #     [ 2, "EXIT" ], ...
    #   }
    for EDGE_TYPE_ID, EDGE_TYPE in sde_bb_json['EDGE_TYPES'][1:]:
        sde_edge_types[EDGE_TYPE_ID] = EDGE_TYPE


def parse_PROCESSES(sde_bb_json=None, sde_procs=None):
    assert(isinstance(sde_bb_json, dict) and isinstance(sde_procs, dict))
    ###########################################################################
    # "PROCESSES" :
    #   [ [ "PROCESS_ID", "PROCESS_DATA" ],
    #     [ 22814, dict{...} ],
    #     [ 958, dict{...} ], ...
    #   ]
    for PROCESS_ID, PROCESS_DATA in sde_bb_json['PROCESSES'][1:]:
        # "PROCESS_DATA" :
        #   { "INSTR_COUNT" : 2134576,
        #     "INSTR_COUNT_PER_THREAD" : [ 1997750, 51676, 19794,18381, ...],
        #     "IMAGES" : [...],
        #     "EDGES" : [...]
        #   }
        # ==> "IMAGES" :
        #       [ [ "IMAGE_ID", "LOAD_ADDR", "SIZE", "IMAGE_DATA" ],
        #         [ 1, "0x400000", 2102216, {...} ],
        #         [ 2, "0x2aaaaaaab000", 1166728, {...} ]
        #       ]
        assert(PROCESS_ID not in sde_procs)
        sde_procs[PROCESS_ID] = {'FileIDs': {}, 'Edges': {}}
        for _, LOAD_ADDR, _, IMAGE_DATA in PROCESS_DATA['IMAGES'][1:]:
            FILE_NAME_ID, symbols, _, bb, _ = parse_IMAGE_DATA(IMAGE_DATA)
            sde_procs[PROCESS_ID]['FileIDs'][FILE_NAME_ID] = {
                'FileID': FILE_NAME_ID, 'Offset': LOAD_ADDR, 'FNs': symbols,
                'Blocks': bb}
        parse_EDGES(PROCESS_DATA['EDGES'], sde_procs[PROCESS_ID]['Edges'])


def parse_IMAGE_DATA(sde_image_data=None):
    assert(isinstance(sde_image_data, dict))
    # "IMAGE_DATA" :
    #   { "FILE_NAME_ID" : 4,
    #     "SYMBOLS" : [...],
    #     "SOURCE_DATA" : [...],
    #     "BASIC_BLOCKS" : [...],
    #     "ROUTINES" : [...]
    #   } -> FIXME: SOURCE_DATA seems optional
    symbols, src_data, bb, routines = None, None, None, None
    if 'FILE_NAME_ID' in sde_image_data:
        FILE_NAME_ID = sde_image_data['FILE_NAME_ID']
    else:
        exit('ERR: no FILE_NAME_ID found in IMAGE_DATA')
    #
    if 'SYMBOLS' in sde_image_data:
        # "SYMBOLS" :
        #   [ [ "NAME", "ADDR_OFFSET", "SIZE" ],
        #     [ "free", "0x127f0", 60],
        #     [ "malloc", "0x12940", 13], ...
        #   ]
        symbols = {}
        for NAME, ADDR_OFFSET, SIZE in sde_image_data['SYMBOLS'][1:]:
            symbols[int(ADDR_OFFSET, 16)] = {
                'Func': NAME, 'Offset': ADDR_OFFSET, 'Size': SIZE}
    #
    if 'SOURCE_DATA' in sde_image_data:
        # "SOURCE_DATA" :
        #   [ [ "FILE_NAME_ID", "LINE_NUM", "ADDR_OFFSET","SIZE","NUM_INSTRS" ],
        #     [ 8, 25, "0x7a8", 4, 1 ],
        #     [ 8, 26, "0x7ac", 5, 1 ], ...
        #   ] -> FIXME: do we need this? maybe if we insert start/stop SDE cmd?
        src_data = {}
    #
    if 'BASIC_BLOCKS' in sde_image_data:
        # "BASIC_BLOCKS" :
        #   [ [ "NODE_ID", "ADDR_OFFSET", "SIZE", "NUM_INSTRS","LAST_INSTR_OFFSET", "COUNT" ],
        #     [4, "0x7a8", 9, 2, 4, 1 ],
        #     [ 5, "0x7b1", 5, 1, 0, 9], ...
        #   ]
        # NODE_ID : is unique across the entire process
        # COUNT : total nr of times this block was executed across all threads
        bb = {}
        for NODE_ID, ADDR_OFFSET, _, _, _, COUNT in sde_image_data['BASIC_BLOCKS'][1:]:
            bb[NODE_ID] = {
                'SymbolID': get_symbol_from_block_offset(
                    symbols,
                    ADDR_OFFSET),
                'Offset': ADDR_OFFSET,
                'ExecCnt': COUNT}
            if bb[NODE_ID]['SymbolID'] is None:
                print('WRN: Cannot find block w/ offset %s in %s' %
                      (ADDR_OFFSET, sde_files[sde_image_data['FILE_NAME_ID']]))
                # can happen sometimes, SDE's SYMBOLS list seems incomplete, so
                # we should record the offset just in case we get ID/name later
    #
    if 'ROUTINES' in sde_image_data:
        # "ROUTINES" :
        #   [ [ "ENTRY_NODE_ID", "EXIT_NODE_IDS", "NODES", "LOOPS" ],
        #     [ 4, [ 4, 5, 6, 7 ],
        #       [ [ "NODE_ID", "IDOM_NODE_ID" ],
        #         [ 4, 4 ], ...
        #       ]]
        #     [ 63, [ 65, 67 ], ... ]
        #   ] -> FIXME: do we need this? not sure yet
        routines = {}
    #
    return (FILE_NAME_ID, symbols, src_data, bb, routines)


def get_symbol_from_block_offset(symbols=None, bb_offset=None):
    assert(isinstance(symbols, dict) and isinstance(bb_offset, str))
    for sym in sorted(symbols.keys(), reverse=False):
        if int(bb_offset, 16) >= sym and int(
                bb_offset, 16) < sym + symbols[sym]['Size']:
            return sym
    else:
        return None
        #exit('ERR: symbol not found basic block')


def parse_EDGES(sde_edge_data=None, edges=None):
    assert(isinstance(sde_edge_data, list) and isinstance(edges, dict))
    # "EDGES" :
    #   [ ["EDGE_ID", "SOURCE_NODE_ID", "TARGET_NODE_ID", "EDGE_TYPE_ID", "COUNT_PER_THREAD" ],
    #     [ 4810, 1683, 3002, 16, [ 0, 0, 1, 1, 1, 0 ] ],
    #     [ 3460, 1597, 1598, 18, [ 1, 0, 0, 0, 0, 0 ] ], ...
    #   ]
    for EDGE_ID, SOURCE_NODE_ID, TARGET_NODE_ID, EDGE_TYPE_ID, COUNT_PER_THREAD in sde_edge_data:
        if SOURCE_NODE_ID in edges:
            assert(EDGE_ID not in edges[SOURCE_NODE_ID])
            edges[SOURCE_NODE_ID][EDGE_ID] = {
                'TgtNode': TARGET_NODE_ID,
                'EdgeTypeID': EDGE_TYPE_ID,
                'ThreadExecCnts': COUNT_PER_THREAD}
        else:
            edges[SOURCE_NODE_ID] = {
                EDGE_ID: {
                    'TgtNode': TARGET_NODE_ID,
                    'EdgeTypeID': EDGE_TYPE_ID,
                    'ThreadExecCnts': COUNT_PER_THREAD}}


sde_files = {}
sde_edge_types = {}
sde_special_nodes = {}  # no clue if needed
sde_procs = {}  # FIXME: are there more than 1?

sde_bb_json = jsonload(bz2open(path.realpath(argv[1]), 'rt'))
if sde_bb_json is not None:
    # get all binary and library data, including assembly from objdump
    parse_FILE_NAMES(sde_bb_json, sde_files, objdp_blocks)
    #
    parse_EDGE_TYPES(sde_bb_json, sde_edge_types)
    #
    parse_SPECIAL_NODES(sde_bb_json, sde_special_nodes)
    #
    parse_PROCESSES(sde_bb_json, sde_procs)

exit()

with bz2open(path.realpath(argv[2]), 'rt') as sde_bb_txt:
    # block number, program counter, #instruction (inst per block * #exec),
    # #exec of this block, #bytes of block, fn the block belongs to, contained
    # in binary/lib
    block_hdr = compile(
        r'^BLOCK:\s+(\d+)\s+PC:\s+(\w+)\s+ICOUNT:\s+(\d+)\s+' +
        r'EXECUTIONS:\s+(\d+)\s+#BYTES:\s+(\d+)(\s+FN:)?((?(6).*|\s+))' +
        r'IMG:\s+(.+)$')
    # program counter, ?, instruction/assembly
    block_part = compile(r'^XDIS\s+(\w+):\s+(\w+)\s+(\w+)\s+(.*)$')
    no_block = compile('^$')
    curr_bb = None
    for line in sde_bb_txt:
        if block_hdr.match(line):
            bb = block_hdr.match(line)
            bb_id, bb_pc, bb_ic, bb_ex, bb_by, bb_fi = \
                bb.group(1).strip(), bb.group(2).strip(), \
                bb.group(3).strip(), bb.group(4).strip(), \
                bb.group(5).strip(), bb.group(8).strip()
            bb_fn = bb.group(6)
            if bb_fn is not None:
                bb_fn = bb.group(7).strip()
            curr_bb = int(bb_id)
            basic_blocks[curr_bb] = \
                {'BlockID': curr_bb, 'ProgCnt': bb_pc,
                 'NumExec': bb_ex, 'Func': bb_fn,
                 'File': bb_fi, 'Inst': []}
        elif block_part.match(line):
            bb = block_part.match(line)
            bb_it = bb.group(2).strip()  # type, like AVX, FMA, etc
            bb_in = bb.group(4).strip()  # assembly as seen by SDE
            basic_blocks[curr_bb]['Inst'].append(bb_in)
        elif no_block.match(line):
            continue
        else:
            exit('ERR: unknown line: %s' % line)
