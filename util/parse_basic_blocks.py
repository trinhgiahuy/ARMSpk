#!/usr/bin/env python3

# usage:
# <script> -j ./dcfg-out.dcfg.json.bz2 -b ./dcfg-out.bb.txt.bz2
#  ( create files: ./dcfg-out.dcfg.json.bz2 and ./dcfg-out.bb.txt.bz2 via
#      `sde64 -dcfg 1 -dcfg:write_bb 1 ...` )

from os import path
from sys import exit
from bz2 import open as bz2open
from json import dumps as jsondumps     # TODO: take this out
from re import compile


def _get_OBJDUMP_ASSEMBLY(sde_files=None):
    assert(isinstance(sde_files, dict))

    import subprocess as subp

    objdp_asm = {}

    # offset, fn name
    fn_hdr = compile(r'^(\w+)\s+<(.+)>\s+\(File Offset:\s+(\w+)\):$')
    # offset, instruction (+comment)
    fn_asm_part = compile(r'^\s+(\w+):\s+([\w\(\.].*)$')
    __trash__ = compile(
        r'^$|^.*:\s+file format elf.*$|^Disassembly of section.*$')

    for fid in sde_files:
        FILE_NAME = sde_files[fid]

        objdp_asm[fid] = {'File': FILE_NAME, 'FileID': fid, 'FNs': {}}

        # ignore kernel lib
        if FILE_NAME == '[vdso]' or not path.exists(FILE_NAME):
            print('WRN 01: skipping missing file %s' % FILE_NAME)
            continue

        # get objdump for each file (main binary and all shared libs)
        p = subp.run(['objdump',
                      '--disassemble',
                      '--disassembler-options=att-mnemonic',
                      '--no-show-raw-insn',
                      '--wide',
                      '--disassemble-zeroes',
                      '--file-offsets',
                      FILE_NAME],
                     stdout=subp.PIPE)

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

                curr_fn_os = int(fn_os, 16) - load_os
                os_str = '0x' + format(curr_fn_os, 'x')

                objdp_asm[fid]['FNs'][curr_fn_os] = {'Offset': os_str,
                                                     'Func': fn_name,
                                                     'ASM': []}

            elif fn_asm_part.match(line):
                assert(curr_fn_os is not None)

                fn = fn_asm_part.match(line)
                fn_asm_os, fn_asm_in = \
                    fn.group(1).strip(), \
                    fn.group(2).strip()  # FIXME: strip comments at end?
                fn_asm_os = '0x' + format(int(fn_asm_os, 16) - load_os, 'x')

                objdp_asm[fid]['FNs'][curr_fn_os]['ASM'].append([fn_asm_os,
                                                                 fn_asm_in])

            elif __trash__.match(line):
                continue
            else:
                exit(
                    'ERR: unknown line (%s) in objdump (%s)' %
                    (line, FILE_NAME))

    return objdp_asm


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
                continue

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
    assert(1 == len(sde_procs.keys()))

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
            symbols[OFFSET] = {'Func': objdp_asm[OFFSET]['Func'],
                               'Offset': '0x' + format(OFFSET, 'x'),
                               'Size': SIZE}
        else:
            # some bigger functions seem to be split in 200k chunks and
            # the initial size is incorrect, but a sum of the .text chunks
            if first:
                symbols[OFFSET]['Size'] = 0
            symbols[OFFSET]['Size'] += SIZE


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
    for NODE_ID, ADDR_OFFSET, _, NUM_INSTRS, _, COUNT in sde_image_data[
            'BASIC_BLOCKS'][1:]:
        bb[NODE_ID] = {
            'Func': _get_fn_from_block_offset(objdp_asm, int(ADDR_OFFSET, 16)),
            'Offset': '0x' + format(int(ADDR_OFFSET, 16), 'x'),
            'NumInst': NUM_INSTRS,
            'ExecCnt': COUNT}
        if bb[NODE_ID]['Func'] is None:
            print('WRN 02: Cannot find block w/ offset %s in %s' %
                  (ADDR_OFFSET, sde_files[sde_image_data['FILE_NAME_ID']]))
            # can happen sometimes, SDE's SYMBOLS list seems incomplete, so
            # we should record the offset just in case we get ID/name later


def _get_fn_from_block_offset(asm=None, bb_offset=None):
    assert(isinstance(asm, dict) and isinstance(bb_offset, int))

    for fn_offset in sorted(asm.keys(), reverse=False):
        if bb_offset >= fn_offset and \
                bb_offset <= int(asm[fn_offset]['ASM'][-1][0], 16):
            return asm[fn_offset]['Func']
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
        routines[ENTRY_NODE_ID] = []
        for nid, dom in NODES[1:]:
            routines[ENTRY_NODE_ID].append(nid)
        routines[ENTRY_NODE_ID].sort()


def _parse_EDGES(sde_edge_data=None, edges=None):
    assert(isinstance(sde_edge_data, list) and isinstance(edges, dict))

    # "EDGES" :
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
    assert(sde_bb_json is not None)

    sde_data['Files'] = _parse_FILE_NAMES(sde_bb_json)
    # first need to get the assembly
    sde_data['ObjDumpAsm'] = _get_OBJDUMP_ASSEMBLY(sde_data['Files'])
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


def convert_sde_data_to_something_usable(sde_data=None):
    assert(isinstance(sde_data, dict))

    data = {}

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
    for file_id in fid2fn:
        fn2fid[fid2fn[file_id]] = file_id
    print('\n', fid2fn, '\n', fn2fid)
    #       special block ID -to- special block name
    sbid2sbn, sbn2sbid = sde_data['SpecialBlockIDs'], {}
    for spec_blk_id in sbid2sbn:
        sbn2sbid[sbid2sbn[spec_blk_id]] = spec_blk_id
    print('\n', sbid2sbn, '\n', sbn2sbid)
    #       edge type ID -to- edge type name
    etid2etn, etn2etid = sde_data['EdgeTypes'], {}
    for et_id in etid2etn:
        etn2etid[etid2etn[et_id]] = et_id
    print('\n', etid2etn, '\n', etn2etid)

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
    for pid in sde_data['Processes']:
        pid_data = sde_data['Processes'][pid]
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
    print('\n', bbid2pid, '\n', bbid2fid)
    #       basic block ID -to- "routine" primary basic block
    bbid2rpbb = {}
    for pid in sde_data['Processes']:
        pid_data = sde_data['Processes'][pid]
        for fid in pid_data['FileIDs']:
            for bid in pid_data['FileIDs'][fid]['Blocks']:
                for rpbbid in pid_data['FileIDs'][fid]['Routines']:
                    if bid in pid_data['FileIDs'][fid]['Routines'][rpbbid]:
                        assert(bid not in bbid2rpbb)
                        bbid2rpbb[bid] = rpbbid
                        break
                else:
                    print('WRN 03: BB ID (%s) missing from "routines" of %s' %
                          (bid, sde_data['Files'][fid]))
    print('\n', bbid2rpbb)
    #       "routine" primary basic block -to- function primary basic block
    rpbb2fpbb = {}
    for pid in sde_data['Processes']:
        pid_data = sde_data['Processes'][pid]
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
                          (rpbbid, sde_data['Files'][fid]))
    print('\n', rpbb2fpbb)
    #       function/file -to- func primary block  // XXX: fn names not unique
    func2fpbb, fpbb2func = {}, {}
    for pid in sde_data['Processes']:
        pid_data = sde_data['Processes'][pid]
        for fid in pid_data['FileIDs']:
            func2fpbb[fid] = {}
            for rpbbid in pid_data['FileIDs'][fid]['Routines']:
                fpbb = rpbb2fpbb[rpbbid]
                # diff. routine BB IDs can point to same function primary BB ID
                if fpbb not in fpbb2func:
                    fpbb2func[fpbb] = \
                        pid_data['FileIDs'][fid]['Blocks'][fpbb]['Func']
                # stupid WRN 04 situation can lead to instance where two (or
                # more rpbbIDs belong to same function but the function itself
                # has no primary block in our list (eg. __cpu_indicator_init)
                if fpbb2func[fpbb] not in func2fpbb[fid]:
                    func2fpbb[fid][fpbb2func[fpbb]] = [fpbb]
                elif fpbb not in func2fpbb[fid][fpbb2func[fpbb]]:
                    func2fpbb[fid][fpbb2func[fpbb]].append(fpbb)
    print('\n', func2fpbb, '\n', fpbb2func)

#        for fid in sde_data['Processes'][pid]['FileIDs']:
#            print('== %s' % sde_data['Files'][fid])
#
#            for sid in sde_data['Processes'][pid]['FileIDs'][fid]['FNs']:
#                prt_hdr = True
#
#                for nid in sde_data['Processes'][pid]['FileIDs'][fid]['Blocks']:
#                    if sde_data['Processes'][pid]['FileIDs'][fid]['Blocks'][nid]['SymbolID'] == sid:
#                        if prt_hdr:
#                            print(
#                                '==== %s' %
#                                sde_data['Processes'][pid]['FileIDs'][fid]['FNs'][sid]['Func'])
#                            prt_hdr = False
#                        print(
#                            '====== %s (#exec: %s)' %
#                            (nid, sde_data['Processes'][pid]['FileIDs'][fid]['Blocks'][nid]['ExecCnt']))
#                        for src_nid in sde_data['Processes'][pid]['Edges']:
#                            for edge_id in sde_data['Processes'][pid]['Edges'][src_nid]:
#                                #print(src_nid, edge_id, sde_data['Processes'][pid]['Edges'][src_nid][edge_id]['TgtNode'], nid)
#                                if sde_data['Processes'][pid]['Edges'][src_nid][edge_id]['TgtNode'] == nid:
#                                    print('======== %s->%s (max per thread: %s)' %
#                                          (src_nid, nid, ','.join([str(x) for x in sde_data['Processes'][pid]['Edges'][src_nid][edge_id]['ThreadExecCnts'] if x > 0])))
#                        for edge_id in sde_data['Processes'][pid]['Edges'][nid]:
#                            print('======== %s->%s (max per thread: %s)' %
#                                  (nid, sde_data['Processes'][pid]['Edges'][nid][edge_id]['TgtNode'], ','.join([str(x) for x in sde_data['Processes'][pid]['Edges'][nid][edge_id]['ThreadExecCnts'] if x > 0])))
#                        # exit()
#
    return data


def _id_in_range(interval, testid):
    if interval[0] > 0 and testid < interval[0]:
        return False
    if interval[1] > 0 and testid > interval[1]:
        return False
    return True


def bb_graph(data, block_id_range, fn_name):
    assert(isinstance(data, dict)
           and isinstance(block_id_range, list)
           and isinstance(fn_name, str))

    import networkx as nx
    import plotly.graph_objs as go

    G = nx.DiGraph()

    # visualize basic blocks
    # for pid in data['Processes']:
    #    for src_nid in data['Processes'][pid]['Edges']:
    #        for edge_id in data['Processes'][pid]['Edges'][src_nid]:
    #            tgt_nid = data['Processes'][pid]['Edges'][src_nid][edge_id]['TgtNode']
    #            src_nid, tgt_nid = int(src_nid), int(tgt_nid)
    #            if sum(data['Processes'][pid]['Edges'][src_nid]
    #                   [edge_id]['ThreadExecCnts']) < 1:
    #                continue
    #            if _id_in_range(block_id_range, src_nid) or _id_in_range(
    #                    block_id_range, tgt_nid):
    #                G.add_edge(src_nid, tgt_nid)
    # visualize functions
    for pid in data['Processes']:
        for src_nid in data['Processes'][pid]['Edges']:
            if src_nid < 4:
                continue
            print('src_nid', src_nid)
            src_fid = [fid for fid in data['Processes'][pid]['FileIDs']
                       if src_nid in data['Processes'][pid]['FileIDs'][fid]['Blocks']]
            if len(src_fid) > 0:
                src_fid = src_fid[0]
            else:
                continue
            print('src_fid', src_fid)
            src_sid = data['Processes'][pid]['FileIDs'][src_fid]['Blocks'][src_nid]['SymbolID']
            if src_sid is None:
                continue
            print(src_fid, src_sid)
            for edge_id in data['Processes'][pid]['Edges'][src_nid]:
                print(edge_id, data['Processes'][pid]
                      ['Edges'][src_nid][edge_id])
                tgt_nid = data['Processes'][pid]['Edges'][src_nid][edge_id]['TgtNode']
                if tgt_nid > 3:
                    tgt_fid = [fid for fid in data['Processes'][pid]['FileIDs']
                               if tgt_nid in data['Processes'][pid]['FileIDs'][fid]['Blocks']]
                    if len(tgt_fid) > 0:
                        tgt_fid = tgt_fid[0]
                    else:
                        continue
                    tgt_sid = data['Processes'][pid]['FileIDs'][tgt_fid]['Blocks'][tgt_nid]['SymbolID']
                else:
                    tgt_fid = -1
                    tgt_sid = -1
                print(tgt_nid, tgt_fid, tgt_sid)
                if src_sid == tgt_sid:
                    continue
                if tgt_fid > 0 and tgt_sid is not None:
                    G.add_edge(data['Processes'][pid]['FileIDs'][src_fid]['FNs'][src_sid]['Func'],
                               data['Processes'][pid]['FileIDs'][tgt_fid]['FNs'][tgt_sid]['Func'])
                else:
                    G.add_edge(data['Processes'][pid]['FileIDs']
                               [src_fid]['FNs'][src_sid]['Func'], 'UNKOWN')
                #print(data['Processes'][pid]['FileIDs'][src_fid]['FNs'][src_sid]['Func'], data['Processes'][pid]['FileIDs'][tgt_fid]['FNs'][tgt_sid]['Func'])
                #    print(data['Processes'][pid]['FileIDs'][fid]['Blocks'][src_nid]['SymbolID'])
                #    if data['Processes'][pid]['FileIDs'][fid]['Blocks'][src_nid]['SymbolID'] != src_sid:
                #        continue
                #    for edge_id in data['Processes'][pid]['Edges'][src_nid]:
                #        print(data['Processes'][pid]['Edges'][src_nid][edge_id]['TgtNode'])
                #        tgt_nid = data['Processes'][pid]['Edges'][src_nid][edge_id]['TgtNode']
                #        print(data['Processes'][pid]['FileIDs'][fid]['Blocks'][tgt_nid]['SymbolID'])
                #        if data['Processes'][pid]['FileIDs'][fid]['Blocks'][tgt_nid]['SymbolID'] == src_sid:
                #            continue
                #        G.add_edge(src_sid, data['Processes'][pid]['FileIDs'][fid]['Blocks'][tgt_nid]['SymbolID'])

    pos = nx.kamada_kawai_layout(G)
    # pos = nx.spectral_layout(G)
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


def init_dash_for_vis(data=None):
    assert(isinstance(data, dict))

    import dash
    import dash_core_components as dcc
    import dash_html_components as html

    block_id_range = [-1, -1]
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
                                    id='my-block-from',
                                    type='text',
                                    value='',
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
                                    id='my-block-to',
                                    type='text',
                                    value='',
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
                                    id='my-fn-input',
                                    type='text',
                                    value='',
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
                                        figure=bb_graph(data, block_id_range,
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
        nonlocal data, block_id_range, search_fn_name
        block_id_range = [int(input1), int(input2)]
        search_fn_name = input2
        return bb_graph(data, [int(input1), int(input2)], input3)

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
    args = vars(arg_parser.parse_args())

    assert(args.get('__sde_json_f__') is not None
           and args.get('__sde_block_f__') is not None)
    assert(path.isfile(path.realpath(args.get('__sde_json_f__')))
           and path.isfile(path.realpath(args.get('__sde_block_f__'))))

    sde_data = {'BasicBlocks': {}}
    parse_SDE_JSON(args, sde_data)
    parse_SDE_BB(args, sde_data['BasicBlocks'])

    data = convert_sde_data_to_something_usable(sde_data)
    del sde_data

    if args.get('__visualize__'):
        dash_app = init_dash_for_vis(data)
        dash_app.run_server(debug=False)


if __name__ == '__main__':
    main()
