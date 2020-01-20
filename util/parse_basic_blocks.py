#!/usr/bin/env python3

# usage:
# <script> ./dcfg-out.dcfg.json.bz2 ./dcfg-out.bb.txt.bz2
#	files:
#		./dcfg-out.dcfg.json.bz2 & ./dcfg-out.bb.txt.bz2 :: sde64 -dcfg 1 -dcfg:write_bb 1 ...

from os import path
from sys import argv, exit
from bz2 import open as bz2open
from json import load as jsonload
from re import compile

if not (len(argv) == 3 and path.isfile(path.realpath(argv[1])) and path.isfile(path.realpath(argv[2]))):
	exit("ERROR: usage: <script> ./dcfg-out.dcfg.json.bz2 ./dcfg-out.bb.txt.bz2")

#with jsonload(bz2open(path.realpath(argv[1]), 'rt')) as sde_bb_json:
#	for FILE_NAME_ID, FILE_NAME in sde_bb_json['FILE_NAMES'][1:]:
#		print (FILE_NAME_ID, FILE_NAME)
#
#		# ignore kernel lib
#		if FILE_NAME == '[vdso]':
#			continue
#
#		# get objdump for each file (main binary and all shared libs)
#
with bz2open(path.realpath(argv[2]), 'rt') as sde_bb_txt:
	# block number, program counter, #instruction (inst per block * #exec), #exec of this block, #bytes of block, fn the block belongs to, contained in binary/lib
	block_hdr  = compile('^BLOCK:\s+(\d+)\s+PC:\s+(\w+)\s+ICOUNT:\s+(\d+)\s+EXECUTIONS:\s+(\d+)\s+#BYTES:\s+(\d+)(\s+FN:)?((?(6).*|\s+))IMG:\s+(.+)$')
	# program counter, ?, instruction/assembly
	block_part = compile('^XDIS\s+(\w+):\s+(\w+)\s+(\w+)\s+(.*)$')
	no_block   = compile('^$')
	curr_bb = None
	basic_blocks = {}
	for line in sde_bb_txt:
		if block_hdr.match(line):
			bb = block_hdr.match(line)
			bb_id, bb_pc, bb_ic, bb_ex, bb_by, bb_fi = \
				bb.group(1).strip(), bb.group(2).strip(), \
				bb.group(3).strip(), bb.group(4).strip(), \
				bb.group(5).strip(), bb.group(8).strip()
			bb_fn = bb.group(6)
			if None != bb_fn:
				bb_fn = bb.group(7).strip()
			curr_bb = bb_id
			basic_blocks[curr_bb] = \
				{'BlockID': curr_bb, 'ProgCnt': bb_pc, \
				 'NumExec': bb_ex, 'Func': bb_fn, \
				 'File': bb_fi, 'Inst': []}
		elif block_part.match(line):
			bb = block_part.match(line)
			bb_it = bb.group(2).strip()	# type, like AVX, FMA, etc
			bb_in = bb.group(4).strip()	# assembly as seen by SDE
			basic_blocks[curr_bb]['Inst'].append(bb_in)
		elif no_block.match(line):
			continue
		else:
			exit('ERROR: unknown line: %s' % line)

