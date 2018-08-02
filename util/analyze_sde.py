#!/usr/bin/env python

from os import walk, path, chdir, getcwd, linesep
from sys import argv, exit
from re import compile

if len(argv) == 3 and (path.isdir(path.realpath(argv[1])) or path.isfile(path.realpath(argv[1]))) and path.isfile(path.realpath(argv[2])):
	if path.isdir(path.realpath(argv[1])):
		must_match = None
		sdeout_dir = path.realpath(argv[1])
	elif path.isfile(path.realpath(argv[1])):
		must_match = path.basename(path.realpath(argv[1]))
		sdeout_dir = path.dirname(path.realpath(argv[1]))
	bestbm_log = path.realpath(argv[2])
else:
	exit("ERROR: Incorrect input directory!" + linesep + linesep + "Usage: %s <folder-with-SDE-output> <log-of-best-run>" % __file__)

chdir(path.dirname(path.realpath(__file__)))
currdir = getcwd()

memr_re = compile('^\*mem-read-(\d+)\s+(\d+)')
memw_re = compile('^\*mem-write-(\d+)\s+(\d+)')
ixfr_re = compile('^\*dataxfer_i(\d+)_(\d+)\s+(\d+)')
fxfr_re = compile('^\*dataxfer_fp_(\w+)_(\d+)\s+(\d+)')
real_re = compile('^\*elements_fp_single_(\d+)\s+(\d+)')
dble_re = compile('^\*elements_fp_double_(\d+)\s+(\d+)')
inte_re = compile('^\*elements_i(\d+)_(\d+)\s+(\d+)')
bran_re = compile('^\*cond-branch-\w+\s+(\d+)')
vfms_re = compile('^VFM.*(SS|SD|PS).*(XMM|YMM|ZMM).*\s+(\d+)')
mult = {'XMM': 2, 'YMM': 4, 'ZMM':8}

mread_in_byte = 0
mwrite_in_byte = 0
num_dataxfer_int = 0
num_dataxfer_flo = 0
num_ops_real = 0
num_ops_dble = 0
num_ops_inte = 0
num_branches = 0

# ignore subfolders
for _, _, files in walk(sdeout_dir):
	for fname in files:
		if must_match and must_match != fname:
			continue
		with open(path.join(sdeout_dir, fname), 'r') as sdeout:
			for line in sdeout:
				if memr_re.match(line):
					m = memr_re.match(line)
					mread_in_byte += int(m.group(1)) * int(m.group(2))
					continue
				if memw_re.match(line):
					m = memw_re.match(line)
					mwrite_in_byte += int(m.group(1)) * int(m.group(2))
					continue
				if ixfr_re.match(line):
					m = ixfr_re.match(line)
					num_dataxfer_int += int(m.group(3))
					continue
				if fxfr_re.match(line):
					m = fxfr_re.match(line)
					num_dataxfer_flo += int(m.group(3))
					continue
				if real_re.match(line):
					m = real_re.match(line)
					num_ops_real += int(m.group(1)) * int(m.group(2))
					continue
				if dble_re.match(line):
					m = dble_re.match(line)
					num_ops_dble += int(m.group(1)) * int(m.group(2))
					continue
				if inte_re.match(line):
					m = inte_re.match(line)
					num_ops_inte += int(m.group(2)) * int(m.group(3))
					continue
				if bran_re.match(line):
					m = bran_re.match(line)
					num_branches += int(m.group(1))
					continue
				if vfms_re.match(line):
					m = vfms_re.match(line)
					if 'SS' in m.group(1):
						num_ops_real += int(m.group(3))
					elif 'SD' in m.group(1):
						num_ops_dble += int(m.group(3))
					elif 'PS' in m.group(1):
						num_ops_real += 2 * mult[m.group(2)] * int(m.group(3))
					elif 'PD' in m.group(1):
						num_ops_dble += mult[m.group(2)] * int(m.group(3))
	break

total_rtime_re = compile('^Total running time:\s+([-+]?\d*\.\d+|\d+|-)')
kernel_rtime_re = compile('^Walltime of the main kernel:\s+([-+]?\d*\.\d+|\d+|-)')

total_rtime, kernel_rtime = 2**64, 2**64
with open(bestbm_log, 'r') as bmlog:
	for line in bmlog:
		if total_rtime_re.match(line):
			m = total_rtime_re.match(line)
			new_total_rtime = float(m.group(1))
			if new_total_rtime < total_rtime:
				total_rtime = new_total_rtime
			continue
		if kernel_rtime_re.match(line):
			m = kernel_rtime_re.match(line)
			new_kernel_rtime = float(m.group(1))
			if new_kernel_rtime < kernel_rtime:
				kernel_rtime = new_kernel_rtime
			continue

if kernel_rtime == 2**64:
	exit("ERROR: could not find walltime of kernel in %s" % bestbm_log)

print("Total BM runtime [in s]:", total_rtime)
print("Walltime of kernel [in s]:", kernel_rtime)

print("GBytes read:", mread_in_byte/(1000.0*1000*1000))
print("GBytes written:", mwrite_in_byte/(1000.0*1000*1000))
print("Number INT data transfers:", num_dataxfer_int)
print("Number FLOAT data transfers:", num_dataxfer_flo)
print("Total #(DP)GFLOP:", num_ops_dble/(1000.0*1000*1000))
print("Total #(SP)GFLOP:", num_ops_real/(1000.0*1000*1000))
print("Total #GINTOP:", num_ops_inte/(1000.0*1000*1000))
print("Total #Gbranches:", num_branches/(1000.0*1000*1000))

print("GFLOP/s (SP):", num_ops_real/(1000.0*1000*1000*kernel_rtime))
print("GFLOP/s (DP):", num_ops_dble/(1000.0*1000*1000*kernel_rtime))
print("GIOP/s:", num_ops_inte/(1000.0*1000*1000*kernel_rtime))
