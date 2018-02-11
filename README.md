# Links
- Spreadsheet: https://docs.google.com/spreadsheets/d/1un0TIi31LXI9yURmwobPkCOatNXTvVPofXbEu_W-SvA/

# Settings
```sh
source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64 >/dev/null
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort
ulimit -s unlimited
ulimit -n 4096
```

- `OMP_NUM_THREADS` **must** be several.
- The number of MPI processes **should** be one.


# MEMORY THROUGHPUTS : LYON0
```sh
lyon0 % perf stat -e l2_requests.miss sleep 1 2>&1 >/dev/null | grep l2_requests | tr -d ',' | awk '{printf ("%d B\n", $1 * 64)}'
487488 B
```

- Note: In order to calculate bandwidth, you must divide bytes by elapsed time outputted by `perf`.
- Each tile has LLC(L2) of 1MB.
- List of counters : https://github.com/TomTheBear/perfmondb/blob/master/KNL/KnightsLanding_core_V6.tsv

### MPI
```sh
lyon0% mkdir p
lyon0% mpiexec -n 8 bash -c 'perf stat -e mem_load_uops_retired.l3_miss sleep 1 >/dev/null 2>p/"$MPI_LOCALRANKID".txt'
lyon0% { for i in p/*.txt; do cat $i | egrep 'sec|miss' | tr -d ','  | sed -e 's/\s\+/ /g' | cut -d ' ' -f 2 | tr '\n' ' '; echo; done } | awk '{ s += $1 / $2 } END { printf ("%f GB/sec\n", s * 64 / (1000 ** 3)) }' 
0.000300 GB/sec
```

# MEMORY THROUGHPUTS : mill\[0-1\] (KNM)
```sh
mill0 % perf stat -e cache-misses sleep 1 2>&1 >/dev/null | grep cache-misses | tr -d ',' | awk '{printf ("%d B\n", $1 * 64)}'
366528 B
```

- The `cache-misses` above is the same as the raw counter `r412e` of KNL's `l2_requests.miss`.

# MEMORY THROUGHPUTS :  KIEV0
```sh
kiev0 % perf stat -e mem_load_uops_retired.l3_miss sleep 1 2>&1 >/dev/null | grep mem_load | tr -d ',' | awk '{printf ("%d B\n", $1 * 64)}'
10880 B
```

- Note: In order to calculate bandwidth, you must divide bytes by elapsed time outputted by `perf`.
- LLC is L3 of 30MB. 
- This is following https://github.com/RRZE-HPC/likwid/blob/master/groups/broadwell/L3CACHE.txt

According to ["Detecting Memory-Boundedness with Hardware Performance Counters" Daniel Molka et al.]( http://www.readex.eu/wp-content/uploads/2017/06/ICPE2017_authors_version.pdf ) :

```
Therefore, the proportion of main memory accesses is severely underestimated by the OFFCORE_RESPONSE events.
However, the sum of the L3 hit and L3 miss events is very close to the number of L1 misses in both cases,
so the number of cache line transfers from the uncore to each core can be measured quite accurately.
```

# FLOP
- Note: `-knl` option must be replaced by `-bdw` on KIEV0.
- sde means [Intel SDE](https://software.intel.com/en-us/articles/intel-software-development-emulator)
- _Considering `FMA FLOP` and `MASKED FLOP` is future work._
  - https://software.intel.com/en-us/articles/calculating-flop-using-intel-software-development-emulator-intel-sde

( https://matsulab.slack.com/files/U755Q4FC0/F8XL55459/calculate.py )
```py
#!/bin/python


##################################
#
#	First, get result.txt:
#	mpirun -np 1 bash -c '../../sde-external-8.12.0-2017-10-23-lin/sde64 -knl -iform 1 -omix tmp/"$MPI_LOCALRANKID".txt -- ./exe'
#	for i in tmp/*.txt; do cat $i | egrep '\*total|elements' | sort -t ' ' -k1,1 -k 2rn | uniq -w 22; done >> result.txt
#
#	Second, get time.txt
#	(time mpirun -n 1 ./exe 2>/dev/null;)2>time.txt
#
##################################

f=open("result.txt","r")
f2=open("time.txt","r")
tmp=f.readline()
time=f2.readline()
time=f2.readline()
time=time.split('\t')
timem=time[1].split('m')[0]
times=time[1].split('m')[1]
times=times.split('s')[0]
#print(timem,times)
timem=float(timem)
times=float(times)
time=times+timem*60
element={}
while  tmp!= '':
    key=tmp.split(' ')[0]
    number=tmp.split(' ')[-1]
    number=int(number)
    if element.has_key(key):
        element[key]+=number
    else:
        element.setdefault(key, 0)
        element[key]+=number
    tmp=f.readline()

print('read file success!')

result={'single':0,'double':0,'int':0,'total':0}
for i in element:
    temp=i.split('_')
    if temp[-1] == 'masked':
        temp[-1] = temp[-2]
    if temp[0]=='*total':
        continue
    if temp[1][0]=='i':
        len = int(temp[1][1:]) * int(temp[-1])
        if len <= 64:
            result['int']+=element[i]
        else :
            result['int']+=element[i] * len/64
    else:
        if temp[2]=='single':
            result['single']+=element[i]*int(temp[-1])
        else: 
            result['double']+=element[i]*int(temp[-1])

result['total']=result['single']+result['double']+result['int']
print(result)
print('Percentage of FP64: %2.2f'%(result['double']*100.0/result['total']))
print('Percentage of FP32: %2.2f'%(result['single']*100.0/result['total']))
print('Percentage of INT:  %2.2f'%(result['int']*100.0/result['total']))
print('total time is: %fs'%(time))
print('Performance (sp):        %f GFLOPS' % (result['single']/time/1000/1000/1000))
print('Performance (dp):        %f GFLOPS' % (result['double']/time/1000/1000/1000))
print('TOTAL GFLOPS (sp):       %f GFLOP' % (result['single']/1.0/1000/1000/1000))
print('TOTAL GFLOPS (dp):       %f GFLOP' % (result['double']/1.0/1000/1000/1000))
f.close()
f2.close()
```