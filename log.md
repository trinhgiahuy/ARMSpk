# LYON0MEMORY_THROUGHPUTS
( https://matsulab.slack.com/archives/C8NENQMC3/p1516867766000196 )

For the time being, let's calculate throughputs on lyon0 from (l2_requests.miss * cache-line-length[64 bytes] / seconds). This throughputs are nearly in proportion to `offcore_response.any_request.ddr` [unit unknown]. Other counters which seem good (e.g. `OFFCORE_RESPONSE.ANY_DATA_RD.DDR`, `OFFCORE_RESPONSE.ANY_DATA_RD.MCDRAM`) are not supported actually.

```
lyon0 % perf stat -e l2_requests.miss -a ls >/dev/null
```

# KIEV0 MEMORY_THROUGHPUTS
( https://matsulab.slack.com/archives/C8NENQMC3/p1516793575000060 )

This is for Xeon.

```
keiv0% perf stat -e uncore_imc_0/cas_count_write/,uncore_imc_0/cas_count_read/,uncore_imc_1/cas_count_write/,uncore_imc_1/cas_count_read/,uncore_imc_4/cas_count_write/,uncore_imc_4/cas_count_read/,uncore_imc_5/cas_count_write/,uncore_imc_5/cas_count_read/ -a sleep 1 2>&1 >/dev/null  | sed -e 's/MiB/@\n/' | grep @ | awk '{s += $1} END {printf ("%.2f MB\n", s * 1.024 * 1.024)}'
53.02 MB
```


# FLOP
( https://matsulab.slack.com/files/U755Q4FC0/F8XL55459/calculate.py )
```
#!/bin/python


##################################
#
#	First, get result.txt:
#	mpirun -np 1 bash -c 
'../../sde-external-8.12.0-2017-10-23-lin/sde64 -knl -iform 1 -omix 
tmp/"$MPI_LOCALRANKID".txt -- ./exe'
#	for i in tmp/*.txt; do cat $i | egrep '\*total|elements' | sort 
-t ' ' -k1,1 -k 2rn | uniq -w 22; done >> result.txt
#
#	Second, get time.txt
#	(time mpirun -n 1 ./exe;)2>time.txt
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
	if temp[0]=='*total':
		continue
	if temp[1][0]=='i':
		result['int']+=element[i]*int(temp[-1])
	else:
		if temp[2]=='single':
			result['single']+=element[i]*int(temp[-1])
		else: 
			result['double']+=element[i]*int(temp[-1])

result['total']=result['single']+result['double']+result['int']
print(result)
print('Percentage of FP64: 
%2.2f'%(result['double']*100.0/result['total']))
print('Percentage of FP32: 
%2.2f'%(result['single']*100.0/result['total']))
print('Percentage of INT:  %2.2f'%(result['int']*100.0/result['total']))
print('total time is: %fs'%(time))
print('Performance (sp):	
%f'%(result['single']/time/1024/1024/1024))
print('Performance (dp):	
%f'%(result['double']/time/1024/1024/1024))
```
