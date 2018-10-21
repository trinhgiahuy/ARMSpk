# Proxy-App Study

## Cloning the repo
```
git clone --recurse-submodules https://gitlab.com/domke/PAstudy.git
```


## Preparation of the environment
```
cd ./PAstudy/
vim ./conf/intel.cfg	# configure correct path to Intel's Parallel Studio XE
vim ./conf/host.cfg	# specify hostnames and CPU freqs (replace kiev,etc.)
./inst/_init.sh		# follow instructions if printed (eps. freq. settings)
```
(Note: root access and setcap will be required; setcap usually does not work
on remote file systems, so the repo should be cloned onto a local disk)


## Examplified installation of a proxy-app
```
./inst/amg.sh
```


## Examplified modification of proxy-app configurations
* allows to specify binary names, changes in input settings, number of benchmark
  trials, maximum running time of the benchmark, which performance tools is
  executed, and MPI/OMP sweep configurations per host system, etc
```
vim ./conf/amg.sh
```


## Running one of the proxy-apps (e.g. AMG)
```
./run/amg/test.sh	# run the MPI/OMP sweep test to find best config
vim ./conf/amg.sh	# replace BESTCONF with config found in last step
./run/amg/best.sh	# execute the presumably best MPI/OMP config many times
./run/amg/freq.sh	# run freqency scaling tests
./run/amg/prof.sh	# use performance tools to extract raw metrics
```
* all logs (results+raw data) will be placed in ./log/\<hostname\>/\<test\>/\<proxy\>


## Author's previously colleted results for reference
* untar the files with our raw data
```
cd ./PAstudy/
tar xjf paper/data/kiev3.tar.bz2	# Broadwell
tar xjf paper/data/lyon1.tar.bz2	# KNL
tar xjf paper/data/mill1.tar.bz2	# KNM
```
* afterwards, raw data will be in ./log/\<hostname\>


## Additional notes:
* this repo does not include 3rd party and proprietary libraries
* if a dependency is not met, the scripts usually point it out and ask the
  user to install the missing lib/tool/etc.
* the license (see ./LICENSE) applies ONLY to all source code, logs, etc.
  written/produced by the authors of this repo (for everything else: please,
  refer to individual licenses of the submodules (proxy-apps) and other
  libraries/dependencies after pulling those)
