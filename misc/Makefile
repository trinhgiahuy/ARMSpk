m.PHONY : libssc.so testC.exe testF.exe

testC.exe : test.c
	mpicc -O3 -fopenmp -o testC.exe test.c

testF.exe : test.f90
	icc -c ssc.c -o ssc.o
	mpif90 -O3 -fopenmp -o testF.exe test.f90 ssc.o

libssc.so : ssc.c
	icc -c -fPIC ssc.c -o ssc.o
	icc -shared -Wl,-soname,libssc.so -o libssc.so ssc.o

clean :
	-rm -vf testC.exe testF.exe libssc.so ssc.o ssc.pyc

run : libssc.so testC.exe testF.exe
	mpirun -n 2 ./testC.exe
	mpirun -n 2 ./testF.exe
	python ./test.py
