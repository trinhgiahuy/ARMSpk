m.PHONY : libssc.so testC.exe testF.exe

testC.exe : test.c
	mpicc -O3 -I${ADVISOR_2018_DIR}/include -fopenmp -o testC.exe test.c -L${ADVISOR_2018_DIR}/lib64 -littnotify

testF.exe : test.f90
	icc -c -I${ADVISOR_2018_DIR}/include ssc.c -o ssc.o
	mpif90 -O3 -fopenmp -o testF.exe test.f90 ssc.o -L${ADVISOR_2018_DIR}/lib64 -littnotify

libssc.so : ssc.c
	icc -c -fPIC -I${ADVISOR_2018_DIR}/include ssc.c -o ssc.o
	icc -shared -Wl,-soname,libssc.so -o libssc.so ssc.o ${ADVISOR_2018_DIR}/lib64/libittnotify.a

clean :
	-rm -vf testC.exe testF.exe libssc.so ssc.o ssc.pyc

run : libssc.so testC.exe testF.exe
	mpirun -n 2 ./testC.exe
	mpirun -n 2 ./testF.exe
	python ./test.py
