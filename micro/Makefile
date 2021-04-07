CC=fccpx
CXX=FCCpx
FC=frtpx
COMMON=common/src/report.o
LDLIBS=-Kopenmp -Bstatic

FFLAGS = -Bstatic $(FOPTIMIZE) 
CFLAGS = -Bstatic $(COPTIMIZE)
CXXFLAGS = -Bstatic $(COPTIMIZE)

%.o : %.F90
	$(FC) -c $(FFLAGS) $< -o $@

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

OBJFILE = 

# main: $(OBJFILE) $(COMMON)
#	$(FC) $^ $(LDLIBS) -o $@
