include ../SLMakes/SLmake.inc.intel

F90 = $(F77) 

PMLIBS = $(PSMMALIB) $(TOOLSLIB) $(BLACSLIB) $(BLASLIB)

#TESTDRV2 = psmma_testcauchy.o
TESTDRV2 = psmma_testtoeplitz.o
#TESTDRV2 = psmma_testhankel.o
#TESTDRV2 = pzfox_fft2.o

all: test2

test2: $(TESTDRV2) $(OBJS) 
	$(F90) $(F77FLAGS) -o $@ $(TESTDRV2) $(OBJS) $(PMLIBS)

%.o:%.f90
	$(F90) -c $(F77FLAGS) $<

%.o:%.f
	$(F90) -c $(F77FLAGS) $<

clean :
	rm -f *.o *.mod
	rm test2
