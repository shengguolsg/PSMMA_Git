#include ../SLMakes/SLmake.inc
include ../SLMakes/SLmake.inc.intel
#include ../SLMakes/SLmake.inc.th2

F90 = $(F77)

DPOBJ = Cauchylowrank.o aux_hss.o BasicMM.o CauchyHssEig.o \
        Scal_Aux.o Auxil_2p5D.o pscauchy_compute.o psmma_cauchy.o \
	pscauchyeig_compute.o psmma_cauchyeig.o \
        psmma_toeplitz.o pstoeplitz_compute.o psmma_hankel.o pshankel_compute.o	

ZDFT =  zzlascl.o zzlascl2.o pzmdftL.o pzmdftL1.o 

lib : all
all : double

double: $(DPOBJ) $(ZDFT)
	$(ARCH) $(ARCHFLAGS) $(PSMMALIB) $(DPOBJ) $(ZDFT)
	$(RANLIB) $(PSMMALIB)

%.o:%.f90
	$(F90) -c $(F77FLAGS) $<

%.o:%.f
	$(F90) -c $(F77FLAGS) $<


clean :
	rm -f $(DPOBJ) $(ZDFT) *.mod

veryclean:
	rm -f ../lib/libpsmma.a
