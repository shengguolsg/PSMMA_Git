F77=mpif90
#-fc=ifort
#F77FLAGS=-g -fopenmp -fcheck=all,no-array-temps
#F77FLAGS=-O3 -qopenmp -heap-arrays -fp-model precise -rdynamic -mavx
#F77FLAGS=-O3 -heap-arrays -fp-model precise -rdynamic -mavx
#F77FLAGS=-O3 -qopenmp -heap-arrays -fp-model precise -rdynamic -CB -CU -mavx
#F77FLAGS=-O0 -g -heap-arrays -fp-model precise -rdynamic -traceback
#F77FLAGS=-O0 -g -heap-arrays -fp-model precise -rdynamic -traceback -CU -CB
#F77FLAGS=-O0 -g -rdynamic -funwind-tables -fcheck=all -fno-backtrace
F77FLAGS=-O0 -g -fbounds-check -Wall -fbacktrace

F77LOADER=mpif90
#-fc=ifort
#F77LOADFLAGS=-O2 -fopenmp
#F77FLAGS=-O3 -qopenmp -heap-arrays -fp-model precise -rdynamic -mavx
#F77FLAGS=-O3 -heap-arrays -fp-model precise -rdynamic -mavx
#F77FLAGS=-O3 -qopenmp -heap-arrays -fp-model precise -rdynamic -CB -CU -mavx
#F77FLAGS=-O0 -g -heap-arrays -fp-model precise -rdynamic -traceback -CU -CB
# -traceback
#F77LOADFLAGS= -g -O0  -rdynamic -funwind-tables -fcheck=all -traceback
#-fno-backtrace
F77LOADFLAGS=-O2 -g -fbounds-check -Wall -fbacktrace

ARCH      = ar
ARCHFLAGS = cr
RANLIB    = ranlib
#CPPFLAG   = -fpp
MACRO_D_PREFIX = -D
INCFLAG   = -I./
#MODFLAG   = -module ./

INSTALLROOT=/home/lsg/MyCodes/THDLA/lib
PSMMALIB=$(INSTALLROOT)/libpsmma.a

#PUMMALIB=/home/lsg/MyProject/PUMMA_OLD_Correct/SPUMMA_Release/libpumma.a

BLASLIB= -L/home/lsg/local/Openblas/lib/ -lscalapack -lopenblas -lgfortran
#BLASLIB= -L/home/lsg/local/lib/ -lscalapack -llapack -lblas

#BLASLIB=-Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm
#BLASLIB=-Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a \
#       $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm
#BLACSLIB=$(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a

#BLASLIB=$(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a \
#        $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm -lgfortran

#BLASLIB=-mkl=parallel
#BLACSLIB=-lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lmkl_blacs_intelmpi_lp64 -lpthread -lm
