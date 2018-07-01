DIR=..
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common

#INCLUDE=-I/share/apps/hdf5-mpich2-i12/include 
#LIB=-L/share/apps/hdf5-mpich2-i12/lib 
#LIBLINK= /share/apps/hdf5-mpich2-i12/lib/libhdf5_fortran.a -lm -lstdc++ -lhdf5   -lz 

#LIB=-L/nasa/ae/tecplot360/lib -lstdc++ -L/nasa/ae/tecplot360/lib/tecio64.a
srcs=modFileIO_Compatible.f90 modTecplotIO.f90 modFFT.f90 modBuffer.f90 modMetrics.f90 modDerivative.f90\
 	interp1d.f90 spectprocess.f90 fftutil.f90 nonlinearsolver.f90 random.f90 tridiagonal.f90\
	modAverage.f90 modGrid1d.f90 modReynoldsBudgets.f90 modTKEBudgets.f90 modRWHDF5.f90 modTecbin.f90 modFlow.f90\
	modNumRecipe.f90 gridgen.f90 orthogonal.f90 unstructured.f90 interp2d.f90
srcs77=rtbis.for
objs=$(srcs:.f90=.o)
objs77=$(srcs77:.for=.o)
.SUFFIXES: .f90 .o
.SUFFIXES: .for .o

all: $(objs77) $(objs)
.f90.o:
#	$(FC) -c $(FFLAG) $< $(LIB)
	$(FC) -c $(FFLAG) $(INCLUDE) $< 
.for.o:
	$(FC) -c $(FFLAG) $(INCLUDE) $<
clean:
	/bin/rm -f *.o *.mod *.pc *.il *.d

#depencency base on module usage
spectprocess.o: modFFT.o
modTecbin.o: modTecplotIO.o
modNumRecipe.o: rtbis.o
orthogonal.o: gridgen.o interp1d.o modMetrics.o modTecbin.o
gridgen.o: tridiagonal.o modMetrics.o
