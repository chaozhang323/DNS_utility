include ../makecomm.inc

srcs=modFileIO_Compatible.f90 modTecplotIO.f90 modFFT.f90 modBuffer.f90 modMetrics.f90 modDerivative.f90\
 	interp1d.f90 spectprocess.f90 fftutil.f90 nonlinearsolver.f90 random.f90 tridiagonal.f90\
	modAverage.f90 modGrid1d.f90 modReynoldsBudgets.f90 modTKEBudgets.f90 modRWHDF5.f90 modTecbin.f90\
        modPRWHDF5.f90
objs=$(srcs:.f90=.o)
.SUFFIXES: .f90 .o

all:$(srcs) $(objs)
.f90.o:
#	$(FC) -c $(FFLAG) $< $(LIB)
	$(FC) -c $(FFLAG) $(INCLUDE) $< 
clean:
	/bin/rm -f *.o *.mod *.pc *.il *.d

#depencency base on module usage
spectprocess.o: modFFT.o
