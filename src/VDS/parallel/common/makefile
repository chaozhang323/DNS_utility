include ../makecomm.inc

srcs=modPRWHDF5.f90
objs=$(srcs:.f90=.o)
.SUFFIXES: .f90 .o

all:$(srcs) $(objs)
.f90.o:
#	$(FC) -c $(FFLAG) $< $(LIB)
	$(FC) -c $(FFLAG) $(INCLUDE) $< 
clean:
	/bin/rm -f *.o *.mod *.pc *.il *.d

#depencency base on module usage
 
