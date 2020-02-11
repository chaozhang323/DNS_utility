DIR=../..
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common

srcs=modReadBin.f90
objs=$(srcs:.f90=.o) 
fortran2hdf5=bin2hdf5
.SUFFIXES: .f90 .o

all: $(srcs) $(fortran2hdf5)

.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) $<
$(fortran2hdf5): $(objs) bin2hdf5.o
	$(FC) -o $@ $(LINKFLAGS)  $(objs) bin2hdf5.o $(LIB) $(LIBLINK)
clean:
	/bin/rm -f $(fortran2hdf5) *.o *.mod *.pc *.il *.d

#depencency base on module usage
bin2hdf5.o: modReadBin.o  

