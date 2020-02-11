DIR=../../..
include $(DIR)/makecomm_2d.inc
COMMDIR=$(DIR)/common

srcs=modParallel.f90 modPHDF5.f90 
objs=$(srcs:.f90=.o) 
Combine=Combine
.SUFFIXES: .f90 .o

all: $(srcs) $(Combine)

.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) $<
$(Combine): $(objs) Combine.o
	$(FC) -o $@ $(LINKFLAGS)  Combine.o $(objs) $(LIB) $(LIBLINK)
clean:
	/bin/rm -f $(Combine) *.o *.mod *.pc *.il *.d

#depencency base on module usage
Combine.o: modParallel.o modPHDF5.o 

