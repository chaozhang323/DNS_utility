DIR=../..
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common

commonobjs=$(COMMDIR)/modBuffer.o $(COMMDIR)/interp1d.o $(COMMDIR)/tridiagonal.o\
           $(COMMDIR)/modMetrics.o $(COMMDIR)/modDerivative.o $(COMMDIR)/modPRWHDF5.o
srcs=InitGridFlow_HDF5.f90
objs=$(srcs:.f90=.o) 
exe=InitGridFlow_HDF5
.SUFFIXES: .f90 .o

all: comm $(srcs) $(exe)
comm:
	cd $(COMMDIR) && make
.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) -I$(COMMDIR) $<
$(exe): $(objs) $(commonobjs)
	$(FC) -o $@ $(LINKFLAG) $(objs) $(commonobjs) $(LIB) $(LIBLINK)
clean:
	/bin/rm -f $(exe) *.o *.mod *.pc *.il *.d

#depencency base on module usage
InitGridFlow_HDF5.o: $(COMMDIR)/modBuffer.o $(COMMDIR)/interp1d.o $(COMMDIR)/tridiagonal.o\
                     $(COMMDIR)/modDerivative.o $(COMMDIR)/modPRWHDF5.o 
