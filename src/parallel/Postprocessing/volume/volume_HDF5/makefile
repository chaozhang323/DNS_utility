DIR=../../..
include $(DIR)/makecomm_2d.inc
COMMDIR=$(DIR)/common

commonobjs=$(COMMDIR)/modFileIO_Compatible.o $(COMMDIR)/modMetrics.o\
           $(COMMDIR)/modDerivative.o $(COMMDIR)/modPRWHDF5.o
srcs= modVolume.f90
objs=$(srcs:.f90=.o) 
exe=volume
.SUFFIXES: .f90 .o

all: comm $(srcs) $(exe)
comm:
	cd $(COMMDIR) && make
.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) $<
$(exe): $(objs) $(commonobjs) volume.o
	$(FC) -o $@ $(LINKFLAG) $(objs) $(commonobjs) volume.o $(LIB) $(LIBLINK)
clean:
	/bin/rm -f $(exe) *.o *.mod *.pc *.il *.d

#depencency base on module usage
volume.o: $(COMMDIR)/modFileIO_Compatible.o $(COMMDIR)/modMetrics.o\
          $(COMMDIR)/modDerivative.o $(COMMDIR)/modPRWHDF5.o modVolume.o

