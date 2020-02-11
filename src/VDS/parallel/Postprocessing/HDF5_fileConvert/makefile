DIR=../..
2decomp_dir=../../../../2decomp_fft
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common

commonobjs=$(COMMDIR)/modPRWHDF5.o 
srcs=
objs=$(srcs:.f90=.o)
exe=HDF5_fileConvert
.SUFFIXES: .f90 .o

all: 2decomp comm $(srcs) $(exe)
2decomp:
	cd $(2decomp_dir)/src && make
comm:
	cd $(COMMDIR) && make

.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) -I$(COMMDIR) $<
$(exe): $(objs) $(commonobjs) HDF5_fileConvert.o
	$(FC) -o $@ $(LINKFLAG) $(objs) $(commonobjs) HDF5_fileConvert.o $(LIB) $(LIBLINK)
clean:
	/bin/rm -f $(exe) *.o *.mod *.pc *.il *.d

#depencency base on module usage
HDF5_fileConvert.o:$(COMMDIR)/modPRWHDF5.o   

