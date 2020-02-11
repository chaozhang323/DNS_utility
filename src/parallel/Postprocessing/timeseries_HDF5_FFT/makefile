DIR=../..
2decomp_dir=../../../2decomp_fft
include ./makecomm.inc
COMMDIR=$(DIR)/common

commonobjs= $(COMMDIR)/modPRWHDF5.o $(COMMDIR)/modFFT.o $(COMMDIR)/fftutil.o

srcs=modSpect.f90
objs=$(srcs:.f90=.o) 
coh=Calcoh_HDF5
.SUFFIXES: .f90 .o

all: 2decomp comm $(srcs) $(coh) 
2decomp:
	cd $(2decomp_dir)/src && make
comm:
	cd $(COMMDIR) && make
.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) -I$(COMMDIR) $<
$(coh): $(objs) $(commonobjs) Calcoh_HDF5.o
	$(FC) -o $@ $(LINKFLAG) $(objs) $(commonobjs) Calcoh_HDF5.o $(LIB) $(LIBLINK)

clean:
	/bin/rm -f $(coh) *.o *.mod *.pc *.il *.d

#depencency base on module usage
Calcoh_HDF5.o:$(COMMDIR)/modPRWHDF5.o $(COMMDIR)/modFFT.o $(COMMDIR)/fftutil.o modSpect.o


