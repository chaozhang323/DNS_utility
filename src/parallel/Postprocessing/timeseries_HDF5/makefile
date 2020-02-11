DIR=../..
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common

commonobjs=$(COMMDIR)/modFFT.o $(COMMDIR)/fftutil.o $(COMMDIR)/interp1d.o\
           $(COMMDIR)/tridiagonal.o $(COMMDIR)/modMetrics.o $(COMMDIR)/modDerivative.o

srcs= modTSData_HDF5.f90 modSpect.f90 modSpecty.f90
objs=$(srcs:.f90=.o) 
spect_HDF5=spectcorr_HDF5
crosscorr=Calcrosscorr_HDF5



.SUFFIXES: .f90 .o

all:  comm $(srcs) $(spect_HDF5) $(crosscorr)
comm:
	cd $(COMMDIR) && make
.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) $<
$(spect_HDF5): $(objs) $(commonobjs) spectcorr_HDF5.o
	$(FC) -o $@ $(LINKFLAG) spectcorr_HDF5.o $(objs) $(commonobjs) $(LIB) $(LIBLINK)
	
$(crosscorr): $(objs) $(commonobjs) Calcrosscorr_HDF5.o
	$(FC) -o $@ $(LINKFLAG) Calcrosscorr_HDF5.o $(objs) $(commonobjs) $(LIB) $(LIBLINK)

clean:
	/bin/rm -f $(spect_HDF5) $(crosscorr) *.o *.mod *.pc *.il *.d

#depencency base on module usage
spectcorr_HDF5.o: $(COMMDIR)/modFFT.o modTSData_HDF5.o modSpect.o modSpecty.o
Calcrosscorr_HDF5.o: $(COMMDIR)/modFFT.o modTSData_HDF5.o modSpect.o 
