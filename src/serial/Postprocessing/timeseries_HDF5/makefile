DIR=../..
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common

commonobjs=$(COMMDIR)/modFFT.o $(COMMDIR)/fftutil.o $(COMMDIR)/interp1d.o\
           $(COMMDIR)/tridiagonal.o $(COMMDIR)/modMetrics.o $(COMMDIR)/modDerivative.o $(COMMDIR)/modTecplotIO.o $(COMMDIR)/modTecbin.o \
            $(COMMDIR)/modSystemTime.o

srcs= modTSData_HDF5.f90 modStat_kplane.f90 modStat_iplane.f90 modStat_jplane.f90 modAcousticSource.f90 modAngle.f90 modSpect.f90 modSpecty.f90
objs=$(srcs:.f90=.o) 
Calstat_HDF5=Calstat_HDF5
CalAngle=CalAngle
spect_HDF5=spectcorr_HDF5
crosscorr=Calcrosscorr_HDF5
coh=Calcoh_HDF5
intpcorr=Interpcorr
coh_fft3d=Calcoh_HDF5_FFT3D
.SUFFIXES: .f90 .o

all:  comm $(srcs) $(Calstat_HDF5) $(CalAngle) $(spect_HDF5) $(crosscorr) $(coh) $(intpcorr) $(coh_fft3d)
comm:
	cd $(COMMDIR) && make
.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) $<
$(Calstat_HDF5): $(objs) Calstat_HDF5.o
	$(FC) -o $@ $(LINKFLAG) Calstat_HDF5.o $(objs) $(commonobjs) $(LIB) $(LIBLINK)
$(CalAngle): $(objs) CalAngle.o
	$(FC) -o $@ $(LINKFLAG) CalAngle.o $(objs) $(commonobjs) $(LIB) $(LIBLINK)	
$(spect_HDF5): $(objs) $(commonobjs) spectcorr_HDF5.o
	$(FC) -o $@ $(LINKFLAG) spectcorr_HDF5.o $(objs) $(commonobjs) $(LIB) $(LIBLINK)
$(crosscorr): $(objs) $(commonobjs) Calcrosscorr_HDF5.o
	$(FC) -o $@ $(LINKFLAG) Calcrosscorr_HDF5.o $(objs) $(commonobjs) $(LIB) $(LIBLINK)
$(coh): $(objs) $(commonobjs) Calcoh_HDF5.o
	$(FC) -o $@ $(LINKFLAG) Calcoh_HDF5.o $(objs) $(commonobjs) $(LIB) $(LIBLINK)
$(intpcorr): $(objs) $(commonobjs) Interpcorr.o
	$(FC) -o $@ $(LINKFLAG) Interpcorr.o $(objs) $(commonobjs) $(LIB) $(LIBLINK)
$(coh_fft3d): $(objs) $(commonobjs) Calcoh_HDF5_FFT3D.o
	$(FC) -o $@ $(LINKFLAG) Calcoh_HDF5_FFT3D.o $(objs) $(commonobjs) $(LIB) $(LIBLINK)

clean:
	/bin/rm -f $(Calstat_HDF5) $(CalAngle) $(spect_HDF5) $(crosscorr) $(coh) $(intpcorr) $(coh_fft3d) *.o *.mod *.pc *.il *.d

#depencency base on module usage
Calstat_HDF5.o: modStat_kplane.o modStat_iplane.o modStat_jplane.o modAcousticSource.o modTSData_HDF5.o
CalAngle.o: modAngle.o modTSData_HDF5.o
spectcorr_HDF5.o: $(COMMDIR)/modFFT.o modTSData_HDF5.o modSpect.o modSpecty.o
Calcrosscorr_HDF5.o: $(COMMDIR)/modFFT.o modTSData_HDF5.o modSpect.o 
Calcoh_HDF5.o: $(COMMDIR)/modFFT.o modTSData_HDF5.o modSpect.o 
Interpcorr.o: $(COMMDIR)/interp1d.o $(COMMDIR)/tridiagonal.o
Calcoh_HDF5_FFT3D.o: $(commonobjs) modTSData_HDF5.o modSpect.o 

