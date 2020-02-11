DIR=../..
2decomp_dir=../../../2decomp_fft
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common

commonobjs=$(COMMDIR)/modFileIO_Compatible.o $(COMMDIR)/modPRWHDF5.o $(COMMDIR)/interp1d.o $(COMMDIR)/tridiagonal.o  
srcs=modFFTW.f90  
objs=$(srcs:.f90=.o)
exe=FFT_4D
.SUFFIXES: .f90 .o

all: 2decomp comm $(srcs) $(exe)
2decomp:
	cd $(2decomp_dir)/src && make
comm:
	cd $(COMMDIR) && make

.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) -I$(COMMDIR) $<
$(exe): $(objs) $(commonobjs)  FFT_4D.o
	$(FC) -o $@ $(LINKFLAG) $(objs) $(commonobjs) FFT_4D.o $(LIB) $(LIBLINK)
clean:
	/bin/rm -f $(exe) *.o *.mod *.pc *.il *.d

#depencency base on module usage
FFT_4D.o: $(COMMDIR)/modFileIO_Compatible.o $(COMMDIR)/modPRWHDF5.o modFFTW.o $(COMMDIR)/interp1d.o $(COMMDIR)/tridiagonal.o  

