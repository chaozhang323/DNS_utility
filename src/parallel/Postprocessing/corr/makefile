DIR=../..
2decomp_dir=../../../2decomp_fft
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common
#decompdir=$(2decomp_dir)/src

commonobjs=$(COMMDIR)/modFileIO_Compatible.o $(COMMDIR)/modPRWHDF5.o $(COMMDIR)/interp1d.o $(COMMDIR)/tridiagonal.o
srcs=
objs=$(srcs:.f90=.o)
exe=correlation
.SUFFIXES: .f90 .o

all: 2decomp comm $(srcs) $(exe)
2decomp:
	cd $(2decomp_dir)/src && make
comm:
	cd $(COMMDIR) && make

.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) -I$(COMMDIR) $<
$(exe): $(objs) $(commonobjs) modFFTW.o correlation.o
	$(FC) -o $@ $(LINKFLAG) $(objs) $(commonobjs) modFFTW.o correlation.o $(LIB) $(LIBLINK)
clean:
	/bin/rm -f $(exe) *.o *.mod *.pc *.il *.d

#depencency base on module usage
correlation.o: $(COMMDIR)/modFileIO_Compatible.o $(COMMDIR)/modPRWHDF5.o modFFTW.o $(COMMDIR)/interp1d.o $(COMMDIR)/tridiagonal.o

