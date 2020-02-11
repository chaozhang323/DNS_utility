DIR=../..
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common

commonobjs=$(COMMDIR)/modFFT.o $(COMMDIR)/fftutil.o \
           $(COMMDIR)/modRWHDF5.o $(COMMDIR)/modFileIO_Compatible.o \
           $(COMMDIR)/interp1d.o $(COMMDIR)/tridiagonal.o\
           $(COMMDIR)/modMetrics.o $(COMMDIR)/modDerivative.o \
           $(COMMDIR)/modTecplotIO.o $(COMMDIR)/modTecbin.o

srcs= modHelmholtzDecomp.f90
objs=$(srcs:.f90=.o) 
planeinlet=Planeinlet
.SUFFIXES: .f90 .o

all: comm $(srcs) $(planeinlet)
comm:
	cd $(COMMDIR) && make
.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) $<
$(planeinlet): $(objs) Planeinlet.o
	$(FC) -o $@ $(LINKFLAG) Planeinlet.o $(objs) $(commonobjs) $(LIB) $(LIBLINK)

clean:
	/bin/rm -f $(planeinlet) *.o *.mod *.pc *.il *.d

#depencency base on module usage
Planeinlet.o: $(COMMDIR)/modRWHDF5.o $(COMMDIR)/interp1d.o $(COMMDIR)/tridiagonal.o\
              $(COMMDIR)/modMetrics.o $(COMMDIR)/modDerivative.o \
              $(COMMDIR)/modTecplotIO.o $(COMMDIR)/modTecbin.o modHelmholtzDecomp.o 

