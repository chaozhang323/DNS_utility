DIR=../..
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common
INCLUDE=-I $(COMMDIR)  -I/share/apps/hdf5-mpich2-i12/include -I/home/czb58/tecplot360/include
LIB=-L/share/apps/hdf5-mpich2-i12/lib -L/share/duan/czb58/dnsmst_utilities/src/fftpack5.1/lib -L/home/czb58/tecplot360/lib -lfftpack
LIBLINK=  /share/apps/hdf5-mpich2-i12/lib/libhdf5_fortran.a /home/czb58/tecplot360/lib/libtecio.a  -lhdf5   -lz  -lm -lstdc++

commonobjs=$(COMMDIR)/modFileIO_Compatible.o $(COMMDIR)/modMetrics.o\
           $(COMMDIR)/modDerivative.o $(COMMDIR)/modTecplotIO.o
srcs=modVolume.f90   modRWHDF5.f90 modTecbin.f90
objs=$(srcs:.f90=.o) 
exe=volume
tsjplane=tsjplane
.SUFFIXES: .f90 .o

all: comm $(srcs) $(exe) $(tsjplane) 
comm:
	cd $(COMMDIR) && make
.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) $<
$(exe): $(objs) $(commonobjs) volume.o
	$(FC) -o $@ $(LINKFLAG) volume.o $(objs) $(commonobjs) $(LIB) $(LIBLINK)
$(tsjplane): $(objs) $(commonobjs) tsjplane.o
	$(FC) -o $@ $(LINKFLAG) tsjplane.o  $(objs) $(commonobjs) $(LIB) $(LIBLINK)
clean:
	/bin/rm -f $(exe) $(tsjplane) *.o *.mod *.pc *.il *.d

#depencency base on module usage
volume.o: $(COMMDIR)/modFileIO_Compatible.o $(COMMDIR)/modMetrics.o\
          $(COMMDIR)/modDerivative.o modRWHDF5.o
tsjplane.o: $(COMMDIR)/modFileIO_Compatible.o $(COMMDIR)/modMetrics.o\
          $(COMMDIR)/modDerivative.o modRWHDF5.o modTecbin.o $(COMMDIR)/modTecplotIO.o
