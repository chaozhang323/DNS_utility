DIR=../..
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common
#INCLUDE=-I $(COMMDIR)  -I/share/apps/hdf5-mpich2-i12/include -I/home/czb58/tecplot360/include
#LIB=-L/share/apps/hdf5-mpich2-i12/lib -L/share/duan/czb58/dnsmst_utilities/src/fftpack5.1/lib -L/home/czb58/tecplot360/lib 
#LIBLINK=  /share/apps/hdf5-mpich2-i12/lib/libhdf5_fortran.a /home/czb58/tecplot360/lib/libtecio.a -lhdf5   -lz  -lm -lstdc++

commonobjs=$(COMMDIR)/modFileIO_Compatible.o $(COMMDIR)/modMetrics.o\
           $(COMMDIR)/modDerivative.o $(COMMDIR)/modRWHDF5.o\
           $(COMMDIR)/modTecplotIO.o $(COMMDIR)/modTecbin.o
srcs=modVolume.f90  
objs=$(srcs:.f90=.o) 
exe=volume
exe2=volume_HDF5
.SUFFIXES: .f90 .o

all: comm $(srcs) $(exe) $(exe2)
comm:
	cd $(COMMDIR) && make
.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) $<
$(exe): $(objs) volume.o
	$(FC) -o $@ $(LINKFLAG) volume.o $(objs) $(commonobjs) $(LIB) $(LIBLINK)
$(exe2): $(objs) volume_HDF5.o
	$(FC) -o $@ $(LINKFLAG) volume_HDF5.o $(objs) $(commonobjs) $(LIB) $(LIBLINK)
clean:
	/bin/rm -f $(exe) $(exe2) *.o *.mod *.pc *.il *.d

#depencency base on module usage
volume.o: $(COMMDIR)/modFileIO_Compatible.o $(COMMDIR)/modMetrics.o\
          $(COMMDIR)/modDerivative.o $(COMMDIR)/modRWHDF5.o\
          $(COMMDIR)/modTecbin.o $(COMMDIR)/modTecplotIO.o
volume_HDF5.o: $(COMMDIR)/modFileIO_Compatible.o $(COMMDIR)/modMetrics.o\
          $(COMMDIR)/modDerivative.o $(COMMDIR)/modRWHDF5.o\
          $(COMMDIR)/modTecbin.o $(COMMDIR)/modTecplotIO.o
