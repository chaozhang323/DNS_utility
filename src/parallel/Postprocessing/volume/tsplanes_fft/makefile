DIR=../../..
include $(DIR)/makecomm_2d.inc
COMMDIR=$(DIR)/common
INCLUDE=-I/share/duan/czb58/Package/2decomp_fft/include -I$(COMMDIR)  -I/share/duan/czb58/Package/HDF5_Szip/include -I/home/czb58/tecplot360/include -I/share/apps/fftw-3.3.4/include
LIB=-L/share/duan/czb58/Package/2decomp_fft/lib -L/share/duan/czb58/Package/HDF5_Szip/lib -L/home/czb58/tecplot360/lib 
LIBLINK= /share/duan/czb58/Package/2decomp_fft/lib/lib2decomp_fft.a /share/duan/czb58/Package/HDF5_Szip/lib/libhdf5_fortran.a /home/czb58/tecplot360/lib/libtecio.a /share/apps/fftw-3.3.4/lib/libfftw3.a -lfftw3f -lfftw3 -lhdf5   -lz  -lm -lstdc++

commonobjs= 
#$(COMMDIR)/modMetrics.o\
#$(COMMDIR)/modDerivative.o  
 
srcs=modMetrics.f90 modDerivative.f90  modRWPHDF5.f90 modVolume.f90    
objs=$(srcs:.f90=.o) 
tsplanes=tsplanes
.SUFFIXES: .f90 .o

all: comm $(srcs) $(tsplanes) 
comm:
	cd $(COMMDIR) && make
.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) $<
$(tsplanes): $(objs) $(commonobjs) tsplanes.o
	$(FC) -o $@ $(LINKFLAG) $(objs) $(commonobjs) tsplanes.o $(LIB) $(LIBLINK)
clean:
	/bin/rm -f $(tsplanes) *.o *.mod *.pc *.il *.d

#depencency base on module usage
tsplanes.o: modMetrics.o\
            modDerivative.o  modRWPHDF5.o modVolume.o
 
            

