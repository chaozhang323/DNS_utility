DIR=../..
2decomp_dir=../../../2decomp_fft
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common
#decompdir=$(2decomp_dir)/2decomp_fft/src

commonobjs=$(COMMDIR)/modPRWHDF5.o $(COMMDIR)/modFileIO_Compatible.o

srcs=
objs=$(srcs:.f90=.o) 
p3d_h5=p3d_h5
.SUFFIXES: .f90 .o

all: 2decomp comm $(srcs) $(p3d_h5)
2decomp:
	cd $(2decomp_dir)/src && make
comm:
	cd $(COMMDIR) && make
.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) -I$(COMMDIR) $<
$(p3d_h5): $(objs) p3d_h5.o
	$(FC) -o $@ $(LINKFLAG) p3d_h5.o $(commonobjs) $(LIB) $(LIBLINK)


clean:
	/bin/rm -f $(p3d_h5) *.o *.mod *.pc *.il *.d

#depencency base on module usage
p3d_h5.o:  $(COMMDIR)/modPRWHDF5.o $(COMMDIR)/modFileIO_Compatible.o
