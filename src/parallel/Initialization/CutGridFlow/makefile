DIR=../..
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common
2decomp_dir=../../../2decomp_fft

commonobjs=$(COMMDIR)/modFileIO_Compatible.o $(COMMDIR)/modPRWHDF5.o
exe=CutGridFlow_parallel

all: 2decomp comm $(srcs) $(exe)
2decomp:
	cd $(2decomp_dir)/src && make
comm:
	cd $(COMMDIR) && make
%.o: %.f90
	$(FC) -c $(FFLAG) $(INCLUDE) -I$(COMMDIR) $<
%: %.o
	$(FC) -o $@ $(LINKFLAG) $^ $(commonobjs) $(LIB) $(LIBLINK)
clean:
	/bin/rm -f $(exe) *.o *.mod *.pc *.il *.d

#depencency base on module usage
CutGridFlow_parallel.o: $(COMMDIR)/modFileIO_Compatible.o $(COMMDIR)/modPRWHDF5.o 