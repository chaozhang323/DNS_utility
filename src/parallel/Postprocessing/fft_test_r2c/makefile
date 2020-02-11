DIR=../..
2decomp_dir=../../../2decomp_fft
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common

commonobjs=

srcs=
objs=$(srcs:.f90=.o) 
test=fft_test_r2c_old
mpi_io=mpi_IO
.SUFFIXES: .f90 .o

all:  $(srcs) $(test)  
comm:
	cd $(COMMDIR) && make
.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) $<
$(test): $(objs) fft_test_r2c_old.o
	$(FC) -o $@ $(LINKFLAG) fft_test_r2c_old.o $(commonobjs) $(LIB) $(LIBLINK)
$(mpi_io): $(objs) mpi_IO.o
	$(FC) -o $@ $(LINKFLAG) mpi_IO.o $(commonobjs) $(LIB) $(LIBLINK)

clean:
	/bin/rm -f $(test) $(mpi_io) *.o *.mod *.pc *.il *.d
#depencency base on module usage
#timeave.0: $(COMMDIR)/modTecplotIO.o
