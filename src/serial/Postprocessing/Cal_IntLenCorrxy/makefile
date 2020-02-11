DIR=../..
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common
EXTRA_FFLAG=-openmp

PROGRAMS=Cal_IntLenCorrxy
commonobjs=$(COMMDIR)/modGrid1d.o $(COMMDIR)/interp1d.o $(COMMDIR)/tridiagonal.o \
           $(COMMDIR)/modTecplotIO.o $(COMMDIR)/modTecbin.o $(COMMDIR)/modRWHDF5.o

all: comm $(PROGRAMS)

comm:
	cd $(COMMDIR) && make

%: %.o
	$(FC) -o $@ $(LINKFLAG) $^ $(commonobjs) $(LIB) $(LIBLINK) $(EXTRA_FFLAG)
%.o: %.f90
	$(FC) $(FFLAG) $(EXTRA_FFLAG) $(INCLUDE) -c $<
	
clean:
	/bin/rm -f $(PROGRAMS) *.o *.mod *.pc *.il *.d

# dependencies
Cal_IntLenCorrxy.o: 