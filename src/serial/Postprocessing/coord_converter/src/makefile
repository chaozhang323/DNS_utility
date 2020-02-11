DIR=../../..
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common

commonobjs=$(COMMDIR)/modRWHDF5.o 

srcs=coord_converter.f90 main.f90
objs=$(srcs:.f90=.o) 
main=coord_converter
.SUFFIXES: .f90 .o

all: comm $(srcs) $(main) 
comm:
	cd $(COMMDIR) && make
.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) $<
$(main): $(objs) 
	$(FC) -o $@  $(LINKFLAG) $(commonobjs) $(objs) $(LIB) $(LIBLINK)


clean:
	/bin/rm -f $(main) *.o *.mod *.pc *.il *.d

#depencency base on module usage
main.o:  $(COMMDIR)/modRWHDF5.o 
