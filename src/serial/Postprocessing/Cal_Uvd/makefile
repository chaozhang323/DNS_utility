DIR=../..
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common

commonobjs=$(COMMDIR)/modTecplotIO.o  $(COMMDIR)/modRWHDF5.o
srcs= modCalInt.f90

objs=$(srcs:.f90=.o) 
Cal_Uvd=Cal_Uvd
.SUFFIXES: .f90 .o

all: comm $(srcs)  $(Cal_Uvd)
comm:
	cd $(COMMDIR) && make
.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) $<
$(Cal_Uvd): $(objs) $(commonobjs) Cal_Uvd.o
	$(FC) -o $@ $(LINKFLAG) Cal_Uvd.o $(objs) $(commonobjs) $(LIB) $(LIBLINK)

clean:
	/bin/rm -f $(Cal_Uvd) *.o *.mod *.pc *.il *.d

#depencency base on module usage
Cal_Uvd.o: $(COMMDIR)/modTecplotIO.o modCalInt.o  $(COMMDIR)/modRWHDF5.o
