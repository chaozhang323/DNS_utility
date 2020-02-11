DIR=../..
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common

commonobjs=$(COMMDIR)/modTecplotIO.o  $(COMMDIR)/modRWHDF5.o
srcs= modCalInt.f90

objs=$(srcs:.f90=.o) 
Cal_IntVar=Cal_IntVar
.SUFFIXES: .f90 .o

all: comm $(srcs)  $(Cal_IntVar)
comm:
	cd $(COMMDIR) && make
.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) $<
$(Cal_IntVar): $(objs) $(commonobjs) Cal_IntVar.o
	$(FC) -o $@ $(LINKFLAG) Cal_IntVar.o $(objs) $(commonobjs) $(LIB) $(LIBLINK)

clean:
	/bin/rm -f $(Cal_IntVar) *.o *.mod *.pc *.il *.d

#depencency base on module usage
Cal_IntVar.o: $(COMMDIR)/modTecplotIO.o modCalInt.o  $(COMMDIR)/modRWHDF5.o
