DIR=../..
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common

commonobjs=$(COMMDIR)/modRWHDF5.o

srcs= 
objs=$(srcs:.f90=.o) 
Cal=Cal_LengthScale
coh=Cal_Length
Cal_TL=Cal_TurbLengthScale
.SUFFIXES: .f90 .o

all:  comm $(srcs) $(coh) $(Cal_TL)
comm:
	cd $(COMMDIR) && make
.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) $<
$(coh): $(objs) Cal_Length.o
	$(FC) -o $@ $(LINKFLAG) Cal_Length.o $(objs) $(commonobjs) $(LIB) $(LIBLINK)
$(Cal_TL): $(objs) Cal_TurbLengthScale.o
	$(FC) -o $@ $(LINKFLAG) Cal_TurbLengthScale.o $(objs) $(commonobjs) $(LIB) $(LIBLINK)


clean:
	/bin/rm -f $(coh) *.o *.mod *.pc *.il *.d
#depencency base on module usage
Cal_Length.o:$(COMMDIR)/modRWHDF5.o
Cal_TurbLengthScale.o:$(COMMDIR)/modRWHDF5.o
