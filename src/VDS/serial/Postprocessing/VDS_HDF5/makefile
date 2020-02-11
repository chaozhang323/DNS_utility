DIR=../..
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common

commonobjs=$(COMMDIR)/modRWHDF5.o

srcs=modVDSHDF5.f90
objs=$(srcs:.f90=.o) 
VDS_HDF5=VDS_HDF5
.SUFFIXES: .f90 .o

all: comm $(srcs) $(VDS_HDF5)
comm:
	cd $(COMMDIR) && make
.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) $<
$(VDS_HDF5): $(objs) VDS_HDF5.o
	$(FC) -o $@ $(LINKFLAG) VDS_HDF5.o $(commonobjs) $(objs)  $(LIB) $(LIBLINK)


clean: 
	/bin/rm -f $(VDS_HDF5) *.o *.mod *.pc *.il *.d
 

#depencency base on module usage
VDS_HDF5.o: $(COMMDIR)/modRWHDF5.o modVDSHDF5.o
