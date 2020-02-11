DIR=../..
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common
#INCLUDE=-I/share/apps/hdf5-mpich2-i12/include -I $(COMMDIR) -I/home/czb58/tecplot360/include
#LIB=-L/share/apps/hdf5-mpich2-i12/lib -L/home/czb58/tecplot360/lib
#LIBLINK= /share/apps/hdf5-mpich2-i12/lib/libhdf5_fortran.a /home/czb58/tecplot360/lib/libtecio.a -lm -lstdc++ -lhdf5   -lz  

commonobjs=$(COMMDIR)/modTecplotIO.o $(COMMDIR)/modRWHDF5.o

srcs= 
objs=$(srcs:.f90=.o) 
Cal_Ub=Cal_Ub
.SUFFIXES: .f90 .o

all:  comm $(srcs) $(Cal_Ub)
comm:
	cd $(COMMDIR) && make
.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) $<
$(Cal_Ub): $(objs) Cal_Ub.o
	$(FC) -o $@ $(LINKFLAG) Cal_Ub.o $(objs) $(commonobjs) $(LIB) $(LIBLINK)


clean:
	/bin/rm -f $(Cal_Ub) *.o *.mod *.pc *.il *.d
#depencency base on module usage
Cal_Ub.o: $(COMMDIR)/modTecplotIO.o $(COMMDIR)/modRWHDF5.o
