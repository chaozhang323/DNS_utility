FC=ifort
FFLAG=-O3 -r8 -zero -warn noalign -traceback

INCLUDE=-I/share/apps/hdf5-mpich2-i12/include
LIB=-L/share/apps/hdf5-mpich2-i12/lib
LIBLINK= /share/apps/hdf5-mpich2-i12/lib/libhdf5_fortran.a -lhdf5   -lz  -lm -lstdc++ -lhdf5

srcs=
objs=$(srcs:.f90=.o) 
p3d_h5_sol=p3d_h5_sol
p3d_h5_grid=p3d_h5_grid
.SUFFIXES: .f90 .o

all:  $(srcs) $(p3d_h5_grid) $(p3d_h5_sol)

.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) $<
$(p3d_h5_sol): $(objs) p3d_h5_sol.o
	$(FC) -o $@ $(LINKFLAG) p3d_h5_sol.o $(objs) $(LIB) $(LIBLINK)
$(p3d_h5_grid): $(objs) p3d_h5_grid.o
	$(FC) -o $@ $(LINKFLAG) p3d_h5_grid.o $(objs) $(LIB) $(LIBLINK)


clean:
	/bin/rm -f $(p3d_h5_sol) $(p3d_h5_grid) *.o *.mod *.pc *.il *.d

#depencency base on module usage
#dnsGridConvert.o: $(COMMDIR)/modFileIO_Compatible.o 
