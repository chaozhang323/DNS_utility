FC=ifort
FFLAG=-O3 -axAVX -xSSE4.1 -r8 -fast -zero -warn noalign -traceback

INCLUDE=-I/share/apps/hdf5-mpich2-i12/include
LIB=-L/share/apps/hdf5-mpich2-i12/lib
LIBLINK= /share/apps/hdf5-mpich2-i12/lib/libhdf5_fortran.a -lhdf5   -lz  -lm -lstdc++ 

srcs=
objs=$(srcs:.f90=.o) 
Convert_timeseries=Convert_timeseries
.SUFFIXES: .f90 .o

all:  $(srcs) $(Convert_timeseries)

.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) $<
$(Convert_timeseries): $(objs) Convert_timeseries.o
	$(FC) -o $@ $(LINKFLAG) Convert_timeseries.o $(objs) $(LIB) $(LIBLINK)
clean:
	/bin/rm -f $(Convert_timeseries) *.o *.mod *.pc *.il *.d

#depencency base on module usage
#dnsGridConvert.o: $(COMMDIR)/modFileIO_Compatible.o 
