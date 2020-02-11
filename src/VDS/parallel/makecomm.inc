# intel compiler
FC=mpif90
FFLAG=-O3 -cpp -r8 -g
#FFLAG=-r8 -cpp -check all,noarg_temp_created -align -fpe0 -g -O0 -traceback -heap-arrays -debug all -ftrapuv

# Cray compiler
#FC=ftn
#FFLAG=-O3 -e Fm -s real64

LINKFLAG=

# Libs at Forge
PHDF5=/duan/data/czb58/Package/HDF5-1.10.1_parallel
FFTW=/opt/fftw/3.3.4/intel2015.2.164/openmpi_ib

 

INCLUDE=-I$(PHDF5)/include -I$(2decomp_dir)/include 
LIB=-L $(2decomp_dir)/lib -L$(PHDF5)/lib -L$(FFTW)/lib 
LIBLINK=$(2decomp_dir)/lib/lib2decomp_fft.a  $(PHDF5)/lib/libhdf5_fortran.a  -lm -lhdf5   -lz   -lstdc++ $(FFTW)/lib/libfftw3.a -lfftw3 








