#FC=mpif90
FC=ifort
#FC=ftn
FFLAG=-O3 -cpp -r8 -g -traceback -check bounds -heap-arrays
#FFLAG=-r8 -check all,noarg_temp_created -align -fpe0 -g -O0 -traceback -heap-arrays -debug all -ftrapuv
#FFLAG=-O3 -assume buffered_io -axAVX -xSSE4.1 -r8 -zero -align -traceback -heap-arrays -fpe0 -ftrapuv -mcmodel medium -fp-model precise
#FFLAG=-O3 -axAVX -xSSE4.1 -ip -r8 -zero -traceback -heap-arrays
#FFLAG=-O3 -cpp -h zero -s real64
FFLAG+=-fopenmp
LINKFLAG+=-fopenmp

# For HDF5/MKL Lib at blue waters
#HDF5=/opt/cray/hdf5-parallel/1.8.16/CRAY/8.3
#MKL=/opt/intel/composer_xe_2015.3.187/mkl

# For HDF5 Lib at forge
MKL=/opt/intel/composer_xe_2015.2.164/mkl
#HDF5=/share/duan/czb58/Package/HDF5_serial
## the following is the library for virtual dataset
HDF5=/duan/data/czb58/Package/HDF5-1.10.1_serial


# For HDF5 Lib at NASA Pleiades
#HDF5=/nasa/hdf5/1.8.18_serial
#MKL=/nasa/intel/Compiler/2017.4.196/compilers_and_libraries_2017.4.196/linux/mkl

# For HDF5 Lib at AirForce ORS Utility Server
#HDF5=/u/US_HOME2/czb58/Package/HDF5
# MKL=/opt/intel/composerxe-2011.2.137/mkl

# For HDF5 Lib at AFRL Spirit
#HDF5=/app/wpostool/COST/hdf5-mpi/1.8.15/intel/sgimpt
#MKL=/app/intel/mkl

# For HDF5 Lib at AFRL Lightning
#HDF5=/opt/cray/hdf5-parallel/1.8.14/INTEL/14.0
#MKL=/opt/intel/composer_xe_2015.3.187/mkl


# For HDF5 Lib at Excalibur
#HDF5=/opt/cray/hdf5/1.8.14/intel/14.0
#MKL=/opt/intel/compilers_and_libraries_2017.1.132/linux/mkl

# For HDF5 Lib at AFRL Thunder
#HDF5=/app/COST/hdf5-mpi/1.8.15/intel/sgimpt
#MKL=/app/intel/composer_xe_2015.3.187/mkl

tec360=$(DIR)/../../tecio
fftpack=$(DIR)/../../fftpack5.1

INCLUDE=-I $(HDF5)/include -I $(COMMDIR) -I $(tec360)/include -I $(MKL)/include/intel64/lp64
LIB=-L$(HDF5)/lib -L $(tec360)/lib -L$(MKL)/lib/intel64   -L$(fftpack)/lib/ 
LIBLINK= $(HDF5)/lib/libhdf5_fortran.a $(tec360)/lib/libtecio.a $(MKL)/lib/intel64/libmkl_intel_lp64.a $(fftpack)/lib/libfftpack.a  -lfftpack -lhdf5   -lz  -lm -lstdc++ \
          -Wl,--start-group \
          -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
          -Wl,--end-group -lpthread



#FC=gfortran
#FFLAG=-O3 -fdefault-real-8
#LINKFLAG=
#LINKFLAG=/nasa/ae/tecplot360/lib/tecio64.a
