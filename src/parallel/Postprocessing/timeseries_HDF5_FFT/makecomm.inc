# intel compiler
FC=mpif90
FFLAG=-O3 -cpp -r8 -g
#FFLAG=-r8 -cpp -check all,noarg_temp_created -align -fpe0 -g -O0 -traceback -heap-arrays -debug all -ftrapuv

# Cray compiler
#FC=ftn
#FFLAG=-O3 -e Fm -s real64

LINKFLAG=

# Libs at nic-cluster
PHDF5=/share/duan/czb58/Package/HDF5_Szip
FFTW=/share/apps/fftw-3.3.4

# Libs at copper utility server
#PHDF5=/u/home/czb58/Package/PHDF5_intel_13.0
#FFTW=/app/projects/COST/fftw-3.3.3

# Libs at NASA Pleiades
#PHDF5=/nasa/hdf5/1.8.18_mpt
#FFTW=/home6/lduan/bin/fftw-3.3.4_install

# Libs at Excalibur
#PHDF5=/opt/cray/hdf5-parallel/1.8.14/CRAY/8.3
#FFTW=/app/unsupported/COST/fftw3-mpi/3.3.5/gnu/craympt

fftpack=$(DIR)/../fftpack5.1

INCLUDE=-I$(PHDF5)/include -I$(2decomp_dir)/include  -I$(FFTW)/include  
LIB=-L $(2decomp_dir)/lib -L$(PHDF5)/lib -L$(FFTW)/lib -L$(fftpack)/lib/ 
LIBLINK=$(2decomp_dir)/lib/lib2decomp_fft.a  $(PHDF5)/lib/libhdf5_fortran.a $(fftpack)/lib/libfftpack.a  -lfftpack -lm -lhdf5   -lz   -lstdc++ $(FFTW)/lib/libfftw3.a -lfftw3 








