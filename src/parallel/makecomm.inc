# intel compiler
FC=mpif90
FFLAG=-O3 -cpp -r8 -g
#FFLAG=-r8 -cpp -check all,noarg_temp_created -align -fpe0 -g -O0 -traceback -heap-arrays -debug all -ftrapuv

# Cray compiler
#FC=ftn
#FFLAG=-O3 -e Fm -s real64

LINKFLAG=

# Libs at Forge
PHDF5=/opt/hdf5/intel2015.2.164/openmpi_ib
FFTW=/opt/fftw/3.3.4/intel2015.2.164/openmpi_ib

# Libs at NASA Pleiades
#PHDF5=/nasa/hdf5/1.8.18_mpt
#FFTW=/home6/lduan/bin/fftw-3.3.4_install

# Libs at Thunder
#PHDF5=/app/COST/hdf5-mpi/1.8.15/intel/sgimpt
#FFTW=/app/COST/fftw3-mpi/3.3.4/intel/sgimpt

# Libs at Excalibur
#PHDF5=/opt/cray/hdf5-parallel/1.8.14/CRAY/8.3
#FFTW=/app/unsupported/COST/fftw3-mpi/3.3.5/gnu/craympt

# Libs at blue waters
#PHDF5=/opt/cray/hdf5-parallel/1.8.16/INTEL/15.0
#FFTW=/opt/cray/fftw/3.3.4.7/interlagos


fftpack=$(DIR)/../fftpack5.1

INCLUDE=-I$(PHDF5)/include -I$(2decomp_dir)/include  -I$(FFTW)/include  
LIB=-L $(2decomp_dir)/lib -L$(PHDF5)/lib -L$(FFTW)/lib -L$(fftpack)/lib/ 
LIBLINK=$(2decomp_dir)/lib/lib2decomp_fft.a  $(PHDF5)/lib/libhdf5_fortran.a $(fftpack)/lib/libfftpack.a  -lm -lhdf5   -lz   -lstdc++ $(FFTW)/lib/libfftw3.a -lfftw3 








