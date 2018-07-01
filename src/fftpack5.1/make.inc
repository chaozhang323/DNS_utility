
LIB=../lib/libfftpack.a

UNAMES := $(shell uname -s)

ifeq ($(UNAMES),Linux)

  INTEL := $(shell ifort 2>&1)

  ifeq ($(INTEL),ifort-Warning-No files to process)

    F90 := ifort ../lib -I../lib
#   F90 := ifort -r8 -module ../lib -I../lib
    CPP := ifort -E

  else

    F90 := ifort -O3 -r8 -L ../lib -I../lib 
    CPP := ifort -E 

  endif

  MAKE := gmake
#  AR := /nobackupp6/lduan/fftpack5.1/ar

endif

ifeq ($(UNAMES),AIX)

  F90 := xlf -qmoddir=../lib -I../lib
# F90 := xlf -qautodbl=dbl -qmoddir=../lib -I../lib -lxlfpmt4 -lxlfpmt8
  CPP := xlf -d -qnoobject
  MAKE := gmake
  AR := /usr/bin/ar

endif

ifeq ($(UNAMES),SunOS)

    AR := /usr/ccs/bin/ar
    F90 := /opt/SUNWspro/bin/f90 -moddir=../lib -I../lib
    CPP := /opt/SUNWspro/bin/f90 -F
    MAKE := /fs/local/bin/make

endif

ifeq ($(UNAMES),IRIX64)

    AR := /usr/bin/ar
    F90 := f90 -I../lib
    CPP := f90 -E
    MAKE := /usr/local/bin/gmake

endif

ifeq ($(UNAMES),Darwin)

    AR := /usr/bin/ar
    F90 := gfortran
#   F90 := gfortran -fdefault-real-8
    CPP := gfortran -cpp
    MAKE := /usr/bin/gmake

endif
