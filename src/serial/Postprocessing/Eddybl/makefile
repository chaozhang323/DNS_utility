FC=ifort
F77FFLAGS=-g -debug -r8 -fixed -warn alignments  -zero  -f77 -fpp -warn noalign  
#-extend_source -fpp -f77rtl -intconstant -ftz -pad-source -sox  \
#          -warn alignments  

exe=eddybl
start=eddybl_start

srcs=eddybl.f eddybl_start.f
objs=$(srcs:.f=.o)
.SUFFIXES: .f .o

all: $(exe) $(start)


.f.o:
	$(FC) -c $(F77FFLAGS) $<

$(exe): eddybl.o
	$(FC) -o $@ eddybl.o

$(start): eddybl_start.o
	  $(FC) -o $@ eddybl_start.o


clean:
	/bin/rm -f $(exe) $(start) *.o *.mod *.pc *.il *.d

#InitSetUp.o:  
