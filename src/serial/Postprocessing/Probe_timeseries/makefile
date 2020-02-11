DIR=../..
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common

PROGRAMS=probe_ts
commonobjs=$(COMMDIR)/modTecplotIO.o $(COMMDIR)/modTecbin.o $(COMMDIR)/modRWHDF5.o\
           $(COMMDIR)/modFileIO_Compatible.o

all: comm $(PROGRAMS)

comm:
	cd $(COMMDIR) && make

%: %.o
	$(FC) -o $@ $(LINKFLAG) $^ $(commonobjs) $(LIB) $(LIBLINK)
%.o: %.f90
	$(FC) $(FFLAG) $(INCLUDE) -c $<
	
clean:
	/bin/rm -f $(PROGRAMS) *.o *.mod *.pc *.il *.d
	
probe_ts.o: $(COMMDIR)/modTecplotIO.o $(COMMDIR)/modTecbin.o $(COMMDIR)/modRWHDF5.o\
               $(COMMDIR)/modFileIO_Compatible.o