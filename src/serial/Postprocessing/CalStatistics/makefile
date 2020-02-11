DIR=../..
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common

program=CalStatistics
commonobjs=$(COMMDIR)/modTecplotIO.o $(COMMDIR)/modTecbin.o  $(COMMDIR)/modRWHDF5.o  $(COMMDIR)/modFlow.o

all: comm $(program)
comm:
	cd $(COMMDIR) && make
	
%: %.o
	$(FC) -o $@ $(LINKFLAG) $^ $(commonobjs) $(LIB) $(LIBLINK)
%.o: %.f90
	$(FC) $(FFLAG) $(INCLUDE) -c $<
	
clean:
	/bin/rm -f $(program) *.o *.mod *.pc *.il *.d
	
CalStatistics.o: $(COMMDIR)/modTecplotIO.o $(COMMDIR)/modTecbin.o $(COMMDIR)/modRWHDF5.o $(COMMDIR)/modFlow.o
