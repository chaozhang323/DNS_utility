DIR=../../..
include $(DIR)/makecomm.inc
COMMDIR=$(DIR)/common


commonobjs=$(COMMDIR)/modTecplotIO.o  $(COMMDIR)/modRWHDF5.o
srcs= modCalInt.f90

objs=$(srcs:.f90=.o) 
Cal_kplane=Cal_kplane
Cal_jplane=Cal_jplane
Cal_iplane=Cal_iplane
Validate_Acoustic=Validate_Acoustic
.SUFFIXES: .f90 .o

all: comm $(srcs)  $(Cal_kplane) $(Cal_jplane) $(Cal_iplane) $(Validate_Acoustic)
comm:
	cd $(COMMDIR) && make
.f90.o:
	$(FC) -c $(FFLAG) $(INCLUDE) $<
$(Cal_kplane): $(objs) Cal_kplane.o
	$(FC) -o $@ $(LINKFLAG) Cal_kplane.o $(objs) $(LIB) $(LIBLINK)
$(Cal_jplane): $(objs) Cal_jplane.o
	$(FC) -o $@ $(LINKFLAG) Cal_jplane.o $(objs) $(LIB) $(LIBLINK)
$(Cal_iplane): $(objs) Cal_iplane.o
	$(FC) -o $@ $(LINKFLAG) Cal_iplane.o $(objs) $(LIB) $(LIBLINK)
$(Validate_Acoustic): $(objs) $(commonobjs) Validate_Acoustic.o
	$(FC) -o $@ $(LINKFLAG) Validate_Acoustic.o $(objs) $(commonobjs) $(LIB) $(LIBLINK)

clean:
	/bin/rm -f $(Cal_kplane) $(Cal_jplane) $(Cal_iplane) $(Validate_Acoustic) *.o *.mod *.pc *.il *.d

#depencency base on module usage
#dnsGridConvert.o: $(COMMDIR)/modFileIO_Compatible.o 
Validate_Acoustic.o: $(COMMDIR)/modTecplotIO.o modCalInt.o  $(COMMDIR)/modRWHDF5.o
