# Makefile for dysk_zwm

#-----------------------------------
# for intel compiler
#-----------------------------------
#CC = /opt/intel/64/bin/icpc
#CFLAGS = -O3 -Wall -ipo
#FF = /opt/intel/64/bin/ifort
#FFLAGS = -O3 -ipo
#LD = $(CC) -ipo
#LDFLAGS = -lifcore
#-----------------------------------
#-----------------------------------
#CC = mpicxx
#CFLAGS = -O3 -DMPI_USED
DYSK_EXE = ../exe/

CC = mpic++
CFLAGS = -O2 -DMPI_USED -DHDF5_USED -I/usr/include/hdf5/openmpi/ -Wno-unused-result

LD = $(CC) 
LDFLAGS = -lm -Llib -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi/ -lhdf5 -lz
#-----------------------------------

CP = cp -fp
RM = rm -f


COBJS = dysk_zwm.o steady.o nsolve.o inout.o konfig.o readini.o\
	boundary.o mdothold.o \
       lubksb.o ludcmp.o mnewt.o nrutil.o random.o spline.o


all: dysk_zwm

.c.o:
	$(CC) $(CFLAGS) -c -o $@ $<

.f.o:
	$(FF) $(FFLAGS) -c -o $@ $<


dysk_zwm: $(COBJS) $(FOBJS)
	$(LD) -o ${DYSK_EXE}dysk_zwm $(COBJS) $(FOBJS) $(LDFLAGS)


clean:
	-$(RM) *.o ${DYSK_EXE}dysk_zwm
