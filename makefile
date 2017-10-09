MKLROOT = /opt/intel/composer_xe_2013.4.183/mkl

all: aout

aout: mod_rate.o rate.o
	ifort -o aout mod_rate.o rate.o -O2 -ipo -xHost /data/home/amber/lib/libcommon.so /data/home/amber/gsl_lib/libgsl.a 

%.o: %.f90
	ifort -c -ipo $<

clean:
	rm *.o aout
