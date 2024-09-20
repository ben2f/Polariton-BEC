ifndef $(compiler)
    compiler = intel
endif

ifeq ($(compiler), gnu)
   FC = gfortran
   FFLAGS = -O3 -Wall #-Wnoexcept
   FFTW_INCS = -I$(MKLROOT)/include/fftw
   MKL_LIB_PATH = -L$(MKLROOT)/lib
   FLIBS = -lmkl_sequential -lmkl_core -lpthread -lmkl_intel_lp64
endif

ifeq ($(compiler), intel)
   FC = ifort
   FFLAGS = -Ofast -xHOST -fp-model strict #-check all
   FFTW_INCS = -I$(MKLROOT)/include/fftw
   MKL_LIB_PATH = -L$(MKLROOT)/lib
   FLIBS = -lmkl_sequential -lmkl_core -lpthread -lmkl_intel_lp64
endif

all: four1d-pbec

real1d-pbec: clean
	$(FC) $(FFLAGS) -c -o real1d-pbec.o src/real1d-pbec.f90
	$(FC) $(FFLAGS) real1d-pbec.o -o real1d-pbec
	rm -rf real1d-pbec.o

four1d-pbec: futil
	$(FC) $(FFLAGS) -c src/four1d-pbec.f90 -o four1d-pbec.o
	$(FC) $(FFLAGS) four1d-pbec.o futil.o -o four1d-pbec $(FLIBS)
	rm -rf four1d-pbec.o futil.o

futil:
	$(FC) $(FFLAGS) $(FFTW_INCS) -c src/futil.f90 -o futil.o

clean:
	rm -rf *.o *.mod real1d-pbec four1d-pbec
