FC = gfortran -t -v
RM = rm
RMFLAGS = -f
MAKE = make
MAKEFLAGS = b
WARN= -Wall -Wextra -Wconversion -pedantic -fbounds-check
DEBUG= -g
OPENMP= -fopenmp
OPTIMIZE= -O3 -fexternal-blas
INCLUDE= -I/sw/include   # directory(s) to search for library include files
LIBDIR= -L/sw/lib        # directory(s) to search for library archives
LIBFITSIO= -lcfitsio     # link to cfitsio library
LIBFFTW3= -lfftw3        # link to FFTW3 library
LIBVEC= -lblas -llapack  # link to BLAS and LAPACK libraries

all: imcom

imcom: imcom_data.o imcom_proc.o imcom_io.o imcom_matrix.o imcom_bisect.o imcom.f90
	$(FC) $(WARN) $(OPTIMIZE) $(OPENMP) imcom_data.o imcom_proc.o imcom_io.o imcom_matrix.o imcom_bisect.o imcom.f90 -o $@ $(LIBDIR) $(INCLUDE) $(LIBVEC) $(LIBFITSIO) $(LIBFFTW3)
	cp ./imcom ../bin/
	cp ./imcom ../example/

imcom_bisect.o: imcom_data.o imcom_proc.o imcom_io.o imcom_matrix.o imcom_bisect.f90
	$(FC) $(WARN) $(OPTIMIZE) $(OPENMP) -c imcom_bisect.f90 -o $@

imcom_matrix.o: imcom_data.o imcom_proc.o imcom_io.o imcom_matrix.f90
	$(FC) $(WARN) $(OPTIMIZE) $(OPENMP) -c imcom_matrix.f90 -o $@

imcom_io.o : imcom_data.o imcom_proc.o imcom_io.f90
	$(FC) $(WARN) $(OPTIMIZE) $(OPENMP) -c imcom_io.f90 -o $@

imcom_proc.o : imcom_data.o imcom_proc.f90
	$(FC) $(WARN) $(OPTIMIZE) $(OPENMP) -c imcom_proc.f90 -o $@

imcom_data.o : imcom_data.f90
	$(FC) $(WARN) $(OPTIMIZE) $(OPENMP) -c imcom_data.f90 -o $@

clean :
	$(RM) $(RMFLAGS) *.o *.mod imcom
