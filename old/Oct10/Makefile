FC = gfortran
AR = ar
AROPTIONS = rvus
RM = rm
RMFLAGS = -f
MAKE = make
MAKEFLAGS = b
DOUBLE= -fdefault-real-8 -fdefault-integer-8
ERRTRACE=#-fbacktrace
WARN= -Wall -Wextra -Wconversion -pedantic -fbounds-check
DEBUG= -g
OPENMP= -fopenmp
WATCH= -watch
OPTIMIZE= -O3 -fexternal-blas
INCLUDE= -I/sw/include/
LIBFITSIO= -L/sw/lib/ -lcfitsio
LIBFFTW3= -L/sw/lib/ -lfftw3
LIBVEC= -lblas -llapack -framework vecLib

all: imcom

imcom: imcom_data.o imcom_proc.o imcom_io.o imcom_matrix.o imcom.f90
	$(FC) $(WARN) $(ERRTRACE) $(OPTIMIZE) $(OPENMP) imcom_data.o imcom_io.o imcom_proc.o imcom_matrix.o imcom.f90 -o $@ $(INCLUDE) $(LIBFITSIO) $(LIBFFTW3) $(LIBVEC) 
	cp ./imcom ../bin/
	cp ./imcom ../jdem/
	cp ./imcom ~/fortran/bin/

imcom_matrix.o: imcom_data.o imcom_proc.o imcom_io.o imcom_matrix.f90
	$(FC) $(WARN) $(ERRTRACE) $(OPTIMIZE) $(OPENMP) -c imcom_matrix.f90 -o $@

imcom_io.o : imcom_data.o imcom_proc.o imcom_io.f90
	$(FC) $(WARN) $(ERRTRACE) $(OPTIMIZE) $(OPENMP) -c imcom_io.f90 -o $@

imcom_proc.o : imcom_data.o imcom_proc.f90
	$(FC) $(WARN) $(ERRTRACE) $(OPTIMIZE) $(OPENMP) -c imcom_proc.f90 -o $@

imcom_data.o : imcom_data.f90
	$(FC) $(WARN) $(ERRTRACE) $(OPTIMIZE) $(OPENMP) -c imcom_data.f90 -o $@

clean :
	$(RM) $(RMFLAGS) *.o *.mod imcom
