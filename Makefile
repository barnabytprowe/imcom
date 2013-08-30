FC = gfortran -t -v
RM = rm
RMFLAGS = -f
MAKE = make
MAKEFLAGS = b
ERRTRACE=#-fbacktrace
WARN= -Wall -Wextra -Wconversion -pedantic -fbounds-check
DEBUG= -g
OPENMP= -fopenmp
OPTIMIZE= -O3 -fexternal-blas
INCLUDE= -I/sw/include/
LIBDIR= -L/sw/lib# -L/Users/barnabyrowe/fortran/lib
LIBFITSIO= -lcfitsio
LIBFFTW3= -lfftw3
LIBVEC= -Wl,-framework -Wl,Accelerate
#LIBVEC= -lblas -llapack_ref

all: imcom

imcom: imcom.f90 imcom_data.o imcom_proc.o imcom_io.o imcom_matrix.o imcom_bisect.o
	$(FC) $(WARN) $(ERRTRACE) $(OPTIMIZE) $(OPENMP) $^ -o $@ $(LIBDIR) $(INCLUDE) $(LIBVEC) $(LIBFITSIO) $(LIBFFTW3)
	cp ./imcom ../bin/
	cp ./imcom ../jdem/
	cp ./imcom ~/fortran/bin/

imcom_bisect.o: imcom_bisect.f90 imcom_data.o imcom_proc.o imcom_io.o imcom_matrix.o
	$(FC) $(WARN) $(ERRTRACE) $(OPTIMIZE) $(OPENMP) -c $< -o $@

imcom_matrix.o: imcom_matrix.f90 imcom_data.o imcom_proc.o imcom_io.o
	$(FC) $(WARN) $(ERRTRACE) $(OPTIMIZE) $(OPENMP) -c $< -o $@

imcom_io.o : imcom_io.f90 imcom_data.o imcom_proc.o
	$(FC) $(WARN) $(ERRTRACE) $(OPTIMIZE) $(OPENMP) -c $< -o $@

imcom_proc.o : imcom_proc.f90 imcom_data.o
	$(FC) $(WARN) $(ERRTRACE) $(OPTIMIZE) $(OPENMP) -c $< -o $@

imcom_data.o : imcom_data.f90
	$(FC) $(WARN) $(ERRTRACE) $(OPTIMIZE) $(OPENMP) -c $< -o $@

clean :
	$(RM) $(RMFLAGS) *.o *.mod imcom
