!    IMCOM_MATRIX.F90 - Fortran 95 module that handles matrix/vector 
!                       operations for the prototype IMCOM package
!    Copyright (C) 2011 Barnaby Rowe
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    DEPENDENCIES - IMCOM_DATA.F90
!                   IMCOM_PROC.F90
!                   IMCOM_IO.F90
!                   BLAS LIBRARY
!                   LAPACK LIBRARY
!
!    CONTAINS
!    Routines for IMCOM for matrix building; lookup table building; 
!    lookup table interpolation; matrix solvers; matrix/vector 
!    multiplication etc.
!
MODULE imcom_matrix

USE imcom_data
USE imcom_io
USE imcom_proc

implicit none

CONTAINS

!---

SUBROUTINE imcom_get_A(npoly, npad, saveA, forcebuild)
implicit none
integer, intent(IN) :: npoly, npad, saveA, forcebuild
logical :: Aexists
integer :: n1_tmp, n2_tmp, bit_tmp
integer :: alstat

inquire(FILE=trim(Afile), EXIST=Aexists)
if (Aexists.and.(forcebuild.eq.0)) then
  call imcom_sizefits(trim(Afile), n1_tmp, n2_tmp, bit_tmp)
  if ((n1_tmp.ne.n).or.(n2_tmp.ne.n)) then
    write(*, FMT='(A)') "IMCOM ERROR: Dimensions of Matrix "//trim(Afile)//" do not match <config_file> image dimensions"
    write(*, FMT='(A,I9)') "IMCOM ERROR: Should be square of side: ",n
    stop
  endif
  allocate(A_aij(n1_tmp, n2_tmp), STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for A matrix -- too large? Sparse Matrix formalism not yet coded."
    stop
  endif
!$omp parallel workshare
  A_aij = 0.d0
!$omp end parallel workshare
  write(*, FMT='(A)') "IMCOM: Reading A matrix from "//trim(Afile)
  call imcom_readfits(trim(Afile), n1_tmp, n2_tmp, A_aij)
else
  call imcom_build_Alookup(npad)
  call imcom_lookup_A(npoly, npad, saveA)
end if
END SUBROUTINE imcom_get_A

!---

SUBROUTINE imcom_get_B(npoly, npad, saveB, forcebuild)
implicit none
integer, intent(IN) :: npoly, npad, saveB, forcebuild
logical :: Bexists
integer :: n1_tmp, n2_tmp, bit_tmp
integer :: alstat

inquire(FILE=trim(Bfile), EXIST=Bexists)
if (Bexists.and.(forcebuild.eq.0)) then
  call imcom_sizefits(trim(Bfile), n1_tmp, n2_tmp, bit_tmp)
  if ((n1_tmp.ne.n).or.(n2_tmp.ne.m)) then
    write(*, FMT='(A)') "IMCOM ERROR: Dimensions of Matrix "//trim(Bfile)//" do not match <config_file> image dimensions"
    write(*, FMT='(A,2I9)') "IMCOM ERROR: Should be rectangle of sides: ",n,m
    stop
  endif
  allocate(B_ia(n1_tmp, n2_tmp), STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for B matrix -- too large? Sparse Matrix formalism not yet coded."
    stop
  endif
  B_ia = 0.d0
  write(*, FMT='(A)') "IMCOM: Reading B matrix from "//trim(Bfile)
  call imcom_readfits(trim(Bfile), n1_tmp, n2_tmp, B_ia)
else
  call imcom_build_Blookup(npad)
  call imcom_lookup_B(npoly, npad, saveB)
end if
END SUBROUTINE imcom_get_B

!---

SUBROUTINE imcom_get_P(saveP, forcebuild)
implicit none
integer, intent(IN) :: saveP, forcebuild
logical :: Pexists
integer :: n1_tmp, n2_tmp, bit_tmp
integer :: alstat

inquire(FILE=trim(Pfile), EXIST=Pexists)
if (Pexists.and.(forcebuild.eq.0)) then
  call imcom_sizefits(trim(Pfile), n1_tmp, n2_tmp, bit_tmp)
  if ((n1_tmp.ne.n).or.(n2_tmp.ne.m)) then
    write(*, FMT='(A)') "IMCOM ERROR: Dimensions of Matrix "//trim(Pfile)//" do not match <config_file> image dimensions"
    write(*, FMT='(A, 2I9)') "IMCOM ERROR: Should be rectangle of sides: ",n,m
    stop
  endif
  allocate(P_ia(n1_tmp, n2_tmp), STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for P matrix -- too large? Sparse Matrix formalism not yet coded."
    stop
  endif
!$omp parallel workshare
  P_ia = 0.d0
!$omp end parallel workshare
  write(*, FMT='(A)') "IMCOM: Reading P matrix from "//trim(Pfile)
  call imcom_readfits(trim(Pfile), n1_tmp, n2_tmp, P_ia)
  allocate(P2_ia(n1_tmp, n2_tmp), STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for P^2 matrix -- too large? Sparse Matrix formalism not yet coded."
    stop
  endif
!$omp parallel workshare
  P2_ia = P_ia * P_ia
!$omp end parallel workshare
else
 call imcom_build_P(saveP)
end if
END SUBROUTINE imcom_get_P

!---

SUBROUTINE imcom_get_QL(saveQL, forcebuild)
implicit none
integer, intent(IN) :: saveQL, forcebuild
logical :: Qexists, Lexists
integer :: n1_tmp, n2_tmp, bit_tmp
integer :: alstat

inquire(FILE=trim(Qfile), EXIST=Qexists)
inquire(FILE=trim(Lfile), EXIST=Lexists)
if (Qexists.and.Lexists.and.(forcebuild.eq.0)) then
  call imcom_sizefits(trim(Qfile), n1_tmp, n2_tmp, bit_tmp)
  if ((n1_tmp.ne.n).or.(n2_tmp.ne.n)) then
    write(*, FMT='(A)') "IMCOM ERROR: Dimensions of Matrix "//trim(Qfile)//" do not match <config_file> image dimensions"
    write(*, FMT='(A, 2I9)') "IMCOM ERROR: Should be square of side: ",n
    stop
  endif
  allocate(Q_ij(n1_tmp, n2_tmp), STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for Q matrix -- too large? Sparse Matrix formalism not yet coded."
    stop
  endif
!$omp parallel workshare
  Q_ij = 0.d0
!$omp end parallel workshare
  write(*, FMT='(A)') "IMCOM: Reading Q matrix from "//trim(Qfile)
  call imcom_readfits(trim(Qfile), n1_tmp, n2_tmp, Q_ij)
  call imcom_sizefits(trim(Lfile), n1_tmp, n2_tmp, bit_tmp)
  if ((n1_tmp.ne.n).or.(n2_tmp.ne.1)) then
    write(*, FMT='(A)') "IMCOM ERROR: Dimensions of vector "//trim(Lfile)//" do not match <config_file> image dimensions"
    write(*, FMT='(A, 2I9)') "IMCOM ERROR: Should be length: ",n
    stop
  endif
  allocate(L_i(n1_tmp), STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for L vector -- too large? Sparse Matrix formalism not yet coded."
    stop
  endif
!$omp parallel workshare
  L_i = 0.d0
!$omp end parallel workshare
  write(*, FMT='(A)') "IMCOM: Reading L eigenvalues from "//trim(Lfile)
  call imcom_readfits(trim(Lfile), n1_tmp, n2_tmp, L_i)
else
 call imcom_eigen_A(saveQL)
end if
END SUBROUTINE imcom_get_QL

!---

SUBROUTINE imcom_lookup_A(npoly, npad, saveA)
! Calculate the A matrix using the lookup table and polynomial interpolation at degree npoly
implicit none
integer, intent(IN) :: npoly, npad, saveA
real(KIND=8) :: Delx, Dely
integer :: i, j, alstat, n1al, n2al, ial

write(*, FMT='(A)') "IMCOM: Building A matrix from lookup table"
allocate(A_aij(n, n), STAT=alstat)
! Initialize, get correct lookup table size
n1al = n1psf * npad
n2al = n2psf * npad
!$omp parallel
!$omp workshare
A_aij = 0.d0
!$omp end workshare
!$omp do schedule(dynamic, 16) private(Delx, Dely, ial)
do i=1, n

  do j=i, n

!    ial = exp_i(i) + exp_i(j) * (exp_i(j) - 1) / 2
! CODE ABOVE IS CORRECT, CODE BELOW IS TO SAVE MEMORY FOR PROJ TESTS
    ial = 1
    Delx = real(npad, 8) * (x_i(j) - x_i(i)) / psfxscale
    Dely = real(npad, 8) * (y_i(j) - y_i(i)) / psfyscale
    A_aij(i, j) = imcom_interp_lookup(n1al, n2al, Alookup(:, :, ial), Delx, &
                                      Dely, npoly)

  end do
  if (i.eq.int(0.2 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 20% complete"
  if (i.eq.int(0.4 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 40% complete"
  if (i.eq.int(0.6 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 60% complete"
  if (i.eq.int(0.8 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 80% complete"

end do
!$omp end do
!$omp do
do i=2, n

  forall(j=1: i-1) A_aij(i, j) = A_aij(j, i)

end do
!$omp end do
!$omp end parallel
write(*, FMT='(A)') "IMCOM: A matrix complete"
if (saveA.eq.1) then 
  write(*, FMT='(A)') "IMCOM: Saving A matrix to "//trim(Afile)
  call imcom_writefits(trim(Afile), n, n, real(A_aij, 8))
end if
deallocate(Alookup, STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot deallocate memory for A lookup tables"
  stop
endif
END SUBROUTINE imcom_lookup_A

!---

SUBROUTINE imcom_lookup_B(npoly, npad, saveB)
! Calculate the B matrix using the lookup table and polynomial interpolation at degree npoly
implicit none
integer, intent(IN) :: npoly, npad, saveB
real(KIND=8) :: Delx, Dely
integer :: i, a, alstat, n1bl, n2bl

write(*, FMT='(A)') "IMCOM: Building B vectors from lookup table"
allocate(B_ia(n, m), STAT=alstat)
! Initialize
!$omp parallel workshare
B_ia = 0.d0
!$omp end parallel workshare
n1bl = n1psf * npad
n2bl = n2psf * npad
!$omp parallel
!$omp do schedule(dynamic, 64) private(Delx, Dely)
do i=1, n

  do a=1, m

    Delx = real(npad, 8) * (x_i(i) - X_a(a)) / psfxscale
    Dely = real(npad, 8) * (y_i(i) - Y_a(a)) / psfyscale 
!    B_ia(i, a) = imcom_interp_lookup(n1bl, n2bl, Blookup(:, :, exp_i(i)), &
!                                     Delx, Dely, npoly)
! CODE ABOVE IS CORRECT, CODE BELOW IS TO SAVE MEMORY FOR PROJ TESTS
    B_ia(i, a) = imcom_interp_lookup(n1bl, n2bl, Blookup(:, :, 1), &
                                     Delx, Dely, npoly)

  end do
  if (i.eq.int(0.2 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 20% complete"
  if (i.eq.int(0.4 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 40% complete"
  if (i.eq.int(0.6 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 60% complete"
  if (i.eq.int(0.8 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 80% complete"

end do
!$omp end do
!$omp end parallel
write(*, FMT='(A)') "IMCOM: B vectors complete"
if (saveB.eq.1) then 
  write(*, FMT='(A)') "IMCOM: Saving B vectors to "//trim(Bfile)
  call imcom_writefits(trim(Bfile), n, m, real(B_ia, 8))
end if
deallocate(Blookup, STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot deallocate memory for B lookup tables"
  stop
endif
END SUBROUTINE imcom_lookup_B

!---

SUBROUTINE imcom_build_Alookup(npad)
! Build the A matrix lookup table using Fourier transforms
!
implicit none
integer, intent(IN) :: npad
complex(KIND=8), dimension(:, :), allocatable ::  A_tmp, ufunc
integer :: alstat, iexp, jexp, ial
integer :: n1big, n2big, n1min, n2min, n1max, n2max, ntotal
integer(KIND=8) :: ftplan
!character(LEN=256) :: outfile

n1big = n1psf * npad
n2big = n2psf * npad
n1min = (n1big - n1psf) / 2 + 1
n2min = (n2big - n2psf) / 2 + 1
n1max = (n1big + n1psf) / 2
n2max = (n2big + n2psf) / 2
ntotal = n1big * n2big * npad * npad * nexp * nexp
write(*, FMT='(A)') "IMCOM: Fourier transforming A matrix lookup table"
!allocate(Alookup(n1big, n2big, nexp * (nexp + 1) / 2), A_tmp(n1big, n2big), &
!         ufunc(n1big, n2big), STAT=alstat)
! CODE ABOVE IS CORRECT, CODE BELOW IS TO SAVE MEMORY FOR PROJ TESTS
allocate(Alookup(n1big, n2big, 1), A_tmp(n1big, n2big), &
         ufunc(n1big, n2big), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for A lookup tables"
  stop
endif
call imcom_plan_invft_c2c(n1big, n2big, 0, ftplan)
!$omp parallel
!$omp workshare
Alookup = 0.d0
!$omp end workshare
!$omp do private(ufunc, A_tmp, ial) schedule(dynamic, 1)
!do iexp=1, nexp
! CODE ABOVE IS CORRECT, CODE BELOW IS TO SAVE MEMORY FOR PROJ TESTS
do iexp = 1, 1

!  do jexp=iexp, nexp
! CODE ABOVE IS CORRECT, CODE BELOW IS TO SAVE MEMORY FOR PROJ TESTS
   do jexp=iexp, 1

! A matrix lookup tables per exposure stored in a packed upper triangular matrix
    ial = iexp + jexp * (jexp - 1) / 2
    ufunc = dcmplx(0.d0, 0.d0)
!  Fix all the padding with circular shifts... Probably not optimal in 
!  N operations, but it saves me making stupid mistakes with counting...
    ufunc(n1min:n1max, n2min:n2max) = cshift(cshift(dconjg(Gt_rot(:, :,  &
                                                                  iexp)) & 
                                    * Gt_rot(:, :, jexp)                 &
                                    * imcom_upsilon2D(ux, uy, 1),        &
                                      n2psf / 2, 2), n1psf / 2, 1)
    ufunc = cshift(cshift(ufunc, n2big / 2, 2), n1big / 2, 1)
    call imcom_invft_c2c(n1big, n2big, ftplan, ufunc, A_tmp)
    Alookup(:, :, ial) = cshift(cshift(real(A_tmp, 8), n2big / 2, 2), &
                                                       n1big / 2, 1)
!    write(outfile, FMT='(A,I1,A,I1,A)') "Alookup.",iexp,"_",jexp,".fits"
!    call imcom_writefits(trim(outfile), n1big, n2big, real(Alookup(:, :, iexp, jexp), 8))
! ...Origin is thus located at centre of pixel (n1big / 2 + 1, n2big / 2 + 1)

  end do

end do
!$omp end do
!$omp end parallel
deallocate(A_tmp, ufunc, STAT=alstat)
call imcom_destroy_plan(ftplan)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot deallocate temporary memory for A lookup tables"
  stop
endif
!$omp parallel workshare
Alookup = Alookup / real(n1psf, 8) / real(n2psf, 8)
!$omp end parallel workshare
END SUBROUTINE imcom_build_Alookup

!---

SUBROUTINE imcom_build_Blookup(npad)
! Build the B matrix lookup table using Fourier transforms
!
implicit none
integer, intent(IN) :: npad
complex(KIND=8), dimension(:, :), allocatable ::  B_tmp, ufunc
integer :: alstat, iexp
integer :: n1min, n2min, n1max, n2max, n1big, n2big
integer(KIND=8) :: ftplan

n1big = n1psf * npad
n2big = n2psf * npad
n1min = (n1big - n1psf) / 2 + 1
n2min = (n2big - n2psf) / 2 + 1
n1max = (n1big + n1psf) / 2
n2max = (n2big + n2psf) / 2
write(*, FMT='(A)') "IMCOM: Fourier transforming B matrix lookup tables"
!allocate(Blookup(n1big, n2big, nexp), B_tmp(n1big, n2big), &
!         ufunc(n1big, n2big), STAT=alstat)
! CODE ABOVE IS CORRECT, CODE BELOW IS TO SAVE MEMORY FOR PROJ TESTS
allocate(Blookup(n1big, n2big, 1), B_tmp(n1big, n2big), &
         ufunc(n1big, n2big), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot allocate memory for B lookup tables"
  stop
endif
call imcom_plan_invft_c2c(n1big, n2big, 0, ftplan)
!$omp parallel
!$omp workshare
Blookup = 0.d0
!$omp end workshare
!$omp do private(ufunc, B_tmp)
!do iexp=1, nexp
! CODE ABOVE IS CORRECT, CODE BELOW IS TO SAVE MEMORY FOR PROJ TESTS
do iexp=1, 1

  ufunc = dcmplx(0.d0, 0.d0)
  ufunc(n1min:n1max, n2min:n2max) = cshift(cshift(dconjg(Gammat(:, :)) &
                                  * Gt_rot(:, :, iexp)                 & 
                                  * imcom_upsilon2D(ux, uy, 1),        &
                                    n2psf / 2, 2), n1psf / 2, 1)
  ufunc = cshift(cshift(ufunc, n2big / 2, 2), n1big / 2, 1)
  call imcom_invft_c2c(n1big, n2big, ftplan, ufunc, B_tmp)
  Blookup(:, :, iexp) = -2.d0 * cshift(cshift(real(B_tmp, 8), n2big / 2, 2), &
                                                              n1big / 2, 1)
! ...Origin thus located at centre of pixel (n1big / 2 + 1, n2big / 2 + 1)

end do
!$omp end do
!$omp end parallel
deallocate(B_tmp, ufunc, STAT=alstat)
call imcom_destroy_plan(ftplan)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot deallocate temporary memory for B lookup tables"
  stop
endif
!$omp parallel workshare
Blookup = Blookup / real(n1psf, 8) / real(n2psf, 8)
!$omp end parallel workshare
END SUBROUTINE imcom_build_Blookup

!---

FUNCTION imcom_interp_lookup(n1, n2, lookup, xarg, yarg, npoly_in)
implicit none
integer, intent(IN) :: n1, n2, npoly_in
real(KIND=8), dimension(n1, n2), intent(IN) :: lookup
real(KIND=8), intent(IN) :: xarg, yarg  ! x, y in lookup table pixel units (not arcsec!)
real(KIND=8) :: imcom_interp_lookup
real(KIND=8) :: icent, jcent, iarg, jarg, f, df
integer :: npoly, npm
integer :: ia_floor, ia_ceiling, ja_floor, ja_ceiling
integer :: imin, imax, jmin, jmax, ii, jj

!
! THIS SHOULD BE REWRITTEN TO BE:
! a) PROPER IN IT'S DEFINITION OF WHAT NPOLY IS!!!!!
! b) MORE EFFICIENT, CURRENTLY IT IS WHOLE FACTOR OF N TOO SLOW!


npoly = npoly_in
if (mod(npoly, 2).eq.0) then
  write(*, FMT='(A)') "IMCOM WARNING: Odd-order polynomial specified for IMCOM_INTERP_LOOKUP..."
  write(*, FMT='(A)') "IMCOM WARNING: Using input NPOLY-1 for interpolation (e.g. quadratic, quartic...)"
  npoly = npoly - 1
end if
npm = npoly / 2
icent = real(n1 / 2 + 1, 8)
jcent = real(n2 / 2 + 1, 8)
iarg = icent + xarg
jarg = jcent + yarg
ia_floor = floor(iarg)
ja_floor = floor(jarg)
ia_ceiling = ceiling(iarg)
ja_ceiling = ceiling(jarg)
if ((ia_floor.ge.1).and.(ja_floor.ge.1).and.(ia_ceiling.le.n1) &
                                       .and.(ja_ceiling.le.n2)) then
 imin = max(1, ia_floor - npm)
 jmin = max(1, ja_floor - npm)
 imax = min(n1, ia_ceiling + npm)
 jmax = min(n2, ja_ceiling + npm)
 call imcom_polin2((/ (real(ii, 8), ii=imin, imax) /), &
                   (/ (real(jj, 8), jj=jmin, jmax) /), &
                   lookup(imin: imax, jmin: jmax), iarg, jarg, f, df)
 imcom_interp_lookup = f
else
 imcom_interp_lookup = 0.d0
end if
END FUNCTION imcom_interp_lookup

!---

SUBROUTINE imcom_calc_C
implicit none
! Calculates C using the Gammat MTFs... easy!
C_a = sum(real(Gammat * dconjg(Gammat) * imcom_upsilon2D(ux, uy, 1), 8)) &
    / real(n1psf, 8) / real(n2psf, 8)
write(*, FMT='(A, E24.16)') "IMCOM: C_a = ", C_a
END SUBROUTINE imcom_calc_C

!---

SUBROUTINE imcom_build_Ndiag
! Make the (diagonals only) Nii trace vector
implicit none
integer :: i, n_i, n_f, alstat

allocate(Ndiag(n), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate space for noise variance vector"
  stop
end if
Ndiag = 0.d0
!$omp parallel
!$omp do schedule(dynamic, 1) private(n_i, n_f)
do i=1, nexp

  n_i = 1 + (i - 1) * (n / nexp)
  if (i.ne.nexp) then
    n_f = i * (n / nexp)
  else
    n_f = n
  end if
  Ndiag(n_i:n_f) = noise(i)

end do
!$omp end do
!$omp end parallel
END SUBROUTINE imcom_build_Ndiag

!---

SUBROUTINE imcom_eigen_A(saveQL)
! Performs eigendecomposition of A matrix using LAPACK v3.X's relatively 
! robust algorithm DSYEVR (other algorithms are available)
! 
! Not threaded with OpenMP, but some LAPACK implementations allow multi-thread
! support for vector operations
!
implicit none
integer, intent(IN) :: saveQL
character(LEN=1), parameter :: jobz = "V", range = "A", uplo = "U"
real(KIND=8), dimension(:, :), allocatable :: Atmp
real(KIND=8), dimension(:), allocatable :: work
integer, dimension(:), allocatable :: isupport, iwork
real(KIND=8), dimension(1) :: worktest
integer, dimension(1) :: iworktest
real(KIND=8) :: abstol, vl, vu
real(KIND=8), external :: DLAMCH
integer :: alstat, il, iu, n_out, info, lwork, liwork

allocate(Q_ij(n, n), Atmp(n, n), L_i(n), isupport(2 * n), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for Eigendecomposition outputs"
  stop
end if
!$omp parallel workshare
L_i = 0.d0
Q_ij = 0.d0
Atmp = A_aij
!$omp end parallel workshare
abstol = DLAMCH('Safe minimum')  ! Suggested in LAPACK notes for DSYEVR
! First query to get range of work matrices
call DSYEVR(jobz, range, uplo, n, Atmp, n, vl, vu, il, iu, abstol, n_out, &
            L_i, Q_ij, n, isupport, worktest, -1, iworktest, -1, info)
lwork = nint(worktest(1))
liwork = iworktest(1)
allocate(work(lwork), iwork(liwork), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for Eigendecomposition work arrays"
  stop
end if
! Re-initialize in case DSYEVR destroyed previous inputs
!$omp parallel workshare
vl = 0.d0
vu = 0.d0
il = 0
iu = 0
isupport = 0
work = 0.d0
iwork = 0
Atmp = A_aij
L_i = 0.d0
Q_ij = 0.d0
n_out = n
!$omp end parallel workshare
! Then solve eigen vectors and values
write(*, FMT='(A)') "IMCOM: Performing eigendecomposition of A matrix"
call DSYEVR(jobz, range, uplo, n, Atmp, n, vl, vu, il, iu, abstol, n_out, &
            L_i, Q_ij, n, isupport, work, lwork, iwork, liwork, info)
! Test
if (info.ne.0) then
  if (info.lt.0) then
    write(*, FMT='(A, I7)') "IMCOM ERROR: Illegal value in argument i = ",-info
    stop
  else
    write(*, FMT='(A)') "IMCOM ERROR: Internal error in DSYEVR"
    stop
  end if
end if
! Reverse results into descending eigenvalue order 
! (important for minimizing roundoff when summing inverse)
!$omp parallel
!$omp workshare
L_i = L_i(n: 1: -1)
!$omp end workshare
!$omp workshare
Q_ij = Q_ij(:, n: 1: -1)
!$omp end workshare
!$omp end parallel
condition = maxval(L_i) / minval(L_i)
if (condition.lt.0.d0) then
  write(*, FMT='(A)') "IMCOM WARNING: Negative eigenvalues!"
  write(*, FMT='(A, E16.8)') "IMCOM WARNING: L_max = ", maxval(L_i)
  write(*, FMT='(A, E16.8)') "IMCOM WARNING: L_min = ", minval(L_i)
else
  write(*, FMT='(A, E24.16)') "IMCOM: A matrix condition number = ", condition
end if
if (saveQL.ne.0) then
  write(*, FMT='(A)') "IMCOM: Writing eigenvectors to "//trim(Qfile)
  call imcom_writefits(trim(Qfile), n, n, real(Q_ij, 8))
  write(*, FMT='(A)') "IMCOM: Writing eigenvalues to "//trim(Lfile)
  call imcom_writefits(trim(Lfile), n, 1, real(L_i, 8))
end if
END SUBROUTINE imcom_eigen_A

!---

SUBROUTINE imcom_build_P(saveP)
! Builds the projection matrix P_ia using the eigenvectors of A_aij and B_ia
implicit none
integer, intent(IN) :: saveP
character(LEN=1) :: transQ = "T"
character(LEN=1) :: transB = "N"
integer :: alstat

write(*, FMT='(A)') "IMCOM: Building projection matrix"
allocate(P_ia(n, m), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate space for projection matrix"
  stop
end if
!$omp parallel workshare
P_ia = 0.d0
!$omp end parallel workshare
call DGEMM(transQ, transB, n, m, n, 1.d0, Q_ij, n, B_ia, n, 0.d0, P_ia, n)
if (saveP.eq.1) then
  write(*, FMT='(A)') "IMCOM: Writing P matrix to "//trim(Pfile)
  call imcom_writefits(trim(Pfile), n, m, real(P_ia, 8))
end if
allocate(P2_ia(n, m), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate space for squared projection matrix"
  stop
end if
!$omp parallel workshare
P2_ia = P_ia * P_ia
!$omp end parallel workshare
END SUBROUTINE imcom_build_P

!---

SUBROUTINE imcom_build_T
implicit none
!character, parameter :: trans = 'T'
character, parameter :: notrans = 'N'
real(KIND=8), parameter :: alpha = 1.d0, beta = 0.d0
integer, parameter :: incx = 1, incy = 1
integer :: a, alstat

write(*, FMT='(A)') "IMCOM: Building T matrix"
allocate(T_ia(n, m), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM: Cannot allocate space for T matrix"
  stop
end if
!$omp parallel
!$omp do schedule(dynamic, 32)
do a=1, m

  call DGEMV(notrans, n, n, alpha, Q_ij, n, P_ia(:, a) / (L_i(:) + K_a(a)), &
             incx, beta, T_ia(:, a), incy)
  if (a.eq.int(0.2 * real(m, 4))) write(*, FMT='(A)') "IMCOM: 20% complete"
  if (a.eq.int(0.4 * real(m, 4))) write(*, FMT='(A)') "IMCOM: 40% complete"
  if (a.eq.int(0.6 * real(m, 4))) write(*, FMT='(A)') "IMCOM: 60% complete"
  if (a.eq.int(0.8 * real(m, 4))) write(*, FMT='(A)') "IMCOM: 80% complete"

end do
!$omp end do
!$omp workshare
T_ia = -T_ia / 2.d0
!$omp end workshare
!$omp end parallel
write(*, FMT='(A)') "IMCOM: T matrix complete"
write(*, FMT='(A)') "IMCOM: Writing T matrix to "//trim(Tfile)
call imcom_writefits(trim(Tfile), n, m, real(T_ia, 8))
END SUBROUTINE imcom_build_T

!---

SUBROUTINE imcom_build_image
! Uses the T matrices and the input images (stored in Im) to get the output
! H_a image
implicit none
character, parameter :: trans = 'T'
!character, parameter :: notrans = 'U'
real(KIND=8), parameter :: alpha = 1.d0, beta = 0.d0
integer, parameter :: incx = 1, incy = 1
integer :: alstat

write(*, FMT='(A)') "IMCOM: Building output image H"
allocate(H_a(m), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate space for image output vectors -- insufficient memory?"
  stop
end if
H_a = 0.d0
call DGEMV(trans, n, m, alpha, T_ia, n, I_i, incx, beta, H_a, incy)
deallocate(I_i, STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot deallocate space for image output vector I_i"
  stop
end if
allocate(H(n1out, n2out), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate space for output image matrix H -- insufficient memory?"
  stop
end if
H = reshape(H_a, (/ n1out, n2out /))
deallocate(H_a, STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot deallocate space for output image matrix H_a"
  stop
end if
write(*, FMT='(A)') "IMCOM: Saving H vector to "//trim(Hfile)
call imcom_writefits(trim(hfile), n1out, n2out, H)
END SUBROUTINE imcom_build_image

!---


SUBROUTINE imcom_build_U
! Uses the T matrices and the input images (stored in Im) to get the output
! Sigma_a noise image
implicit none
!real(KIND=8), dimension(m) :: ATTplusBT
integer :: a !i, a_i, a_f, nthreads
integer :: alstat

! THIS SHOULD BE REWRITTEN TO NOT DEPEND ON A or B BUT ONLY ON Qij, Pia or Tia.
! ...THAT WOULD ALLOW A & B TO BE DEALLOCATED AFTER THE CALCULATION OF Pia.
!

write(*, FMT='(A)') "IMCOM: Building output leakage map U"
allocate(U_a(m), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for leakage matrix U"
  write(*, FMT='(A)') "IMCOM ERROR: -- too large? Sparse Matrix formalism not yet coded."
  stop
end if
U_a = 0.d0
!$omp parallel do
do a=1, m

  U_a(a) = imcom_U_from_kappa(a, n, K_a(a))

end do
!$omp end parallel do
allocate(U(n1out, n2out), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for output leakage matrix U -- insufficient memory?"
  stop
endif
U = reshape(U_a, (/ n1out, n2out /))
deallocate(U_a, STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot deallocate memory for leakage matrix U_a"
  stop
endif
write(*, FMT='(A)') "IMCOM: Saving U vector to "//trim(Ufile)
call imcom_writefits(trim(Ufile), n1out, n2out, U)
END SUBROUTINE imcom_build_U

!---

SUBROUTINE imcom_build_S
! Uses the T matrices and to get the output S image
implicit none
real(KIND=8), dimension(n, m) :: TN
integer :: alstat, i, a

write(*, FMT='(A)') "IMCOM: Building output noise map S"
! Begin loop, do OPEN MP stuff
!$omp parallel do
do i=1, n

  TN(i, :) = T_ia(i, :) * Ndiag(i)

end do
!$omp end parallel do
allocate(S_a(m), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate space for noise map output matrices -- insufficient memory?"
  stop
end if
S_a = 0.d0
! Not sure whether BLAS or parallel routines are quicker... test?
!$omp parallel workshare
forall(a=1:m) S_a(a) = sum(T_ia(:, a) * TN(:, a))
!$omp end parallel workshare
allocate(S(n1out, n2out), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate space for noise map output matrix -- insufficient memory?"
  stop
end if
S = reshape(S_a, (/ n1out, n2out /))
deallocate(S_a, STAT=alstat)
write(*, FMT='(A)') "IMCOM: Saving S vector to "//trim(Sfile)
call imcom_writefits(trim(Sfile), n1out, n2out, S)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot deallocate space for noise map output matrix -- insufficient memory?"
  stop
end if
END SUBROUTINE imcom_build_S

!---

FUNCTION imcom_A_plus_kNdiag(n_dum, A_in, kNdiag_in)
! Adds (kappa * N_ii) to the diagonals A_ii of an input A matrix
implicit none
integer, intent(IN) :: n_dum
real(KIND=8), dimension(n_dum, n_dum), intent(IN) :: A_in
real(KIND=8), dimension(n_dum), intent(IN) :: kNdiag_in
real(KIND=8), dimension(n_dum, n_dum) :: imcom_A_plus_kNdiag  ! output
integer :: i

! Begin OPEN MP stuff
!$omp parallel
!$omp workshare
imcom_A_plus_kNdiag = A_in
!$omp end workshare
!$omp workshare
forall(i=1: n_dum) imcom_A_plus_kNdiag(i, i) = imcom_A_plus_kNdiag(i, i) &
                                             + kNdiag_in(i)
!$omp end workshare
!$omp end parallel
END FUNCTION imcom_A_plus_kNdiag

!---

FUNCTION imcom_U_from_kappa(a_dum, n_dum, k_dum)
! Thread higher up the food chain.
implicit none
integer, intent(IN) :: a_dum, n_dum
real(KIND=8), intent(IN) :: k_dum
real(KIND=8) :: imcom_U_from_kappa
real(KIND=8), dimension(n_dum) :: vec1s
real(KIND=8), dimension(n_dum) :: denom
real(KIND=8) :: Usum
real(KIND=8), external :: DDOT

! Note that summing vector elements with BLAS requires a dot product with 
! a vector of ones... Is this faster than SUM()?  Test...
vec1s = 1.d0
denom = (L_i(1:n_dum) + k_dum) * (L_i(1:n_dum) + k_dum)
Usum = DDOT(n_dum, (L_i(1:n_dum) + 2.d0 * k_dum) * P2_ia(1:n_dum, a_dum) &
            / denom, 1, vec1s, 1)
imcom_U_from_kappa = (C_a - Usum / 4.d0) / C_a
END FUNCTION imcom_U_from_kappa

!---

FUNCTION imcom_U_from_kappa_test(a_dum, n_dum, k_dum)
! Thread higher up the food chain.
implicit none
integer, intent(IN) :: a_dum, n_dum
real(KIND=8), intent(IN) :: k_dum
real(KIND=8) :: imcom_U_from_kappa_test
real(KIND=8), dimension(n_dum) :: vec1s
real(KIND=8), dimension(n_dum) :: denom, denom2
real(KIND=8) :: Usum, Usum2
real(KIND=8), external :: DDOT

! Note that summing vector elements with BLAS requires a dot product with 
! a vector of ones... Is this faster than SUM()?  Test...
vec1s = 1.d0
denom = (L_i(1:n_dum) + k_dum) * (L_i(1:n_dum) + k_dum)
denom2 = (L_i(n_dum: 1: -1) + k_dum) * (L_i(n_dum: 1: -1) + k_dum)
Usum = DDOT(n_dum, (L_i(1:n_dum) + 2.d0 * k_dum) * P2_ia(1:n_dum, a_dum) &
     / denom, 1, vec1s, 1)
Usum2 = DDOT(n_dum, (L_i(n_dum: 1: -1) + 2.d0 * k_dum) &
      * P2_ia(n_dum: 1: -1, a_dum) / denom, 1, vec1s, 1)
imcom_U_from_kappa_test = (C_a - Usum / 4.d0) / C_a
write(*, FMT='(5E24.16)') C_a, Usum / 4.d0, (C_a - Usum / 4.d0) / C_a, &
                          1.d0 - Usum / 4.d0 / C_a, Usum / 4.d0 / C_a
write(*, FMT='(5E24.16)') C_a, Usum2 / 4.d0, (C_a - Usum2 / 4.d0) / C_a, &
                          1.d0 - Usum2 / 4.d0 / C_a, Usum2 / 4.d0 / C_a
END FUNCTION imcom_U_from_kappa_test

!---

FUNCTION imcom_S_from_kappa(a_dum, n_dum, k_dum)
! Thread higher up the food chain.
implicit none
integer, intent(IN) :: a_dum, n_dum
real(KIND=8), intent(IN) :: k_dum
real(KIND=8) :: imcom_S_from_kappa
real(KIND=8), dimension(n_dum) :: vec1s
real(KIND=8), dimension(n_dum) :: denom
real(KIND=8) :: Ssum
real(KIND=8), external :: DDOT

vec1s = 1.d0
denom = (L_i(1:n_dum) + k_dum) * (L_i(1:n_dum) + k_dum)
Ssum = DDOT(n_dum, Ndiag(1:n_dum) * P2_ia(1:n_dum, a_dum) / denom, 1, vec1s, 1)
imcom_S_from_kappa = Ssum / 4.d0
END FUNCTION imcom_S_from_kappa

!---

END MODULE imcom_matrix
