!
! Module for matrix solvers, building(?) routines etc. used in the IMCOM
! package - calls LAPACK library
!
!---------------------------------------------------

MODULE imcom_matrix

USE imcom_data
USE imcom_io
USE imcom_proc
USE omp_lib

implicit none

CONTAINS

!---

SUBROUTINE imcom_get_A(npoly, saveA, forcebuild)
implicit none
integer, intent(IN) :: npoly, saveA, forcebuild
logical :: Aexists
integer :: n1_tmp, n2_tmp, bit_tmp
integer :: alstat

inquire(FILE=trim(Afile), EXIST=Aexists)
if (Aexists.and.(forcebuild.eq.0)) then
  call imcom_sizefits(trim(Afile), n1_tmp, n2_tmp, bit_tmp)
  if ((n1_tmp.ne.n).or.(n2_tmp.ne.n)) then
    write(*, FMT='(A)') "IMCOM ERROR: Dimensions of Matrix "//trim(Afile)//" do not match <config_file> image dimensions"
    write(*, FMT='(A,I9)') "IMCOM ERROR: Shoud be square of side: ",n
    stop
  endif
  allocate(A_aij(n1_tmp, n2_tmp), STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for A matrix -- too large? Sparse Matrix formalism not yet coded."
    stop
  endif
  A_aij = 0.d0
  write(*, FMT='(A)') "IMCOM: Reading A matrix from "//trim(Afile)
  call imcom_readfits(trim(Afile), n1_tmp, n2_tmp, A_aij)
else
  call imcom_build_Alookup
  call imcom_lookup_A(npoly, saveA)
end if
END SUBROUTINE imcom_get_A

!---

SUBROUTINE imcom_get_B(npoly, saveB, forcebuild)
implicit none
integer, intent(IN) :: npoly, saveB, forcebuild
logical :: Bexists
integer :: n1_tmp, n2_tmp, bit_tmp
integer :: alstat

inquire(FILE=trim(Bfile), EXIST=Bexists)
if (Bexists.and.(forcebuild.eq.0)) then
  call imcom_sizefits(trim(Bfile), n1_tmp, n2_tmp, bit_tmp)
  if ((n1_tmp.ne.n).or.(n2_tmp.ne.m)) then
    write(*, FMT='(A)') "IMCOM ERROR: Dimensions of Matrix "//trim(Bfile)//" do not match <config_file> image dimensions"
    write(*, FMT='(A,2I9)') "IMCOM ERROR: Shoud be rectangle of sides: ",n,m
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
  call imcom_build_Blookup
  call imcom_lookup_B(npoly, saveB)
end if
END SUBROUTINE imcom_get_B

!---

SUBROUTINE imcom_get_T(saveT, forcebuild)
implicit none
integer, intent(IN) :: saveT, forcebuild
character(LEN=1) :: equed
logical :: Texists
integer :: n1_tmp, n2_tmp, bit_tmp
integer :: alstat

inquire(FILE=trim(Tfile), EXIST=Texists)
if (Texists.and.(forcebuild.eq.0)) then
  call imcom_sizefits(trim(Tfile), n1_tmp, n2_tmp, bit_tmp)
  if ((n1_tmp.ne.n).or.(n2_tmp.ne.m)) then
    write(*, FMT='(A)') "IMCOM ERROR: Dimensions of Matrix "//trim(Tfile)//" do not match <config_file> image dimensions"
    write(*, FMT='(A,2I9)') "IMCOM ERROR: Shoud be rectangle of sides: ",n,m
    stop
  endif
  allocate(T_ia(n1_tmp, n2_tmp), STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for T matrix -- too large? Sparse Matrix formalism not yet coded."
    stop
  endif
  T_ia = 0.d0
  write(*, FMT='(A)') "IMCOM: Reading T matrix from "//trim(Tfile)
  call imcom_readfits(trim(Tfile), n1_tmp, n2_tmp, T_ia)
else
! Put the A matrix in packed format, equilibrate if necesssary, calculate
! the Cholesky decomposition and solve!
  call imcom_pack_A_add_kN
  call imcom_equilibrate_A(equed)
  call imcom_chols_Apequ
  call imcom_build_T(saveT, equed)
end if
END SUBROUTINE imcom_get_T

!---

SUBROUTINE imcom_build_A(saveA)
! Build the A matrix from combo.pdf
!
implicit none
integer, intent(IN) :: saveA
real(KIND=8) :: Tx, Ty
real(KIND=8) :: dux, duy  !, d2u
real(KIND=8) :: Delx, Dely, sincfac
integer :: i, j
integer :: iexp, jexp, nperexp
integer :: alstat
! Internal work arrays for storing the contents of the integral
integer :: nmask, imask
real(KIND=8) :: maxfreq2
real(KIND=8), dimension(:), allocatable :: ux_m, uy_m, arg
COMPLEX(KIND=8), dimension(:, :, :), allocatable :: integrand
real(KIND=8), external :: ZDOTU

write(*, FMT='(A)') "IMCOM: Building A matrix"
allocate(A_aij(n, n), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for A matrix -- too large? Sparse Matrix formalism not yet coded."
  write(*, FMT='(A)') "IMCOM ERROR: Or, ensure x_i and y_i are allocated?"
  stop
endif
Tx = psfxscale * real(n1psf, 8)  ! Linear extent of the PSF map, one period
Ty = psfyscale * real(n2psf, 8)  ! in arcsec
dux = 1.d0 / Tx                 ! in cycles (not radians) per arcsec
duy = 1.d0 / Ty                 !
!d2u = dux * duy
! Build masked (much smaller) ux, uy, Gt arrays for the non-zero regions 
! of the MTF, those for |u| < maxfreq...
maxfreq2 = maxfreq * maxfreq
nmask = count(logical((ux * ux + uy * uy).le.maxfreq2, 1))
allocate(ux_m(nmask), uy_m(nmask), arg(nmask), &
         integrand(nmask, nexp, nexp), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for masked work arrays in IMCOM_BUILD_A"
  stop
endif
ux_m = 0.d0
uy_m = 0.d0
integrand = dcmplx(0.d0, 0.d0)
imask = 0
! Define masked ux and uy arrays
do i=1, n1psf

  do j=1, n2psf

    if ((ux(i, j) * ux(i, j) + uy(i, j) * uy(i, j)).le.maxfreq2) then
      imask = imask + 1
      ux_m(imask) = ux(i, j)
      uy_m(imask) = uy(i, j)
    endif

  end do

end do
! Begin with OPEN MP stuff (note private specifications v. important)
! ...note that only the first loop is split up, so each thread takes a 
! different exposure...
!$omp parallel
!$omp do private(imask, iexp, jexp, i, j)
do iexp=1, nexp

  do jexp=1, nexp

    imask = 0
    do i=1, n1psf

      do j=1, n2psf

        if ((ux(i, j) * ux(i, j) + uy(i, j) * uy(i, j)).le.maxfreq2) then
          imask = imask + 1
          integrand(imask, iexp, jexp) = imcom_upsilon(ux(i, j), uy(i, j), 1) &
                                       * dconjg(Gt_rot(i, j, iexp))           &
                                       * Gt_rot(i, j, jexp)
        endif

      end do

    end do

  end do

end do
!$omp end do
!$omp end parallel
! Initialize
A_aij = 0.d0
arg = 0.d0
nperexp = n / nexp
iexp = 0
jexp = 0
!write(*, FMT='(A)') "IMCOM: Calculating A elements"
! Begin loop, do OPEN MP stuff (note private specifications v. important)
!$omp parallel
!$omp do schedule(dynamic, 32) private(iexp, jexp, Delx, Dely, sincfac, arg)
do i=1, n

  iexp = 1 + (i - 1) / nperexp
  do j=i, n

    jexp = 1 + (j - 1) / nperexp 
    Delx = x_i(j) - x_i(i)
    Dely = y_i(j) - y_i(i)
    sincfac = 1.d0 !imcom_dsinc2(dux * Delx) * imcom_dsinc2(duy * Dely)
    arg = Delx * ux_m + Dely * uy_m
    A_aij(i, j) = sincfac * sum(real(zexp(dcmplx(0.d0, 2.d0 * pi  * arg)) &
                          * integrand(:, iexp, jexp), 8))

  end do
  if (i.eq.int(0.2 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 20% complete"
  if (i.eq.int(0.4 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 40% complete"
  if (i.eq.int(0.6 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 60% complete"
  if (i.eq.int(0.8 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 80% complete"

end do
!$omp end do
!$omp do schedule(dynamic, 16)
do i=2, n

  forall(j=1:i-1) A_aij(i, j) = A_aij(j, i)

end do
!$omp end do
! Then do a normalization "correctly", including Ns in the FFT
!$omp workshare
A_aij = A_aij / real(n1psf, 8) / real(n2psf, 8)
!$omp end workshare
!$omp end parallel
write(*, FMT='(A)') "IMCOM: A matrix complete"
if (saveA.eq.1) then 
  write(*, FMT='(A)') "IMCOM: Saving A matrix to "//trim(Afile)
  call imcom_writefits(trim(Afile), n, n, real(A_aij, 8))
end if
deallocate(ux_m, uy_m, arg, integrand, STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot deallocate memory for masked work arrays in IMCOM_BUILD_A"
  stop
endif
END SUBROUTINE imcom_build_A

!---

SUBROUTINE imcom_lookup_A(npoly, saveA)
! Calculate the A matrix using the lookup table and polynomial interpolation at degree npoly
implicit none
integer, intent(IN) :: npoly, saveA
real(KIND=8) :: Delx, Dely
integer :: iexp, jexp, i, j, alstat, nperexp

write(*, FMT='(A)') "IMCOM: Building A matrix from lookup table"
allocate(A_aij(n, n), STAT=alstat)
! Initialize
A_aij = 0.d0
nperexp = n / nexp
iexp = 0
jexp = 0
!$omp parallel
!$omp do schedule(dynamic, 32) private(iexp, jexp, Delx, Dely)
do i=1, n

  iexp = 1 + (i - 1) / nperexp
  do j=i, n

    jexp = 1 + (j - 1) / nperexp 
    Delx = (x_i(j) - x_i(i)) / psfxscale
    Dely = (y_i(j) - y_i(i)) / psfyscale
    A_aij(i, j) = imcom_interp_lookup(Alookup(:, :, iexp, jexp), Delx, Dely, &
                                      npoly)

  end do
  if (i.eq.int(0.2 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 20% complete"
  if (i.eq.int(0.4 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 40% complete"
  if (i.eq.int(0.6 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 60% complete"
  if (i.eq.int(0.8 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 80% complete"

end do
!$omp end do
!$omp do schedule(dynamic, 16)
do i=2, n

  forall(j=1:i-1) A_aij(i, j) = A_aij(j, i)

end do
!$omp end do
!$omp end parallel
write(*, FMT='(A)') "IMCOM: A matrix complete"
if (saveA.eq.1) then 
  write(*, FMT='(A)') "IMCOM: Saving A matrix to "//trim(Afile)
  call imcom_writefits(trim(Afile), n, n, real(A_aij, 8))
end if
END SUBROUTINE imcom_lookup_A

!---

SUBROUTINE imcom_lookup_B(npoly, saveB)
! Calculate the B matrix using the lookup table and polynomial interpolation at degree npoly
implicit none
integer, intent(IN) :: npoly, saveB
real(KIND=8) :: Delx, Dely
integer :: iexp, i, a, alstat, nperexp

write(*, FMT='(A)') "IMCOM: Building B vectors from lookup table"
allocate(B_ia(n, m), STAT=alstat)
! Initialize
B_ia = 0.d0
nperexp = n / nexp
iexp = 0
!$omp parallel
!$omp do schedule(dynamic, 32) private(iexp, Delx, Dely)
do i=1, n

  iexp = 1 + (i - 1) / nperexp
  do a=1, m

    Delx = (x_i(i) - X_a(a)) / psfxscale
    Dely = (y_i(i) - Y_a(a)) / psfyscale
    B_ia(i, a) = imcom_interp_lookup(Blookup(:, :, iexp), Delx, Dely, npoly)

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
END SUBROUTINE imcom_lookup_B

!---

SUBROUTINE imcom_build_Alookup
! Build the A matrix lookup table using Fourier transforms
!
implicit none
complex(KIND=8), dimension(:, :), allocatable ::  A_tmp, ufunc
integer :: alstat, iexp, jexp

write(*, FMT='(A)') "IMCOM: Fourier transforming A matrix lookup table"
allocate(Alookup(n1psf, n2psf, nexp, nexp), A_tmp(n1psf, n2psf), &
         ufunc(n1psf, n2psf), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot allocate memory for A lookup tables"
  stop
endif
do iexp=1, nexp

  do jexp=1, nexp

    ufunc = dconjg(Gt_rot(:, :, iexp)) * Gt_rot(:, :, jexp)         &
                                       * imcom_upsilon2D(ux, uy, 1)
    call imcom_invft_c2c(n1psf, n2psf, ufunc, A_tmp, 0)
    Alookup(:, :, iexp, jexp) = cshift(cshift(real(A_tmp, 8), n2psf / 2, 2), &
                                       n1psf / 2, 1)
! Origin thus located at centre of pixel (n1psf / 2 + 1, n2psf / 2 + 1)

  end do

end do
Alookup = Alookup / real(n1psf, 8) / real(n2psf, 8)
deallocate(A_tmp, ufunc, STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot deallocate temporary memory for A lookup tables"
  stop
endif
END SUBROUTINE imcom_build_Alookup

!---

SUBROUTINE imcom_build_Blookup
! Build the B matrix lookup table using Fourier transforms
!
implicit none
complex(KIND=8), dimension(:, :), allocatable ::  B_tmp, ufunc
integer :: alstat, iexp

write(*, FMT='(A)') "IMCOM: Fourier transforming B matrix lookup tables"
allocate(Blookup(n1psf, n2psf, nexp), B_tmp(n1psf, n2psf), &
         ufunc(n1psf, n2psf), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot allocate memory for B lookup tables"
  stop
endif
do iexp=1, nexp

  ufunc = dconjg(Gammat(:, :)) * Gt_rot(:, :, iexp) * imcom_upsilon2D(ux, uy, 1)
  call imcom_invft_c2c(n1psf, n2psf, ufunc, B_tmp, 0)
  Blookup(:, :, iexp) = -2.d0 * cshift(cshift(real(B_tmp, 8), n2psf / 2, 2), &
                                       n1psf / 2, 1)
! Origin thus located at centre of pixel (n1psf / 2 + 1, n2psf / 2 + 1)

end do
Blookup = Blookup / real(n1psf, 8) / real(n2psf, 8)
deallocate(B_tmp, ufunc, STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot deallocate temporary memory for B lookup tables"
  stop
endif
END SUBROUTINE imcom_build_Blookup

!---

FUNCTION imcom_interp_lookup(lookup, xarg, yarg, npoly_in)
implicit none
real(KIND=8), dimension(:, :), intent(IN) :: lookup
real(KIND=8), intent(IN) :: xarg, yarg  ! x, y in lookup table pixel units (not arcsec!)
integer, intent(IN) :: npoly_in
real(KIND=8) :: imcom_interp_lookup
real(KIND=8) :: icent, jcent, iarg, jarg, f, df
integer :: npoly, npm
integer :: ia_floor, ia_ceiling, ja_floor, ja_ceiling
integer :: imin, imax, jmin, jmax, ii, jj

npoly = npoly_in
if (mod(npoly, 2).eq.0) then
  write(*, FMT='(A)') "IMCOM WARNING: Even-order polynomial specified for IMCOM_INTERP_LOOKUP..."
  write(*, FMT='(A)') "IMCOM WARNING: Using NPOLY-1 for interpolation (e.g. linear, cubic, quintic...)"
  npoly = npoly - 1
end if
npm = npoly / 2
icent = real(n1psf / 2 + 1, 8)
jcent = real(n2psf / 2 + 1, 8)
iarg = icent + xarg
jarg = jcent + yarg
ia_floor = floor(iarg)
ja_floor = floor(jarg)
ia_ceiling = ceiling(iarg)
ja_ceiling = ceiling(jarg)
if ((ia_floor.ge.1).and.(ja_floor.ge.1).and.(ia_ceiling.le.n1psf) &
                                       .and.(ja_ceiling.le.n2psf)) then
 imin = max(1, ia_floor - npm)
 jmin = max(1, ja_floor - npm)
 imax = min(n1psf, ia_ceiling + npm)
 jmax = min(n2psf, ja_ceiling + npm)
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

C_a = sum(real(Gammat * dconjg(Gammat) * imcom_upsilon2D(ux, uy, 1), 8)) &
    / real(n1psf, 8) / real(n2psf, 8)
write(*, FMT='(A, E21.15)') "IMCOM: C_a = ", C_a
END SUBROUTINE imcom_calc_C

!---

SUBROUTINE imcom_build_kN
! Make the (diagonals only) noise and kappa * N trace vectors
implicit none
integer :: i, n_i, n_f, alstat

allocate(kN(n), Ndiag(n), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for kN matrix"
  stop
endif
kN = 0.d0
Ndiag = 0.d0
do i=1, nexp

  n_i = 1 + (i - 1) * (n / nexp)
  if (i.ne.nexp) then
    n_f = i * (n / nexp)
  else
    n_f = n
  end if
  Ndiag(n_i:n_f) = noise(i)
  kN(n_i:n_f) = kappa * noise(i)

end do
END SUBROUTINE imcom_build_kN

!---

SUBROUTINE imcom_build_B(saveB)
! Build the B matrix from combo.pdf
!
implicit none
integer, intent(IN) :: saveB
real(KIND=8) :: Tx, Ty
!real(KIND=8) :: Txcut, Tycut
real(KIND=8) :: dux, duy  !, d2u
real(KIND=8) :: Delx, Dely, sincfac
integer :: a, i, j
integer :: iexp, nperexp
integer :: alstat
! Internal work arrays for storing the contents of the integral
integer :: nmask, imask
real(KIND=8) :: maxfreq2
real(KIND=8), dimension(:), allocatable :: ux_m, uy_m
real(KIND=8), dimension(:), allocatable :: arg
COMPLEX(KIND=8), dimension(:, :), allocatable :: integrand

allocate(B_ia(n, m), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for B matrix -- too large? Sparse Matrix formalism not yet coded."
  stop
endif
write(*, FMT='(A)') "IMCOM: Building B vectors"
Tx = psfxscale * real(n1psf, 8)  ! Linear extent of the PSF map, one period
Ty = psfyscale * real(n2psf, 8)  ! in arcsec
dux = 1.d0 / Tx                 ! in cycles (not radians) per arcsec
duy = 1.d0 / Ty                 !
!d2u = dux * duy
! Build masked (much smaller) ux, uy, Gt arrays for the non-zero regions 
! of the MTF, those for |u| < maxfreq...
maxfreq2 = maxfreq * maxfreq
nmask = count(logical((ux * ux + uy * uy).le.maxfreq2, 1))
allocate(ux_m(nmask), uy_m(nmask), arg(nmask), integrand(nmask, nexp), &
         STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for masked work arrays in IMCOM_BUILD_B"
  stop
endif
ux_m = 0.d0
uy_m = 0.d0
integrand = dcmplx(0.d0, 0.d0)
! Define masked ux and uy arrays
imask = 0
do i=1, n1psf

  do j=1, n2psf

    if ((ux(i, j) * ux(i, j) + uy(i, j) * uy(i, j)).le.maxfreq2) then
      imask = imask + 1
      ux_m(imask) = ux(i, j)
      uy_m(imask) = uy(i, j)
    endif

  end do

end do
! Begin loop, do OPEN MP stuff (note private specifications v. important)
! ...note that only the first loop is split up, so each thread takes a 
! different exposure...
!$omp parallel
!$omp do private(imask)
do iexp=1,nexp

  imask = 0
  do i=1, n1psf

    do j=1, n2psf

      if ((ux(i, j) * ux(i, j) + uy(i, j) * uy(i, j)).le.maxfreq2) then
        imask = imask + 1
        integrand(imask, iexp) = imcom_upsilon(ux(i, j), uy(i, j), 1) &
                               * Gt_rot(i, j, iexp) * dconjg(Gammat(i, j))
      endif

    end do

  end do

end do
!$omp end do
!$omp end parallel
! Initialize
B_ia = 0.d0
arg = 0.d0
nperexp = n / nexp
iexp = 0
! Begin loop, do OPEN MP stuff
!$omp parallel
!$omp do schedule(dynamic, 32) private(iexp, arg, Delx, Dely, sincfac)
do i=1, n

  iexp = 1 + (i - 1) / nperexp
  do a=1, m

    Delx = x_i(i) - X_a(a)
    Dely = y_i(i) - Y_a(a)
    sincfac = 1.d0 !imcom_dsinc2(dux * Delx) * imcom_dsinc2(duy * Dely)
    arg = (ux_m * Delx + uy_m * Dely)
    B_ia(i, a) = sincfac * &
                 sum(real(zexp(dcmplx(0.d0, 2.d0 * pi * arg)) &
                        * integrand(:, iexp), 8))

  end do
  if (i.eq.int(0.2 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 20% complete"
  if (i.eq.int(0.4 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 40% complete"
  if (i.eq.int(0.6 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 60% complete"
  if (i.eq.int(0.8 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 80% complete"

end do
!$omp end do
! Then do a normalization "correctly", including Ns in the FFT
!$omp workshare
B_ia = -B_ia * 2.d0 / real(n1psf, 8) / real(n2psf, 8)
!$omp end workshare
!$omp end parallel
write(*, FMT='(A)') "IMCOM: All B vectors complete"
if (saveB.eq.1) then 
  write(*, FMT='(A)') "IMCOM: Saving B matrix to "//trim(Bfile)
  call imcom_writefits(trim(Bfile), n, m, real(B_ia, 8))
end if
deallocate(ux_m, uy_m, arg, integrand, STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot deallocate memory for masked work arrays in IMCOM_BUILD_B"
  stop
endif
END SUBROUTINE imcom_build_B

!---

SUBROUTINE imcom_pack_A
! Pack the symmetric positive definite matrix Aij + kNij in upper 
! triangular form, column-wise, as required by LAPACK
implicit none
integer :: i, j, ip, alstat

allocate(A_pack(n * (n + 1) / 2), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate space for packed A matrix -- insufficient memory?"
  stop
end if
ip = 0   ! starting location in A_packed
! Begin loop, do OPEN MP stuff
!$omp parallel
!$omp do schedule(dynamic, 32) private(ip)
do j=1, n

  do i=1, j

    ip = (j - 1) * j / 2 + i
    A_pack(ip) = A_aij(i, j)

  end do

end do
!$omp end do
!$omp end parallel
END SUBROUTINE imcom_pack_A

!---

SUBROUTINE imcom_pack_A_add_kN
! Pack the symmetric positive definite matrix Aij + kNij in upper 
! triangular form, column-wise, as required by LAPACK
implicit none
integer :: i, j, ip, alstat

allocate(ApkN_aij(n, n), A_pack(n * (n + 1) / 2), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate space for packed A matrices -- insufficient memory?"
  stop
end if
ApkN_aij = A_aij   ! kN on diagnonal added in jsut one second...
A_pack = 0.d0
ip = 0   ! starting location in A_packed
! Begin loop, do OPEN MP stuff
!$omp parallel
!$omp workshare
forall(i = 1: n) ApkN_aij(i, i) = A_aij(i, i) + kN(i)
!$omp end workshare
!$omp do schedule(dynamic, 32) private(ip)
do j=1, n

  do i=1, j

    ip = (j - 1) * j / 2 + i
    A_pack(ip) = ApkN_aij(i, j)

  end do

end do
!$omp end do
!$omp end parallel

!deallocate(A_aij, STAT=dealstat)
!if (dealstat.ne.0) then
!  write(*, FMT='(A)') "IMCOM ERROR: Cannot deallocate space for A_aij matrix"
!  stop
!end if
END SUBROUTINE imcom_pack_A_add_kN

!---

SUBROUTINE imcom_equilibrate_A(equed_out)
! Calls LAPACK routines to equilibrate A (if necessary, whether action 
! was eprformed is specified by EQUED_OUT
implicit none
character(LEN=1), intent(OUT) :: equed_out
character(LEN=1), parameter :: uplo = 'U'
real(KIND=8) :: Amax, scond
integer :: info, alstat
!real(KIND=8), dimension(:, :), allocatable :: A_test

allocate(A_pequ(n * (n + 1) / 2), scaling(n), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate space for equilibrating A matrix -- insufficient memory?"
  stop
end if
scaling = 0.d0
scond = 0.d0
A_pequ = A_pack
call DPPEQU(uplo, n, A_pequ, scaling, scond, Amax, info)
if (info.ne.0) then
  if (info.lt.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Equilibration failed!"
    write(*, FMT='(A, I12)') "IMCOM ERROR: Illegal value in A_pack argument i = ",-info
    stop
  else
    write(*, FMT='(A)') "IMCOM ERROR: Equilibration failed!"
    write(*, FMT='(A, I12)')"IMCOM ERROR: Non-positive i-th diagonal element in A_pack: i = ",info
    stop
  end if
endif
! Use these computed scaling factors to build A_pequ
call DLAQSP(uplo, n, A_pequ, scaling, scond, Amax, equed_out)
END SUBROUTINE imcom_equilibrate_A

!---

SUBROUTINE imcom_chols_Apequ
! Calls LAPACK routines to compute the Cholesky decomposition of the
! packed-storage (upper triangular form, see imcom_pack_A), equilibrated, 
! real, symmetric positive definite matrix A
implicit none
character(LEN=1), parameter :: uplo = 'U' ! Upper triangular packing
integer :: info, alstat

allocate(A_chls(n * (n + 1) / 2), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate space for Cholesky-factorized, packed A matrix -- insufficient memory?"
  stop
end if
A_chls = A_pequ
write(*, FMT='(A)') "IMCOM: Factorizing A matrix"
call DPPTRF(uplo, n, A_chls, info)
if (info.ne.0) then
  if (info.lt.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cholesky factorization failed!"
    write(*, FMT='(A, I12)') "IMCOM ERROR: Illegal value in A_packed argument i = ",-info
    stop
  else
    write(*, FMT='(A)') "IMCOM ERROR: Cholesky factorization failed!"
    write(*, FMT='(A, I12)') "IMCOM ERROR: Non-positive definite A_packed in leading minor order i = ", info
    stop
  end if
endif
END SUBROUTINE imcom_chols_Apequ

!---

SUBROUTINE imcom_build_T(saveT, equed_A_dum)
! Solves the system of combo.pdf eq 16) using multithreaded LAPACK solving routines
implicit none
integer, intent(IN) :: saveT
character(LEN=1), intent(IN) :: equed_A_dum
integer :: nthmax, a_i, a_f, i
integer :: alstat

! Ready the solution matrix
allocate(T_ia(n, m), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for solution matrix T_ia"
  write(*, FMT='(A)') "IMCOM ERROR: -- too large? Sparse Matrix formalism not yet coded."
  stop
end if
nthmax = 2!OMP_GET_MAX_THREADS()
T_ia = 0.d0
write(*, FMT='(A)') "IMCOM: Solving T matrix"
! Begin loop, do OPEN MP stuff
!$omp parallel
!$omp do schedule(dynamic, 1) private(a_i, a_f)
do i=1, nthmax

  a_i = 1 + (i - 1) * (m / nthmax)
  if (i.ne.nthmax) then
    a_f = i * (m / nthmax)
  else
    a_f = m
  end if
  call imcom_solve_ATB(n, (1 + a_f - a_i), A_pack, A_chls, equed_A_dum, &
                       scaling, -0.5d0 * B_ia(:, a_i:a_f), T_ia(:, a_i:a_f))

end do
!$omp end do
!$omp end parallel
if (saveT.eq.1) then
  write(*, FMT='(A)') "IMCOM: Saving T matrix to "//trim(Tfile)
  call imcom_writefits(trim(Tfile), n, m, T_ia)
endif
END SUBROUTINE imcom_build_T

!---

SUBROUTINE imcom_solve_ATB(n_dum, m_dum, Ap_dum, Apch_dum, equed_dum, s_dum, &
                           B_dum, T_dum)
!  Solves the AT = B equation for supplied matrices: note these can be 
! sub-matrices of the full thing, to take advantage of parallel processing...
!  Uses the LAPACK solver for positive definite matrices, DPPSVX
!
implicit none
integer, intent(IN) :: n_dum, m_dum
real(KIND=8), dimension(n_dum * (n_dum + 1) / 2), intent(IN) :: Ap_dum
real(KIND=8), dimension(n_dum * (n_dum + 1) / 2), intent(IN) :: Apch_dum
real(KIND=8), dimension(n_dum), intent(IN) :: s_dum
character(LEN=1), intent(IN) :: equed_dum ! Equilibration done ?
real(KIND=8), dimension(n_dum, m_dum), intent(IN) :: B_dum
real(KIND=8), dimension(n_dum, m_dum), intent(OUT) :: T_dum
character(LEN=1), parameter :: fact = 'F'  ! Factorization already done
character(LEN=1), parameter :: uplo = 'U'  ! Upper triangular form
real(KIND=8) :: rcond
real(KIND=8), dimension(m_dum) :: ferr, berr
real(KIND=8), dimension(3 * n_dum) :: workspace
integer, dimension(n_dum) :: iworkspace
integer :: info

ferr = 0.d0
berr = 0.d0
workspace = 0.d0
iworkspace = 0
call DPPSVX(fact, uplo, n_dum, m_dum, Ap_dum, Apch_dum, equed_dum, s_dum, &
            B_dum, n_dum, T_dum, n_dum, rcond, ferr, berr, workspace,     & 
            iworkspace, info)
if (info.ne.0) then
  if (info.lt.0) then
    write(*, FMT='(A)') "IMCOM ERROR: LAPACK solver failed!"
    write(*, FMT='(A, I12)') "IMCOM ERROR: Illegal value in A_packed argument i = ",-info
    stop
  else
    if (info.le.n_dum) then
      write(*, FMT='(A)') "IMCOM ERROR: LAPACK solver failed!"
      write(*, FMT='(A)') "IMCOM ERROR: Non-positive definite A_packed in leading minor order i = ",info
      stop
    else if (info.eq.(n_dum + 1)) then
      write(*, FMT='(A, I12)') "IMCOM WARNING: LAPACK found Ap_dum to be singular to within working precision"
    endif
    write(*, FMT='(A, I12)') "IMCOM ERROR: LAPACK solver failed but produced weird INFO code = ", info
    stop
  end if
endif
END SUBROUTINE imcom_solve_ATB

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
allocate(H_a(m), I_i(n), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate space for image output vectors -- insufficient memory?"
  stop
end if
I_i = 0.d0
I_i = reshape(Im, (/ n /))
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

SUBROUTINE imcom_build_noise
! Uses the T matrices and the input images (stored in Im) to get the output
! H_a image
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
END SUBROUTINE imcom_build_noise

!---

SUBROUTINE imcom_build_U
! Uses the T matrices and the input images (stored in Im) to get the output
! Sigma_a noise image
implicit none
real(KIND=8), dimension(m) :: ATTplusBT
integer :: nthmax, i, a_i, a_f
integer :: alstat

write(*, FMT='(A)') "IMCOM: Building output leakage map U"
! C must already be built!
! Then the ATT + BT part
nthmax = 1! OMP_GET_MAX_THREADS()   ! Note that, while testing, the routine imcom_calc_ATTpBT uses I/O so should not be multi-threaded
ATTplusBT = 0.d0
! Ready the solution matrix
allocate(U_a(m), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for leakage matrix U"
  write(*, FMT='(A)') "IMCOM ERROR: -- too large? Sparse Matrix formalism not yet coded."
  stop
end if
U_a = 0.d0
! Begin loop, do OPEN MP stuff
!$omp parallel
!$omp do schedule(dynamic, 1) private(a_i, a_f)
do i=1, nthmax

  a_i = 1 + (i - 1) * (m / nthmax)
  if (i.ne.nthmax) then
    a_f = i * (m / nthmax)
  else
    a_f = m
  end if
  call imcom_calc_ATTpBT(n, (1 + a_f - a_i), A_aij, B_ia(:, a_i:a_f), &
                         T_ia(:, a_i:a_f), U_a(a_i:a_f))

end do
!$omp end do
!$omp end parallel
U_a = U_a + C_a
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
call imcom_writefits(trim(Ufile), n1out, n2out, U / C_a)
END SUBROUTINE imcom_build_U

!---

SUBROUTINE imcom_calc_U
! Uses the T matrices and the input images (stored in Im) to get the output
! Sigma_a noise image
implicit none

real(KIND=8) :: Tx, Ty
!real(KIND=8) :: Txcut, Tycut
real(KIND=8) :: dux, duy !, d2u
real(KIND=8) :: sincfac
integer :: a, i, j
integer :: iexp, nperexp
integer :: alstat
! Internal work arrays for storing the contents of the integral
integer :: nmask, imask
real(KIND=8) :: maxfreq2
real(KIND=8), dimension(:), allocatable :: ux_m, uy_m
real(KIND=8), dimension(:), allocatable :: argGam, argG
COMPLEX(KIND=8), dimension(:, :), allocatable :: integrandG
COMPLEX(KIND=8), dimension(:), allocatable :: integrandGam
real(KIND=8), dimension(:), allocatable :: integrand
COMPLEX(KIND=8), dimension(:), allocatable :: Gterm, Gamterm

write(*, FMT='(A)') "IMCOM: Calculating U image by force"
Tx = psfxscale * real(n1psf, 8)  ! Linear extent of the PSF map, one period
Ty = psfyscale * real(n2psf, 8)  ! in arcsec
dux = 1.d0 / Tx                 ! in cycles (not radians) per arcsec
duy = 1.d0 / Ty                 !
!d2u = dux * duy
! Build masked (much smaller) ux, uy, Gt arrays for the non-zero regions 
! of the MTF, those for |u| < maxfreq...
maxfreq2 = maxfreq * maxfreq
nmask = count(logical((ux * ux + uy * uy).le.maxfreq2, 1))
allocate(ux_m(nmask), uy_m(nmask), argG(nmask), argGam(nmask),           &
         integrand(nmask), integrandG(nmask, nexp), integrandGam(nmask), &
         Gterm(nmask), Gamterm(nmask), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for masked work arrays in IMCOM_BUILD_B"
  stop
endif
i = 0
j = 0
ux_m = 0.d0
uy_m = 0.d0
integrand = 0.d0
!print *, real(size(integrandG), 8) * 8.d0 ! size in bytes
integrandG = complex(0.d0, 0.d0)
integrandGam = complex(0.d0, 0.d0)
! Begin loop, do ux,uy then do OPEN MP stuff (note private specifications v. important)
! ...note that only the first loop is split up, so each thread takes a 
! different exposure...
imask = 0
do i=1, n1psf

  do j=1, n2psf

    if ((ux(i, j) * ux(i, j) + uy(i, j) * uy(i, j)).le.maxfreq2) then
      imask = imask + 1
      ux_m(imask) = ux(i, j)
      uy_m(imask) = uy(i, j)
      integrandGam(imask) = dconjg(Gammat(i, j))
    endif

  end do

end do
nperexp = n / nexp
iexp = 0
!$omp parallel
!$omp do private(imask, iexp, i, j)
do iexp=1, nexp

  imask = 0
!  iexp = 1 + (in - 1) / nperexp
  do i=1, n1psf

    do j=1, n2psf

      if ((ux(i, j) * ux(i, j) + uy(i, j) * uy(i, j)).le.maxfreq2) then
        imask = imask + 1
        integrandG(imask, iexp) = dconjg(Gt_rot(i, j, iexp))
      endif

    end do

  end do

end do
!$omp end do
!$omp end parallel
! Ready the solution matrix
allocate(U_calc_a(m), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for leakage matrix U"
  write(*, FMT='(A)') "IMCOM ERROR: -- too large? Sparse Matrix formalism not yet coded."
  stop
end if
! Initialize stuff
U_calc_a = 0.d0
argG = 0.d0
argGam = 0.d0
nperexp = n / nexp
iexp = 0
! Begin loop, do OPEN MP stuff
!$omp parallel
!$omp do schedule(dynamic, 32) private(iexp, i, argG, argGam, Gterm, Gamterm)
do a=1, m

  sincfac = 1.d0
  argGam = -(ux_m * X_a(a) + uy_m * Y_a(a))
  Gamterm = zexp(dcmplx(0.d0, 2.d0 * pi * argGam)) * integrandGam
  Gterm = dcmplx(0.d0, 0.d0)
  do i=1, n

    iexp = 1 + (i - 1) / nperexp
    argG = -(ux_m * x_i(i) + uy_m * y_i(i))
    Gterm = Gterm + dcmplx(T_ia(i, a), 0.d0)             &
                  * zexp(dcmplx(0.d0, 2.d0 * pi * argG)) &
                  * integrandG(:, iexp)

  end do
  integrand = real(dconjg(Gterm - Gamterm) * (Gterm - Gamterm) &
                  * imcom_upsilon1D(ux_m, uy_m, 1), 8)
  U_calc_a(a) =  sincfac * sum(integrand)
  if (a.eq.int(0.2 * real(m, 4))) write(*, FMT='(A)') "IMCOM: 20% complete"
  if (a.eq.int(0.4 * real(m, 4))) write(*, FMT='(A)') "IMCOM: 40% complete"
  if (a.eq.int(0.6 * real(m, 4))) write(*, FMT='(A)') "IMCOM: 60% complete"
  if (a.eq.int(0.8 * real(m, 4))) write(*, FMT='(A)') "IMCOM: 80% complete"

end do
!$omp end do
! Then do a normalization "correctly", including Ns in the FFT
!$omp workshare
U_calc_a = U_calc_a / real(n1psf, 8) / real(n2psf, 8)
!$omp end workshare
!$omp end parallel
write(*, FMT='(A)') "IMCOM: U calculation complete"
allocate(U_calc(n1out, n2out), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for output leakage matrix U -- insufficient memory?"
  stop
endif
U_calc = reshape(U_calc_a, (/ n1out, n2out /))
deallocate(U_calc_a, STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot deallocate memory for leakage matrix U_a"
  stop
endif
write(*, FMT='(A)') "IMCOM: Saving U vector to "//trim(Ucalcfile)
call imcom_writefits(trim(Ucalcfile), n1out, n2out, U_calc / C_a)
END SUBROUTINE imcom_calc_U

!---

SUBROUTINE imcom_calc_ATTpBT(n_dum, m_dum, A_dum, B_dum, T_dum, U_dum)

implicit none
integer, intent(IN) :: n_dum, m_dum
real(KIND=8), dimension(n_dum, n_dum), intent(IN) :: A_dum
real(KIND=8), dimension(n_dum, m_dum), intent(IN) :: B_dum, T_dum
real(KIND=8), dimension(m_dum), intent(OUT) :: U_dum
real(KIND=8), dimension(m_dum) :: ATT, BT
character, parameter :: side= 'L', uplo = 'U'
real(KIND=8), parameter :: alpha = 1.d0, beta = 0.d0
real(KIND=8), dimension(n_dum, m_dum) :: AT_ia
integer :: a

AT_ia = 0.d0
call DSYMM(side, uplo, n_dum, m_dum, alpha, A_dum, n_dum, T_dum, n_dum, &
           beta, AT_ia, n_dum)
forall(a=1:m_dum) ATT(a) = sum(T_dum(:, a) * AT_ia(:, a))
forall(a=1:m_dum) BT(a) = sum(T_dum(:, a) * B_dum(:, a))
!call imcom_writefits('ATT_test.fits', n1out, n2out, reshape(ATT, (/n1out, n2out/)))
!call imcom_writefits('BT_test.fits', n1out, n2out, reshape(BT, (/n1out, n2out/)))
forall(a=1:m_dum) U_dum(a) = sum(T_dum(:, a) * (AT_ia(:, a) + B_dum(:, a)))
!call imcom_writefits('ATTpBT_test.fits', n1out, n2out, reshape(U_dum, (/n1out, n2out/)))
END SUBROUTINE imcom_calc_ATTpBT

!---

SUBROUTINE imcom_calc_C_maxfreq
! Calculate C_a as in combo.pdf eq 15) 
implicit none
! Internal work arrays for storing the contents of the integral
integer :: nmask, imask
real(KIND=8) :: maxfreq2
real(KIND=8), dimension(:), allocatable :: ux_m, uy_m
COMPLEX(KIND=8), dimension(:), allocatable :: integrand
real(KIND=8) :: Tx, Ty, dux, duy
integer :: alstat, i, j

Tx = psfxscale * real(n1psf, 8)  ! Linear extent of the PSF map, one period
Ty = psfyscale * real(n2psf, 8)  ! in arcsec
dux = 1.d0 / Tx                  ! in cycles (not radians) per arcsec
duy = 1.d0 / Ty                  !
!d2u = dux * duy

maxfreq2 = maxfreq * maxfreq
nmask = count(logical((ux * ux + uy * uy).le.maxfreq2 * 100.d0, 1))
allocate(ux_m(nmask), uy_m(nmask), integrand(nmask), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for masked work arrays in IMCOM_CALC_C"
  stop
endif
imask = 0
do i=1, n1psf

  do j=1, n2psf

    if ((ux(i, j) * ux(i, j) + uy(i, j) * uy(i, j)).le.maxfreq2) then
      imask = imask + 1
      ux_m(imask) = ux(i, j)
      uy_m(imask) = uy(i, j)
      integrand(imask) = imcom_upsilon(ux(i, j), uy(i, j), 1) &
                       * Gammat(i, j) * dconjg(Gammat(i, j))
     endif

  end do

end do
C_a = sum(real(integrand, 8)) / real(n1psf, 8) &
                              / real(n2psf, 8)
write(*, FMT='(A, E21.15)') "IMCOM: C_a = ", C_a
END SUBROUTINE imcom_calc_C_maxfreq

!---

END MODULE imcom_matrix
