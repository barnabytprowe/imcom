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
    write(*, FMT='(A,I9)') "IMCOM ERROR: Shoud be square of side: ",n
    stop
  endif
  allocate(A0_aij(n1_tmp, n2_tmp), STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for A matrix -- too large? Sparse Matrix formalism not yet coded."
    stop
  endif
  A0_aij = 0.d0
  write(*, FMT='(A)') "IMCOM: Reading A matrix from "//trim(Afile)
  call imcom_readfits(trim(Afile), n1_tmp, n2_tmp, A0_aij)
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
  call imcom_build_Blookup(npad)
  call imcom_lookup_B(npoly, npad, saveB)
end if
END SUBROUTINE imcom_get_B

!---

SUBROUTINE imcom_lookup_A(npoly, npad, saveA)
! Calculate the A matrix using the lookup table and polynomial interpolation at degree npoly
implicit none
integer, intent(IN) :: npoly, npad, saveA
real(KIND=8) :: Delx, Dely
integer :: iexp, jexp, i, j, alstat, nperexp

write(*, FMT='(A)') "IMCOM: Building A matrix from lookup table"
allocate(A0_aij(n, n), STAT=alstat)
! Initialize
A0_aij = 0.d0
nperexp = n / nexp
iexp = 0
jexp = 0
!$omp parallel
!$omp do schedule(dynamic, 32) private(iexp, jexp, Delx, Dely)
do i=1, n

  iexp = 1 + (i - 1) / nperexp
  do j=i, n

    jexp = 1 + (j - 1) / nperexp 
    Delx = real(npad, 8) * (x_i(j) - x_i(i)) / psfxscale
    Dely = real(npad, 8) * (y_i(j) - y_i(i)) / psfyscale
    A0_aij(i, j) = imcom_interp_lookup(Alookup(:, :, iexp, jexp), Delx, Dely, &
                                       npoly, npad)

  end do
  if (i.eq.int(0.2 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 20% complete"
  if (i.eq.int(0.4 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 40% complete"
  if (i.eq.int(0.6 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 60% complete"
  if (i.eq.int(0.8 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 80% complete"

end do
!$omp end do
!$omp do schedule(dynamic, 16)
do i=2, n

  forall(j=1:i-1) A0_aij(i, j) = A0_aij(j, i)

end do
!$omp end do
!$omp end parallel
write(*, FMT='(A)') "IMCOM: A matrix complete"
if (saveA.eq.1) then 
  write(*, FMT='(A)') "IMCOM: Saving A matrix to "//trim(Afile)
  call imcom_writefits(trim(Afile), n, n, real(A0_aij, 8))
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

    Delx = (x_i(i) - X_a(a)) / psfxscale / real(npad, 8)
    Dely = (y_i(i) - Y_a(a)) / psfyscale / real(npad, 8)
    B_ia(i, a) = imcom_interp_lookup(Blookup(:, :, iexp), Delx, Dely, npoly, &
                                     npad)

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
integer :: alstat, iexp, jexp
integer :: n1big, n2big, n1min, n2min, n1max, n2max

n1big = n1psf * npad
n2big = n2psf * npad
n1min = (n1big - n1psf) / 2 + 1
n2min = (n2big - n2psf) / 2 + 1
n1max = (n1big + n1psf) / 2
n2max = (n2big + n2psf) / 2
write(*, FMT='(A)') "IMCOM: Fourier transforming A matrix lookup table"
allocate(Alookup(n1big, n2big, nexp, nexp), A_tmp(n1big, n2big), &
         ufunc(n1big, n2big), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot allocate memory for A lookup tables"
  stop
endif

!$omp parallel
!$omp do private(iexp, jexp, ufunc, A_tmp)
do iexp=1, nexp

  do jexp=1, nexp

    ufunc = dcmplx(0.d0, 0.d0)
!  Fix all the padding with circular shifts... Probably not optimal in 
!  N operations, but it saves me making stupid mistakes with counting...
    ufunc(n1min: n1max, n2min: n2max) = cshift(cshift(dconjg(Gt_rot(:, :,    &
                                                             iexp))          & 
                                                     * Gt_rot(:, :, jexp)    &
                                                     * imcom_upsilon2D(ux,   &
                                                      uy, 1), n2psf / 2, 2), &
                                                      n1psf / 2, 1)
    ufunc = cshift(cshift(ufunc, n2big / 2, 2), n1big / 2, 1)
    call imcom_invft_c2c(n1big, n2big, ufunc, A_tmp, 0)
    Alookup(:, :, iexp, jexp) = cshift(cshift(real(A_tmp, 8), n2big / 2, 2), &
                                                              n1big / 2, 1)

! Origin thus located at centre of pixel (n1big / 2 + 1, n2big / 2 + 1)

  end do

end do
!$omp end do
!$omp end parallel
deallocate(A_tmp, ufunc, STAT=alstat)
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
integer :: n1min, n2min, n1max, n2max

n1min = (n1psf * npad - n1psf) / 2 + 1
n2min = (n2psf * npad - n2psf) / 2 + 1
n1max = (n1psf * npad + n1psf) / 2
n2max = (n2psf * npad + n2psf) / 2
write(*, FMT='(A)') "IMCOM: Fourier transforming B matrix lookup tables"
allocate(Blookup(n1psf * npad, n2psf * npad, nexp), &
         B_tmp(n1psf * npad, n2psf * npad),         &
         ufunc(n1psf * npad, n2psf * npad), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot allocate memory for B lookup tables"
  stop
endif
ufunc = dcmplx(0.d0, 0.d0)
!$omp parallel
!$omp do private(iexp, ufunc, B_tmp)
do iexp=1, nexp

  ufunc(n1min: n1max, n2min: n2max) = dconjg(Gammat(:, :)) &
                                    * Gt_rot(:, :, iexp)   &
                                    * imcom_upsilon2D(ux, uy, 1)
  call imcom_invft_c2c(n1psf * npad, n2psf * npad, ufunc, B_tmp, 0)
  Blookup(:, :, iexp) = -2.d0 * cshift(cshift(real(B_tmp, 8),       &
                                              n2psf * npad / 2, 2), &
                                       n1psf * npad / 2, 1)
! Origin thus located at centre of pixel (n1psf * npad / 2 + 1, n2psf * npad / 2 + 1)

end do
!$omp end do
!$omp end parallel
deallocate(B_tmp, ufunc, STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot deallocate temporary memory for B lookup tables"
  stop
endif
!$omp parallel workshare
Blookup = Blookup / real(n1psf, 8) / real(n2psf, 8)
!$omp end parallel workshare
END SUBROUTINE imcom_build_Blookup

!---

FUNCTION imcom_interp_lookup(lookup, xarg, yarg, npoly_in, npad_in)
implicit none
real(KIND=8), dimension(:, :), intent(IN) :: lookup
real(KIND=8), intent(IN) :: xarg, yarg  ! x, y in lookup table pixel units (not arcsec!)
integer, intent(IN) :: npoly_in, npad_in
real(KIND=8) :: imcom_interp_lookup
real(KIND=8) :: icent, jcent, iarg, jarg, f, df
integer :: npoly, npm, n1big, n2big
integer :: ia_floor, ia_ceiling, ja_floor, ja_ceiling
integer :: imin, imax, jmin, jmax, ii, jj

npoly = npoly_in
if (mod(npoly, 2).eq.0) then
  write(*, FMT='(A)') "IMCOM WARNING: Even-order polynomial specified for IMCOM_INTERP_LOOKUP..."
  write(*, FMT='(A)') "IMCOM WARNING: Using NPOLY-1 for interpolation (e.g. linear, cubic, quintic...)"
  npoly = npoly - 1
end if
npm = npoly / 2
n1big = n1psf * npad_in
n2big = n2psf * npad_in
icent = real(n1big / 2 + 1, 8)
jcent = real(n2big / 2 + 1, 8)
iarg = icent + xarg
jarg = jcent + yarg
ia_floor = floor(iarg)
ja_floor = floor(jarg)
ia_ceiling = ceiling(iarg)
ja_ceiling = ceiling(jarg)
if ((ia_floor.ge.1).and.(ja_floor.ge.1).and.(ia_ceiling.le.n1big) &
                                       .and.(ja_ceiling.le.n2big)) then
 imin = max(1, ia_floor - npm)
 jmin = max(1, ja_floor - npm)
 imax = min(n1big, ia_ceiling + npm)
 jmax = min(n2big, ja_ceiling + npm)
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
write(*, FMT='(A, E21.15)') "IMCOM: C_a = ", C_a
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
  call imcom_calc_ATTpBT(n, (1 + a_f - a_i), A0_aij, B_ia(:, a_i:a_f), &
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

END MODULE imcom_matrix
