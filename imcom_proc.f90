!    IMCOM_PROC.F90 - Fortran 95 module that contains basic procedures for
!                     the prototype IMCOM package
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
!                   FFTW3 LIBRARIES
!
MODULE imcom_proc

USE imcom_data

implicit none

CONTAINS

FUNCTION imcom_xarr(n1, n2)
! Outputs an array of x values at the centre of pixels, starting at (0., 0.)
! etc. (all of the values in any particular column will be identical)
implicit none
integer, intent(IN) :: n1, n2              ! x and y dimensions of image
real(KIND=8), dimension(n1, n2) :: imcom_xarr ! output
integer :: i
real(KIND=8), dimension(n1) :: x 

x = (/ (real(i, 8) - 1.d0, i=1, n1) /)
imcom_xarr = spread(x, 2, n2)
END FUNCTION imcom_xarr

!---

FUNCTION imcom_yarr(n1, n2)
! Outputs an array of y values at the centre of pixels, starting at (0., 0.)
! etc. (all of the values in any particular row will be identical)
implicit none
integer, intent(IN) :: n1, n2      ! x and y dimensions of image
real(KIND=8), dimension(n1, n2) :: imcom_yarr  ! output

imcom_yarr = transpose(imcom_xarr(n2, n1))
END FUNCTION imcom_yarr

!---

SUBROUTINE imcom_build_xy
! Make the x, y, x_i & x_i (including dimension check) arrays for the
! galaxy images
implicit none
integer :: alstat, dealstat, i

allocate(x_unmasked(n1gim, n2gim, nexp), y_unmasked(n1gim, n2gim, nexp), &
         STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot allocate memory for x, y"
  stop
endif
do i=1, nexp

  x_unmasked(:, :, i) = imcom_xarr(n1gim, n2gim) * gimxscale(i) + dither(i, 1)
  y_unmasked(:, :, i) = imcom_yarr(n1gim, n2gim) * gimyscale(i) + dither(i, 2)

end do
allocate(x_unmasked_i(n_unmasked), y_unmasked_i(n_unmasked), x_i(n), y_i(n), &
         STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot allocate memory for x_i, y_i"
  stop
endif
x_unmasked_i = reshape(x_unmasked, (/ n_unmasked /))
y_unmasked_i = reshape(y_unmasked, (/ n_unmasked /))
forall(i=1: n) x_i(i) = x_unmasked_i(mask_subs_i(i))
forall(i=1: n) y_i(i) = y_unmasked_i(mask_subs_i(i))
deallocate(x_unmasked, y_unmasked, x_unmasked_i, y_unmasked_i, STAT=dealstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot deallocate memory for x, y"
  stop
endif
END SUBROUTINE imcom_build_xy

!---

SUBROUTINE imcom_build_XaYa
! Make the x, y arrays for the galaxy images
implicit none
integer :: alstat, dealstat

allocate(Xa(n1out, n2out), Ya(n1out, n2out), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot allocate memory for Xa, Ya"
  stop
endif
Xa(:, :) = imcom_xarr(n1out, n2out) * outxscale + outpos1
Ya(:, :) = imcom_yarr(n1out, n2out) * outyscale + outpos2
allocate(X_a(m), Y_a(m), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot allocate memory for X_a & Y_a"
  stop
endif
X_a = 0.d0
Y_a = 0.d0
X_a = reshape(Xa, (/ m /))
Y_a = reshape(Ya, (/ m /))
deallocate(Xa, Ya, STAT=dealstat)
if (dealstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot deallocate memory for Xa & Ya"
  stop
endif
END SUBROUTINE imcom_build_XaYa

!---

SUBROUTINE imcom_build_uxuy
! Make the ux, uy arrays for the PSF images...
implicit none
integer :: alstat, j

allocate(ux(n1psf, n2psf), uy(n1psf, n2psf), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot allocate memory for ux, uy"
  stop
endif
ux = 0.d0
uy = 0.d0
if (mod(n1psf, 2).eq.0) then
  ux(1 : 1 + n1psf / 2, :) = imcom_xarr(1 + n1psf / 2, n2psf) &
                           / psfxscale / real(n1psf, 8)
  forall(j=1 : n1psf / 2 - 1) ux(1 + j + n1psf / 2, :) &
                           = -ux(1 - j + n1psf / 2, :)
else
  ux(1 : 1 + n1psf / 2, :) = imcom_xarr(1 + n1psf / 2, n2psf) &
                           / psfxscale / real(n1psf, 8)
  forall(j=1 : n1psf / 2) ux(1 + j + n1psf / 2, :) &
                       = -ux(2 - j + n1psf / 2, :)
end if
if (mod(n2psf, 2).eq.0) then
  uy(:, 1 : 1 + n2psf / 2) = imcom_yarr(n1psf, 1 + n2psf / 2) &
                           / psfyscale / real(n2psf, 8)
  forall(j=1 : n2psf / 2 - 1) uy(:, 1 + j + n2psf / 2) &
                           = -uy(:, 1 - j + n2psf / 2)
else
  uy(:, 1 : 1 + n2psf / 2) = imcom_yarr(n1psf, 1 + n2psf / 2) &
                           / psfyscale / real(n2psf, 8)
  forall(j=1 : n2psf / 2) uy(:, 1 + j + n2psf / 2) &
                       = -uy(:, 2 - j + n2psf / 2)
end if
END SUBROUTINE imcom_build_uxuy

!---

SUBROUTINE imcom_rotate_xy(x_dum, y_dum, theta_deg, xp_dum, yp_dum)
! Linearly rotate arrays of x and y to give a new x' and y'
!
implicit none
real(KIND=8), dimension(:, :), intent(IN) :: x_dum, y_dum
real(KIND=8), intent(IN) :: theta_deg
real(KIND=8), dimension(size(x_dum, 1), size(x_dum, 2)), intent(OUT) :: xp_dum, yp_dum
real(KIND=8) :: theta_rad

if (size(x_dum, 1).ne.size(y_dum, 1)) then
  write(*, FMT='(A)') "IMCOM ERROR: X and Y input to IMCOM_ROTATE_XY of unequal dimensions"
  stop
endif
if (size(x_dum, 2).ne.size(y_dum, 2)) then
  write(*, FMT='(A)') "IMCOM ERROR: X and Y input to IMCOM_ROTATE_XY of unequal dimensions"
  stop
endif
theta_rad = mod(theta_deg, 360.d0) * pi / 180.d0
xp_dum =  x_dum * cos(theta_rad) - y_dum * sin(theta_rad)
yp_dum =  x_dum * sin(theta_rad) + y_dum * cos(theta_rad)
END SUBROUTINE imcom_rotate_xy

!---

SUBROUTINE imcom_rotate_image(im_in, theta_deg, npoly_in, conserve_flux, im_out)
! Linearly rotate an image about its centre using polynomial interpolation
!
implicit none
real(KIND=8), dimension(:, :), intent(IN) :: im_in
real(KIND=8), intent(IN) :: theta_deg
integer, intent(IN) :: npoly_in, conserve_flux
real(KIND=8), dimension(size(im_in, 1), size(im_in, 2)), intent(OUT) :: im_out
real(KIND=8) :: theta_rad, cost, sint
integer :: nim1, nim2
real(KIND=8) :: xs, ys, x_cent, y_cent, f, df
integer :: is_floor, is_ceiling, js_floor, js_ceiling
real(KIND=8), dimension(size(im_in, 1), size(im_in, 2)) :: xin, yin
real(KIND=8) :: xout, yout
integer :: npoly, npm, iout, jout, imin, jmin, imax, jmax, ii, jj
!real(KIND=8) :: int1, int2, delx, dely ! delx, dely < 1

npoly = npoly_in
if (mod(npoly, 2).eq.0) then
  write(*, FMT='(A)') "IMCOM WARNING: Odd-order polynomial interpolation specified..."
  write(*, FMT='(A)') "IMCOM WARNING: Using NPOLY-1 for interpolation (e.g. quadratic...)"
  npoly = npoly - 1
end if
npm = npoly / 2    ! integer division rounds down, i.e. 1/2 = 0, 3/2 = 1
theta_rad = mod(theta_deg, 360.d0) * pi / 180.d0
cost = cos(theta_rad)
sint = sin(theta_rad)
nim1 = size(im_in, 1)
nim2 = size(im_in, 2)
x_cent = real(nim1 + 1, 8) / 2.d0
y_cent = real(nim2 + 1, 8) / 2.d0
xin = imcom_xarr(size(im_in, 1), size(im_in, 2)) + 1.d0
yin = imcom_yarr(size(im_in, 1), size(im_in, 2)) + 1.d0
im_out = 0.d0
!$omp parallel do private(xout, yout, xs, ys, is_floor, js_floor, is_ceiling, js_ceiling, imin, jmin, imax, jmax, ii, jj, f, df)
do iout=1, nim1

  xout = real(iout, 8) - x_cent
  do jout=1, nim2

    yout = real(jout, 8) - y_cent
    xs =  xout * cost + yout * sint + x_cent
    ys = -xout * sint + yout * cost + y_cent
    is_floor = floor(xs)
    js_floor = floor(ys)
    is_ceiling = ceiling(xs)
    js_ceiling = ceiling(ys)
    if ((is_floor.ge.1).and.(js_floor.ge.1)      &
                       .and.(is_ceiling.le.nim1) &
                       .and.(js_ceiling.le.nim2)) then
      imin = max(1, is_floor - npm)
      jmin = max(1, js_floor - npm)
      imax = min(nim1, is_ceiling + npm)
      jmax = min(nim2, js_ceiling + npm)
      call imcom_polin2((/ (real(ii, 8), ii=imin, imax) /), &
                        (/ (real(jj, 8), jj=jmin, jmax) /), &
                        im_in(imin: imax, jmin: jmax), xs, ys, f, df)
      im_out(iout, jout) = f
    else
      im_out(iout, jout) = 0.d0
    end if

  end do

end do
!$omp end parallel do
if (conserve_flux.ne.0) im_out = sum(im_in) * im_out / sum(im_out)
END SUBROUTINE imcom_rotate_image

!---

SUBROUTINE imcom_rotate_cimage(im_in, theta_deg, npoly, conserve_flux, im_out)
! Linearly rotate a complex image about its centre using blinear interpolation 
! (...code up bicubic?) 
!
implicit none
complex(KIND=8), dimension(:, :), intent(IN) :: im_in
real(KIND=8), intent(IN) :: theta_deg
integer, intent(IN) :: npoly, conserve_flux
complex(KIND=8), dimension(size(im_in, 1), size(im_in, 2)), intent(OUT) :: im_out
real(KIND=8), dimension(size(im_in, 1), size(im_in, 2)) :: rin, rout
real(KIND=8), dimension(size(im_in, 1), size(im_in, 2)) :: iin, iout
real(KIND=8) :: flux_in, flux_out

im_out = dcmplx(0.d0, 0.d0)
rin = real(im_in, 8)
iin = dimag(im_in)
call imcom_rotate_image(rin, theta_deg, npoly, 0, rout)
call imcom_rotate_image(iin, theta_deg, npoly, 0, iout)
im_out = dcmplx(rout, iout)
if (conserve_flux.ne.0) then
  flux_in = sum(abs(im_in))
  flux_out = sum(abs(im_out))
  im_out = dcmplx(flux_in, 0.d0) * im_out / dcmplx(flux_out, 0.d0)
endif
END SUBROUTINE imcom_rotate_cimage

!---

FUNCTION imcom_upsilon(ux_dum, uy_dum, option)
! Real symmetric positive definite kernel (p2, combo.pdf)
!
implicit none
real(KIND=8), intent(IN) :: ux_dum
real(KIND=8), intent(IN) :: uy_dum
integer, intent(IN) :: option
complex(KIND=8) :: imcom_upsilon  ! output var

if (option.eq.1) then
  imcom_upsilon = dcmplx(1.d0, 0.d0)
else if (option.eq.2) then
  imcom_upsilon = zexp(-dcmplx((ux_dum * ux_dum + uy_dum * uy_dum), 0.d0))
else
  imcom_upsilon = dcmplx(0.d0, 0.d0)
  write(*, FMT='(A)') "IMCOM ERROR: Upsilon function must take OPTION = 1 or 2"
  stop
end if
END FUNCTION imcom_upsilon

!---

FUNCTION imcom_upsilon1D(ux_dum, uy_dum, option)
! Real symmetric positive definite kernel (p2, combo.pdf)
implicit none
real(KIND=8), dimension(:), intent(IN) :: ux_dum
real(KIND=8), dimension(:), intent(IN) :: uy_dum
integer, intent(IN) :: option
complex(KIND=8), dimension(size(ux_dum)) :: imcom_upsilon1D 
! output var

if (size(ux_dum).ne.size(uy_dum)) then
  write(*, FMT='(A)') "IMCOM ERROR: UX and UY arguments to UPSILON must be equal shape"
  stop
endif
if (option.eq.1) then
  imcom_upsilon1D = dcmplx(1.d0, 0.d0)
else if (option.eq.2) then
  imcom_upsilon1D = zexp(-dcmplx((ux_dum * ux_dum + uy_dum * uy_dum), 0.d0))
else
  imcom_upsilon1D = dcmplx(0.d0, 0.d0)
  write(*, FMT='(A)') "IMCOM ERROR: Upsilon2D function must take OPTION = 1 or 2"
  stop
end if
END FUNCTION imcom_upsilon1D

!---

FUNCTION imcom_upsilon2D(ux_dum, uy_dum, option)
! Real symmetric positive definite kernel (p2, combo.pdf)
implicit none
real(KIND=8), dimension(:, :), intent(IN) :: ux_dum
real(KIND=8), dimension(:, :), intent(IN) :: uy_dum
integer, intent(IN) :: option
complex(KIND=8), dimension(size(ux_dum, 1), size(uy_dum, 2)) :: imcom_upsilon2D 
! output var

if (size(ux_dum, 1).ne.size(uy_dum, 1).or.size(ux_dum, 2).ne.size(uy_dum, 2)) then
!  write(*, FMT='(A)') "IMCOM ERROR: UX and UY arguments to UPSILON must be equal shape"
!  stop
endif
if (option.eq.1) then
  imcom_upsilon2D = dcmplx(1.d0, 0.d0)
else if (option.eq.2) then
  imcom_upsilon2D = zexp(-dcmplx((ux_dum * ux_dum + uy_dum * uy_dum), 0.d0))
else
  imcom_upsilon2D = dcmplx(0.d0, 0.d0)
  write(*, FMT='(A)') "IMCOM ERROR: Upsilon2D function must take OPTION = 1 or 2"
  stop
end if
END FUNCTION imcom_upsilon2D

!---

PURE FUNCTION imcom_dsinc(x)
implicit none
real(KIND=8), intent(IN) :: x
real(KIND=8) :: imcom_dsinc     ! Output variable
real(KIND=8) :: xx

xx = dabs(x)
if (xx.lt.1.d-15) then
  imcom_dsinc = 1.d0
else
  imcom_dsinc = dsin(pi * x) / x / pi
end if
END FUNCTION imcom_dsinc

!---

PURE FUNCTION imcom_dsinc2(x)
implicit none
real(KIND=8), intent(IN) :: x
real(KIND=8) :: imcom_dsinc2    ! Output variable
real(KIND=8) :: xx, sintmp

xx = dabs(x)
if (xx.lt.1.d-15) then
  imcom_dsinc2 = 1.d0
else
  sintmp = dsin(pi * x)
  imcom_dsinc2 = sintmp * sintmp / x / x / pi / pi
end if
END FUNCTION imcom_dsinc2

!---

FUNCTION imcom_iminloc(arr)
implicit none
REAL(KIND=8), DIMENSION(:), INTENT(IN) :: arr
INTEGER, DIMENSION(1) :: imin
INTEGER :: imcom_iminloc
imin = minloc(arr(:))
imcom_iminloc = imin(1)
END FUNCTION imcom_iminloc

!---

SUBROUTINE imcom_polint(xa, ya, x, y, dy, eps)
! Adapted from Numerical recipes polint.f90
implicit none
REAL(KIND=8), DIMENSION(:), INTENT(IN) :: xa, ya
REAL(KIND=8), INTENT(IN) :: x, eps
REAL(KIND=8), INTENT(OUT) :: y, dy
INTEGER :: m, n, ns
REAL(KIND=8), DIMENSION(size(xa)) :: c, d, den, ho

n = size(xa)
if (size(ya).ne.n) then
  write(*, FMT='(A)') "IMCOM ERROR: x and y input must be same size for IMCOM_POLINT"
  stop
end if
c = ya
d = ya
ho = xa - x
ns = imcom_iminloc(abs(x - xa))
y = ya(ns)
ns = ns-1
do m=1, n-1

  den(1: n - m) = ho(1: n - m) - ho(1 + m: n)
  if (any(abs(den(1:n-m)).le.eps)) then
    write(*, FMT='(A)') "IMCOM ERROR: Two or more identical x input locations to IMCOM_POLINT"
    stop
  end if
  den(1: n - m) = (c(2: n - m + 1) - d(1: n - m)) / den(1: n - m)
  d(1: n - m) = ho(1 + m: n) * den(1: n - m)
  c(1: n - m) = ho(1: n - m) * den(1: n - m)
  if (2 * ns < n - m) then
    dy = c(ns + 1)
  else
    dy = d(ns)
    ns = ns - 1
  end if
  y = y + dy

end do
END SUBROUTINE imcom_polint

!---

SUBROUTINE imcom_polin2(x1a, x2a, ya, x1, x2, y, dy)
! Adapted from Numerical Recipes polin2.f90
implicit none
REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x1a, x2a
REAL(KIND=8), DIMENSION(:, :), INTENT(IN) :: ya
REAL(KIND=8), INTENT(IN) :: x1, x2
REAL(KIND=8), INTENT(OUT) :: y, dy
INTEGER :: j, m, ndum
REAL(KIND=8), DIMENSION(size(x1a)) :: ymtmp
REAL(KIND=8), DIMENSION(size(x2a)) :: yntmp
REAL(KIND=8) :: eps
REAL(KIND=8), EXTERNAL :: DLAMCH

m = size(x1a) ! assert_eq(size(x1a),size(ya,1),'polin2: m')
if (size(ya, 1).ne.m) then
  write(*, FMT='(A)') "IMCOM ERROR: x1a and ya(:, 1) input must be same size for IMCOM_POLIN2"
  stop
end if
ndum = size(x2a) ! assert_eq(size(x2a),size(ya,2),'polin2: ndum')
if (size(ya, 2).ne.ndum) then
  write(*, FMT='(A)') "IMCOM ERROR: x2a and ya(1, :) input must be same size for IMCOM_POLIN2"
  stop
end if
eps = DLAMCH('Eps')   ! From LAPACK, machine epsilon, used in polint
do j=1, m

  yntmp = ya(j, :)
  call imcom_polint(x2a, yntmp, x2, ymtmp(j), dy, eps)

end do
call imcom_polint(x1a, ymtmp, x1, y, dy, eps)
END SUBROUTINE imcom_polin2

!---

SUBROUTINE imcom_ft_r2c(n1, n2, plan, F_dum, Ft_dum)
! Perform the real-to-complex Discrete Fourier Transform using the linked
! FFTW3 library routines... NOTE this currently only ESTIMATES the best plan 
! for performing the DFT and doesn't MEASURE it (see FFTW docs)
implicit none
integer, intent(IN) :: n1, n2
integer(KIND=8), intent(IN) :: plan
real(KIND=8), dimension(0 : n1 - 1, 0 : n2 - 1), intent(IN) :: F_dum
complex(KIND=8), dimension(0 : n1 - 1, 0 : n2 - 1), intent(OUT) :: Ft_dum
real(KIND=8), dimension(0: n1 - 1, 0: n2 - 1) :: F_tmp
complex(KIND=8), dimension(0: n1 / 2, 0 : n2 - 1) :: Ft_tmp
complex(KIND=8), dimension(0: n1 / 2, 0 : n2 - 1) :: Ft_conj
integer :: i

F_tmp = real(F_dum, 8)     ! Fill working array with data before FFT
call DFFTW_EXECUTE_DFT_R2C(plan, F_tmp, Ft_tmp)
! From FFTW3 documentation:
!In many practical applications, the input data in[i] are purely real numbers, in which case the DFT output satisfies the Hermitian redundancy: out[i] is the conjugate of out[n-i]. It is possible to take advantage of these circumstances in order to achieve roughly a factor of two improvement in both speed and memory usage.
!
! Fix this...
!$omp parallel
!$omp workshare
Ft_conj = dconjg(Ft_tmp)
Ft_dum(0 :n1 / 2, :) = Ft_tmp(0 :n1 / 2, :)
!$omp end workshare
Ft_dum(n1 / 2 + 1: n1 - 1, 0) = Ft_conj((n1 - 1) / 2: 1: -1, 0)
if (mod(n1, 2).eq.0) then
!$omp workshare
  forall(i = 1:(n1 - 1) / 2) Ft_dum(n1 / 2 + i, 1: n2 - 1)    &
                           = Ft_conj(n1 / 2 - i, n2 - 1: 1: -1)
!$omp end workshare
else if (mod(n1, 2).eq.1) then
!$omp workshare
  forall(i = 1:(n1 - 1) / 2) Ft_dum(n1 / 2 + i, 1: n2 - 1)    &
                           = Ft_conj(n1 / 2 + 1 - i, n2 - 1: 1: -1)
!$omp end workshare
end if
!$omp end parallel
END SUBROUTINE imcom_ft_r2c

!---

SUBROUTINE imcom_ft_c2c(n1, n2, plan, F_dum, Ft_dum)
! Perform the complex-to-complex Discrete Fourier Transform using the linked
! FFTW3 library routines... NOTE this currently only ESTIMATES the best plan 
! for performing the DFT and doesn't MEASURE it (see FFTW docs)
implicit none
integer, intent(IN) :: n1, n2
integer(KIND=8), intent(IN) :: plan
complex(KIND=8), dimension(n1, n1), intent(IN) :: F_dum
complex(KIND=8), dimension(n1, n2), intent(OUT) :: Ft_dum
complex(KIND=8), dimension(n1, n2) :: F_tmp
complex(KIND=8), dimension(n1, n2) :: Ft_tmp

F_tmp = F_dum     ! Fill working array with data before FFT
call DFFTW_EXECUTE_DFT(plan, F_tmp, Ft_tmp)
Ft_dum = Ft_tmp   ! Fill output array with FFT 
END SUBROUTINE imcom_ft_c2c

!---

SUBROUTINE imcom_invft_c2c(n1, n2, plan, Ft_dum, F_dum)
! Perform the complex-to-complex inverse Discrete Fourier Transform using 
! the linked FFTW3 library routines... NOTE this currently only ESTIMATES the 
! best plan for performing the DFT and doesn't MEASURE it (see FFTW docs)
implicit none
integer, intent(IN) :: n1, n2
integer(KIND=8), intent(IN) :: plan
complex(KIND=8), dimension(n1, n1), intent(IN) :: Ft_dum
complex(KIND=8), dimension(n1, n2), intent(OUT) :: F_dum
complex(KIND=8), dimension(n1, n2) :: F_tmp
complex(KIND=8), dimension(n1, n2) :: Ft_tmp

Ft_tmp = Ft_dum     ! Fill working array with data before FFT
call DFFTW_EXECUTE_DFT(plan, Ft_tmp, F_tmp)
F_dum = F_tmp       ! Fill output array with inverse FFT 
END SUBROUTINE imcom_invft_c2c

!---

SUBROUTINE imcom_plan_ft_r2c(n1, n2, info, plan)
implicit none
integer, intent(IN) :: n1, n2, info
integer(KIND=8), intent(OUT) :: plan
real(KIND=8), dimension(0: n1 - 1, 0: n2 - 1) :: F_tmp
complex(KIND=8), dimension(0: n1 / 2, 0 : n2 - 1) :: Ft_tmp

call DFFTW_PLAN_DFT_R2C_2D(plan, n1, n2, F_tmp, Ft_tmp, &
                           FFTW_ESTIMATE + FFTW_PRESERVE_INPUT)
if (info.ne.0) then 
  write(*, FMT='(A)') "IMCOM: Outputing FFTW3 plan output"
  call DFFTW_PRINT_PLAN(plan)
endif
END SUBROUTINE imcom_plan_ft_r2c

!---

SUBROUTINE imcom_plan_ft_c2c(n1, n2, info, plan)
implicit none
integer, intent(IN) :: n1, n2, info
integer(KIND=8), intent(OUT) :: plan
complex(KIND=8), dimension(n1, n2) :: F_tmp
complex(KIND=8), dimension(n1, n2) :: Ft_tmp

call DFFTW_PLAN_DFT_2D(plan, n1, n2, F_tmp, Ft_tmp, FFTW_FORWARD, &
                       FFTW_ESTIMATE + FFTW_PRESERVE_INPUT)
if (info.ne.0) then 
  write(*, FMT='(A)') "IMCOM: Outputing FFTW3 plan output"
  call DFFTW_PRINT_PLAN(plan)
endif
END SUBROUTINE imcom_plan_ft_c2c

!---

SUBROUTINE imcom_plan_invft_c2c(n1, n2, info, plan)
implicit none
integer, intent(IN) :: n1, n2, info
integer(KIND=8), intent(OUT) :: plan
complex(KIND=8), dimension(n1, n2) :: F_tmp
complex(KIND=8), dimension(n1, n2) :: Ft_tmp

call DFFTW_PLAN_DFT_2D(plan, n1, n2, Ft_tmp, F_tmp, FFTW_BACKWARD, &
                       FFTW_ESTIMATE + FFTW_PRESERVE_INPUT)
if (info.ne.0) then 
  write(*, FMT='(A)') "IMCOM: Outputing FFTW3 plan output"
  call DFFTW_PRINT_PLAN(plan)
endif
END SUBROUTINE imcom_plan_invft_c2c

!---

SUBROUTINE imcom_destroy_plan(plan)
implicit none
integer(KIND=8), intent(IN) :: plan

call DFFTW_DESTROY_PLAN(plan)
END SUBROUTINE imcom_destroy_plan

!---

FUNCTION imcom_test_psfconst()
implicit none
integer :: imcom_test_psfconst
integer :: i

imcom_test_psfconst = 1
do i=2, nexp

  if ((psffile(i).ne.psffile(1)).or.(rotangdeg(i).ne.rotangdeg(1))) then
    imcom_test_psfconst = 0
    exit
  end if

end do
END FUNCTION imcom_test_psfconst

!---

END MODULE imcom_proc
