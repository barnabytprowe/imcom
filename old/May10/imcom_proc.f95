!
! Contains basic procedures for the imcom package
!
!--------------------------------------------------------------------

MODULE imcom_proc

USE imcom_data

implicit none

CONTAINS

FUNCTION imcom_xarr(n1, n2)
! Outputs an array of x values at the centre of pixels, starting at (0., 0.)
! etc. (all of the values in any particular column will be identical)
implicit none
real(KIND=8), dimension(n1, n2) :: imcom_xarr ! output
integer, intent(IN) :: n1, n2              ! x and y dimensions of image
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
real(KIND=8), dimension(n1, n2) :: imcom_yarr  ! output
integer, intent(IN) :: n1, n2      ! x and y dimensions of image

imcom_yarr = transpose(imcom_xarr(n2, n1))
END FUNCTION imcom_yarr

!---

SUBROUTINE imcom_build_xy
! Make the x, y, x_i & x_i (including dimension check) arrays for the
! galaxy images
implicit none
integer :: alstat, dealstat, i

allocate(x(n1gim, n2gim, nexp), y(n1gim, n2gim, nexp), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot allocate memory for x, y"
  stop
endif
do i=1, nexp

  x(:, :, i) = imcom_xarr(n1gim, n2gim) * gimxscale(i) + dither(i, 1)
  y(:, :, i) = imcom_yarr(n1gim, n2gim) * gimyscale(i) + dither(i, 2)

end do
if (size(x).ne.n) then
  write(*, FMT='(A)') "IMCOM ERROR: Dimensions of x do not match x_i!?"
  stop
endif
if (size(y).ne.n) then
  write(*, FMT='(A)') "IMCOM ERROR: Dimensions of y do not match y_i!?"
  stop
endif
allocate(x_i(n), y_i(n), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot allocate memory for x_i, y_i"
  stop
endif
x_i = 0.d0
y_i = 0.d0
x_i = reshape(x, (/ n /))
y_i = reshape(y, (/ n /))
deallocate(x, y, STAT=dealstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot deallocate memory for x, y"
  stop
endif
END SUBROUTINE imcom_build_xy

!---

SUBROUTINE imcom_build_XaYa
! Make the x, y arrays for the galaxy images
implicit none
integer :: alstat, dealstat, i

allocate(Xa(n1out, n2out), Ya(n1out, n2out), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot allocate memory for Xa, Ya"
  stop
endif
Xa(:, :) = imcom_xarr(n1out, n2out) * outxscale + outpos1
Ya(:, :) = imcom_yarr(n1out, n2out) * outyscale + outpos2
allocate(X_a(m), Y_a(m), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot allocate memory for X_a, Y_a"
  stop
endif
X_a = 0.d0
Y_a = 0.d0
X_a = reshape(Xa, (/ m /))
Y_a = reshape(Ya, (/ m /))
deallocate(Xa, Ya, STAT=dealstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot deallocate memory for Xa, Ya"
  stop
endif
END SUBROUTINE imcom_build_XaYa

!---

SUBROUTINE imcom_build_uxuy
! Make the ux, uy arrays for the PSF images...
implicit none
integer :: alstat, i, j

allocate(ux(n1psf, n2psf), uy(n1psf, n2psf), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot allocate memory for ux, uy"
  stop
endif
ux = 0.d0
uy = 0.d0
if (modulo(n1psf, 2).eq.0) then
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
if (modulo(n2psf, 2).eq.0) then
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

PURE FUNCTION imcom_upsilon(ux_dum, uy_dum, option)
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
  imcom_upsilon = zexp( - (ux_dum * ux_dum + uy_dum * uy_dum) &
                * dcmplx(1.d0, 0.d0))
end if
END FUNCTION imcom_upsilon

!---

PURE FUNCTION imcom_upsilon2D(ux_dum, uy_dum, option)
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
  imcom_upsilon2D = zexp( - (ux_dum * ux_dum + uy_dum * uy_dum) &
                  * dcmplx(1.d0, 0.d0))
end if
END FUNCTION imcom_upsilon2D

!---

PURE FUNCTION imcom_dsinc(x)
implicit none
real(KIND=8), intent(IN) :: x
real(KIND=8) :: imcom_dsinc     ! Output variable
real(KIND=8) :: xx

xx = dabs(x)
if (xx.lt.1.e-15) then
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
if (xx.lt.1.e-15) then
  imcom_dsinc2 = 1.d0
else
  sintmp = dsin(pi * x)
  imcom_dsinc2 = sintmp * sintmp / x / x / pi / pi
end if
END FUNCTION imcom_dsinc2

!---

SUBROUTINE imcom_ft(n1, n2, F_dum, Ft_dum, info)
! Perform the Discrete Fourier Transform using the linked FFTW3 library routines
implicit none
integer, intent(IN) :: n1, n2, info
real(KIND=8), dimension(n1, n1), intent(IN) :: F_dum
complex(KIND=8), dimension(n1, n2), intent(OUT) :: Ft_dum
real(KIND=8), dimension(n1, n2) :: F_tmp
complex(KIND=8), dimension(n1 / 2 + 1, n2) :: Ft_tmp
complex(KIND=8), dimension(n1 / 2 + 1, n2) :: Ft_conj
integer(KIND=8) :: plan
integer :: i

call DFFTW_PLAN_DFT_R2C_2D(plan, n1, n2, F_tmp, Ft_tmp, &
                           FFTW_MEASURE + FFTW_PRESERVE_INPUT)
if (info.ne.0) then 
  call DFFTW_PRINT_PLAN(plan)
  print *, F_dum(1:info, 1)
endif
F_tmp = real(F_dum, 8)     ! Fill working array with data before FFT
call DFFTW_EXECUTE_DFT_R2C(plan, F_tmp, Ft_tmp)
call DFFTW_DESTROY_PLAN(plan)
! From FFTW3 documentation:
!In many practical applications, the input data in[i] are purely real numbers, in which case the DFT output satisfies the Hermitian redundancy: out[i] is the conjugate of out[n-i]. It is possible to take advantage of these circumstances in order to achieve roughly a factor of two improvement in both speed and memory usage.
!
! Fix this...
Ft_conj = conjg(Ft_tmp)
Ft_dum(1 : n1 / 2 + 1, :) = Ft_tmp
forall(i =  1 : n1 / 2) Ft_dum(n1 / 2 + 1 + i, :) = Ft_conj(n1 / 2 + 2 - i, :)
END SUBROUTINE imcom_ft

!---

END MODULE imcom_proc
