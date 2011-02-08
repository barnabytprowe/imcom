PROGRAM imcom

USE imcom_data
USE imcom_io
USE imcom_proc
USE imcom_matrix
USE imcom_bisect
USE omp_lib

implicit none
integer :: nargs
character(LEN=256) :: config
character(LEN=1) :: US_string
character(LEN=256) :: buffer
integer :: userxy, i, npoly, npad
integer(KIND=8) :: ftplan
real(KIND=8) :: eps
real(KIND=8), external :: DLAMCH

npoly = 5    ! Make these user specifiable!
npad = 3     !

nargs = iargc()
if (nargs.eq.7) then
  call getarg(1, config)
  call getarg(2, US_string)
  if((US_string.eq."U").or.(US_string.eq."u")) then
    USB = .TRUE.    ! Solve for T_ia with U_a < US_min
  else if((US_string.eq."S").or.(US_string.eq."s")) then
    USB = .FALSE.   ! Solve for T_ia with S_a < US_min
  else
    call imcom_usage_message
    stop
  endif
  call getarg(3, buffer)
  read(buffer, FMT=*) US_max
  call getarg(4, buffer)
  read(buffer, FMT=*) US_tol
  call getarg(5, buffer)
  read(buffer, FMT=*) kappa_min
  call getarg(6, buffer)
  read(buffer, FMT=*) kappa_max
  call getarg(7, buffer)
  read(buffer, FMT=*) maxNbis
else
  call imcom_usage_message
  stop
endif

if (kappa_min.ge.kappa_max) write(*, FMT='(A)') "IMCOM WARNING: kappa_min >= kappa_max!"



! Welcome banner
write(*, FMT='(A)') "IMCOM: IMage COMbination v0.1 (B. Rowe & C. Hirata 2010)"
! Read the config files
call imcom_read_config(config, userxy)

! Read and allocate the input images
call imcom_alloc_gims(n1gim, n2gim)

if (userxy.eq.1) then
! Read the x-y arrays
  call imcom_alloc_xy(n1gim, n2gim)
  call imcom_alloc_XaYa(n1out, n2out)
else
! Make the x-y arrays
  call imcom_build_xy
  call imcom_build_XaYa
end if
! Read and allocate the PSFs and Gamma PSF
call imcom_alloc_psfs(n1psf, n2psf)
! Make the ux, uy FT arrays
call imcom_build_uxuy
! Rotate PSFs and calculate the Fourier transform Gt(ux, uy) of each
write(*, FMT='(A)') "IMCOM: Rotating and Fourier transforming PSF images"
eps = DLAMCH('Epsilon')
call imcom_plan_ft_r2c(n1psf, n2psf, 0, ftplan)
do i=1, nexp

  if (abs(rotangdeg(i)).gt.(2.d0 * pi * eps)) then
    call imcom_rotate_image(G_unrot(:, :, i), rotangdeg(i), npoly, 1, &
                            G_rot(:, :, i))
  else
    G_rot(:, :, i) = G_unrot(:, :, i)
  end if
  call imcom_ft_r2c(n1psf, n2psf, ftplan, G_rot(:, :, i), Gt_rot(:, :, i))

end do
! Calculate the Fourier transform Gammat(ux, uy) of Gamma(y, y)
call imcom_ft_r2c(n1psf, n2psf, ftplan, Gamma(:, :), Gammat(:, :))
call imcom_destroy_plan(ftplan)

! Output some checks
!call imcom_writefits('./fits/ux_check.fits', n1psf, n2psf, ux)
!call imcom_writefits('./fits/uy_check.fits', n1psf, n2psf, uy)
!call imcom_writefits('./fits/mtf1.real.fits', n1psf, n2psf, real(Gt_rot(:, :, 1), 8))
!call imcom_writefits('./fits/mtf1.imag.fits', n1psf, n2psf, dimag(Gt_rot(:, :, 1)))
!call imcom_writefits('./fits/mtf1.abs.fits', n1psf, n2psf, abs(Gt_rot(:, :, 1)))

! Calculate the N (diagonal only, currently) noise array
call imcom_build_Ndiag

! Calulate C (really quick) and build the (currently diagonal only) noise cov.
call imcom_calc_C

! Read / build the A array
call imcom_get_A(npoly, npad, 1, 0)     ! arg#3 = save, arg#4 = force build
Qfile = Qdeffile
Lfile = Ldeffile
call imcom_eigen_A

! Read / build the B array
call imcom_get_B(npoly, npad, 1, 0)

! Construct the projection matrices
Pfile = Pdeffile
call imcom_build_P

! Read / build the T array
call imcom_get_T(0)

! Build the output products using the galaxy images, A, B, C and our calculated T transformation
call imcom_build_image

END PROGRAM imcom
