PROGRAM imcom

USE imcom_data
USE imcom_io
USE imcom_proc
USE imcom_matrix
USE omp_lib

implicit none
integer :: nargs
character(LEN=256) :: config
integer :: userxy
integer :: i, npoly

npoly = 5

nargs = iargc()
if (nargs.eq.1) then
  call getarg(1, config)
else
  write(*, '(A)') "IMCOM [usage]: "
  write(*, '(A)') " ./imcom <config_file>"
  stop
endif

! Welcome banner
write(*, FMT='(A)') "IMCOM: IMage COMbination v1.0"! (B. Rowe & C. Hirata 2010)"
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
do i=1, nexp

   call imcom_rotate_image(G_unrot(:, :, i), rotangdeg(i), npoly, 1, &
                           G_rot(:, :, i))
   call imcom_ft_r2c(n1psf, n2psf, G_rot(:, :, i), Gt_rot(:, :, i), 0)

end do
! Calculate the Fourier transform Gammat(ux, uy) of each Gamma(y, y)
call imcom_ft_r2c(n1psf, n2psf, Gamma(:, :), Gammat(:, :), 0)

! Output some checks
call imcom_writefits('./fits/ux_check.fits', n1psf, n2psf, ux)
call imcom_writefits('./fits/uy_check.fits', n1psf, n2psf, uy)
call imcom_writefits('./fits/mtf1.real.fits', n1psf, n2psf, real(Gt_rot(:, :, 1), 8))
call imcom_writefits('./fits/mtf1.imag.fits', n1psf, n2psf, dimag(Gt_rot(:, :, 1)))
call imcom_writefits('./fits/mtf1.abs.fits', n1psf, n2psf, abs(Gt_rot(:, :, 1)))

! Make the (diagonal only!) noise matrix
call imcom_build_kN

! Read / build the A array
call imcom_get_A(npoly, 1, 0)        ! arg#2 = save, arg#3 = force build
! Read / build the B array
call imcom_get_B(npoly, 1, 0)
! Calulate C (really quick)
call imcom_calc_C
! Read / build the T array
call imcom_get_T(1, 0)

! Build the output products using the galaxy images, A, B, C and our calculated T transofrmation
call imcom_build_image
call imcom_build_noise
call imcom_build_U  ! The leakage function map

END PROGRAM imcom
