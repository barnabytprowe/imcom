PROGRAM imcom

USE imcom_data
USE imcom_io
USE imcom_proc
USE imcom_matrix
USE omp_lib

implicit none
integer :: nargs
character(LEN=256) :: config
character(LEN=256) :: rdbuffer
integer :: alstat, dealstat, i
integer :: userxy, info

nargs = iargc()
if (nargs.eq.1) then
  call getarg(1, config)
else
  write(*, '(A)') "IMCOM [usage]: "
  write(*, '(A)') " ./imcom <config_file>"
  stop
endif

! Welcome banner
write(*, FMT='(A)') "IMCOM: IMage COMbination v0.1"! (B. Rowe & C. Hirata 2010)"
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
! Read and allocate the PSFs
call imcom_alloc_psfs(n1psf, n2psf)
! Make the ux, uy FT arrays
call imcom_build_uxuy
! Calculate the Fourier transform Gt(ux, uy) of each G(x, y)
write(*, FMT='(A)') "IMCOM: Fourier transforming PSF images"
do i=1, nexp

  call imcom_ft(n1psf, n2psf, G(:, :, i), Gt(:, :, i), 0)

end do
! Read and allocate the Gamma PSF
call imcom_alloc_gamma(n1gam, n2gam)
! Calculate the Fourier transform Gammat(ux, uy) of each Gamma(y, y)
call imcom_ft(n1gam, n2gam, Gamma(:, :), Gammat(:, :), 0)
! Make the X-Y arrays
call imcom_writefits('ux.fits', n1psf, n2psf, ux)
call imcom_writefits('uy.fits', n1psf, n2psf, uy)
call imcom_writefits('u.fits', n1psf, n2psf, sqrt(ux * ux + uy * uy))
call imcom_writefits('mtf.real.fits', n1psf, n2psf, real(Gt(:, :, 1), 8))
call imcom_writefits('mtf.imag.fits', n1psf, n2psf, real(dcmplx(0., -1.) * Gt(:, :, 1), 8))
call imcom_writefits('mtf.abs.fits', n1psf, n2psf, real(Gt(:, :, 1) * dconjg(Gt(:, :, 1)), 8))

! Make the (diagonal only!) noise matrix
call imcom_build_kN
! Read / build the A array
call imcom_get_A(1, 0)        ! arg#1 = save, arg#2 = force build
! Read / build the B array
call imcom_get_B(1, 0)
! Read / build the T array
call imcom_get_T(1, 0)

! Build the output products using the galaxy images, A, B, C and our calculated T transofrmation
call imcom_build_image
call imcom_build_noise
call imcom_build_U  ! The leakage function map

!deallocate(G, Gt, x_i, y_i, STAT=dealstat)

END PROGRAM imcom
