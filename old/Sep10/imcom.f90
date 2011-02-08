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
integer :: alstat, i, n1gam, n2gam, npoly

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
write(*, FMT='(A)') "IMCOM: IMage COMbination v0.1"! (B. Rowe & C. Hirata 2010)"
! Read the config files
call imcom_read_config(config, userxy)

! Read and allocate the input images
call imcom_alloc_gims(n1gim, n2gim)
if (userxy.eq.1) then
! Read the x-y arrays
  call imcom_alloc_xy(n1gim, n2gim)
  call imcom_alloc_XaYa(n1out, n2out)
  write(*, FMT='(A)') "IMCOM: HACK - creating GIMXSCALE and GIMYSCALE as 10xPSFSCALE"
  allocate(gimxscale(nexp), gimyscale(nexp), STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate space for GIMXSCALE and GIMYSCALE arrays"
    stop
  endif
  gimxscale = psfxscale * 10.d0
  gimyscale = psfyscale * 10.d0

else
! Make the x-y arrays
  call imcom_build_xy
  call imcom_build_XaYa
end if

! Read and allocate the PSFs
call imcom_alloc_psfs(n1psf, n2psf)
! Read and allocate the Gamma PSF
call imcom_alloc_gamma(n1gam, n2gam)
if ((n1psf.ne.n1gam).or.(n2psf.ne.n2gam)) then
  write(*, FMT='(A)') "IMCOM ERROR: PSF and GAMMA images must be same size"
  stop
end if
! Make the ux, uy FT arrays
call imcom_build_uxuy
! Rotate PSFs and calculate the Fourier transform Gt(ux, uy) of each
write(*, FMT='(A)') "IMCOM: Rotating and Fourier transforming PSF images"
do i=1, nexp

   call imcom_rotate_image(G_unrot(:, :, i), rotangdeg(i), npoly, 1, G_rot(:, :, i))
   call imcom_ft_r2c(n1psf, n2psf, G_rot(:, :, i), Gt_rot(:, :, i), 0)
!  call imcom_ft_r2c(n1psf, n2psf, G_unrot(:, :, i), Gt_unrot(:, :, i), 0)
!  call imcom_rotate_cimage(cshift(cshift(Gt_unrot(:, :, i), n2psf / 2, 2), n1psf / 2, 1), rotangdeg(i), Gt_rot(:, :, i), 1)
!  Gt_rot(:, :, i) = cshift(cshift(Gt_rot(:, :, i), -n2psf / 2, 2), -n1psf / 2, 1)

end do
! Calculate the Fourier transform Gammat(ux, uy) of each Gamma(y, y)
call imcom_ft_r2c(n1psf, n2psf, Gamma(:, :), Gammat(:, :), 0)

! Output some checks
call imcom_writefits('./fits/ux_check.fits', n1psf, n2psf, ux)
call imcom_writefits('./fits/uy_check.fits', n1psf, n2psf, uy)
call imcom_writefits('./fits/mtf1.real.fits', n1psf, n2psf, real(Gt_rot(:, :, 1), 8))
call imcom_writefits('./fits/mtf1.imag.fits', n1psf, n2psf, dimag(Gt_rot(:, :, 1)))
call imcom_writefits('./fits/mtf1.abs.fits', n1psf, n2psf, abs(Gt_rot(:, :, 1)))

call imcom_writefits('./fits/mtf2.real.fits', n1psf, n2psf, real(Gt_rot(:, :, 2), 8))
call imcom_writefits('./fits/mtf2.imag.fits', n1psf, n2psf, real(dcmplx(0.d0, -1.d0) * Gt_rot(:, :, 2), 8))
call imcom_writefits('./fits/mtf2.abs.fits', n1psf, n2psf, real(Gt_rot(:, :, 2) * dconjg(Gt_rot(:, :, 2)), 8))

call imcom_writefits('./fits/mtf3.real.fits', n1psf, n2psf, real(Gt_rot(:, :, 3), 8))
call imcom_writefits('./fits/mtf3.imag.fits', n1psf, n2psf, real(dcmplx(0.d0, -1.d0) * Gt_rot(:, :, 3), 8))
call imcom_writefits('./fits/mtf3.abs.fits', n1psf, n2psf, real(Gt_rot(:, :, 3) * dconjg(Gt_rot(:, :, 3)), 8))

call imcom_writefits('./fits/mtfgam.real.fits', n1psf, n2psf, real(Gammat(:, :), 8))
call imcom_writefits('./fits/mtfgam.imag.fits', n1psf, n2psf, real(dcmplx(0.d0, -1.d0) * Gammat(:, :), 8))
call imcom_writefits('./fits/mtfgam.abs.fits', n1psf, n2psf, sqrt(real(Gammat(:, :) * dconjg(Gammat(:, :)), 8)))


! Make the (diagonal only!) noise matrix
call imcom_build_kN
call imcom_get_Al(npoly, 1)
call imcom_calc_C

stop
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


call imcom_calc_U   ! make (brute force) leakage function

!deallocate(G, Gt, x_i, y_i, STAT=dealstat)

END PROGRAM imcom
