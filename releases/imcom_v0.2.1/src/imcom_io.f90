!    IMCOM_IO.F90 - Fortran 95 module that handles basic I/O for the
!                   prototype IMCOM package
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
!                   FITSIO/CFITSIO LIBRARIES
!
!    CONTAINS
!    Routines for IMCOM that include FITSIO wrapper routines; usage error 
!    handling; reading of input config files; input, allocation and
!    bad-pixel masking of image and x-y arrays; command line parsing etc.
!
MODULE imcom_io

USE imcom_data
USE imcom_proc

CONTAINS

!---

SUBROUTINE imcom_alloc_gims
! Explore & read in the galaxy images, then allocate the Im arrays in which
! to store them
implicit none
integer, dimension(:), allocatable :: n1, n2, bitpix
integer :: i, imask, alstat, dealstat

allocate(n1(nexp), n2(nexp), bitpix(nexp), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for image property arrays"
  stop
end if
do i=1, nexp

  call imcom_sizefits(trim(gimfile(i)), n1(i), n2(i), bitpix(i))
  if (i.ne.1) then
    if (n1(i).ne.(n1(i-1)).or.n2(i).ne.(n2(i-1))) then
      write(*, FMT='(A)') "IMCOM ERROR: Variable size input images not supported in this version"
      stop
    end if
  end if

end do
! Then set the global variables for the image size
n1gim = n1(1)
n2gim = n2(1)
! Then allocate the array for the images
allocate(Im_unmasked(n1gim, n2gim, nexp), exp_unmasked(n1gim, n2gim, nexp), &
         STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for image arrays"
  stop
end if
Im_unmasked = 0.d0
do i=1, nexp

  write(*, FMT='(A)') "IMCOM: Reading image from "//trim(gimfile(i))
  call imcom_readfits(trim(gimfile(i)), n1gim, n2gim, Im_unmasked(:, :, i))
  exp_unmasked(:, :, i) = i

end do
n_unmasked = n1gim * n2gim * nexp
allocate(I_unmasked_i(n_unmasked), exp_unmasked_i(n_unmasked), &
         mask_i(n_unmasked), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for image arrays"
  stop
end if
I_unmasked_i = reshape(Im_unmasked, (/ n_unmasked /))
exp_unmasked_i = reshape(exp_unmasked, (/ n_unmasked /))
deallocate(Im_unmasked, exp_unmasked, STAT=dealstat)
if (dealstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot deallocate memory for image arrays"
  stop
end if
! Define a 1D mask wherever the image is less than the saturation value
where(I_unmasked_i.lt.saturation)
  mask_i = 1
elsewhere
  mask_i = 0
end where
! Define the n value:
n = sum(mask_i)
! Allocate 1D mask and mask subscript arrays
allocate(mask_subs_i(n), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for image mask arrays"
  stop
end if
imask = 0   ! start mask counter
! Build the subscript array - is there a quicker way to do this?
do i=1, n_unmasked

  if (mask_i(i).eq.1) then
    imask = imask + 1
    mask_subs_i(imask) = i
  end if

end do
! Allocate the 1D masked image array...
allocate(I_i(n), exp_i(n), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for masked image array"
  stop
end if
! Then construct it..
forall(i=1: n) I_i(i) = I_unmasked_i(mask_subs_i(i))
forall(i=1: n) exp_i(i) = exp_unmasked_i(mask_subs_i(i))
deallocate(n1, n2, bitpix, I_unmasked_i, exp_unmasked_i, STAT=dealstat)
if (dealstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot de allocate memory for image property arrays"
  stop
end if
END SUBROUTINE imcom_alloc_gims

!---

SUBROUTINE imcom_alloc_xy
! Explore & read in the galaxy images, then allocate the Im arrays in which
! to store them
implicit none
integer, dimension(:), allocatable :: n1x, n2x, n1y, n2y, bitpix
integer :: i, alstat, dealstat

allocate(n1x(nexp), n2x(nexp), n1y(nexp), n2y(nexp), bitpix(nexp), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for image property arrays"
  stop
end if
do i=1, nexp

  call imcom_sizefits(trim(gimxfile(i)), n1x(i), n2x(i), bitpix(i))
  call imcom_sizefits(trim(gimyfile(i)), n1y(i), n2y(i), bitpix(i))
  if (i.ne.1) then
    if ((n1x(i).ne.n1gim).or.(n2x(i).ne.n2gim)) then
      write(*, FMT='(A)') "IMCOM ERROR: X coord image not matching input image"
      stop
    end if
    if ((n1y(i).ne.n1gim).or.(n2y(i).ne.n2gim)) then
      write(*, FMT='(A)') "IMCOM ERROR: Y coord image not matching input image"
      stop
    end if
  end if

end do
! Then allocate the initial storage arrays
allocate(x_unmasked(n1gim, n2gim, nexp), y_unmasked(n1gim, n2gim, nexp), &
         STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for image X & Y arrays"
  stop
end if
do i=1, nexp

  write(*, FMT='(A)') "IMCOM: Reading X coordinates from "//trim(gimxfile(i))
  call imcom_readfits(trim(gimxfile(i)), n1gim, n2gim, x_unmasked(:, :, i))
  write(*, FMT='(A)') "IMCOM: Reading Y coordinates from "//trim(gimyfile(i))
  call imcom_readfits(trim(gimyfile(i)), n1gim, n2gim, y_unmasked(:, :, i))

end do
deallocate(n1x, n2x, n1y, n2y, bitpix, STAT=dealstat)
if (dealstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot de allocate memory for image x_i & y_i property arrays"
  stop
end if
allocate(x_unmasked_i(n_unmasked), y_unmasked_i(n_unmasked), &
         x_i(n), y_i(n), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot allocate memory for x_i, y_i"
  stop
endif
x_unmasked_i = reshape(x_unmasked, (/ n_unmasked /))
y_unmasked_i = reshape(y_unmasked, (/ n_unmasked /))
! Then construct masked x and y..
forall(i=1: n) x_i(i) = x_unmasked_i(mask_subs_i(i))
forall(i=1: n) y_i(i) = y_unmasked_i(mask_subs_i(i))
deallocate(x_unmasked, y_unmasked, x_unmasked_i, y_unmasked_i, STAT=dealstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot deallocate memory for x, y"
  stop
endif
END SUBROUTINE imcom_alloc_xy

!---

SUBROUTINE imcom_alloc_XaYa
! Explore & read in the output XY arrays
implicit none
integer :: n1x, n2x, n1y, n2y, bitpix
integer :: alstat, dealstat

call imcom_sizefits(trim(outxfile), n1x, n2x, bitpix)
call imcom_sizefits(trim(outyfile), n1y, n2y, bitpix)
if ((n1x.ne.n1y).or.(n2x.ne.n2y)) then
  write(*, FMT='(A)') "IMCOM ERROR: Size of X_a array "//trim(outxfile)//&
                      " does not match Y_a array "//trim(outyfile)
  stop
endif
n1out = n1x
n2out = n2x
! Then allocate the 2D storage arrays
allocate(Xa(n1out, n2out), Ya(n1out, n2out), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for image X & Y arrays"
  stop
end if
write(*, FMT='(A)') "IMCOM: Reading X coordinates from "//trim(outxfile)
call imcom_readfits(trim(outxfile), n1out, n2out, Xa(:, :))
write(*, FMT='(A)') "IMCOM: Reading Y coordinates from "//trim(outyfile)
call imcom_readfits(trim(outyfile), n1out, n2out, Ya(:, :))
! From the information available define total size of output array
m = n1out * n2out
allocate(X_a(m), Y_a(m), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot allocate memory for X_a & Y_a"
  stop
endif
X_a = reshape(Xa, (/ m /))
Y_a = reshape(Ya, (/ m /))
deallocate(Xa, Ya, STAT=dealstat)
if (alstat.ne.0) then
  write(*, FMT='(A)')"IMCOM ERROR: Cannot deallocate memory for Xa, Ya"
  stop
endif
END SUBROUTINE imcom_alloc_XaYa

!---

SUBROUTINE imcom_alloc_psfs(n1pad, n2pad)
! Explore & read in the PSF files, pad them to sufficient size, 
! then allocate the G and Gt arrays to store them
implicit none
integer, intent(OUT) :: n1pad, n2pad
integer, dimension(:), allocatable :: n1, n2, bitpix
real(KIND=8), dimension(:, :), allocatable :: Gtemp
integer :: i, alstat, dealstat
integer :: n1min, n2min, n1max, n2max

allocate(n1(nexp + 1), n2(nexp + 1), bitpix(nexp + 1), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for PSF image property arrays"
  stop
end if
do i=1, nexp

  call imcom_sizefits(trim(psffile(i)), n1(i), n2(i), bitpix(i))
  if (i.ne.1) then
    if (n1(i).ne.(n1(i-1)).or.n2(i).ne.(n2(i-1))) then
      write(*, FMT='(A)') "IMCOM ERROR: Variable size PSF images not supported in this version"
      stop
    end if
  end if

end do
call imcom_sizefits(trim(gamfile), n1(nexp + 1), n2(nexp + 1), bitpix(nexp + 1))
if (n1(nexp + 1).ne.(n1(nexp)).or.n2(nexp + 1).ne.(n2(nexp))) then
  write(*, FMT='(A)') "IMCOM ERROR: Variable size PSF images not supported in this version"
  write(*, FMT='(A)') "IMCOM ERROR: Gamma image different size from input PSF images"
  stop
end if
!
! Zero-pad arrays to at least twice the size of the input image, or double current size, whichever is larger. Set n1pad = n2pad
n1pad = 2 * nint((max(maxval(x_i), maxval(X_a)) - min(minval(x_i), minval(X_a))) / psfxscale)
n2pad = 2 * nint((max(maxval(y_i), maxval(Y_a)) - min(minval(y_i), minval(Y_a))) / psfyscale)
n1pad = max(n1pad, n2pad)
n2pad = n1pad
if (n1pad.lt.(n1(1) * 2).or.(n2pad.lt.(n2(1) * 2))) then
  n1pad = max(n1(1) * 2, n2(1) * 2)
  n2pad = n1pad
end if
! Then allocate the final storage array for the PSF and its transform
allocate(G_unrot(n1pad, n2pad, nexp), G_rot(n1pad, n2pad, nexp),   &
         Gt_unrot(n1pad, n2pad, nexp), Gt_rot(n1pad, n2pad, nexp), &
         Gamma(n1pad, n2pad), Gammat(n1pad, n2pad), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for PSF arrays"
  stop
end if
G_unrot = 0.d0
G_rot = 0.d0
Gamma = 0.d0
Gt_unrot = dcmplx(0.d0, 0.d0)
Gt_rot = dcmplx(0.d0, 0.d0)
Gammat = dcmplx(0.d0, 0.d0)
allocate(Gtemp(n1(1), n2(1)), STAT=alstat)  ! all n(i) the same...
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for temporary PSF image array"
  stop
end if
do i=1, nexp

  Gtemp = 0.d0
  n1min = (n1pad - n1(i)) / 2 + 1
  n2min = (n2pad - n2(i)) / 2 + 1
  n1max = (n1pad + n1(i)) / 2
  n2max = (n2pad + n2(i)) / 2
  write(*, FMT='(A)') "IMCOM: Reading PSF from "//trim(psffile(i))
  call imcom_readfits(trim(psffile(i)), n1(i), n2(i), Gtemp)
  G_unrot(n1min:n1max, n2min:n2max, i) = Gtemp / sum(Gtemp)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot deallocate memory for temporary PSF image array"
    stop
  end if

end do
Gtemp = 0.d0
n1min = (n1pad - n1(nexp +1)) / 2 + 1
n2min = (n2pad - n2(nexp +1)) / 2 + 1
n1max = (n1pad + n1(nexp +1)) / 2
n2max = (n2pad + n2(nexp +1)) / 2
write(*, FMT='(A)') "IMCOM: Reading ouput PSF from "//trim(gamfile)
call imcom_readfits(trim(gamfile), n1(nexp +1), n2(nexp +1), Gtemp)
Gamma(n1min:n1max, n2min:n2max) = Gtemp / sum(Gtemp)
deallocate(Gtemp, n1, n2, bitpix, STAT=dealstat)
if (dealstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot deallocate memory for PSF image property arrays"
  stop
end if
END SUBROUTINE imcom_alloc_psfs

!---

SUBROUTINE imcom_sizefits(filename, naxis1, naxis2, bitpix)
! Gets important information from a fits file without reading it
!
implicit none
character(LEN=*), intent(IN) :: filename
integer, intent(OUT) :: naxis1, naxis2, bitpix
integer :: unit, stat, naxis
integer, parameter :: rwmode = 0    ! Read only for this subroutine
integer, parameter :: maxdim = 2    ! Maximum dimension for images, 2D only
integer, dimension(maxdim) :: naxes

stat = 0                                             ! Initialize the status
call FTGIOU(unit, stat)                              ! Get a unit...
call FTIOPN(unit, trim(filename), rwmode, stat)      ! Open the .fits file 
                                                     ! (image specific call) 
call FTGIPR(unit, maxdim, bitpix, naxis, naxes, stat) ! Get image size etc.
naxis1 = naxes(1)
naxis2 = naxes(2)
call FTCLOS(unit, stat)       ! Close the file
call FTFIOU(unit, stat)       ! Free the I/O unit
call FTRPRT('STDERR', stat)   ! Error report
END SUBROUTINE imcom_sizefits

!---

SUBROUTINE imcom_readfits(filename, naxis1, naxis2, image)
! Reads in a fits file from filename, outputs the array dimensions and
! a pointer to the location of the image array (1D) in memory
implicit none
character(LEN=*), intent(IN) :: filename
integer, intent(IN) :: naxis1, naxis2
real(KIND=8), dimension(naxis1, naxis2), intent(OUT) :: image
integer, parameter :: rwmode = 0    ! Read only for this subroutine
integer, parameter :: maxdim = 2    ! Maximum dimension for images, 2D only
real(KIND=8), parameter :: nullvar = -66.D6 ! Assigned to undefined pixel values
real(KIND=8), dimension(naxis1 * naxis2) :: imtemp
integer :: unit, stat, bitpix, naxis
integer, dimension(maxdim) :: naxes
logical :: anyf
integer :: group, firstpix

stat = 0                                         ! Initialize the status
call FTGIOU(unit, stat)                          ! Get a unit...
call FTIOPN(unit, trim(filename), rwmode, stat)  ! Open the .fits file 
                                                 ! (image specific call) 
call FTGIPR(unit, maxdim, bitpix, naxis, naxes, stat) ! Get image size etc.
if ((bitpix.ne.-64)) then
  write(*, '(A)') "IMCOM WARNING: Input fits file not 64-bit floating point data... will probably end badly!"
end if
if (naxis1.ne.naxes(1)) then
  write(*, FMT='(A)') "IMCOM ERROR: Specified NAXIS1 does not match up to image in "//trim(filename) 
  stop
endif
if (naxis2.ne.naxes(2)) then
  write(*, FMT='(A)') "IMCOM ERROR: Specified NAXIS2 does not match up to image in "//trim(filename) 
  stop
endif
group = 1
firstpix = 1
imtemp = 0.d0
! Read in image from primary data unit
call FTGPVD(unit, group, firstpix, naxis1 * naxis2, nullvar, &
            imtemp(1:naxis1 * naxis2), anyf, stat)
image = reshape(imtemp, (/ naxis1, naxis2 /))
call FTCLOS(unit, stat)       ! Close the file
call FTFIOU(unit, stat)       ! Free the I/O unit
call FTRPRT('STDERR', stat)   ! Error report
END SUBROUTINE imcom_readfits

!---

SUBROUTINE imcom_writefits(filename, naxis1, naxis2, image)
! Wrapper for writing fits files
!
implicit none
character(LEN=*), intent(IN) :: filename
integer, intent(IN) :: naxis1, naxis2
real(KIND=8), dimension(naxis1, naxis2), intent(IN) :: image
real(KIND=8), dimension(:), allocatable :: imvec
integer, parameter :: rwmode = 1    ! Write mode for this subroutine
integer, parameter :: maxdim = 2    ! Maximum dimension for images, 2D only
integer :: unit, stat, alstat
integer, parameter :: naxis = maxdim
integer, parameter :: blocksize = 1
integer, parameter :: bitpix = -64  ! DOUBLE precision
integer, dimension(maxdim) :: naxes
integer :: group, firstpix, exists

stat = 0                                         ! Initialize the status
naxes = (/ naxis1, naxis2 /)                     ! Define naxes vector
call FTEXIST(trim(filename), exists, stat)
if (exists.eq.1) then
  call FTGIOU(unit, stat)                          ! Get a unit...
  call FTIOPN(unit, trim(filename), rwmode, stat)  ! Open the .fits file 
                                                   ! (image specific call)
  call FTDELT(unit, stat) ! Ruthlessly delete the old file!
endif
call FTGIOU(unit, stat)                            ! Get a new unit...
call FTINIT(unit, trim(filename), blocksize, stat) ! Open the .fits file 
                                                   ! (image specific call) 
call FTIIMG(unit, bitpix, naxis, naxes, stat)
group = 1
firstpix = 1
allocate(imvec(naxis1 * naxis2), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate space for FITS output array, too large?"
  write(*, FMT='(A)') "IMCOM ERROR: "
  write(*, FMT='(A, I4)') "IMCOM ERROR: ALSTAT = ", alstat
  stop
end if
imvec = reshape(image, (/naxis1 * naxis2/))
call FTPPRD(unit, group, firstpix, naxis1 * naxis2, imvec, stat)
deallocate(imvec)
call FTCLOS(unit, stat)       ! Close the file
call FTFIOU(unit, stat)       ! Free the I/O unit
call FTRPRT('STDERR', stat)   ! Error report
END SUBROUTINE imcom_writefits

!---

SUBROUTINE imcom_writefitscube(filename, naxis1, naxis2, naxis3, dcube)
! Wrapper for writing fits files
!
implicit none
character(LEN=*), intent(IN) :: filename
integer, intent(IN) :: naxis1, naxis2, naxis3
real(KIND=8), dimension(naxis1, naxis2, naxis3), intent(IN) :: dcube
integer, parameter :: rwmode = 1    ! Write mode for this subroutine
integer, parameter :: maxdim = 3    ! Maximum dimension for images, 3D only
integer :: unit, stat
integer, parameter :: naxis = maxdim
integer, parameter :: blocksize = 1
integer, parameter :: bitpix = -64  ! DOUBLE precision
integer, dimension(maxdim) :: naxes
integer :: group, firstpix, exists

stat = 0                                         ! Initialize the status
naxes = (/ naxis1, naxis2, naxis3 /)                     ! Define naxes vector
call FTEXIST(trim(filename), exists, stat)
if (exists.eq.1) then
  call FTGIOU(unit, stat)                          ! Get a unit...
  call FTIOPN(unit, trim(filename), rwmode, stat)  ! Open the .fits file 
                                                   ! (image specific call)
  call FTDELT(unit, stat) ! Ruthlessly delete the old file!
endif
call FTGIOU(unit, stat)                            ! Get a new unit...
call FTINIT(unit, trim(filename), blocksize, stat) ! Open the .fits file 
                                                   ! (image specific call) 
call FTIIMG(unit, bitpix, naxis, naxes, stat)
group = 1
firstpix = 1
call FTPPRD(unit, group, firstpix, naxis1 * naxis2 * naxis3, &
            reshape(dcube, (/ naxis1 * naxis2 * naxis3 /)), stat)
call FTCLOS(unit, stat)       ! Close the file
call FTFIOU(unit, stat)       ! Free the I/O unit
call FTRPRT('STDERR', stat)   ! Error report
END SUBROUTINE imcom_writefitscube

!---

SUBROUTINE imcom_read_config(filename, userxy)
! Reads the main config file, outputs and uses the USERXY parameter
implicit none
character(LEN=*), intent(IN) :: filename
integer, intent(OUT) :: userxy
character(LEN=256) :: scratch
integer :: iostat, alstat, i

open(UNIT=10, FILE=trim(filename), STATUS="OLD", FORM="FORMATTED", &
     ACTION="READ", IOSTAT=iostat)
if (iostat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot open "//trim(filename)
  write(*, FMT='(A)') "IMCOM: Check filename and path"
  stop
end if
read(UNIT=10, FMT=*, IOSTAT=iostat) scratch, psfxscale
if (iostat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(filename)
  write(*, FMT='(A)') "IMCOM: Expecting: PSFXSCALE  <psfxscale>"
  write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
  write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
  write(*, FMT='(A)') "IMCOM: Check <psfxscale> floating point"
  stop
end if
read(UNIT=10, FMT=*, IOSTAT=iostat) scratch, psfyscale
if (iostat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(filename)
  write(*, FMT='(A)') "IMCOM: Expecting: PSFYSCALE  <psfyscale>"
  write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
  write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
  write(*, FMT='(A)') "IMCOM: Check <psfyscale> floating point"
  stop
end if
read(UNIT=10, FMT=*, IOSTAT=iostat) scratch, nexp
if (iostat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Incorrect first-line format for "//trim(filename)
  write(*, FMT='(A)') "IMCOM: Expecting: NEXP      <integer>"
  write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
  stop
end if
read(UNIT=10, FMT=*, IOSTAT=iostat) scratch, userxy
if (iostat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Incorrect second-line format for "//trim(filename)
  write(*, FMT='(A)') "IMCOM: Expecting: USERXY    <integer>"
  write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
  stop
end if
allocate(inconfig(nexp), STAT=alstat)
if (alstat.ne.0) then
   write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for input image filenames"
   stop
end if 
do i=1, nexp

  read(UNIT=10, FMT='(A12, A)', IOSTAT=iostat) scratch, inconfig(i)
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(filename)
    write(*, FMT='(A)') "IMCOM: Expecting: INCONFIG  <filename>"
    write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
    write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
    write(*, FMT='(A)') "IMCOM: Check <filename>: must be 256 characters or less"
    stop
  end if  

end do
! Then read in each of these config files... (these all read from IO UNIT=20)
call imcom_read_inconfigs(userxy)

! Then continuing reading this main config file...
read(UNIT=10, FMT='(A12, A)', IOSTAT=iostat) scratch, outconfig
if (iostat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(filename)
  write(*, FMT='(A)') "IMCOM: Expecting: OUTCONFIG <filename>"
  write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
  write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
  write(*, FMT='(A)') "IMCOM: Check <filename>: must be 256 characters or less"
  stop
end if

! Then read in the outconfig file (read from IO UNIT=30)
call imcom_read_outconfig(userxy)

! Optional read-in of A, B, filenames (sets to default if EOF encountered)
read(UNIT=10, FMT='(A12, A)', IOSTAT=iostat) scratch, Afile
if (iostat.lt.0) then
  Afile = Adeffile
  Bfile = Bdeffile
  Qfile = Qdeffile
  Lfile = Ldeffile
  Pfile = Pdeffile
  close(UNIT=10, IOSTAT=iostat)
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Unable to close "//trim(filename)
    stop
  end if
  return
else if (iostat.gt.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(filename)
  write(*, FMT='(A)') "IMCOM: Expecting: AFILE     <filename>"
  write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
  write(*, FMT='(A)') "IMCOM: Check <filename>: must be 256 characters or less"
  stop
end if

read(UNIT=10, FMT='(A12, A)', IOSTAT=iostat) scratch, Bfile
if (iostat.lt.0) then
  Bfile = Bdeffile
  Qfile = Qdeffile
  Lfile = Ldeffile
  Pfile = Pdeffile
  close(UNIT=10, IOSTAT=iostat)
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Unable to close "//trim(filename)
    stop
  end if
  return
else if (iostat.gt.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(filename)
  write(*, FMT='(A)') "IMCOM: Expecting: BFILE     <filename>"
  write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
  write(*, FMT='(A)') "IMCOM: Check <filename>: must be 256 characters or less"
  stop
end if

read(UNIT=10, FMT='(A12, A)', IOSTAT=iostat) scratch, Qfile
if (iostat.lt.0) then
  Qfile = Qdeffile
  Lfile = Ldeffile
  Pfile = Pdeffile
  close(UNIT=10, IOSTAT=iostat)
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Unable to close "//trim(filename)
    stop
  end if
  return
else if (iostat.gt.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(filename)
  write(*, FMT='(A)') "IMCOM: Expecting: QFILE     <filename>"
  write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
  write(*, FMT='(A)') "IMCOM: Check <filename>: must be 256 characters or less"
  stop
end if

read(UNIT=10, FMT='(A12, A)', IOSTAT=iostat) scratch, Lfile
if (iostat.lt.0) then
  Lfile = Ldeffile
  Pfile = Pdeffile
  close(UNIT=10, IOSTAT=iostat)
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Unable to close "//trim(filename)
    stop
  end if
  return
else if (iostat.gt.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(filename)
  write(*, FMT='(A)') "IMCOM: Expecting: LFILE     <filename>"
  write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
  write(*, FMT='(A)') "IMCOM: Check <filename>: must be 256 characters or less"
  stop
end if

read(UNIT=10, FMT='(A12, A)', IOSTAT=iostat) scratch, Pfile
if (iostat.lt.0) then
  Pfile = Pdeffile
  close(UNIT=10, IOSTAT=iostat)
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Unable to close "//trim(filename)
    stop
  end if
  return
else if (iostat.gt.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(filename)
  write(*, FMT='(A)') "IMCOM: Expecting: PFILE     <filename>"
  write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
  write(*, FMT='(A)') "IMCOM: Check <filename>: must be 256 characters or less"
  stop
end if

END SUBROUTINE imcom_read_config

!---

SUBROUTINE imcom_read_inconfigs(userxy)
! Read the configuration files for the input images and PSFs
implicit none
integer, intent(IN) :: userxy
character(LEN=256) :: scratch
integer :: iostat, alstat, i

allocate(gimfile(nexp), psffile(nexp), noise(nexp), rotangdeg(nexp), &
         STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Unable to allocate space for image filename arrays"
  stop
end if
! Slightly different setup for user/non-user defined X and Y arrays
if (userxy.eq.1) then
  allocate(gimxfile(nexp), gimyfile(nexp), STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Unable to allocate space for image X & Y coord filename arrays"
    stop
  endif
else
  allocate(gimxscale(nexp), gimyscale(nexp), dither(nexp, 2), STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Unable to allocate space for image X & Y coord filename arrays"
    stop
  endif
end if

! Then read in the info
do i=1, nexp

  open(UNIT=20, FILE=trim(adjustl(inconfig(i))), STATUS="OLD", &
       FORM="FORMATTED", ACTION="READ", IOSTAT=iostat)
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot open "//trim(inconfig(i))
    write(*, FMT='(A)') "IMCOM: Check filename and path"
    stop
  end if
  read(UNIT=20, FMT='(A12, A)', IOSTAT=iostat) scratch, psffile(i)
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(inconfig(i))
    write(*, FMT='(A)') "IMCOM: Expecting: PSFFILE   <filename>"
    write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
    write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
    write(*, FMT='(A)') "IMCOM: Check <filename>: must be 256 characters or less"
    stop
  end if    
  read(UNIT=20, FMT='(A12, A)', IOSTAT=iostat) scratch, gimfile(i)
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(inconfig(i))
    write(*, FMT='(A)') "IMCOM: Expecting: GIMFILE   <filename>"
    write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
    write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
    write(*, FMT='(A)') "IMCOM: Check <filename>: must be 256 characters or less"
    stop
  end if
  read(UNIT=20, FMT=*, IOSTAT=iostat) scratch, rotangdeg(i)
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(inconfig(i))
    write(*, FMT='(A)') "IMCOM: Expecting: ROTANGDEG <psfrotdeg>"
    write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
    write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
    write(*, FMT='(A)') "IMCOM: Check <rotangdeg> floating point"
    stop
  end if
  read(UNIT=20, FMT=*, IOSTAT=iostat) scratch, noise(i)
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(inconfig(i))
    write(*, FMT='(A)') "IMCOM: Expecting: NOISE   <noise>"
    write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
    write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
    write(*, FMT='(A)') "IMCOM: Check <noise> floating point"
    stop
  end if  
!
! Then read in galaxy image exposure information: format depends on setting 
! of USERXY = 0/1 in main config file
!
  if (userxy.eq.1) then
    read(UNIT=20, FMT='(A12, A)', IOSTAT=iostat) scratch, gimxfile(i)
    if (iostat.ne.0) then
      write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(inconfig(i))
      write(*, FMT='(A)') "IMCOM: Expecting: GIMXFILE  <filename>"
      write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
      write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
      write(*, FMT='(A)') "IMCOM: Check USERXY: must be zero if not supplying X & Y coordinate arrays"
      write(*, FMT='(A)') "IMCOM: Check <filename>: must be 256 characters or less"
      stop
    end if
    read(UNIT=20, FMT='(A12, A)', IOSTAT=iostat) scratch, gimyfile(i)
    if (iostat.ne.0) then
      write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(inconfig(i))
      write(*, FMT='(A)') "IMCOM: Expecting: GIMYFILE  <filename>"
      write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
      write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
      write(*, FMT='(A)') "IMCOM: Check USERXY: must be zero if not supplying X & Y coordinate arrays"
      write(*, FMT='(A)') "IMCOM: Check <filename>: must be 256 characters or less"
      stop
    end if
  else  ! Otherwise just read in xscale, yscale & dithers
    read(UNIT=20, FMT=*, IOSTAT=iostat) scratch, gimxscale(i)
    if (iostat.ne.0) then
      write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(inconfig(i))
      write(*, FMT='(A)') "IMCOM: Expecting: GIMXSCALE  <gimxscale>"
      write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
      write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
      write(*, FMT='(A)') "IMCOM: Check USERXY: must be 1 if supplying X & Y coordinate arrays"
      write(*, FMT='(A)') "IMCOM: Check <gimscale> floating point"
      stop
    end if  
    read(UNIT=20, FMT=*, IOSTAT=iostat) scratch, gimyscale(i)
    if (iostat.ne.0) then
      write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(inconfig(i))
      write(*, FMT='(A)') "IMCOM: Expecting: GIMYSCALE  <gimyscale>"
      write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
      write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
      write(*, FMT='(A)') "IMCOM: Check USERXY: must be 1 if supplying X & Y coordinate arrays"
      write(*, FMT='(A)') "IMCOM: Check <gimyscale> floating point"
      stop
    end if  
    read(UNIT=20, FMT=*, IOSTAT=iostat) scratch, dither(i, 1), dither(i, 2)
    if (iostat.ne.0) then
      write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(inconfig(i))
      write(*, FMT='(A)') "IMCOM: Expecting: DITHER   <dither_x> <dither_y>"
      write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
      write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
      write(*, FMT='(A)') "IMCOM: Check USERXY: must be 1 if supplying X & Y coordinate arrays"
      write(*, FMT='(A)') "IMCOM: Check <dither_x>, <dither_y> floating point"
      stop
    end if
  end if
  close(UNIT=20, IOSTAT=iostat)
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Unable to close "//trim(inconfig(i))
    stop
  end if

end do
END SUBROUTINE imcom_read_inconfigs

!---

SUBROUTINE imcom_read_outconfig(userxy)
! Read the configuration files for the output images and Gamma PSFs
implicit none
integer, intent(IN) :: userxy
character(LEN=256) :: scratch
integer :: iostat

open(UNIT=30, FILE=trim(outconfig), STATUS="OLD", FORM="FORMATTED", &
     ACTION="READ", IOSTAT=iostat)
if (iostat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot open "//trim(outconfig)
  write(*, FMT='(A)') "IMCOM: Check filename and path"
  stop
end if
read(UNIT=30, FMT='(A12, A)', IOSTAT=iostat) scratch, gamfile
if (iostat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(outconfig)
  write(*, FMT='(A)') "IMCOM: Expecting: GAMFILE   <filename>"
  write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
  write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
  write(*, FMT='(A)') "IMCOM: Check <filename>: must be 256 characters or less"
  stop
end if
read(UNIT=30, FMT='(A12, A)', IOSTAT=iostat) scratch, Hfile
if (iostat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(outconfig)
  write(*, FMT='(A)') "IMCOM: Expecting: HFILE     <filename>"
  write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
  write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
  write(*, FMT='(A)') "IMCOM: Check <filename>: must be 256 characters or less"
  stop
end if
read(UNIT=30, FMT='(A12, A)', IOSTAT=iostat) scratch, Kfile
if (iostat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(outconfig)
  write(*, FMT='(A)') "IMCOM: Expecting: KFILE     <filename>"
  write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
  write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
  write(*, FMT='(A)') "IMCOM: Check <filename>: must be 256 characters or less"
  stop
end if
read(UNIT=30, FMT='(A12, A)', IOSTAT=iostat) scratch, Tfile
if (iostat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(outconfig)
  write(*, FMT='(A)') "IMCOM: Expecting: TFILE     <filename>"
  write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
  write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
  write(*, FMT='(A)') "IMCOM: Check <filename>: must be 256 characters or less"
  stop
end if
read(UNIT=30, FMT='(A12, A)', IOSTAT=iostat) scratch, Sfile
if (iostat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(outconfig)
  write(*, FMT='(A)') "IMCOM: Expecting: SFILE     <filename>"
  write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
  write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
  write(*, FMT='(A)') "IMCOM: Check <filename>: must be 256 characters or less"
  stop
end if
read(UNIT=30, FMT='(A12, A)', IOSTAT=iostat) scratch, Ufile
if (iostat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(outconfig)
  write(*, FMT='(A)') "IMCOM: Expecting: UFILE     <filename>"
  write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
  write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
  write(*, FMT='(A)') "IMCOM: Check <filename>: must be 256 characters or less"
  stop
end if
!
!  Then read in user-supplied X,Y, or simply the scales and dither (for regulat rectilinear grids)
!
if (userxy.eq.1) then
  read(UNIT=30, FMT='(A12, A)', IOSTAT=iostat) scratch, outxfile
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(outconfig)
    write(*, FMT='(A)') "IMCOM: Expecting: OUTXFILE  <filename>"
    write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
    write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
    write(*, FMT='(A)') "IMCOM: Check USERXY: must be zero if not supplying X & Y coordinate arrays"
    write(*, FMT='(A)') "IMCOM: Check <filename>: must be 256 characters or less"
    stop
  end if
  read(UNIT=30, FMT='(A12, A)', IOSTAT=iostat) scratch, outyfile
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(outconfig)
    write(*, FMT='(A)') "IMCOM: Expecting: OUTYFILE  <filename>"
    write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
    write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
    write(*, FMT='(A)') "IMCOM: Check USERXY: must be zero if not supplying X & Y coordinate arrays"
    write(*, FMT='(A)') "IMCOM: Check <filename>: must be 256 characters or less"
    stop
  end if
else
  read(UNIT=30, FMT=*, IOSTAT=iostat) scratch, outxscale
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(outconfig)
    write(*, FMT='(A)') "IMCOM: Expecting: OUTXSCALE  <outxscale>"
    write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
    write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
    write(*, FMT='(A)') "IMCOM: Check USERXY: must be 1 if supplying X & Y coordinate arrays"
    write(*, FMT='(A)') "IMCOM: Check <gimscale> floating point"
    stop
  end if  
  read(UNIT=30, FMT=*, IOSTAT=iostat) scratch, outyscale
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(outconfig)
    write(*, FMT='(A)') "IMCOM: Expecting: OUTYSCALE  <outyscale>"
    write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
    write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
    write(*, FMT='(A)') "IMCOM: Check USERXY: must be 1 if supplying X & Y coordinate arrays"
    write(*, FMT='(A)') "IMCOM: Check <gimscale> floating point"
    stop
  end if  
  read(UNIT=30, FMT=*, IOSTAT=iostat) scratch, outpos1, outpos2
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(outconfig)
    write(*, FMT='(A)') "IMCOM: Expecting: OUTPOS    <outpos_x> <outpos_y>"
    write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
    write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
    write(*, FMT='(A)') "IMCOM: Check <outpos_x>, <outpos_y> floating point"
    stop
  end if
  read(UNIT=30, FMT=*, IOSTAT=iostat) scratch, n1out, n2out
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(outconfig)
    write(*, FMT='(A)') "IMCOM: Expecting: NOUT      <n1out> <n2out>"
    write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
    write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)" 
    write(*, FMT='(A)') "IMCOM: Check <n1out>, <n2out> integer"
    stop
  end if
! use this information to calculate M
  m = n1out * n2out
end if
close(UNIT=30, IOSTAT=iostat)
if (iostat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Unable to close "//trim(outconfig)
  stop
end if

END SUBROUTINE imcom_read_outconfig

!---

SUBROUTINE imcom_welcome
implicit none
! Welcome banner
write(*, FMT='(A)') "IMCOM: IMage COMbination v0.2 (B. Rowe & C. Hirata 2010-2011)"
END SUBROUTINE

!---

SUBROUTINE imcom_usage_message
implicit none
call imcom_welcome
write(*, '(A)') "IMCOM: [usage]"
write(*, '(A)') " ./imcom <config_file> <U/S> <U/S_max> <U/S_tol> [...<forceT> <forceSys> <saturation>]"
write(*, '(A)') " "
write(*, '(A)') "  <config_file> [string]  : Global config file containing image locations etc."
!write(*, '(A)') " "
write(*, '(A)') "  <U/S> [string]          : U/u or S/s to specify minimization of U or S in output image"
!write(*, '(A)') " "
write(*, '(A)') "  <U/S_max> [dbl]         : Required maximum U or S"
!write(*, '(A)') " "
write(*, '(A)') "  <U/S_tol> [dbl]         : Absolute tolerance on U or S for interval bisection"
!write(*, '(A)') " "
write(*, '(A)') "  [ ...<forceT> [int]     : 1 = force build T matrix , 0 = not (default)                 ]"
write(*, '(A)') "  [ ...<forceSys> [int]   : 1 = force build system matrices A, B etc., 0 = not (default) ]"
write(*, '(A)') "  [ ...<saturation> [dbl] : Image saturation/bad pixel value (default = Float overflow)  ]"
write(*, '(A)') " "
END SUBROUTINE imcom_usage_message

!---

SUBROUTINE imcom_get_cmdline
implicit none
integer :: nargs
character(LEN=1) :: US_string
character(LEN=256) :: buffer
real(KIND=8), external :: DLAMCH
! Gets the command line arguments or prints the usage message
nargs = iargc()
if (nargs.ge.4) then
  call getarg(1, config)
  call getarg(2, US_string)
  if((US_string.eq."U").or.(US_string.eq."u")) then
    USB = .TRUE.    ! Solve for T_ia with U_a < US_min
  else if((US_string.eq."S").or.(US_string.eq."s")) then
    USB = .FALSE.   ! Solve for T_ia with S_a < US_min
  else
    call imcom_usage_message
    stop
  end if
  call getarg(3, buffer)
  read(buffer, FMT=*) US_max
  call getarg(4, buffer)
  read(buffer, FMT=*) US_tol
  ! Then set defaults for optional params
  forceT = 0
  forceSys = 0
  saturation = DLAMCH('O')   ! Largest number available to double 
                             ! precision...
  if (nargs.ge.5) then
    call getarg(5, buffer)
    read(buffer, FMT=*) forceT
    if (nargs.ge.6) then
      call getarg(6, buffer)
      read(buffer, FMT=*) forceSys
      if (nargs.eq.7) then
        call getarg(7, buffer)
        read(buffer, FMT=*) saturation
      end if
    end if
  end if
else
  call imcom_usage_message
  stop
endif
END SUBROUTINE imcom_get_cmdline

!---

END MODULE imcom_io
