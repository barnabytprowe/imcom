!
! Routines for imcom.f95 that handle fitsio routines and error handling, 
! reading of config files, outputing plotting images etc.
!
!-----------------------------------------------------------------------

MODULE imcom_io

USE imcom_data

CONTAINS

!---

SUBROUTINE imcom_alloc_gims(n1dum, n2dum)
! Explore & read in the galaxy images, then allocate the Im arrays in which
! to store them
implicit none
integer, intent(OUT) :: n1dum, n2dum
integer, dimension(:), allocatable :: n1, n2, bitpix
integer :: i, alstat, dealstat

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
n1dum = n1(1)
n2dum = n2(1)
! Then allocate the final storage arrays for the PSF and its transform
allocate(Im(n1dum, n2dum, nexp), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for image arrays"
  stop
end if
Im = 0.
do i=1, nexp

  write(*, FMT='(A)') "IMCOM: Reading image from "//trim(gimfile(i))
  call imcom_readfits(trim(gimfile(i)), n1dum, n2dum, Im(:, :, i))

end do
! Define total size N of input array
n = n1dum * n2dum * nexp
deallocate(n1, n2, bitpix, STAT=dealstat)
if (dealstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot de allocate memory for image property arrays"
  stop
end if
END SUBROUTINE imcom_alloc_gims

!---

SUBROUTINE imcom_alloc_xy(n1dum, n2dum)
! Explore & read in the galaxy images, then allocate the Im arrays in which
! to store them
implicit none
integer, intent(IN) :: n1dum, n2dum ! note that these must now be INPUT and equal to n1, n2 as output by imcom_alloc_gims
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
    if (n1x(i).ne.(n1x(i-1)).or.n2x(i).ne.(n2x(i-1))) then
      write(*, FMT='(A)') "IMCOM ERROR: Variable size X coord images: not supported in this version"
      stop
    end if
    if (n1y(i).ne.(n1y(i-1)).or.n2y(i).ne.(n2y(i-1))) then
      write(*, FMT='(A)') "IMCOM ERROR: Variable size Y coord images: not supported in this version"
      stop
    end if
  end if

end do
if ((n1x(1).ne.n1dum).or.(n2x(1).ne.n2dum)) then
  write(*, FMT='(A)') "IMCOM ERROR: Size of image X array does not match input image size"
  stop
endif
if ((n1y(1).ne.n1dum).or.(n2y(1).ne.n2dum)) then
  write(*, FMT='(A)') "IMCOM ERROR: Size of image Y array does not match input image size"
  stop
endif
! Then allocate the final storage arrays
allocate(x(n1dum, n2dum, nexp), y(n1dum, n2dum, nexp), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for image X & Y arrays"
  stop
end if
x = 0.
y = 0.
do i=1, nexp

  write(*, FMT='(A)') "IMCOM: Reading X coordinates from "//trim(gimxfile(i))
  call imcom_readfits(trim(gimxfile(i)), n1dum, n2dum, x(:, :, i))
  write(*, FMT='(A)') "IMCOM: Reading Y coordinates from "//trim(gimyfile(i))
  call imcom_readfits(trim(gimyfile(i)), n1dum, n2dum, y(:, :, i))

end do
! Define total size N of input array
deallocate(n1x, n2x, n1y, n2y, bitpix, STAT=dealstat)
if (dealstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot de allocate memory for image X & Y property arrays"
  stop
end if
END SUBROUTINE imcom_alloc_xy

!---

SUBROUTINE imcom_alloc_XaYa(n1dum, n2dum)
! Explore & read in the galaxy images, then allocate the Im arrays in which
! to store them
implicit none
integer, intent(IN) :: n1dum, n2dum ! note that these must now be INPUT and equal to n1, n2 as output by imcom_alloc_gims
integer, dimension(:), allocatable :: n1x, n2x, n1y, n2y, bitpix
integer :: i, alstat, dealstat

allocate(n1x(nexp), n2x(nexp), n1y(nexp), n2y(nexp), bitpix(nexp), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for image property arrays"
  stop
end if
call imcom_sizefits(trim(gimxfile(i)), n1x(i), n2x(i), bitpix(i))
call imcom_sizefits(trim(gimyfile(i)), n1y(i), n2y(i), bitpix(i))
if ((n1x(1).ne.n1dum).or.(n2x(1).ne.n2dum)) then
  write(*, FMT='(A)') "IMCOM ERROR: Size of output Xa array does not match specified output image size"
  stop
endif
if ((n1y(1).ne.n1dum).or.(n2y(1).ne.n2dum)) then
  write(*, FMT='(A)') "IMCOM ERROR: Size of output Ya array does not match specified output image size"
  stop
endif

! Then allocate the final storage arrays for the PSF and its transform
allocate(Xa(n1dum, n2dum), Ya(n1dum, n2dum), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for image X & Y arrays"
  stop
end if
x = 0.
y = 0.
do i=1, nexp

  write(*, FMT='(A)') "IMCOM: Reading X coordinates from "//trim(gimxfile(i))
  call imcom_readfits(trim(gimxfile(i)), n1dum, n2dum, x(:, :, i))
  write(*, FMT='(A)') "IMCOM: Reading Y coordinates from "//trim(gimyfile(i))
  call imcom_readfits(trim(gimyfile(i)), n1dum, n2dum, y(:, :, i))

end do
! Define total size N of input array
deallocate(n1x, n1y, n2x, n2y, bitpix, STAT=dealstat)
if (dealstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot de allocate memory for image X & Y property arrays"
  stop
end if
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

allocate(n1(nexp), n2(nexp), bitpix(nexp), STAT=alstat)
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
n1pad = 2 * n1gim * anint(maxval(gimxscale / psfxscale))
n2pad = 2 * n2gim * anint(maxval(gimyscale / psfyscale))
if (n1pad.lt.n1(1)) n1pad = n1(1)
if (n2pad.lt.n2(1)) n2pad = n2(1)
! Then allocate the final storage arrays for the PSF and its transform
allocate(G(n1pad, n2pad, nexp), Gt(n1pad, n2pad, nexp), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for PSF arrays"
  stop
end if
G = 0.d0
Gt = 0.d0
do i=1, nexp

  if ((n1(i).gt.n1pad).or.(n2(i).gt.n2pad)) then
    write(*, FMT='(A)') "IMCOM ERROR: PSF input image is larger than galaxy input image: FAIL! Sorry, needs recoding!"
    stop
  end if
  n1min = (n1pad - n1(i)) / 2 + 1
  n2min = (n2pad - n2(i)) / 2 + 1
  n1max = (n1pad + n1(i)) / 2
  n2max = (n2pad + n2(i)) / 2
  allocate(Gtemp(n1(i), n2(i)), STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for temporary PSF image array"
    stop
  end if
  write(*, FMT='(A)') "IMCOM: Reading PSF from "//trim(psffile(i))
  call imcom_readfits(trim(psffile(i)), n1(i), n2(i), Gtemp)
  G(n1min:n1max, n2min:n2max, i) = Gtemp
!  print *, n1pad, n2pad
!  print *, shape(G(n1min:n1max, n2min:n2max, i))
!  print *, n1min, n1max, n2min, n2max
  deallocate(Gtemp, STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot deallocate memory for temporary PSF image array"
    stop
  end if

end do
call imcom_writefits('pad1.fits', n1pad, n2pad, G(:, :, 1))
deallocate(n1, n2, bitpix, STAT=dealstat)
if (dealstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot deallocate memory for PSF image property arrays"
  stop
end if
END SUBROUTINE imcom_alloc_psfs

!---

SUBROUTINE imcom_alloc_gamma(n1pad, n2pad)
! Explore & read in the PSF files, then allocate the G and Gt arrays to store
! them
implicit none
integer, intent(OUT) :: n1pad, n2pad
integer :: n1, n2
integer :: bitpix
integer :: alstat, dealstat
integer :: n1min, n2min, n1max, n2max
real(KIND=8), dimension(:, :), allocatable :: Gamtemp

call imcom_sizefits(trim(gamfile), n1, n2, bitpix)
if ((abs((gamxscale - psfxscale) / psfxscale).gt.1.e-7).or.&
    (abs((gamyscale - psfyscale) / psfyscale).gt.1.e-7)) then
  write(*, FMT='(A)') "IMCOM ERROR: Input PSFs and output Gamma PSF must have same pixel scale in this version, sorry!"
  write(*, FMT='(A, 2F14.5)') "IMCOM: PSFXSCALE, PSFYSCALE = ", psfxscale, psfyscale
  write(*, FMT='(A, 2F14.5)') "IMCOM: GAMXSCALE, GAMYSCALE = ", gamxscale, gamyscale
  print *,(abs(gamxscale-psfxscale)),(abs(gamyscale-psfyscale))
  stop
endif
n1pad = 2 * n1gim * anint(maxval(gimxscale / gamxscale))
n2pad = 2 * n2gim * anint(maxval(gimyscale / gamyscale))
if (n1pad.lt.n1) n1pad = n1
if (n2pad.lt.n2) n2pad = n2
allocate(Gamma(n1pad, n2pad), Gammat(n1pad, n2pad), Gamtemp(n1, n2), &
         STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for Gamma arrays"
  stop
end if
Gamma = 0.d0
Gammat = 0.d0
write(*, FMT='(A)') "IMCOM: Reading ouput PSF from "//trim(gamfile)
call imcom_readfits(trim(gamfile), n1, n2, Gamtemp)
n1min = (n1pad - n1) / 2 + 1
n2min = (n2pad - n2) / 2 + 1
n1max = (n1pad + n1) / 2
n2max = (n2pad + n2) / 2
Gamma(n1min:n1max, n2min:n2max) = Gamtemp
call imcom_writefits('gam1.fits', n1pad, n2pad, Gamma(:, :))
END SUBROUTINE imcom_alloc_gamma

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
integer :: group, firstpix, alstat

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
integer, parameter :: rwmode = 1    ! Write mode for this subroutine
integer, parameter :: maxdim = 2    ! Maximum dimension for images, 2D only
real(KIND=8), parameter :: nullvar = -66.D6 ! Assigned to undefined pixel values
integer :: unit, stat
integer, parameter :: naxis = maxdim
integer, parameter :: blocksize = 1
integer, parameter :: bitpix = -64  ! DOUBLE precision
integer, dimension(maxdim) :: naxes
logical :: anyf
integer :: group, firstpix, exists, alstat

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
call FTPPRD(unit, group, firstpix, naxis1 * naxis2, &
            reshape(image, (/ naxis1 * naxis2 /)), stat)
call FTCLOS(unit, stat)       ! Close the file
call FTFIOU(unit, stat)       ! Free the I/O unit
call FTRPRT('STDERR', stat)   ! Error report
END SUBROUTINE imcom_writefits

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
read(UNIT=10, FMT=*, IOSTAT=iostat) scratch, kappa
if (iostat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(filename)
  write(*, FMT='(A)') "IMCOM: Expecting: KAPPA     <kappa>"
  write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
  write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
  write(*, FMT='(A)') "IMCOM: Check <kappa> floating point"
  stop
end if
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

! Optional read-in of A, B, T, filenames (sets to default if EOF encountered)
read(UNIT=10, FMT='(A12, A)', IOSTAT=iostat) scratch, Afile
if (iostat.lt.0) then
  Afile = Adeffile
  Bfile = Bdeffile
  Tfile = Tdeffile
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
  Tfile = Tdeffile
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
read(UNIT=10, FMT='(A12, A)', IOSTAT=iostat) scratch, Tfile
if (iostat.lt.0) then
  Tfile = Tdeffile
  close(UNIT=10, IOSTAT=iostat)
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Unable to close "//trim(filename)
    stop
  end if
  return
else if (iostat.gt.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(filename)
  write(*, FMT='(A)') "IMCOM: Expecting: TFILE     <filename>"
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
real, dimension(:), allocatable :: psfxscale_tmp, psfyscale_tmp
real, dimension(:), allocatable :: maxfreq_tmp

allocate(gimfile(nexp), psffile(nexp), noise(nexp), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Unable to allocate space for image filename arrays"
  stop
end if
! Slightly different setup for user/non-user defined X and Y arrays
if (userxy.eq.1) then
  allocate(psfxscale_tmp(nexp), psfyscale_tmp(nexp), maxfreq_tmp(nexp), &
           gimxfile(nexp), gimyfile(nexp), STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Unable to allocate space for image X & Y coord filename arrays"
    stop
  endif
else
  allocate(psfxscale_tmp(nexp), psfyscale_tmp(nexp), maxfreq_tmp(nexp), &
           gimxscale(nexp), gimyscale(nexp), dither(nexp, 2), STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Unable to allocate space for image X & Y coord filename arrays"
    stop
  endif
end if

! Then read in the info
do i=1, nexp

  open(UNIT=20, FILE=trim(inconfig(i)), STATUS="OLD", FORM="FORMATTED", &
     ACTION="READ", IOSTAT=iostat)
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
  read(UNIT=20, FMT=*, IOSTAT=iostat) scratch, psfxscale_tmp(i)
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(inconfig(i))
    write(*, FMT='(A)') "IMCOM: Expecting: PSFXSCALE <psfxscale>"
    write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
    write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
    write(*, FMT='(A)') "IMCOM: Check <psfxscale> floating point"
    stop
  end if 
  read(UNIT=20, FMT=*, IOSTAT=iostat) scratch, psfyscale_tmp(i)
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(inconfig(i))
    write(*, FMT='(A)') "IMCOM: Expecting: PSFYSCALE <psfxscale>"
    write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
    write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
    write(*, FMT='(A)') "IMCOM: Check <psfyscale> floating point"
    stop
  end if
  read(UNIT=20, FMT=*, IOSTAT=iostat) scratch, maxfreq_tmp(i)
  if (iostat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(inconfig(i))
    write(*, FMT='(A)') "IMCOM: Expecting: MAXFREQ   <maxfreq>"
    write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
    write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
    write(*, FMT='(A)') "IMCOM: Check <maxfreq> floating point [units cycles per arcsec]"
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
! Then read in galaxy image exposure information: format depends on setting of USERXY = 0/1 in main config file
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
! Then check that PSFXSCALE, PSFYSCALE and MAXFREQ are all equally valued, and assign
if (all(psfxscale_tmp.eq.psfxscale_tmp(1))) then
  psfxscale = psfxscale_tmp(1)
else
  write(*, FMT='(A)') "IMCOM ERROR: PSFXSCALE different for different exposures, not supported in this version"
  stop
end if
if (all(psfyscale_tmp.eq.psfyscale_tmp(1))) then
  psfyscale = psfyscale_tmp(1)
else
  write(*, FMT='(A)') "IMCOM ERROR: PSFYSCALE different for different exposures, not supported in this version"
  stop
end if
if (all(maxfreq_tmp.eq.maxfreq_tmp(1))) then
  maxfreq = maxfreq_tmp(1)
else
  write(*, FMT='(A)') "IMCOM ERROR: PSFXSCALE different for different exposures, not supported in this version"
  stop
end if
deallocate(psfxscale_tmp, psfyscale_tmp, maxfreq_tmp, STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot deallocate space for temporary PSF property arrays"
  stop
end if
END SUBROUTINE imcom_read_inconfigs

!---

SUBROUTINE imcom_read_outconfig(userxy)
! Read the configuration files for the output images and Gamma PSFs
implicit none
integer, intent(IN) :: userxy
character(LEN=256) :: scratch
integer :: iostat, alstat, i

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
read(UNIT=30, FMT=*, IOSTAT=iostat) scratch, gamxscale
if (iostat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(outconfig)
  write(*, FMT='(A)') "IMCOM: Expecting: GAMXSCALE  <gamxscale>"
  write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
  write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
  write(*, FMT='(A)') "IMCOM: Check <gamxscale> floating point"
  stop
end if
read(UNIT=30, FMT=*, IOSTAT=iostat) scratch, gamyscale
if (iostat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Incorrect line format for "//trim(outconfig)
  write(*, FMT='(A)') "IMCOM: Expecting: GAMYSCALE  <gamyscale>"
  write(*, FMT='(A,I7)') "IMCOM: IOSTAT = ", iostat
  write(*, FMT='(A)') "IMCOM: Check EOF (if IOSTAT negative)"
  write(*, FMT='(A)') "IMCOM: Check <gamyscale> floating point"
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
end if
close(UNIT=30, IOSTAT=iostat)
if (iostat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Unable to close "//trim(outconfig)
  stop
end if
! From this information define total size of output array
m = n1out * n2out
END SUBROUTINE imcom_read_outconfig

!---

END MODULE imcom_io
