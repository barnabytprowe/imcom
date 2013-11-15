!    IMCOM_DATA.F90 - Fortran 95 module that acts as a shared-data / COMMON 
!                     block for the prototype IMCOM package
!    Copyright (C) 2011-2013  Barnaby Rowe
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
MODULE imcom_data

implicit none
save

! Command line variables:
logical :: USB
character(LEN=256) :: config
real(KIND=8) :: US_max, US_tol, kappa_min, kappa_max, saturation
integer :: maxNbis, forceT, forceSys

! Config file variables:
! i) PSF pixel scales and number of exposures we are combining
integer :: nexp
real(KIND=8) :: psfxscale, psfyscale

! ii) galaxy image & PSF image, plus x,y maps: config file and filenames
character(LEN=256), dimension(:), allocatable :: inconfig
character(LEN=256), dimension(:), allocatable :: gimfile, psffile
character(LEN=256), dimension(:), allocatable :: gimxfile, gimyfile
real(KIND=8), dimension(:), allocatable :: gimxscale, gimyscale
real(KIND=8), dimension(:), allocatable :: rotangdeg

! iii) dither position (coordinates of lower left corner pixel centre in each exposure) & coordinates of lower-left pixel in output image
real(KIND=8), dimension(:), allocatable :: noise
real(KIND=8), dimension(:, :), allocatable :: dither
real(KIND=8) :: outpos1, outpos2

! iv) output image "H" and Gamma info/config file and filenames
character(LEN=256) :: outconfig, Hfile, gamfile
character(LEN=256) :: outxfile, outyfile
real(KIND=8) :: outxscale, outyscale ! input/desired image pixel units

! Switch for whether all PSFs are actually the same file
integer :: psfconst, nexppsf

! Sizes of PSF, galaxy input and output images, including number of masked image pixels
integer :: n1psf, n2psf
integer :: n1gim, n2gim
integer :: n_unmasked
integer :: n1out, n2out

integer :: n, m  ! Dimensions specifiers for A, B, K, T, S, U matrices

! Default A, B and T storage filenames
character(LEN=256), parameter :: Adeffile = "A.fits"
character(LEN=256), parameter :: Bdeffile = "B.fits"
character(LEN=256), parameter :: Qdeffile = "Q.fits"
character(LEN=256), parameter :: Ldeffile = "L.fits"
character(LEN=256), parameter :: Pdeffile = "P.fits"

! A, B and T, U, S storage filenames
character(LEN=256) :: Afile
character(LEN=256) :: Bfile
character(LEN=256) :: Qfile
character(LEN=256) :: Lfile
character(LEN=256) :: Pfile
character(LEN=256) :: Kfile
character(LEN=256) :: Tfile
character(LEN=256) :: Sfile
character(LEN=256) :: Ufile


! Real space inputs and outputs
real(KIND=8), dimension(:, :, :), allocatable :: Im_unmasked
integer, dimension(:, :, :), allocatable :: exp_unmasked
real(KIND=8), dimension(:, :), allocatable :: H
real(KIND=8), dimension(:, :), allocatable :: U
real(KIND=8), dimension(:, :), allocatable :: S
real(KIND=8), dimension(:, :), allocatable :: K  ! kappa now varies across output image
real(KIND=8), dimension(:, :, :), allocatable :: x_unmasked
real(KIND=8), dimension(:, :, :), allocatable :: y_unmasked
real(KIND=8), dimension(:, :), allocatable :: Xa
real(KIND=8), dimension(:, :), allocatable :: Ya

! Similar things but arranged into big long 1D vectors...
real(KIND=8), dimension(:), allocatable :: I_i, I_unmasked_i
integer, dimension(:), allocatable :: mask_i, mask_subs_i, exp_unmasked_i, exp_i
real(KIND=8), dimension(:), allocatable :: x_i, x_unmasked_i
real(KIND=8), dimension(:), allocatable :: y_i, y_unmasked_i
real(KIND=8), dimension(:), allocatable :: H_a
real(KIND=8), dimension(:), allocatable :: X_a
real(KIND=8), dimension(:), allocatable :: Y_a
real(KIND=8), dimension(:), allocatable :: U_a
real(KIND=8), dimension(:), allocatable :: S_a
real(KIND=8), dimension(:), allocatable :: K_a

! Diagonal noise matrix (only diagonal noise implemented)
real(KIND=8), dimension(:), allocatable :: Ndiag

! (PSFs: one per exposure)
real(KIND=8), dimension(:, :, :), allocatable :: G_unrot
real(KIND=8), dimension(:, :, :), allocatable :: G_rot

! (Desired PSF)
real(KIND=8), dimension(:, :), allocatable :: Gamma

! Fourier space counterparts
complex(KIND=8), dimension(:, :, :), allocatable :: Gt_unrot
complex(KIND=8), dimension(:, :, :), allocatable :: Gt_rot
real(KIND=8), dimension(:, :), allocatable :: ux
real(KIND=8), dimension(:, :), allocatable :: uy
complex(KIND=8), dimension(:, :), allocatable :: Gammat
  
! Problem/work matrices
real(KIND=8), dimension(:, :), allocatable :: A_aij
real(KIND=8), dimension(:, :, :), allocatable :: Alookup ! now packed up-triang.
real(KIND=8), dimension(:, :), allocatable :: Q_ij
real(KIND=8), dimension(:), allocatable :: L_i
real(KIND=8), dimension(:, :), allocatable :: B_ia
real(KIND=8), dimension(:, :, :), allocatable :: Blookup
real(KIND=8) :: C_a
real(KIND=8), dimension(:, :), allocatable :: T_ia
real(KIND=8), dimension(:, :), allocatable :: P_ia, P2_ia

! A matrix condition number
real(KIND=8) :: condition

! Constants
real(KIND=8), parameter :: pi = 3.141592653589793d0

! FFTW3 planning constants required by DFT calls (directly lifted from fftw3.f in directory
! containing fftw3.h)
!
      INTEGER FFTW_R2HC
      PARAMETER (FFTW_R2HC=0)
      INTEGER FFTW_HC2R
      PARAMETER (FFTW_HC2R=1)
      INTEGER FFTW_DHT
      PARAMETER (FFTW_DHT=2)
      INTEGER FFTW_REDFT00
      PARAMETER (FFTW_REDFT00=3)
      INTEGER FFTW_REDFT01
      PARAMETER (FFTW_REDFT01=4)
      INTEGER FFTW_REDFT10
      PARAMETER (FFTW_REDFT10=5)
      INTEGER FFTW_REDFT11
      PARAMETER (FFTW_REDFT11=6)
      INTEGER FFTW_RODFT00
      PARAMETER (FFTW_RODFT00=7)
      INTEGER FFTW_RODFT01
      PARAMETER (FFTW_RODFT01=8)
      INTEGER FFTW_RODFT10
      PARAMETER (FFTW_RODFT10=9)
      INTEGER FFTW_RODFT11
      PARAMETER (FFTW_RODFT11=10)
      INTEGER FFTW_FORWARD
      PARAMETER (FFTW_FORWARD=-1)
      INTEGER FFTW_BACKWARD
      PARAMETER (FFTW_BACKWARD=+1)
      INTEGER FFTW_MEASURE
      PARAMETER (FFTW_MEASURE=0)
      INTEGER FFTW_DESTROY_INPUT
      PARAMETER (FFTW_DESTROY_INPUT=1)
      INTEGER FFTW_UNALIGNED
      PARAMETER (FFTW_UNALIGNED=2)
      INTEGER FFTW_CONSERVE_MEMORY
      PARAMETER (FFTW_CONSERVE_MEMORY=4)
      INTEGER FFTW_EXHAUSTIVE
      PARAMETER (FFTW_EXHAUSTIVE=8)
      INTEGER FFTW_PRESERVE_INPUT
      PARAMETER (FFTW_PRESERVE_INPUT=16)
      INTEGER FFTW_PATIENT
      PARAMETER (FFTW_PATIENT=32)
      INTEGER FFTW_ESTIMATE
      PARAMETER (FFTW_ESTIMATE=64)
      INTEGER FFTW_ESTIMATE_PATIENT
      PARAMETER (FFTW_ESTIMATE_PATIENT=128)
      INTEGER FFTW_BELIEVE_PCOST
      PARAMETER (FFTW_BELIEVE_PCOST=256)
      INTEGER FFTW_NO_DFT_R2HC
      PARAMETER (FFTW_NO_DFT_R2HC=512)
      INTEGER FFTW_NO_NONTHREADED
      PARAMETER (FFTW_NO_NONTHREADED=1024)
      INTEGER FFTW_NO_BUFFERING
      PARAMETER (FFTW_NO_BUFFERING=2048)
      INTEGER FFTW_NO_INDIRECT_OP
      PARAMETER (FFTW_NO_INDIRECT_OP=4096)
      INTEGER FFTW_ALLOW_LARGE_GENERIC
      PARAMETER (FFTW_ALLOW_LARGE_GENERIC=8192)
      INTEGER FFTW_NO_RANK_SPLITS
      PARAMETER (FFTW_NO_RANK_SPLITS=16384)
      INTEGER FFTW_NO_VRANK_SPLITS
      PARAMETER (FFTW_NO_VRANK_SPLITS=32768)
      INTEGER FFTW_NO_VRECURSE
      PARAMETER (FFTW_NO_VRECURSE=65536)
      INTEGER FFTW_NO_SIMD
      PARAMETER (FFTW_NO_SIMD=131072)
      INTEGER FFTW_NO_SLOW
      PARAMETER (FFTW_NO_SLOW=262144)
      INTEGER FFTW_NO_FIXED_RADIX_LARGE_N
      PARAMETER (FFTW_NO_FIXED_RADIX_LARGE_N=524288)
      INTEGER FFTW_ALLOW_PRUNING
      PARAMETER (FFTW_ALLOW_PRUNING=1048576)
      INTEGER FFTW_WISDOM_ONLY
      PARAMETER (FFTW_WISDOM_ONLY=2097152)

END MODULE imcom_data
