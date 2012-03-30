!    IMCOM.F90 - Fortran 95 prototype package for the implementation of
!                the Optimal Linear IMage COMbination algorithm described
!                by Rowe, Hirata & Rhodes (2011)
!    Version 0.2 (alpha)
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
!                   IMCOM_IO.F90
!                   IMCOM_MATRIX.F90
!                   IMCOM_BISECT.F90
!
!    DESCRIPTION
!    Please see the file IMCOM_DOCUMENTATION.PDF for details on the
!    installation, configuration and running of the IMCOM package.
!    For a full description of this implementation and its uses, please
!    see Rowe, Hirata & Rhodes (2011).
!
PROGRAM imcom

USE imcom_data
USE imcom_io
USE imcom_proc
USE imcom_matrix
USE imcom_bisect

implicit none
integer :: userxy, i, npoly, npad, alstat
integer(KIND=8) :: ftplan
real(KIND=8) :: eps
real(KIND=8), external :: DLAMCH
logical :: Aexists, Bexists, Texists
integer :: nexppsf
integer :: saveA, saveB, saveQL, saveP

npoly = 7     ! Make these user specifiable?
npad = 3      !
maxNbis = 53       ! This since 2^-53 is the machine precision at real(KIND=8)
! Also, currently kappa_min and kappa_max are now set via C_a * machine_epsilon
! and C_a / machine_epsilon [accessed with DLAMCH('Eps')], rather than
! being user-input... is this a good idea?

call imcom_get_cmdline
call imcom_welcome

! Read the config files
call imcom_read_config(config, userxy)

! Read and allocate the input images
call imcom_alloc_gims

if (userxy.eq.1) then
! Read the x-y arrays
  call imcom_alloc_xy
  call imcom_alloc_XaYa
else
! Make the x-y arrays
  call imcom_build_xy
  call imcom_build_XaYa
end if

psfconst = imcom_test_psfconst() ! Set global switch psfconst
if (psfconst.eq.1) then ! and related nexppsf
  nexppsf = 1
else
  nexppsf = nexp
end if

! Read and allocate the PSFs and Gamma PSF
call imcom_alloc_psfs(n1psf, n2psf)
! Make the ux, uy FT arrays
call imcom_build_uxuy

inquire(FILE=trim(Afile), EXIST=Aexists)
inquire(FILE=trim(Bfile), EXIST=Bexists)
inquire(FILE=trim(Tfile), EXIST=Texists)
if ((((Bexists.eqv..FALSE.).or.(Aexists.eqv..FALSE.)) &
    .and.(Texists.eqv..FALSE.)).or.(forceSys.eq.1)) then
! Rotate PSFs and calculate the Fourier transform Gt(ux, uy) of each
    write(*, FMT='(A)') "IMCOM: Rotating and Fourier transforming PSF images"
    eps = DLAMCH('Epsilon')
    call imcom_plan_ft_r2c(n1psf, n2psf, 0, ftplan)
    do i=1, nexppsf

      if (abs(rotangdeg(i)).gt.(2.d0 * pi * eps)) then
        call imcom_rotate_image(G_unrot(:, :, i), rotangdeg(i), npoly, 1, &
                                G_rot(:, :, i))
      else
        G_rot(:, :, i) = G_unrot(:, :, i)
      end if
      call imcom_ft_r2c(n1psf, n2psf, ftplan, G_rot(:, :, i), Gt_rot(:, :, i))

    end do
    call imcom_destroy_plan(ftplan)
endif
! Calculate the Fourier transform Gammat(ux, uy) of Gamma(y, y)
call imcom_plan_ft_r2c(n1psf, n2psf, 0, ftplan)
call imcom_ft_r2c(n1psf, n2psf, ftplan, Gamma(:, :), Gammat(:, :))
call imcom_destroy_plan(ftplan)

! Calculate the N (diagonal only, currently) noise array
call imcom_build_Ndiag

! Calulate C (really quick) and build the (currently diagonal only) noise cov.
call imcom_calc_C

! Set kappa min / max from machine precision limits
eps = DLAMCH('Epsilon')
kappa_min = C_a * eps
kappa_max = C_a / eps

! New feature: if the specified system matrix filenames are 'None' then don't bother saving them...
if (Afile.eq."None") then
  saveA = 0
else
  saveA = 1
end if
if (Bfile.eq."None") then
  saveB = 0
else
  saveB = 1
end if
if ((Qfile.eq."None").and.(Lfile.eq."None")) then
  saveQL = 0
else
  saveQL = 1
end if
if (Pfile.eq."None") then
  saveP = 0
else
  saveP = 1
end if
! If they can't be loaded in, or if T needs to be rebuilt, build the systems matrices
if (((forceSys.ne.0).or.(forceT.ne.0)).or.(Texists.eqv..False.)) then
! Read / build the A array
  call imcom_get_A(npoly, npad, saveA, forceSys)  ! arg#3 = save, arg#4 = force build
  call imcom_get_QL(saveQL, forceSys)
! Read / build the B array
  call imcom_get_B(npoly, npad, saveB, forceSys)
! Construct the projection matrices
  call imcom_get_P(saveP, forceSys)
  deallocate(A_aij, B_ia, STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot deallocate A & B matrices"
    stop
  end if
end if
! Read / solve the T array
call imcom_get_T(forceT)
! Build the output products using the galaxy images, A, B, C and our calculated T 
! transformation
call imcom_build_image

END PROGRAM imcom
