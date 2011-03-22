!    IMCOM_BISECT.F90 - Fortran 95 module that handles interval bisection  
!                       root-finding for the prototype IMCOM package
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
!
!    CONTAINS
!    Routines for IMCOM that find the pixel-by-pixel solution for
!    kappa that minimizes the noise (S) or leakage objective (U).
!
MODULE imcom_bisect

USE imcom_data
USE imcom_io
USE imcom_proc
USE imcom_matrix

implicit none

CONTAINS

!---

SUBROUTINE imcom_get_T(forcebuild)
implicit none
integer, intent(IN) :: forcebuild
logical :: Texists
integer :: n1_tmp, n2_tmp, bit_tmp
integer :: alstat

inquire(FILE=trim(Tfile), EXIST=Texists)
if (Texists.and.(forcebuild.eq.0)) then
  call imcom_sizefits(trim(Tfile), n1_tmp, n2_tmp, bit_tmp)
  if ((n1_tmp.ne.n).or.(n2_tmp.ne.m)) then
    write(*, FMT='(A)') "IMCOM ERROR: Dimensions of Matrix "//trim(Tfile)//" do not match <config_file> image dimensions"
    write(*, FMT='(A,2I9)') "IMCOM ERROR: Shoud be rectangle of sides: ",n,m
    stop
  endif
  allocate(T_ia(n1_tmp, n2_tmp), STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for T matrix -- too large? Sparse Matrix formalism not yet coded."
    stop
  endif
  T_ia = 0.d0
  write(*, FMT='(A)') "IMCOM: Reading T matrix from "//trim(Tfile)
  call imcom_readfits(trim(Tfile), n1_tmp, n2_tmp, T_ia)
!  call imcom_build_U  don't bother with these... they exist, if T exists.
!  call imcom_build_S
else
  call imcom_bisect_kappa(USB, US_max, kappa_min, kappa_max, US_tol, maxNbis)
  call imcom_build_T
  if (USB.eqv.(.TRUE.)) then
    call imcom_build_S
  else if (USB.eqv.(.FALSE.)) then
    call imcom_build_U
  end if
end if
END SUBROUTINE imcom_get_T

!---

SUBROUTINE imcom_bisect_kappa(UnotS, US, k_min, k_max, tol, maxnits)
!
! Solves AT = -2B
implicit none
logical, intent(IN) :: UnotS
real(KIND=8), intent(IN) :: US, k_min, k_max, tol
integer, intent(IN) :: maxnits
integer :: alstat, a

write(*, FMT='(A)') "IMCOM: Solving for kappa using interval bisection"
if (UnotS.eqv.(.TRUE.)) then
  allocate(K_a(m), U_a(m), STAT=alstat)
  K_a = 0.d0
  U_a = 0.d0
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for kappa and U storage arrays"
    stop
  end if
  !$omp parallel do schedule(dynamic, 64)
  do a=1, m

    call imcom_bisect_by_U(n, a, k_min, k_max, US, tol, maxnits, K_a(a), &
                           U_a(a))
    if (a.eq.int(0.2 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 20% complete"
    if (a.eq.int(0.4 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 40% complete"
    if (a.eq.int(0.6 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 60% complete"
    if (a.eq.int(0.8 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 80% complete"

  end do
  !$omp end parallel do
  write(*, FMT='(A)') "IMCOM: 100% complete"
  allocate(U(n1out, n2out), STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for U output array"
    stop
  end if
  U = reshape(U_a, (/ n1out, n2out /))
  deallocate(U_a, STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot deallocate memory for U output array"
    stop
  end if
  write(*, FMT='(A)') "IMCOM: Saving U leakage map to "//trim(Ufile)
  call imcom_writefits(trim(Ufile), n1out, n2out, real(U, 8))
else if (UnotS.eqv.(.FALSE.)) then
  allocate(K_a(m), S_a(m), STAT=alstat)
  K_a = 0.d0
  S_a = 0.d0
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for T storage arrays"
    stop
  end if
  !$omp parallel do schedule(dynamic, 64)
  do a=1, m

    call imcom_bisect_by_S(n, a, k_min, k_max, US, tol, maxnits, K_a(a), &
                           S_a(a))
    if (a.eq.int(0.2 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 20% complete"
    if (a.eq.int(0.4 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 40% complete"
    if (a.eq.int(0.6 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 60% complete"
    if (a.eq.int(0.8 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 80% complete"

  end do
  !$omp end parallel do
  write(*, FMT='(A)') "IMCOM: 100% complete"
  allocate(S(n1out, n2out), STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for S output array"
    stop
  end if
  S = reshape(S_a, (/ n1out, n2out /))
  deallocate(S_a, STAT=alstat)
  if (alstat.ne.0) then
    write(*, FMT='(A)') "IMCOM ERROR: Cannot deallocate memory for S output array"
    stop
  end if
  write(*, FMT='(A)') "IMCOM: Saving S variance map to "//trim(Sfile)
  call imcom_writefits(trim(Sfile), n1out, n2out, real(S, 8))
end if
allocate(K(n1out, n2out), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for kappa output array"
  stop
end if
K = reshape(K_a, (/ n1out, n2out /))
! Note do not deallocate K_a(:), it's small and highly useful later...
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot deallocate memory for kappa storage array"
  stop
end if
write(*, FMT='(A)') "IMCOM: Saving K solution map to "//trim(Kfile)
call imcom_writefits(trim(Kfile), n1out, n2out, real(K, 8))
END SUBROUTINE imcom_bisect_kappa

!---

SUBROUTINE imcom_bisect_by_U(n_dum, a_dum, k_min, k_max, U_max, tol, maxnits, &
                             k_root, U_out)
! This subroutine solves (via interval bisection) the T_i at a given output 
! pixel a, keeping U < max(U_max, U(kappa = k_min))
!
implicit none
integer, intent(IN) :: n_dum, a_dum, maxnits  ! a must be in range 1,...,m
real(KIND=8), intent(IN) :: k_min, k_max, U_max, tol
real(KIND=8), intent(OUT) :: k_root, U_out
integer :: i
real(KIND=8) :: hi, lo, mid, root, dx, fmid, f, eps  ! log(k)_hi, etc.
real(KIND=8), external :: DLAMCH

eps = DLAMCH('Espilon')
U_out = imcom_U_from_kappa(a_dum, n_dum, k_min)
f = U_out - U_max
if (f.lt.0.d0) then
  lo = dlog10(k_min)
  hi = dlog10(k_max)
else
  k_root = k_min
  return
end if
root = lo
dx = (hi - lo)
do i=1, maxnits

  dx = 0.5d0 * dx
  mid = root + dx
  fmid = imcom_U_from_kappa(a_dum, n_dum, (10.d0)**mid) - U_max
  if (fmid.lt.0.d0) then
    root = mid
    U_out = fmid + U_max
    if (U_out.lt.0.d0) then 
      write(*, FMT='(A, I7)') "IMCOM WARNING: Negative U at alpha = ", a_dum
      write(*, FMT='(A, E14.6, A, E14.6)') "IMCOM WARNING: log(k), U = ", &
                                           root, ", ", U_out
    end if
     if ((abs(fmid).lt.tol)) exit
  end if
  if (abs(fmid).lt.(C_a * eps)) exit

end do
if (U_out.lt.0.d0) then
  write(*, FMT='(A, I4)') "IMCOM WARNING: Negative U for a = ", a_dum
end if
k_root = 10.d0**root
return
END SUBROUTINE imcom_bisect_by_U

!---

SUBROUTINE imcom_bisect_by_S(n_dum, a_dum, k_min, k_max, S_max, tol, maxnits, &
                             k_root, S_out)
! This subroutine solves (via interval bisection) the T_i at a given output 
! pixel a, keeping S < max(U_max, S(kappa = k_max))
!
implicit none
integer, intent(IN) :: n_dum, a_dum, maxnits  ! a must be in range 1,...,m
real(KIND=8), intent(IN) :: k_min, k_max, S_max, tol
real(KIND=8), intent(OUT) :: k_root, S_out
integer :: i
real(KIND=8) :: hi, lo, mid, root, dx, fmid, f, eps  ! log(k)_hi, etc.
real(KIND=8), external :: DLAMCH

eps = DLAMCH('Espilon')
S_out = imcom_S_from_kappa(a_dum, n_dum, k_max)
f = S_out - S_max
if (f.lt.0.d0) then
  lo = dlog10(k_max)
  hi = dlog10(k_min)
else
  k_root = k_max
  return
end if
root = lo
dx = (hi - lo)
do i=1, maxnits

  dx = 0.5d0 * dx
  mid = root + dx
  fmid = imcom_S_from_kappa(a_dum, n_dum, (10.d0)**mid) - S_max
  if (fmid.lt.0.d0) then
     root = mid
     S_out = fmid + S_max
     if (S_out.lt.0.d0) print *, a_dum, root, S_out, fmid, tol, S_max
     if (abs(fmid).lt.tol) exit
  end if
  if (abs(fmid).lt.(C_a * eps)) exit

end do
k_root = 10.d0**root
return
END SUBROUTINE imcom_bisect_by_S

!---

END MODULE imcom_bisect

