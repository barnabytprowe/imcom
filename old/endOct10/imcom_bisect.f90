!
! Module for the routines associated with finding the pixel-by-pixel solution
! for kappa that minimizes the noise (S) or leakage objective (U)
!
!-----------------------------------------------------------

MODULE imcom_bisect

USE imcom_data
USE imcom_io
USE imcom_proc
USE imcom_matrix
USE omp_lib

implicit none

CONTAINS

!---

SUBROUTINE imcom_get_T(saveTKU, forcebuild)
implicit none
integer, intent(IN) :: saveTKU, forcebuild
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
else
  call imcom_eigen_A
  call imcom_bisect_T(USB, US_min, kappa_min, kappa_max, US_tol, maxNbis, saveTKU)
end if
END SUBROUTINE imcom_get_T

!---

SUBROUTINE imcom_bisect_T(UnotS, US_min, k_min, k_max, tol, maxnits, saveTKU)
!
! Solves AT = -2B
implicit none
logical, intent(IN) :: UnotS
real(KIND=8), intent(IN) :: US_min, k_min, k_max, tol
integer, intent(IN) :: maxnits, saveTKU
integer :: alstat, a
real(KIND=8), dimension(:), allocatable :: Tidum

allocate(T_ia(n, m), Tidum(n), K_a(m), U_a(m), STAT=alstat)
T_ia = 0.d0
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for T storage array"
  stop
end if
write(*, FMT='(A)') "IMCOM: Solving for T matrix using interval bisection"
if (UnotS.eqv.(.TRUE.)) then
  !$omp parallel do schedule(dynamic, 64)
  do a=1, m

    call imcom_bisect_TibyU(n, a, k_min, k_max, US_min, tol, maxnits, K_a(a), &
                            T_ia(:, a), U_a(a))
    if (a.eq.int(0.1 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 10% complete"
    if (a.eq.int(0.2 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 20% complete"
    if (a.eq.int(0.3 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 30% complete"
    if (a.eq.int(0.4 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 40% complete"
    if (a.eq.int(0.5 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 50% complete"
    if (a.eq.int(0.6 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 60% complete"
    if (a.eq.int(0.7 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 70% complete"
    if (a.eq.int(0.8 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 80% complete"
    if (a.eq.int(0.9 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 90% complete"
    print *, a, K_a(a), U_a(a)

  end do
  !$omp end parallel do
else if (UnotS.eqv.(.FALSE.)) then
  !$oomp parallel do
  do a=1, m

    call imcom_bisect_TibyS(n, a, k_min, k_max, US_min, tol, maxnits, K_a(a), &
                            T_ia(:, a))
    if (a.eq.int(0.2 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 20% complete"
    if (a.eq.int(0.4 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 40% complete"
    if (a.eq.int(0.6 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 60% complete"
    if (a.eq.int(0.8 * real(n, 4))) write(*, FMT='(A)') "IMCOM: 80% complete"

  end do
  !$oomp end parallel do
end if
write(*, FMT='(A)') "IMCOM: T matrix complete"
allocate(K(n1out, n2out), U(n1out, n2out), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for kappa and U output array"
  stop
end if
K = reshape(K_a, (/ n1out, n2out /))
U = reshape(U_a, (/ n1out, n2out /))
deallocate(K_a, U_a, STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot deallocate memory for kappa storage array"
  stop
end if
if (saveTKU.eq.1) then 
  write(*, FMT='(A)') "IMCOM: Saving K solution map to "//trim(Kfile)
  call imcom_writefits(trim(Kfile), n1out, n2out, real(K, 8))
  write(*, FMT='(A)') "IMCOM: Saving T solution matrix to "//trim(Tfile)
  call imcom_writefits(trim(Tfile), n, m, real(T_ia, 8))
  write(*, FMT='(A)') "IMCOM: Saving U leakage map to "//trim(Ufile)
  call imcom_writefits(trim(Ufile), n1out, n2out, real(U, 8))
end if
END SUBROUTINE imcom_bisect_T

!---

SUBROUTINE imcom_bisect_TibyU(n_dum, a_dum, k_min, k_max, U_min, tol, maxnits, &
                              k_root, Ti_out, U_out)
! This subroutine solves (via interval bisection) the T_i at a given output 
! pixel a, keeping U < max(U_min, U(kappa = k_min))
!
implicit none
integer, intent(IN) :: n_dum, a_dum, maxnits  ! a must be in range 1,...,m
real(KIND=8), intent(IN) :: k_min, k_max, U_min, tol
real(KIND=8), intent(OUT) :: k_root, U_out
real(KIND=8), dimension(n_dum), intent(OUT) :: Ti_out
real(KIND=8), dimension(n_dum) :: Bi
integer :: i
real(KIND=8) :: hi, lo, mid, root, dx, fmid, f, eps  ! log(k)_hi, etc.
real(KIND=8), external :: DLAMCH

eps = DLAMCH('Espilon')
Bi = B_ia(:, a_dum)
Ti_out = 0.d0
lo = dlog10(k_min)
U_out = imcom_bisect_func(.TRUE., lo, n_dum, a_dum, Bi, Ti_out)
f = U_out - U_min
if (f.lt.0.d0) then
  lo = lo
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
  fmid = imcom_bisect_func(.TRUE., mid, n_dum, a_dum, Bi, Ti_out) - U_min
  if (fmid.le.0.d0) then
     root = mid
     U_out = fmid + U_min
     if (abs(fmid).lt.tol) exit
  end if
  if (abs(fmid).lt.(C_a * eps)) exit
  print *, mid, fmid

end do
k_root = 10.d0**root
print *, mid, fmid
return
END SUBROUTINE imcom_bisect_TibyU

!---

SUBROUTINE imcom_bisect_TibyS(n_dum, a_dum, k_min, k_max, S_min, tol, maxnits, &
                              k_root, Ti_out)
! This subroutine solves (via interval bisection) the T_i at a given output 
! pixel a, keeping S < max(S_min, U(kappa = k_max))
!
implicit none
integer, intent(IN) :: n_dum, a_dum, maxnits  ! a must be in range 1,...,m
real(KIND=8), intent(IN) :: k_min, k_max, S_min, tol
real(KIND=8), intent(OUT) :: k_root
real(KIND=8), dimension(n_dum), intent(OUT) :: Ti_out
real(KIND=8), dimension(n_dum) :: Bi
integer :: i
real(KIND=8) :: hi, lo, mid, root, dx  ! log(k)_hi, etc.
real(KIND=8) :: f, fmid, eps
real(KIND=8), external :: DLAMCH

eps = DLAMCH('Espilon')
Bi = B_ia(:, a_dum)
Ti_out = 0.d0
hi = dlog10(k_max)
f = imcom_bisect_func(.FALSE., hi, n_dum, a_dum, Bi, Ti_out) - S_min
if (f.lt.0.d0) then
  lo = hi
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
  fmid = imcom_bisect_func(.FALSE., mid, n_dum, a_dum, Bi, Ti_out) - S_min
  if (fmid.le.0.d0) root = mid
  if ((abs(dx).lt.tol).or.(abs(fmid).lt.(C_a * eps))) exit
  print *, mid, fmid

end do
k_root = 10.d0**root
return
END SUBROUTINE imcom_bisect_TibyS

!---

FUNCTION imcom_bisect_func(doU, lk_dum, n_dum, a_dum, Bi_in, Ti_out)
! Uses the input value of log(kappa) to solve A_aij T_ia = B_ia for a single 
! vector specified by a = a_dum
!
! This T_i is then used to calculate U or S (as specified by UnotS_dum) for 
! this a, and the return value of the function is then U - U_min or S - S_min
! The T_i is also output (it is useful!)
implicit none
logical, intent(IN) :: doU
real(KIND=8), intent(IN) :: lk_dum
integer, intent(IN) :: n_dum, a_dum
! **** Uses saved, common A0_aij  ****
real(KIND=8), dimension(n_dum), intent(IN) :: Bi_in
real(KIND=8), dimension(n_dum), intent(OUT) :: Ti_out
real(KIND=8) :: imcom_bisect_func
real(KIND=8) :: ktrial

! Sanity checks
if (n_dum.ne.n) then
  write(*, FMT='(A)') "IMCOM ERROR: N argument to IMCOM_BISECT_FUNC does not match N dimension of A matrix"
  stop
endif
if ((a_dum.lt.0).or.(a_dum.gt.m)) then
  write(*, FMT='(A)') "IMCOM ERROR: A argument to IMCOM_BISECT_FUNC does fit within M bounds of output matrix"
  stop
endif
! Initialize kappa
ktrial = (10.d0)**lk_dum
Ti_out = imcom_Ti_from_kappa(a_dum, n_dum, ktrial)
if (doU.eqv.(.TRUE.)) then
  imcom_bisect_func = imcom_Ua(n_dum, Bi_in, Ti_out)
else if (doU.eqv.(.FALSE.)) then
  imcom_bisect_func = imcom_Sa(n_dum, Ti_out)
end if
END FUNCTION imcom_bisect_func

!---

FUNCTION imcom_Ua(n_dum, Bi_dum, Tj_dum)
! Also uses unpacked (for speed), pre calculated A0_aij matrix and C_a scalar
implicit none
integer, intent(IN) :: n_dum
real(KIND=8), dimension(n_dum), intent(IN) :: Bi_dum, Tj_dum
real(KIND=8) :: imcom_Ua  ! output
real(KIND=8), dimension(n_dum) :: ATi  ! internal work vector
character, parameter :: uplo = 'U'
real(KIND=8), parameter :: alpha = 1.d0, beta= 0.d0
real(KIND=8), external :: DDOT

ATi = 0.d0
! Use the BLAS sym. matrix * vector routines and dot product 
call DSYMV(uplo, n_dum, alpha, A0_aij, n_dum, Tj_dum, 1, beta, ATi, 1)
imcom_Ua = (C_a + DDOT(n_dum, ATi + Bi_dum, 1, Tj_dum, 1)) / C_a
END FUNCTION imcom_Ua

!---

FUNCTION imcom_Sa(n_dum, Ti_dum)
! Also uses the diagonals of the noise covariance NN(i=1,n)
implicit none
integer, intent(IN) :: n_dum
real(KIND=8), dimension(n_dum), intent(IN) :: Ti_dum
real(KIND=8) :: imcom_Sa
real(KIND=8), external :: DDOT

imcom_Sa = DDOT(n_dum, Ndiag(1: n_dum) * Ti_dum, 1, Ti_dum, 1)
END FUNCTION imcom_Sa

!---

FUNCTION imcom_A_plus_kNdiag(n_dum, A_in, kNdiag_in)
! Adds (kappa * N_ii) to the diagonals A_ii of an input A matrix
implicit none
integer, intent(IN) :: n_dum
real(KIND=8), dimension(n_dum, n_dum), intent(IN) :: A_in
real(KIND=8), dimension(n_dum), intent(IN) :: kNdiag_in
real(KIND=8), dimension(n_dum, n_dum) :: imcom_A_plus_kNdiag  ! output
integer :: i

! Begin OPEN MP stuff
!$omp parallel
!$omp workshare
imcom_A_plus_kNdiag = A_in
!$omp end workshare
!$omp workshare
forall(i=1: n_dum) imcom_A_plus_kNdiag(i, i) = imcom_A_plus_kNdiag(i, i) &
                                             + kNdiag_in(i)
!$omp end workshare
!$omp end parallel
END FUNCTION imcom_A_plus_kNdiag

!---

SUBROUTINE imcom_eigen_A
! Performs eigen decomposition of A matrix using LAPACK v3.X's relatively 
! robust algorithm DSYEVR (other algorithms are available)
! 
! Not threaded with OpenMP, but some LAPACK implementations allow multi-thread
! support for vector operations
!
implicit none
character(LEN=1), parameter :: jobz = "V", range = "A", uplo = "U"
real(KIND=8), dimension(:, :), allocatable :: A0tmp
real(KIND=8), dimension(:), allocatable :: work
integer, dimension(:), allocatable :: isupport, iwork
real(KIND=8), dimension(1) :: worktest
integer, dimension(1) :: iworktest
real(KIND=8) :: abstol, vl, vu
real(KIND=8), external :: DLAMCH
integer :: alstat, il, iu, n_out, info, lwork, liwork

allocate(R_A(n, n), A0tmp(n, n), L0_A(n), isupport(2 * n), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for Eigendecomposition outputs"
  stop
end if
A0tmp = A0_aij
abstol = DLAMCH('Safe minimum')  ! Suggested in LAPACK notes for DSYEVR
! First query to get range of work matrices
call DSYEVR(jobz, range, uplo, n, A0tmp, n, vl, vu, il, iu, abstol, n_out, &
            L0_A, R_A, n, isupport, worktest, -1, iworktest, -1, info)
lwork = nint(worktest(1))
liwork = iworktest(1)
allocate(work(lwork), iwork(liwork), STAT=alstat)
if (alstat.ne.0) then
  write(*, FMT='(A)') "IMCOM ERROR: Cannot allocate memory for Eigendecomposition work arrays"
  stop
end if
! Re-initialize in case DSYEVR destroyed previous inputs
vl = 0.d0
vu = 0.d0
il = 0
iu = 0
isupport = 0
work = 0.d0
iwork = 0
A0tmp = A0_aij
L0_A = 0.d0
R_A = 0.d0
n_out = n
! Then solve eigen vectors and values
write(*, FMT='(A)') "IMCOM: Performing eigendecomposition of A matrix"
call DSYEVR(jobz, range, uplo, n, A0tmp, n, vl, vu, il, iu, abstol, n_out, &
            L0_A, R_A, n, isupport, work, lwork, iwork, liwork, info)
! Test
if (info.ne.0) then
  if (info.lt.0) then
    write(*, FMT='(A, I7)') "IMCOM ERROR: Illegal value in argument i = ", -info
    stop
  else
    write(*, FMT='(A)') "IMCOM ERROR: Internal error in DSYEVR"
    stop
  end if
end if
! Reverse results into descending eigenvalue order 
! (important for minimizing roundoff when summing inverse)
!$omp parallel
!$omp workshare
L0_A = L0_A(n: 1: -1)
!$omp end workshare
!$omp workshare
R_A = R_A(:, n: 1: -1)
!$omp end workshare
!$omp end parallel
call imcom_writefits("./fits/Rtest.fits", n, n, real(R_A, 8))
call imcom_writefits("./fits/Ltest.fits", n, 1, real(L0_A, 8))
END SUBROUTINE imcom_eigen_A

!---

FUNCTION imcom_Ti_from_kappa(a_dum, n_dum, k_dum)
! Not threaded, for control of T summation order.
! Thread higher up the food chain.
implicit none
integer, intent(IN) :: a_dum, n_dum
real(KIND=8), intent(IN) :: k_dum
real(KIND=8), dimension(n_dum) :: imcom_Ti_from_kappa
real(KIND=8), dimension(n_dum) :: Bi, Ti_tmp
real(KIND=8), external :: DDOT
integer :: iev

Bi = B_ia(1:n_dum, a_dum)
Ti_tmp = 0.d0
do iev=1, n_dum

  Ti_tmp = Ti_tmp &
         + R_A(1:n_dum, iev) * DDOT(n_dum, R_A(1:n_dum, iev), 1, Bi, 1) &
         / (L0_A(iev) + k_dum * Ndiag(iev))

end do
imcom_Ti_from_kappa = -0.5d0 * Ti_tmp
END FUNCTION imcom_Ti_from_kappa

!---

SUBROUTINE imcom_Ti_from_kappa_test(a_dum, n_dum, k_dum, tol, flag)
! Not threaded, for control of T summation order.
! Thread higher up the food chain.
implicit none
integer, intent(IN) :: a_dum, n_dum
real(KIND=8), intent(IN) :: k_dum, tol
real(KIND=8), dimension(n_dum) :: Bi, Ti_tmp, Ti_rev_tmp
integer, intent(OUT) :: flag
real(KIND=8), external :: DDOT
integer :: iev

Bi = B_ia(1:n_dum, a_dum)
flag = 0
Ti_tmp = 0.d0
do iev=1, n_dum

  Ti_tmp = Ti_tmp &
         + R_A(1:n_dum, iev) * DDOT(n_dum, R_A(1:n_dum, iev), 1, Bi, 1) &
         / (L0_A(iev) + k_dum * Ndiag(iev))


end do
Ti_rev_tmp = 0.d0
do iev=n_dum, 1, -1

  Ti_rev_tmp = Ti_rev_tmp &
             + R_A(1:n_dum, iev) * DDOT(n_dum, R_A(1:n_dum, iev), 1, Bi, 1) &
             / (L0_A(iev) + k_dum * Ndiag(iev))

end do
print *, maxval((Ti_tmp - Ti_rev_tmp) / Ti_tmp)
if (maxval(abs((Ti_tmp - Ti_rev_tmp) / Ti_tmp)).gt.tol) flag = 1
END SUBROUTINE imcom_Ti_from_kappa_test


!---

END MODULE imcom_bisect

