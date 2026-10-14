!
! Test the Helmholtz FMM at high frequency, where the boxes at the
! top levels of the tree are many wavelengths in size and the
! multipole-to-local translations at those levels use the diagonal
! form (see h3ddiag.f90).
!
! Sources and targets are placed on the surface of the unit sphere.
! The potential and gradient at the targets due to charges and
! dipoles at the sources are compared with direct evaluation at a
! subset of the targets.
!
program test_hfmm3d_wideband
  use iso_fortran_env
  !$ use omp_lib
  implicit none
  integer(int64) :: ntests, nsuccess, isucc, ik, iprec
  integer(int64) :: ns, nt, ncores
  real(real64) :: eps
  complex(real64) :: zk

  ntests = 0
  nsuccess = 0

  call prini(6,13)

  ns = 200000
  nt = 200000

  ncores = 1
  !$ ncores = omp_get_max_threads()

  write(*,'(a,i0,a,i0,a,i0)') ' ns = ', ns, ', nt = ', nt, &
      ', cores (OpenMP threads) = ', ncores
  write(*,'(a)') ' pts/sec = (ns + nt)/time'
  write(*,*)
  write(*,'(a16,a10,a12,a10,a13,a13)') 'zk (re,im)', 'eps', &
      'rel err', 'time (s)', 'pts/sec', 'pts/sec/core'

  do ik = 1,2
    if (ik .eq. 1) zk = 150.0d0
    if (ik .eq. 2) zk = 300.0d0
    do iprec = 1,3
      if (iprec .eq. 1) eps = 0.51d-3
      if (iprec .eq. 2) eps = 0.51d-6
      if (iprec .eq. 3) eps = 0.51d-9
      call test_wideband(zk, eps, ns, nt, ncores, isucc)
      ntests = ntests + 1
      nsuccess = nsuccess + isucc
    end do
  end do

  open(unit=33, file='print_testres.txt', access='append')
  write(33,'(a,i2,a,i2,a)') 'Successfully completed ', nsuccess, &
      ' out of ', ntests, ' in hfmm3d wideband testing suite'
  close(33)

end program test_hfmm3d_wideband





subroutine test_wideband(zk, eps, ns, nt, ncores, isuccess)
  use iso_fortran_env
  !$ use omp_lib
  implicit none
  complex(real64), intent(in) :: zk
  real(real64), intent(in) :: eps
  integer(int64), intent(in) :: ns, nt, ncores
  integer(int64), intent(out) :: isuccess

  integer(int64) :: ntest, i, ier, nd
  real(real64) :: pi, thet, phi, thresh, err, ra, t0, t1, pps
  real(real64), allocatable :: source(:,:), targ(:,:)
  complex(real64), allocatable :: charge(:), dipvec(:,:)
  complex(real64), allocatable :: pot(:), grad(:,:), pottarg(:), gradtarg(:,:)
  complex(real64), allocatable :: potex(:), gradex(:,:)
  real(real64), external :: hkrand
  complex(real64), parameter :: eye = (0.0d0,1.0d0)

  pi = 4*atan(1.0d0)
  nd = 1

  ntest = 100

  allocate(source(3,ns), charge(ns), dipvec(3,ns))
  allocate(pot(ns), grad(3,ns))
  allocate(targ(3,nt), pottarg(nt), gradtarg(3,nt))
  allocate(potex(ntest), gradex(3,ntest))

  do i = 1,ns
    thet = acos(2*hkrand(0)-1)
    phi = 2*pi*hkrand(0)
    source(1,i) = sin(thet)*cos(phi)
    source(2,i) = sin(thet)*sin(phi)
    source(3,i) = cos(thet)
    charge(i) = hkrand(0) - 0.5d0 + eye*(hkrand(0) - 0.5d0)
    dipvec(1,i) = hkrand(0) - 0.5d0 + eye*(hkrand(0) - 0.5d0)
    dipvec(2,i) = hkrand(0) - 0.5d0 + eye*(hkrand(0) - 0.5d0)
    dipvec(3,i) = hkrand(0) - 0.5d0 + eye*(hkrand(0) - 0.5d0)
  end do

  do i = 1,nt
    thet = acos(2*hkrand(0)-1)
    phi = 2*pi*hkrand(0)
    targ(1,i) = sin(thet)*cos(phi)
    targ(2,i) = sin(thet)*sin(phi)
    targ(3,i) = cos(thet)
  end do

  pot = 0
  grad = 0
  pottarg = 0
  gradtarg = 0

  call cpu_time(t0)
  !$ t0 = omp_get_wtime()
  call hfmm3d_st_cd_g(eps, zk, ns, source, charge, dipvec, pot, grad, &
      nt, targ, pottarg, gradtarg, ier)
  call cpu_time(t1)
  !$ t1 = omp_get_wtime()

  potex = 0
  gradex = 0
  thresh = 2.0d0**(-51)
  call h3ddirectcdg(nd, zk, source, charge, dipvec, ns, targ, ntest, &
      potex, gradex, thresh)

  err = 0
  ra = 0
  do i = 1,ntest
    err = err + abs(pottarg(i)-potex(i))**2
    err = err + sum(abs(gradtarg(:,i)-gradex(:,i))**2)
    ra = ra + abs(potex(i))**2 + sum(abs(gradex(:,i))**2)
  end do
  err = sqrt(err/ra)

  pps = (ns + nt)/(t1 - t0)
  write(*,'(2f8.2,es10.2,es12.4,f10.2,2i13)') zk, eps, err, t1-t0, &
      nint(pps, int64), nint(pps/ncores, int64)

  isuccess = 0
  if (ier .eq. 0 .and. err .lt. eps) isuccess = 1

  return
end subroutine test_wideband
