!
! Test the high-frequency diagonal multipole-to-local translation
! operator in h3ddiag.f90.
!
! Sources are placed in a box centered at c1, targets in a box
! centered at c2, with c2-c1 a list 2 translation vector. The
! potential at the targets is computed via
!
!   (1) direct evaluation
!   (2) h3dformmpc + h3dmploc + h3dtaevalp
!   (3) h3dformmpc + h3dmp2sig + h3ddiagtrans + h3ddiagapply
!       + h3dsig2loc + h3dtaevalp
!
program test_h3ddiag
  implicit none
  integer(8) :: ntests, nsuccess, isucc, ii, iprec, ioff
  integer(8) :: ioffs(3,2)
  real(8) :: pi, eps
  complex(8) :: zk

  pi = 4*atan(1.0d0)

  ! closest list 2 translation and a generic one
  ioffs(:,1) = (/ 2, 0, 0 /)
  ioffs(:,2) = (/ 2, -1, 3 /)

  ntests = 0
  nsuccess = 0

  call prini(6,13)

  do ii = 1,3
    if (ii .eq. 1) zk = 20*pi
    if (ii .eq. 2) zk = 40*pi
    if (ii .eq. 3) zk = 40*pi + dcmplx(0.0d0,0.1d0)
    do iprec = 1,3
      if (iprec .eq. 1) eps = 1.0d-3
      if (iprec .eq. 2) eps = 1.0d-6
      if (iprec .eq. 3) eps = 1.0d-9
      do ioff = 1,2
        call test_diag_m2l(zk, eps, ioffs(1,ioff), isucc)
        ntests = ntests + 1
        nsuccess = nsuccess + isucc
      end do
    end do
  end do

  open(unit=33, file='print_testres.txt', access='append')
  write(33,'(a,i2,a,i2,a)') 'Successfully completed ', nsuccess, &
      ' out of ', ntests, ' in h3ddiag testing suite'
  close(33)

end program test_h3ddiag





subroutine test_diag_m2l(zk, eps, ioff, isuccess)
  implicit none
  complex(8), intent(in) :: zk
  real(8), intent(in) :: eps
  integer(8), intent(in) :: ioff(3)
  integer(8), intent(out) :: isuccess

  integer(8) :: nd, ns, nt, nterms, nlterms, nth, nph, nsig, nquad
  integer(8) :: nlege, lw, lused, i, j, ix, iy, iz, ifwhts, ier
  real(8) :: bs, rscale, radius, thresh, ra, e1, e2
  real(8) :: t0, t1, tpas, tm2s, tapp, ts2l
  real(8) :: c1(3), c2(3), rvec(3)
  real(8), allocatable :: sources(:,:), targs(:,:)
  real(8), allocatable :: xnodes(:), wts(:), ynms(:,:,:)
  real(8), allocatable :: xq(:), wq(:), wlege(:)
  complex(8), allocatable :: charge(:), potex(:), pot1(:), pot2(:)
  complex(8), allocatable :: mpole(:,:), loc1(:,:), loc2(:,:)
  complex(8), allocatable :: ephs(:,:), sigout(:), sigin(:)
  complex(8), allocatable :: trans(:)
  real(8), external :: hkrand

  nd = 1
  bs = 1.0d0
  rscale = 1.0d0

  ns = 200
  nt = 200

  allocate(sources(3,ns), targs(3,nt), charge(ns))
  allocate(potex(nt), pot1(nt), pot2(nt))

  c1 = 0

  ix = ioff(1)
  iy = ioff(2)
  iz = ioff(3)
  c2(1) = c1(1) + ix*bs
  c2(2) = c1(2) + iy*bs
  c2(3) = c1(3) + iz*bs
  rvec = c2 - c1

  do i = 1,ns
    do j = 1,3
      sources(j,i) = c1(j) + (hkrand(0)-0.5d0)*bs
    end do
    charge(i) = dcmplx(hkrand(0)-0.5d0, hkrand(0)-0.5d0)
  end do

  do i = 1,nt
    do j = 1,3
      targs(j,i) = c2(j) + (hkrand(0)-0.5d0)*bs
    end do
  end do

  !
  ! put the first 8 sources and targets on the box corners, which
  ! gives the largest source-to-target displacement relative to the
  ! box centers
  !
  do i = 1,8
    sources(1,i) = c1(1) + (mod(i-1,2_8)-0.5d0)*bs
    sources(2,i) = c1(2) + (mod((i-1)/2,2_8)-0.5d0)*bs
    sources(3,i) = c1(3) + ((i-1)/4-0.5d0)*bs
    targs(1,i) = c2(1) + (mod(i-1,2_8)-0.5d0)*bs
    targs(2,i) = c2(2) + (mod((i-1)/2,2_8)-0.5d0)*bs
    targs(3,i) = c2(3) + ((i-1)/4-0.5d0)*bs
  end do
  potex = 0
  pot1 = 0
  pot2 = 0

  thresh = 1.0d-15
  call h3ddirectcp(nd, zk, sources, charge, ns, targs, nt, potex, &
      thresh)

  call h3dterms(bs, zk, eps, nterms)
  call h3ddiagterms(bs, zk, eps, nterms, nlterms, ier)

  nlege = nterms + 10
  lw = 4*(nlege+1)**2
  allocate(wlege(lw))
  call ylgndrfwini(nlege, wlege, lw, lused)

  allocate(mpole(0:nterms,-nterms:nterms))
  allocate(loc1(0:nterms,-nterms:nterms))
  allocate(loc2(0:nterms,-nterms:nterms))
  mpole = 0
  loc1 = 0
  loc2 = 0

  call h3dformmpc(nd, zk, rscale, sources, charge, ns, c1, nterms, &
      mpole, wlege, nlege)

  !
  ! rotate and shoot
  !
  nquad = max(6, int(2.2*nterms))
  allocate(xq(nquad), wq(nquad))
  ifwhts = 1
  call legewhts(nquad, xq, wq, ifwhts)
  radius = bs/2*sqrt(3.0d0)
  call cpu_time(t0)
  call h3dmploc(nd, zk, rscale, c1, mpole, nterms, rscale, c2, loc1, &
      nterms, radius, xq, wq, nquad)
  call cpu_time(t1)
  tpas = t1-t0

  !
  ! diagonal form
  !
  call h3ddiaggridsize(nterms, nlterms, nth, nph)
  nsig = nth*nph
  allocate(xnodes(nth), wts(nth), ynms(0:nterms,0:nterms,nth))
  allocate(ephs(-nterms:nterms,nph))
  allocate(sigout(nsig), sigin(nsig), trans(nsig))
  call h3ddiaginit(nterms, nth, nph, xnodes, wts, ynms, ephs)
  call h3ddiagtrans(nlterms, zk, rvec, nth, nph, xnodes, trans)

  call cpu_time(t0)
  call h3dmp2sig(nd, nterms, rscale, mpole, nth, nph, ynms, ephs, &
      sigout)
  call cpu_time(t1)
  tm2s = t1-t0

  sigin = 0
  call cpu_time(t0)
  call h3ddiagapply(nd, nsig, trans, sigout, sigin)
  call cpu_time(t1)
  tapp = t1-t0

  call cpu_time(t0)
  call h3dsig2loc(nd, nterms, rscale, sigin, nth, nph, wts, ynms, &
      ephs, loc2)
  call cpu_time(t1)
  ts2l = t1-t0

  call h3dtaevalp(nd, zk, rscale, c2, loc1, nterms, targs, nt, pot1, &
      wlege, nlege)
  call h3dtaevalp(nd, zk, rscale, c2, loc2, nterms, targs, nt, pot2, &
      wlege, nlege)

  ra = sum(abs(potex)**2)
  e1 = sqrt(sum(abs(pot1-potex)**2)/ra)
  e2 = sqrt(sum(abs(pot2-potex)**2)/ra)

  write(*,'(a,2f9.3,a,es8.1,a,3i3,a,i4,a,i4,a,i4,a,i4)') ' zk=', zk, &
      ' eps=', eps, ' off=', ioff, ' nterms=', nterms, ' L=', nlterms, &
      ' nth=', nth, ' nph=', nph
  write(*,'(a,es11.4,a,es11.4)') '   rel err rotate+shoot: ', e1, &
      '   diagonal: ', e2
  write(*,'(a,es10.3,a,3es10.3)') '   time rotate+shoot: ', tpas, &
      '   mp2sig, apply, sig2loc: ', tm2s, tapp, ts2l

  !
  ! if the diagonal form is flagged as unstable at this precision,
  ! the FMM does not use it, so there is nothing to check
  !
  isuccess = 0
  if (ier .ne. 0) then
    write(*,'(a)') '   diagonal form flagged unstable at this eps'
    isuccess = 1
  else if (e2 .lt. eps) then
    isuccess = 1
  end if

  return
end subroutine test_diag_m2l
