!
! This file contains routines for the high-frequency diagonal form
! of the multipole-to-local translation operator (Rokhlin, 1990, 1993)
! for the Helmholtz equation in three dimensions.
!
! These operators are only stable when the boxes are large relative
! to the wavelength. The FMM uses them at levels where the boxes are
! many wavelengths across, and the low-frequency operators elsewhere.
!
! Conventions follow helmrouts3d.f: an outgoing expansion is
!
!   u(x) = sum_{n,m} mpole(n,m) h_n(k r) Y_nm(theta,phi)
!
! and an incoming expansion is
!
!   u(x) = sum_{n,m} local(n,m) j_n(k r) Y_nm(theta,phi),
!
! where Y_nm = sqrt(2n+1) sqrt((n-|m|)!/(n+|m|)!) P_n^|m| e^{i m phi}
! and the stored coefficients are scaled by rscale^(-n) and
! rscale^(n), respectively.
!
! The far-field signature of an outgoing expansion is
!
!   F(s) = sum_{n,m} (-i)^n mpole(n,m) Y_nm(s),   |s| = 1.
!
! For a source box centered at c1 and a target box centered at c2,
! with R = c2 - c1, the potential at c2 + y is
!
!   u(c2+y) = 1/(4 pi) \int_S e^{i k s.y} T_L(s,R) F(s) ds
!
! with the diagonal transfer function
!
!   T_L(s,R) = sum_{l=0}^{L} i^l (2l+1) h_l(k|R|) P_l(s.R/|R|).
!
! The incoming signature G = sum T_L F is converted to local
! coefficients via
!
!   local(n,m) = i^n/(4 pi) \int_S G(s) conj(Y_nm(s)) ds.
!
! The sphere is discretized with nth Gauss-Legendre nodes in
! cos(theta) and nph equispaced nodes in phi. Signatures are stored
! as sig(nd,nph,nth).
!
!-----------------------------------------------------------------------
!
! h3ddiagterms: order of the transfer function for list 2 boxes
!
! h3ddiaggridsize: number of quadrature nodes on the sphere
!
! h3ddiaginit: precompute quadrature and spherical harmonic tables
!
! h3dmp2sig: convert a multipole expansion into its far-field
!            signature
!
! h3dsig2loc: convert an incoming signature into a local expansion
!             and add it
!
! h3ddiagtrans: tabulate the diagonal transfer function for a given
!               translation vector
!
! h3ddiagtransall: tabulate the transfer functions for all list 2
!                  translation vectors at a level
!
! h3ddiagind: index of a list 2 translation vector in the table
!
! h3ddiagapply: multiply a signature by a transfer function and
!               accumulate
!
! h3dmplocdiaglev: all list 2 multipole-to-local translations at one
!                  level of the FMM tree, via the diagonal form
!
!-----------------------------------------------------------------------



subroutine h3ddiagterms(bs, zk, eps, nterms, nlterms, ier)
  !
  ! Determine the order of the diagonal transfer function for list 2
  ! interactions between boxes of size bs, and whether the diagonal
  ! form can attain the precision eps at this box size.
  !
  ! The transfer function is the truncation of the Gegenbauer series
  !
  !   h_0(k|R+d|) = sum_l (2l+1) h_l(k|R|) j_l(k|d|) P_l(.)
  !
  ! where d is the displacement from a source to a target relative
  ! to their box centers, |d| <= sqrt(3) bs, and |R| >= 2 bs is the
  ! distance between the box centers. With z1 = 2 k bs and
  ! z2 = sqrt(3) k bs, the truncation error and the roundoff error in
  ! the transfer function are estimated by
  !
  !   etrunc(L) = ctrunc z1 (2L+1) |h_L(z1) j_L(z2)|
  !   eround(L) = cround z1 epsmach sum_{l<=L} (2l+1) |h_l(z1)|
  !
  ! The constants ctrunc and cround were calibrated against the
  ! measured worst-case error (sources and targets on opposite box
  ! corners) for boxes 10 and 20 wavelengths in size.
  !
  ! The terms of the series decay only once l > |z2|, so the order is
  ! the smallest L > |z2| with etrunc(L) < eps. Since h_l(z1)
  ! grows rapidly for l > z1, eround(L) increases with L, and the
  ! diagonal form is unstable if eround(L) > eps. This is always the
  ! case for small boxes, and for large boxes when eps is small.
  !
  ! Since the multipole and local expansions have order nterms, the
  ! translation is exact (up to quadrature) for L = 2*nterms, and the
  ! order is capped there.
  !
  ! Input:
  !   bs      - box size
  !   zk      - Helmholtz parameter
  !   eps     - precision requested
  !   nterms  - order of the multipole and local expansions
  !
  ! Output:
  !   nlterms - order of the transfer function
  !   ier     - error flag
  !               ier = 0, diagonal form is stable at precision eps
  !               ier = 1, diagonal form is unstable at precision eps
  !
  implicit none
  real(8), intent(in) :: bs, eps
  complex(8), intent(in) :: zk
  integer(8), intent(in) :: nterms
  integer(8), intent(out) :: nlterms, ier

  real(8), parameter :: ctrunc = 1.0d-3, cround = 1.0d-1
  integer(8) :: ifder, l, lmin, lmax
  real(8) :: epsmach, etrunc, eround, hsum
  complex(8) :: z1, z2, fhder(0:1), fjder(0:1)
  complex(8), allocatable :: hfun(:), jfun(:)

  epsmach = epsilon(1.0d0)
  lmax = max(2*nterms, 1_8)
  allocate(hfun(0:lmax), jfun(0:lmax))

  z1 = zk*bs*2
  z2 = zk*bs*sqrt(3.0d0)

  ifder = 0
  call h3dall(lmax, z1, 1.0d0, hfun, ifder, fhder)
  call besseljs3d(lmax, z2, 1.0d0, jfun, ifder, fjder)

  ier = 0
  nlterms = lmax
  lmin = min(ceiling(abs(z2), 8), lmax)
  hsum = 0
  do l = 0,lmax
    hsum = hsum + (2*l+1)*abs(hfun(l))
    if (l .lt. lmin) cycle
    etrunc = ctrunc*abs(z1)*(2*l+1)*abs(hfun(l)*jfun(l))
    if (etrunc .lt. eps) then
      nlterms = l
      exit
    end if
  end do

  eround = cround*abs(z1)*epsmach*hsum
  if (eround .gt. eps) ier = 1

  return
end subroutine h3ddiagterms





subroutine h3ddiaggridsize(nterms, nlterms, nth, nph)
  !
  ! Returns the size of the quadrature grid on the sphere. The grid
  ! integrates exactly the product of a signature of degree nterms,
  ! a transfer function of degree nlterms and a spherical harmonic
  ! of degree nterms.
  !
  ! Input:
  !   nterms  - order of the multipole/local expansions
  !   nlterms - order of the transfer function
  !
  ! Output:
  !   nth     - number of Gauss-Legendre nodes in cos(theta)
  !   nph     - number of equispaced nodes in phi
  !
  implicit none
  integer(8), intent(in) :: nterms, nlterms
  integer(8), intent(out) :: nth, nph

  integer(8) :: ndeg

  ndeg = 2*nterms + nlterms
  nth = ndeg/2 + 1
  nph = ndeg + 1
  if (mod(nph,2_8) .eq. 1) nph = nph + 1

  return
end subroutine h3ddiaggridsize





subroutine h3ddiaginit(nterms, nth, nph, xnodes, wts, ynms, ephs)
  !
  ! Precompute the quadrature nodes and weights on the sphere and
  ! tables of normalized Legendre functions and exponentials.
  !
  ! Input:
  !   nterms - order of the expansions
  !   nth    - number of Gauss-Legendre nodes in cos(theta)
  !   nph    - number of equispaced nodes in phi
  !
  ! Output:
  !   xnodes - Gauss-Legendre nodes in cos(theta)
  !   wts    - quadrature weights, including the factor 2 pi/nph
  !            from the phi integration
  !   ynms   - ynms(n,m,i) = Y_n^m at cos(theta) = xnodes(i)
  !   ephs   - ephs(m,j) = e^{i m phi_j}
  !
  implicit none
  integer(8), intent(in) :: nterms, nth, nph
  real(8), intent(out) :: xnodes(nth), wts(nth)
  real(8), intent(out) :: ynms(0:nterms,0:nterms,nth)
  complex(8), intent(out) :: ephs(-nterms:nterms,nph)

  integer(8) :: i, j, m, ifwhts
  real(8) :: pi, phi

  pi = 4*atan(1.0d0)

  ifwhts = 1
  call legewhts(nth, xnodes, wts, ifwhts)
  do i = 1,nth
    wts(i) = wts(i)*2*pi/nph
    call ylgndr(nterms, xnodes(i), ynms(0,0,i))
  end do

  do j = 1,nph
    phi = 2*pi*(j-1)/nph
    do m = -nterms,nterms
      ephs(m,j) = dcmplx(cos(m*phi), sin(m*phi))
    end do
  end do

  return
end subroutine h3ddiaginit





subroutine h3dmp2sig(nd, nterms, rscale, mpole, nth, nph, ynms, &
    ephs, sig)
  !
  ! Compute the far-field signature of a multipole expansion
  !
  !   sig(s) = sum_{n,m} (-i)^n rscale^n mpole(n,m) Y_nm(s)
  !
  ! on the quadrature grid. The cost is O(nterms^3).
  !
  ! Input:
  !   nd      - number of expansions
  !   nterms  - order of the multipole expansion
  !   rscale  - scaling parameter of the multipole expansion
  !   mpole   - multipole expansion
  !   nth,nph - grid size
  !   ynms    - Legendre table from h3ddiaginit
  !   ephs    - exponential table from h3ddiaginit
  !
  ! Output:
  !   sig     - far-field signature
  !
  implicit none
  integer(8), intent(in) :: nd, nterms, nth, nph
  real(8), intent(in) :: rscale, ynms(0:nterms,0:nterms,nth)
  complex(8), intent(in) :: mpole(nd,0:nterms,-nterms:nterms)
  complex(8), intent(in) :: ephs(-nterms:nterms,nph)
  complex(8), intent(out) :: sig(nd,nph,nth)

  integer(8) :: i, j, m, n
  complex(8), allocatable :: fac(:), g(:,:), scmp(:,:,:)
  complex(8), parameter :: eye = (0.0d0,1.0d0)

  allocate(fac(0:nterms), g(nd,-nterms:nterms))
  allocate(scmp(nd,0:nterms,-nterms:nterms))

  fac(0) = 1
  do n = 1,nterms
    fac(n) = fac(n-1)*(-eye)*rscale
  end do

  do m = -nterms,nterms
    do n = abs(m),nterms
      scmp(:,n,m) = mpole(:,n,m)*fac(n)
    end do
  end do

  do i = 1,nth
    do m = -nterms,nterms
      g(:,m) = 0
      do n = abs(m),nterms
        g(:,m) = g(:,m) + scmp(:,n,m)*ynms(n,abs(m),i)
      end do
    end do

    do j = 1,nph
      sig(:,j,i) = 0
      do m = -nterms,nterms
        sig(:,j,i) = sig(:,j,i) + g(:,m)*ephs(m,j)
      end do
    end do
  end do

  return
end subroutine h3dmp2sig





subroutine h3dsig2loc(nd, nterms, rscale, sig, nth, nph, wts, ynms, &
    ephs, local)
  !
  ! Convert an incoming signature into a local expansion and add it
  ! to local:
  !
  !   local(n,m) = local(n,m) +
  !       rscale^n i^n/(4 pi) \int_S sig(s) conj(Y_nm(s)) ds
  !
  ! The cost is O(nterms^3).
  !
  ! Input:
  !   nd      - number of expansions
  !   nterms  - order of the local expansion
  !   rscale  - scaling parameter of the local expansion
  !   sig     - incoming signature
  !   nth,nph - grid size
  !   wts     - quadrature weights from h3ddiaginit
  !   ynms    - Legendre table from h3ddiaginit
  !   ephs    - exponential table from h3ddiaginit
  !
  ! Input/output:
  !   local   - local expansion, incremented
  !
  implicit none
  integer(8), intent(in) :: nd, nterms, nth, nph
  real(8), intent(in) :: rscale, wts(nth)
  real(8), intent(in) :: ynms(0:nterms,0:nterms,nth)
  complex(8), intent(in) :: sig(nd,nph,nth), ephs(-nterms:nterms,nph)
  complex(8), intent(inout) :: local(nd,0:nterms,-nterms:nterms)

  integer(8) :: i, j, m, n
  complex(8), allocatable :: fac(:), h(:,:), tmp(:,:,:)
  complex(8), parameter :: eye = (0.0d0,1.0d0)
  real(8) :: pi, dtmp

  pi = 4*atan(1.0d0)

  allocate(fac(0:nterms), h(nd,-nterms:nterms))
  allocate(tmp(nd,0:nterms,-nterms:nterms))

  fac(0) = 1/(4*pi)
  do n = 1,nterms
    fac(n) = fac(n-1)*eye*rscale
  end do

  tmp = 0

  do i = 1,nth
    h = 0
    do j = 1,nph
      do m = -nterms,nterms
        h(:,m) = h(:,m) + sig(:,j,i)*dconjg(ephs(m,j))
      end do
    end do

    do m = -nterms,nterms
      do n = abs(m),nterms
        dtmp = wts(i)*ynms(n,abs(m),i)
        tmp(:,n,m) = tmp(:,n,m) + dtmp*h(:,m)
      end do
    end do
  end do

  do m = -nterms,nterms
    do n = abs(m),nterms
      local(:,n,m) = local(:,n,m) + fac(n)*tmp(:,n,m)
    end do
  end do

  return
end subroutine h3dsig2loc





subroutine h3ddiagtrans(nlterms, zk, rvec, nth, nph, xnodes, trans)
  !
  ! Tabulate the diagonal transfer function
  !
  !   trans(s) = sum_{l=0}^{nlterms} i^l (2l+1) h_l(k|R|) P_l(s.R/|R|)
  !
  ! on the quadrature grid.
  !
  ! Input:
  !   nlterms - order of the transfer function
  !   zk      - Helmholtz parameter
  !   rvec    - translation vector R = (target center - source center)
  !   nth,nph - grid size
  !   xnodes  - Gauss-Legendre nodes in cos(theta)
  !
  ! Output:
  !   trans   - transfer function, trans(nph,nth)
  !
  implicit none
  integer(8), intent(in) :: nlterms, nth, nph
  real(8), intent(in) :: rvec(3), xnodes(nth)
  complex(8), intent(in) :: zk
  complex(8), intent(out) :: trans(nph,nth)

  integer(8) :: i, j, l, ifder
  real(8) :: r, pi, phi, ct, st, x, pm1, p0, p1
  complex(8), allocatable :: hfac(:)
  complex(8) :: z, hder(0:1), ztmp
  complex(8), parameter :: eye = (0.0d0,1.0d0)

  pi = 4*atan(1.0d0)

  ! h3dall always computes h_0 and h_1
  allocate(hfac(0:max(nlterms,1_8)))

  r = sqrt(rvec(1)**2 + rvec(2)**2 + rvec(3)**2)
  z = zk*r
  ifder = 0
  call h3dall(nlterms, z, 1.0d0, hfac, ifder, hder)
  ztmp = 1
  do l = 0,nlterms
    hfac(l) = hfac(l)*ztmp*(2*l+1)
    ztmp = ztmp*eye
  end do

  do i = 1,nth
    ct = xnodes(i)
    st = sqrt(1-ct**2)
    do j = 1,nph
      phi = 2*pi*(j-1)/nph
      x = (st*cos(phi)*rvec(1) + st*sin(phi)*rvec(2) + ct*rvec(3))/r
      pm1 = 1
      p0 = x
      trans(j,i) = hfac(0)
      if (nlterms .ge. 1) trans(j,i) = trans(j,i) + hfac(1)*x
      do l = 2,nlterms
        p1 = ((2*l-1)*x*p0 - (l-1)*pm1)/l
        trans(j,i) = trans(j,i) + hfac(l)*p1
        pm1 = p0
        p0 = p1
      end do
    end do
  end do

  return
end subroutine h3ddiagtrans





subroutine h3ddiagtransall(nlterms, zk, bs, nth, nph, xnodes, trans)
  !
  ! Tabulate the diagonal transfer functions for all translation
  ! vectors bs*(ix,iy,iz), -3 <= ix,iy,iz <= 3, that occur in list 2.
  ! Translations between adjacent boxes are skipped.
  !
  ! The translation (ix,iy,iz) is stored in
  !   trans(:,:,h3ddiagind(ix,iy,iz))
  !
  ! Input:
  !   nlterms - order of the transfer function
  !   zk      - Helmholtz parameter
  !   bs      - box size
  !   nth,nph - grid size
  !   xnodes  - Gauss-Legendre nodes in cos(theta)
  !
  ! Output:
  !   trans   - transfer functions, trans(nph,nth,343)
  !
  implicit none
  integer(8), intent(in) :: nlterms, nth, nph
  real(8), intent(in) :: bs, xnodes(nth)
  complex(8), intent(in) :: zk
  complex(8), intent(out) :: trans(nph,nth,343)

  integer(8) :: ix, iy, iz, ind
  integer(8), external :: h3ddiagind
  real(8) :: rvec(3)

  !$omp parallel do default(shared) private(ix,iy,iz,ind,rvec) &
  !$omp collapse(3) schedule(dynamic)
  do iz = -3,3
    do iy = -3,3
      do ix = -3,3
        if (max(abs(ix),abs(iy),abs(iz)) .ge. 2) then
          ind = h3ddiagind(ix, iy, iz)
          rvec(1) = bs*ix
          rvec(2) = bs*iy
          rvec(3) = bs*iz
          call h3ddiagtrans(nlterms, zk, rvec, nth, nph, xnodes, &
              trans(1,1,ind))
        end if
      end do
    end do
  end do
  !$omp end parallel do

  return
end subroutine h3ddiagtransall





function h3ddiagind(ix, iy, iz)
  !
  ! Index of the translation vector (ix,iy,iz), -3 <= ix,iy,iz <= 3,
  ! in the table of transfer functions.
  !
  implicit none
  integer(8) :: h3ddiagind
  integer(8), intent(in) :: ix, iy, iz

  h3ddiagind = (ix+3) + 7*(iy+3) + 49*(iz+3) + 1

  return
end function h3ddiagind





subroutine h3ddiagapply(nd, nsig, trans, sigout, sigin)
  !
  ! Apply the diagonal translation operator:
  !
  !   sigin = sigin + trans .* sigout
  !
  implicit none
  integer(8), intent(in) :: nd, nsig
  complex(8), intent(in) :: trans(nsig), sigout(nd,nsig)
  complex(8), intent(inout) :: sigin(nd,nsig)

  integer(8) :: i

  do i = 1,nsig
    sigin(:,i) = sigin(:,i) + trans(i)*sigout(:,i)
  end do

  return
end subroutine h3ddiagapply





subroutine h3dmplocdiaglev(nd, zk, nterms, nlterms, rscale, bs, &
    ibstart, ibend, nboxes, centers, isrcse, itargse, iexpcse, ifpgh, &
    ifpghtarg, iaddr, rmlexp, mnlist2, nlist2, list2)
  !
  ! Add to the local expansions of the boxes ibstart:ibend the
  ! contributions of the multipole expansions of the boxes in their
  ! list 2, using the diagonal form of the translation operator.
  !
  ! All boxes ibstart:ibend are on the same level, with box size bs,
  ! expansion order nterms and scaling parameter rscale. The caller
  ! is responsible for checking that the diagonal form is stable at
  ! this level (see h3ddiagterms).
  !
  ! Input:
  !   nd        - number of expansions
  !   zk        - Helmholtz parameter
  !   nterms    - order of the multipole and local expansions
  !   nlterms   - order of the transfer function (from h3ddiagterms)
  !   rscale    - scaling parameter of the expansions
  !   bs        - box size
  !   ibstart   - first box on the level
  !   ibend     - last box on the level
  !   nboxes    - number of boxes in the tree
  !   centers   - box centers
  !   isrcse    - source index ranges of boxes
  !   itargse   - target index ranges of boxes
  !   iexpcse   - expansion center index ranges of boxes
  !   ifpgh     - flag for evaluation at sources
  !   ifpghtarg - flag for evaluation at targets
  !   iaddr     - addresses of the expansions in rmlexp, in units of
  !               real(8), as set up by mpalloc
  !   mnlist2   - leading dimension of list2
  !   nlist2    - number of boxes in list 2
  !   list2     - list 2 of each box
  !
  ! Input/output:
  !   rmlexp    - multipole and local expansions; the local
  !               expansions are incremented
  !
  implicit none
  integer(8), intent(in) :: nd, nterms, nlterms, ibstart, ibend
  integer(8), intent(in) :: nboxes, ifpgh, ifpghtarg, mnlist2
  integer(8), intent(in) :: isrcse(2,nboxes), itargse(2,nboxes)
  integer(8), intent(in) :: iexpcse(2,nboxes), iaddr(2,nboxes)
  integer(8), intent(in) :: nlist2(nboxes), list2(mnlist2,nboxes)
  real(8), intent(in) :: rscale, bs, centers(3,nboxes)
  complex(8), intent(in) :: zk
  complex(8), intent(inout) :: rmlexp(*)

  integer(8) :: nth, nph, nsig, nsrcbox, ibox, jbox, i, npts
  integer(8) :: ix, iy, iz, ind, impole, ilocal
  integer(8), external :: h3ddiagind
  integer(8), allocatable :: isig(:)
  real(8), allocatable :: xnodes(:), wts(:), ynms(:,:,:)
  complex(8), allocatable :: ephs(:,:), trans(:,:), sigout(:,:)
  complex(8), allocatable :: sigin(:)

  call h3ddiaggridsize(nterms, nlterms, nth, nph)
  nsig = nth*nph

  allocate(xnodes(nth), wts(nth), ynms(0:nterms,0:nterms,nth))
  allocate(ephs(-nterms:nterms,nph), trans(nsig,343))

  call h3ddiaginit(nterms, nth, nph, xnodes, wts, ynms, ephs)
  call h3ddiagtransall(nlterms, zk, bs, nth, nph, xnodes, trans)

  !
  ! far-field signatures of the boxes with sources
  !
  allocate(isig(ibstart:ibend))
  nsrcbox = 0
  do ibox = ibstart,ibend
    isig(ibox) = 0
    if (isrcse(2,ibox) .ge. isrcse(1,ibox)) then
      nsrcbox = nsrcbox + 1
      isig(ibox) = nsrcbox
    end if
  end do

  allocate(sigout(nd*nsig,nsrcbox))

  !$omp parallel do default(shared) private(ibox,impole) &
  !$omp schedule(dynamic)
  do ibox = ibstart,ibend
    if (isig(ibox) .gt. 0) then
      impole = (iaddr(1,ibox)+1)/2
      call h3dmp2sig(nd, nterms, rscale, rmlexp(impole), nth, nph, &
          ynms, ephs, sigout(1,isig(ibox)))
    end if
  end do
  !$omp end parallel do

  !
  ! diagonal translations and conversion to local expansions
  !
  !$omp parallel default(shared) &
  !$omp private(ibox,jbox,i,npts,ix,iy,iz,ind,ilocal,sigin)
  allocate(sigin(nd*nsig))
  !$omp do schedule(dynamic)
  do ibox = ibstart,ibend
    npts = iexpcse(2,ibox) - iexpcse(1,ibox) + 1
    if (ifpghtarg .gt. 0) then
      npts = npts + itargse(2,ibox) - itargse(1,ibox) + 1
    end if
    if (ifpgh .gt. 0) then
      npts = npts + isrcse(2,ibox) - isrcse(1,ibox) + 1
    end if
    if (npts .le. 0) cycle

    sigin = 0
    npts = 0
    do i = 1,nlist2(ibox)
      jbox = list2(i,ibox)
      if (isig(jbox) .eq. 0) cycle
      ix = nint((centers(1,ibox) - centers(1,jbox))/bs, 8)
      iy = nint((centers(2,ibox) - centers(2,jbox))/bs, 8)
      iz = nint((centers(3,ibox) - centers(3,jbox))/bs, 8)
      ind = h3ddiagind(ix, iy, iz)
      call h3ddiagapply(nd, nsig, trans(1,ind), sigout(1,isig(jbox)), &
          sigin)
      npts = npts + 1
    end do

    if (npts .gt. 0) then
      ilocal = (iaddr(2,ibox)+1)/2
      call h3dsig2loc(nd, nterms, rscale, sigin, nth, nph, wts, ynms, &
          ephs, rmlexp(ilocal))
    end if
  end do
  !$omp end do
  deallocate(sigin)
  !$omp end parallel

  return
end subroutine h3dmplocdiaglev
