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
! h3ddiaghfac: coefficients of the Legendre series of the transfer
!              function
!
! h3ddiagtfun: evaluate the transfer function at one point
!
! h3ddiagtrans: tabulate the diagonal transfer function for a given
!               translation vector, directly (reference version)
!
! h3ddiagtabsize: size of the one-dimensional table of the transfer
!                 function
!
! h3ddiagtab: build the one-dimensional table of the transfer
!             function for a given translation length
!
! h3ddiagtranstab: tabulate the diagonal transfer function for a
!                  given translation vector, by interpolation from
!                  the one-dimensional table
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
  use iso_fortran_env
  implicit none
  real(real64), intent(in) :: bs, eps
  complex(real64), intent(in) :: zk
  integer(int64), intent(in) :: nterms
  integer(int64), intent(out) :: nlterms, ier

  real(real64), parameter :: ctrunc = 1.0d-3, cround = 1.0d-1
  integer(int64) :: ifder, l, lmin, lmax
  real(real64) :: epsmach, etrunc, eround, hsum
  complex(real64) :: z1, z2, fhder(0:1), fjder(0:1)
  complex(real64), allocatable :: hfun(:), jfun(:)

  epsmach = epsilon(1.0d0)
  lmax = max(2*nterms, 1_int64)
  allocate(hfun(0:lmax), jfun(0:lmax))

  z1 = zk*bs*2
  z2 = zk*bs*sqrt(3.0d0)

  ifder = 0
  call h3dall(lmax, z1, 1.0d0, hfun, ifder, fhder)
  call besseljs3d(lmax, z2, 1.0d0, jfun, ifder, fjder)

  ier = 0
  nlterms = lmax
  lmin = min(ceiling(abs(z2), int64), lmax)
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
  use iso_fortran_env
  implicit none
  integer(int64), intent(in) :: nterms, nlterms
  integer(int64), intent(out) :: nth, nph

  integer(int64) :: ndeg

  ndeg = 2*nterms + nlterms
  nth = ndeg/2 + 1
  nph = ndeg + 1
  if (mod(nph,2_int64) .eq. 1) nph = nph + 1

  return
end subroutine h3ddiaggridsize





subroutine h3ddiaginit(nterms, nth, nph, xnodes, wts, ynms, wsave)
  !
  ! Precompute the quadrature nodes and weights on the sphere, the
  ! table of normalized Legendre functions, and the FFT in phi.
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
  !   wsave  - initialized work array for the FFT of length nph,
  !            dimension 4*nph+20
  !
  use iso_fortran_env
  implicit none
  integer(int64), intent(in) :: nterms, nth, nph
  real(real64), intent(out) :: xnodes(nth), wts(nth)
  real(real64), intent(out) :: ynms(0:nterms,0:nterms,nth)
  real(real64), intent(out) :: wsave(4*nph+20)

  integer(int64) :: i, ifwhts
  ! the FFT routines in dfft.f take default integers
  integer :: nfft
  real(real64) :: pi

  pi = 4*atan(1.0d0)

  ifwhts = 1
  call legewhts(nth, xnodes, wts, ifwhts)
  do i = 1,nth
    wts(i) = wts(i)*2*pi/nph
    call ylgndr(nterms, xnodes(i), ynms(0,0,i))
  end do

  nfft = int(nph)
  call zffti(nfft, wsave)

  return
end subroutine h3ddiaginit





subroutine h3dmp2sig(nd, nterms, rscale, mpole, nth, nph, ynms, &
    wsave, sig)
  !
  ! Compute the far-field signature of a multipole expansion
  !
  !   sig(s) = sum_{n,m} (-i)^n rscale^n mpole(n,m) Y_nm(s)
  !
  ! on the quadrature grid. The sum over m is computed with an FFT
  ! at each latitude, and the cost is O(nth nterms^2 +
  ! nth nph log(nph)).
  !
  ! Input:
  !   nd      - number of expansions
  !   nterms  - order of the multipole expansion
  !   rscale  - scaling parameter of the multipole expansion
  !   mpole   - multipole expansion
  !   nth,nph - grid size
  !   ynms    - Legendre table from h3ddiaginit
  !   wsave   - FFT work array from h3ddiaginit
  !
  ! Output:
  !   sig     - far-field signature
  !
  use iso_fortran_env
  implicit none
  integer(int64), intent(in) :: nd, nterms, nth, nph
  real(real64), intent(in) :: rscale, ynms(0:nterms,0:nterms,nth)
  real(real64), intent(in) :: wsave(4*nph+20)
  complex(real64), intent(in) :: mpole(nd,0:nterms,-nterms:nterms)
  complex(real64), intent(out) :: sig(nd,nph,nth)

  integer(int64) :: i, j, m, n, idim
  integer :: nfft
  real(real64), allocatable :: w(:)
  complex(real64), allocatable :: fac(:), g(:,:), scmp(:,:,:), c(:)
  complex(real64), parameter :: eye = (0.0d0,1.0d0)

  allocate(fac(0:nterms), g(nd,-nterms:nterms))
  allocate(scmp(nd,0:nterms,-nterms:nterms))

  ! the FFT overwrites part of its work array, so each call needs
  ! its own copy
  allocate(w(4*nph+20), c(nph))
  w = wsave
  nfft = int(nph)

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

    ! sig(j) = sum_m g(m) e^{i m phi_j}
    do idim = 1,nd
      c = 0
      do m = -nterms,nterms
        j = modulo(m,nph) + 1
        c(j) = c(j) + g(idim,m)
      end do
      call zfftb(nfft, c, w)
      sig(idim,:,i) = c
    end do
  end do

  return
end subroutine h3dmp2sig





subroutine h3dsig2loc(nd, nterms, rscale, sig, nth, nph, wts, ynms, &
    wsave, local)
  !
  ! Convert an incoming signature into a local expansion and add it
  ! to local:
  !
  !   local(n,m) = local(n,m) +
  !       rscale^n i^n/(4 pi) \int_S sig(s) conj(Y_nm(s)) ds
  !
  ! The integral in phi is computed with an FFT at each latitude,
  ! and the cost is O(nth nterms^2 + nth nph log(nph)).
  !
  ! Input:
  !   nd      - number of expansions
  !   nterms  - order of the local expansion
  !   rscale  - scaling parameter of the local expansion
  !   sig     - incoming signature
  !   nth,nph - grid size
  !   wts     - quadrature weights from h3ddiaginit
  !   ynms    - Legendre table from h3ddiaginit
  !   wsave   - FFT work array from h3ddiaginit
  !
  ! Input/output:
  !   local   - local expansion, incremented
  !
  use iso_fortran_env
  implicit none
  integer(int64), intent(in) :: nd, nterms, nth, nph
  real(real64), intent(in) :: rscale, wts(nth)
  real(real64), intent(in) :: ynms(0:nterms,0:nterms,nth)
  real(real64), intent(in) :: wsave(4*nph+20)
  complex(real64), intent(in) :: sig(nd,nph,nth)
  complex(real64), intent(inout) :: local(nd,0:nterms,-nterms:nterms)

  integer(int64) :: i, m, n, idim
  integer :: nfft
  real(real64), allocatable :: w(:)
  complex(real64), allocatable :: fac(:), h(:,:), tmp(:,:,:), c(:)
  complex(real64), parameter :: eye = (0.0d0,1.0d0)
  real(real64) :: pi, dtmp

  pi = 4*atan(1.0d0)

  allocate(fac(0:nterms), h(nd,-nterms:nterms))
  allocate(tmp(nd,0:nterms,-nterms:nterms))

  ! the FFT overwrites part of its work array, so each call needs
  ! its own copy
  allocate(w(4*nph+20), c(nph))
  w = wsave
  nfft = int(nph)

  fac(0) = 1/(4*pi)
  do n = 1,nterms
    fac(n) = fac(n-1)*eye*rscale
  end do

  tmp = 0

  do i = 1,nth
    ! h(m) = sum_j sig(j) e^{-i m phi_j}
    do idim = 1,nd
      c = sig(idim,:,i)
      call zfftf(nfft, c, w)
      do m = -nterms,nterms
        h(idim,m) = c(modulo(m,nph) + 1)
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





subroutine h3ddiaghfac(nlterms, zk, r, hfac)
  !
  ! Coefficients of the Legendre series of the diagonal transfer
  ! function for a translation of length r:
  !
  !   hfac(l) = i^l (2l+1) h_l(k r),   l = 0,...,nlterms
  !
  ! Input:
  !   nlterms - order of the transfer function
  !   zk      - Helmholtz parameter
  !   r       - length of the translation vector
  !
  ! Output:
  !   hfac    - coefficients, dimension 0:max(nlterms,1)
  !
  use iso_fortran_env
  implicit none
  integer(int64), intent(in) :: nlterms
  real(real64), intent(in) :: r
  complex(real64), intent(in) :: zk
  complex(real64), intent(out) :: hfac(0:max(nlterms,1_int64))

  integer(int64) :: l, ifder
  complex(real64) :: z, hder(0:1), ztmp, ima

  ima = (0,1)

  ! h3dall always computes h_0 and h_1
  z = zk*r
  ifder = 0
  call h3dall(nlterms, z, 1.0d0, hfac, ifder, hder)
  ztmp = 1
  do l = 0,nlterms
    hfac(l) = hfac(l)*ztmp*(2*l+1)
    ztmp = ztmp*ima
  end do

  return
end subroutine h3ddiaghfac





function h3ddiagtfun(nlterms, hfac, x)
  !
  ! Evaluate the diagonal transfer function
  !
  !   T(x) = sum_{l=0}^{nlterms} hfac(l) P_l(x)
  !
  ! by the Legendre recurrence. The cost is O(nlterms).
  !
  ! Input:
  !   nlterms - order of the transfer function
  !   hfac    - coefficients from h3ddiaghfac
  !   x       - cosine of the angle between the direction on the
  !             sphere and the translation vector
  !
  use iso_fortran_env
  implicit none
  complex(real64) :: h3ddiagtfun
  integer(int64), intent(in) :: nlterms
  complex(real64), intent(in) :: hfac(0:max(nlterms,1_int64))
  real(real64), intent(in) :: x

  integer(int64) :: l
  real(real64) :: pm1, p0, p1
  complex(real64) :: t

  pm1 = 1
  p0 = x
  t = hfac(0)
  if (nlterms .ge. 1) t = t + hfac(1)*x
  do l = 2,nlterms
    p1 = ((2*l-1)*x*p0 - (l-1)*pm1)/l
    t = t + hfac(l)*p1
    pm1 = p0
    p0 = p1
  end do
  h3ddiagtfun = t

  return
end function h3ddiagtfun





subroutine h3ddiagtrans(nlterms, zk, rvec, nth, nph, xnodes, trans)
  !
  ! Tabulate the diagonal transfer function
  !
  !   trans(s) = sum_{l=0}^{nlterms} i^l (2l+1) h_l(k|R|) P_l(s.R/|R|)
  !
  ! on the quadrature grid, by evaluating the Legendre series at
  ! every grid point. The cost is O(nth nph nlterms). This is the
  ! reference version; the FMM uses h3ddiagtab and h3ddiagtranstab.
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
  use iso_fortran_env
  implicit none
  integer(int64), intent(in) :: nlterms, nth, nph
  real(real64), intent(in) :: rvec(3), xnodes(nth)
  complex(real64), intent(in) :: zk
  complex(real64), intent(out) :: trans(nph,nth)

  integer(int64) :: i, j
  real(real64) :: r, pi, phi, ct, st, x
  complex(real64), allocatable :: hfac(:)
  complex(real64), external :: h3ddiagtfun

  pi = 4*atan(1.0d0)

  allocate(hfac(0:max(nlterms,1_int64)))

  r = sqrt(rvec(1)**2 + rvec(2)**2 + rvec(3)**2)
  call h3ddiaghfac(nlterms, zk, r, hfac)

  do i = 1,nth
    ct = xnodes(i)
    st = sqrt(1-ct**2)
    do j = 1,nph
      phi = 2*pi*(j-1)/nph
      x = (st*cos(phi)*rvec(1) + st*sin(phi)*rvec(2) + ct*rvec(3))/r
      trans(j,i) = h3ddiagtfun(nlterms, hfac, x)
    end do
  end do

  return
end subroutine h3ddiagtrans





subroutine h3ddiagtabsize(nlterms, norder, npan)
  !
  ! Size of the piecewise Chebyshev table of the transfer function
  ! built by h3ddiagtab.
  !
  ! As a function of the angle gamma, x = cos(gamma), the transfer
  ! function is a trigonometric polynomial of degree nlterms. In the
  ! table variable w (see h3ddiagtab) its phase advances by at most
  ! 2 sqrt(2) nlterms per unit w, so with npan = 8*nlterms panels on
  ! [0, 1/sqrt(2)] the phase changes by at most 1/4 across a panel.
  ! The truncation error of the interpolant of order 8 is then about
  ! 2 (1/16)^8/8! = 1e-14 relative to the size of the function, which
  ! is below the roundoff error in the Legendre series itself.
  !
  ! Input:
  !   nlterms - order of the transfer function
  !
  ! Output:
  !   norder  - number of Chebyshev coefficients per panel
  !   npan    - number of panels on each of the two halves of [-1,1]
  !
  use iso_fortran_env
  implicit none
  integer(int64), intent(in) :: nlterms
  integer(int64), intent(out) :: norder, npan

  norder = 8
  npan = 8*max(nlterms, 1_int64)

  return
end subroutine h3ddiagtabsize





subroutine h3ddiagtab(nlterms, zk, r, norder, npan, coefs)
  !
  ! Build a table of the diagonal transfer function
  !
  !   T(x) = sum_{l=0}^{nlterms} i^l (2l+1) h_l(k r) P_l(x)
  !
  ! for a translation of length r, as a function of x on [-1,1].
  ! T depends on the translation vector only through its length, so
  ! one table serves all translation vectors of that length and all
  ! points of the quadrature grid.
  !
  ! The table is a piecewise Chebyshev interpolant in the variable
  !
  !   w = sqrt((1-x)/2),   x >= 0     (iside = 1)
  !   w = sqrt((1+x)/2),   x <  0     (iside = 2)
  !
  ! with w in [0, 1/sqrt(2)] split into npan equal panels. In terms
  ! of the angle gamma, w = sin(gamma/2) on the first half, so T
  ! oscillates at a nearly uniform rate in w, and w is obtained from
  ! x with a square root.
  !
  ! Input:
  !   nlterms - order of the transfer function
  !   zk      - Helmholtz parameter
  !   r       - length of the translation vector
  !   norder  - number of Chebyshev coefficients per panel
  !   npan    - number of panels on each half (see h3ddiagtabsize)
  !
  ! Output:
  !   coefs   - coefs(k,ipan,iside), Chebyshev coefficients of T on
  !             panel ipan of half iside
  !
  use iso_fortran_env
  implicit none
  integer(int64), intent(in) :: nlterms, norder, npan
  real(real64), intent(in) :: r
  complex(real64), intent(in) :: zk
  complex(real64), intent(out) :: coefs(norder,npan,2)

  integer(int64) :: ipan, iside, j, k
  real(real64) :: pi, h, w, x
  real(real64), allocatable :: tnodes(:), umat(:,:)
  complex(real64), allocatable :: hfac(:), f(:)
  complex(real64), external :: h3ddiagtfun

  pi = 4*atan(1.0d0)

  allocate(hfac(0:max(nlterms,1_int64)), f(norder))
  allocate(tnodes(norder), umat(norder,norder))

  call h3ddiaghfac(nlterms, zk, r, hfac)

  ! Chebyshev nodes on [-1,1] and the values-to-coefficients matrix
  do j = 1,norder
    tnodes(j) = cos(pi*(j-0.5d0)/norder)
    do k = 1,norder
      umat(k,j) = 2*cos(pi*(k-1)*(j-0.5d0)/norder)/norder
    end do
    umat(1,j) = umat(1,j)/2
  end do

  h = sqrt(0.5d0)/npan
  do iside = 1,2
    do ipan = 1,npan
      do j = 1,norder
        w = (ipan-1)*h + h*(tnodes(j)+1)/2
        x = 1 - 2*w*w
        if (iside .eq. 2) x = -x
        f(j) = h3ddiagtfun(nlterms, hfac, x)
      end do
      do k = 1,norder
        coefs(k,ipan,iside) = 0
        do j = 1,norder
          coefs(k,ipan,iside) = coefs(k,ipan,iside) + umat(k,j)*f(j)
        end do
      end do
    end do
  end do

  return
end subroutine h3ddiagtab





subroutine h3ddiagtranstab(rvec, nth, nph, xnodes, cphi, sphi, &
    norder, npan, coefs, trans)
  !
  ! Tabulate the diagonal transfer function for the translation
  ! vector rvec on the quadrature grid, by interpolation from the
  ! table built by h3ddiagtab for a translation of length |rvec|.
  ! The cost is O(nth nph norder).
  !
  ! Input:
  !   rvec    - translation vector R = (target center - source center)
  !   nth,nph - grid size
  !   xnodes  - Gauss-Legendre nodes in cos(theta)
  !   cphi    - cos(phi_j), phi_j = 2 pi (j-1)/nph
  !   sphi    - sin(phi_j)
  !   norder  - number of Chebyshev coefficients per panel
  !   npan    - number of panels on each half
  !   coefs   - table from h3ddiagtab for r = |rvec|
  !
  ! Output:
  !   trans   - transfer function, trans(nph,nth)
  !
  use iso_fortran_env
  implicit none
  integer(int64), intent(in) :: nth, nph, norder, npan
  real(real64), intent(in) :: rvec(3), xnodes(nth), cphi(nph), sphi(nph)
  complex(real64), intent(in) :: coefs(norder,npan,2)
  complex(real64), intent(out) :: trans(nph,nth)

  integer(int64) :: i, j, k, ipan, iside
  real(real64) :: r, ct, st, a, b, c, x, w, s, u, hinv
  complex(real64) :: b0, b1, b2

  r = sqrt(rvec(1)**2 + rvec(2)**2 + rvec(3)**2)
  hinv = npan*sqrt(2.0d0)

  do i = 1,nth
    ct = xnodes(i)
    st = sqrt(1-ct**2)
    a = ct*rvec(3)/r
    b = st*rvec(1)/r
    c = st*rvec(2)/r
    do j = 1,nph
      x = a + b*cphi(j) + c*sphi(j)
      if (x .ge. 0) then
        iside = 1
        w = sqrt(max(1-x, 0.0d0)/2)
      else
        iside = 2
        w = sqrt(max(1+x, 0.0d0)/2)
      end if
      s = w*hinv
      ipan = min(int(s, int64), npan-1)
      u = 2*(s - ipan) - 1

      ! Clenshaw recurrence on panel ipan+1
      b1 = 0
      b2 = 0
      do k = norder,2,-1
        b0 = coefs(k,ipan+1,iside) + 2*u*b1 - b2
        b2 = b1
        b1 = b0
      end do
      trans(j,i) = coefs(1,ipan+1,iside) + u*b1 - b2
    end do
  end do

  return
end subroutine h3ddiagtranstab





subroutine h3ddiagtransall(nlterms, zk, bs, nth, nph, xnodes, trans)
  !
  ! Tabulate the diagonal transfer functions for all translation
  ! vectors bs*(ix,iy,iz), -3 <= ix,iy,iz <= 3, that occur in list 2.
  ! Translations between adjacent boxes are skipped.
  !
  ! The translation (ix,iy,iz) is stored in
  !   trans(:,:,h3ddiagind(ix,iy,iz))
  !
  ! The transfer function depends on the translation vector only
  ! through its length and the angle to the direction on the sphere.
  ! A one-dimensional table is built for each distinct length
  ! bs*sqrt(ix^2+iy^2+iz^2), and the transfer functions on the grid
  ! are interpolated from the tables.
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
  use iso_fortran_env
  implicit none
  integer(int64), intent(in) :: nlterms, nth, nph
  real(real64), intent(in) :: bs, xnodes(nth)
  complex(real64), intent(in) :: zk
  complex(real64), intent(out) :: trans(nph,nth,343)

  ! largest value of ix^2+iy^2+iz^2
  integer(int64), parameter :: n2max = 27
  integer(int64) :: ix, iy, iz, ind, n2, j, norder, npan
  integer(int64) :: ifuse(n2max)
  integer(int64), external :: h3ddiagind
  real(real64) :: rvec(3), r, pi, phi
  real(real64), allocatable :: cphi(:), sphi(:)
  complex(real64), allocatable :: coefs(:,:,:,:)

  pi = 4*atan(1.0d0)

  allocate(cphi(nph), sphi(nph))
  do j = 1,nph
    phi = 2*pi*(j-1)/nph
    cphi(j) = cos(phi)
    sphi(j) = sin(phi)
  end do

  ! distinct lengths of the translation vectors
  ifuse = 0
  do iz = -3,3
    do iy = -3,3
      do ix = -3,3
        if (max(abs(ix),abs(iy),abs(iz)) .ge. 2) then
          ifuse(ix*ix + iy*iy + iz*iz) = 1
        end if
      end do
    end do
  end do

  call h3ddiagtabsize(nlterms, norder, npan)
  allocate(coefs(norder,npan,2,n2max))

  !$omp parallel do default(shared) private(n2,r) schedule(dynamic)
  do n2 = 1,n2max
    if (ifuse(n2) .eq. 1) then
      r = bs*sqrt(1.0d0*n2)
      call h3ddiagtab(nlterms, zk, r, norder, npan, coefs(1,1,1,n2))
    end if
  end do
  !$omp end parallel do

  !$omp parallel do default(shared) private(ix,iy,iz,ind,n2,rvec) &
  !$omp collapse(3) schedule(dynamic)
  do iz = -3,3
    do iy = -3,3
      do ix = -3,3
        if (max(abs(ix),abs(iy),abs(iz)) .ge. 2) then
          ind = h3ddiagind(ix, iy, iz)
          n2 = ix*ix + iy*iy + iz*iz
          rvec(1) = bs*ix
          rvec(2) = bs*iy
          rvec(3) = bs*iz
          call h3ddiagtranstab(rvec, nth, nph, xnodes, cphi, sphi, &
              norder, npan, coefs(1,1,1,n2), trans(1,1,ind))
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
  use iso_fortran_env
  implicit none
  integer(int64) :: h3ddiagind
  integer(int64), intent(in) :: ix, iy, iz

  h3ddiagind = (ix+3) + 7*(iy+3) + 49*(iz+3) + 1

  return
end function h3ddiagind





subroutine h3ddiagapply(nd, nsig, trans, sigout, sigin)
  !
  ! Apply the diagonal translation operator:
  !
  !   sigin = sigin + trans .* sigout
  !
  use iso_fortran_env
  implicit none
  integer(int64), intent(in) :: nd, nsig
  complex(real64), intent(in) :: trans(nsig), sigout(nd,nsig)
  complex(real64), intent(inout) :: sigin(nd,nsig)

  integer(int64) :: i

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
  !               real(real64), as set up by mpalloc
  !   mnlist2   - leading dimension of list2
  !   nlist2    - number of boxes in list 2
  !   list2     - list 2 of each box
  !
  ! Input/output:
  !   rmlexp    - multipole and local expansions; the local
  !               expansions are incremented
  !
  use iso_fortran_env
  implicit none
  integer(int64), intent(in) :: nd, nterms, nlterms, ibstart, ibend
  integer(int64), intent(in) :: nboxes, ifpgh, ifpghtarg, mnlist2
  integer(int64), intent(in) :: isrcse(2,nboxes), itargse(2,nboxes)
  integer(int64), intent(in) :: iexpcse(2,nboxes), iaddr(2,nboxes)
  integer(int64), intent(in) :: nlist2(nboxes), list2(mnlist2,nboxes)
  real(real64), intent(in) :: rscale, bs, centers(3,nboxes)
  complex(real64), intent(in) :: zk
  complex(real64), intent(inout) :: rmlexp(*)

  integer(int64) :: nth, nph, nsig, nsrcbox, ibox, jbox, i, npts
  integer(int64) :: ix, iy, iz, ind, impole, ilocal
  integer(int64), external :: h3ddiagind
  integer(int64), allocatable :: isig(:)
  real(real64), allocatable :: xnodes(:), wts(:), ynms(:,:,:)
  real(real64), allocatable :: wsave(:)
  complex(real64), allocatable :: trans(:,:), sigout(:,:)
  complex(real64), allocatable :: sigin(:)

  call h3ddiaggridsize(nterms, nlterms, nth, nph)
  nsig = nth*nph

  allocate(xnodes(nth), wts(nth), ynms(0:nterms,0:nterms,nth))
  allocate(wsave(4*nph+20), trans(nsig,343))

  call h3ddiaginit(nterms, nth, nph, xnodes, wts, ynms, wsave)
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
          ynms, wsave, sigout(1,isig(ibox)))
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
      ix = nint((centers(1,ibox) - centers(1,jbox))/bs, int64)
      iy = nint((centers(2,ibox) - centers(2,jbox))/bs, int64)
      iz = nint((centers(3,ibox) - centers(3,jbox))/bs, int64)
      ind = h3ddiagind(ix, iy, iz)
      call h3ddiagapply(nd, nsig, trans(1,ind), sigout(1,isig(jbox)), &
          sigin)
      npts = npts + 1
    end do

    if (npts .gt. 0) then
      ilocal = (iaddr(2,ibox)+1)/2
      call h3dsig2loc(nd, nterms, rscale, sigin, nth, nph, wts, ynms, &
          wsave, rmlexp(ilocal))
    end if
  end do
  !$omp end do
  deallocate(sigin)
  !$omp end parallel

  return
end subroutine h3dmplocdiaglev
