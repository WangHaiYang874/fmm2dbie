! (C) Haiyang Wang @ Flatiron Institute 2023
! This file aims to solve the transmission boundary value problem
! for the 2D Helmholtz equation with one obstacle in free space.

! PDE:
!   (interior) u1 + k1^2 u1 = 0
!   (exterior) u0 + k0^2 u0 = 0

! The integral representation is
!   (interior) u1 = ep1^2 S_{k1}[\lambda] + ep1*D_{k1}[\rho]
!   (exterior) u0 = ep0^2 S_{k0}[\lambda] + ep0*D_{k0}[\rho]

! The boundary integral equation is
!
!   (ep0+ep1)/2 * rho
!       + (ep0*D_{k0} - ep1*D_{k1})[rho]
!       + (ep0^2*S_{k0} - ep1^2*S_{k1})[sigma]
!       = - uinc
!       = u0 - u1
!
!   (ep0+ep1)/2 * sigma
!       + (D'_{k0} - D'_{k1})[rho]
!       + (ep0*S'_{k0} - ep1*S'_{k1})[sigma]
!       = - uinc/ep0
!       = 1/ep0 u0' - 1/ep1 u1'
!
! The main task is to evaluate the LHS of the two equations above,
! then solve for the density rho and sigma iteratively using GMRES

! essentially there are following subroutines
!  - oncurvequad
!  - lpcomp
!  - solver

! zpars = (omega, ep0, mu0, ep1, mu1)

subroutine getoncurvequad_helm_comb_trans_2d(nch, norders, &
    ixys, iptype, npts, srccoefs, srcvals, adjs, eps, zpars, &
    iquadtype, nnz, row_ptr, col_ind, iquad, nquad, wnear)

! this subroutine generates the on curve quadrature
! for the representation
!
! u1 = ep1^2 S_{k1}[\lambda]+ep1*D_{k1}[\rho] (interior representation)
! u0 = ep0^2 S_{k0}[\lambda]+ep0*D_{k0}[\rho] (exterior representation)
!
! wnear must be of size 4*nquad as 4 different layer potentials are returned
!  * ep0^2 S_{k0} - ep1^2 D_{k1}
!  * ep0 D_{k0}   - ep1 D_{k1}
!  * ep0 S_{k0}'  - ep1 S_{k1}'
!  * D_{k0}'      - D_{k1}'
!
! zpars(5) = omega, ep0, mu0, ep1, mu1
!
! k0 = omega * sqrt(ep0*mu0)
! k1 = omega * sqrt(ep1*mu1)
! --------------------------------
! Input arguments:
! --------------------------------
!    - nch: integer
!     number of chunks
!    - norders: integer(nch)
!        order of discretization on each patch
!    - ixys: integer(nch+1)
!        starting location of data on patch i
!    - iptype: integer(nch)
!        type of patch
!        iptype = 1 -> chunk discretized with Gauss Legendre nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (6,npts)
!        basis coefficients of x,y,dxdt,dydt,dxdt2,dydt2
!        if(iptype.eq.1) then basis = legendre polynomials
!    - srcvals: real *8 (8,npts)
!        x,y,dxdt,dydt,dxdt2,dydt2,rnx,rny at the discretization nodes
!    - adjs: integer (2, nch)
!        adjacency information of each trunk.
!    - eps: real *8
!        precision requested
!    - zpars: complex *16 (5)
!        kernel parameters, as in
!        zpars(1:5) = omega, ep0, mu0, ep1, mu1
!    - iquadtype - integer
!        quadrature type
!        iquadtype = 1, using ggq for on curve quadrature.
!    - nnz: integer
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(ntarg+1)
!        row_ptr(i) is the pointer to col_ind array where list of
!        relevant source patches for target i start
!    - col_ind: integer (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(nnz+1)
!        location in wnear array where quadrature for col_ind(i)
!        starts
!    - nquad: integer
!        number of entries in wnear
!
!  Output
!  -------------------------------
!    - wnear: complex *16(4*nquad)
!        the desired near field quadrature

implicit none

! input
integer nch, norders(nch), npts, nquad, iquadtype, nnz, row_ptr(npts + 1), col_ind(nnz), iquad(nnz + 1)
integer ixys(nch + 1), iptype(nch), adjs(2, nch)
real*8 srccoefs(6, npts), srcvals(8, npts), eps
! srccoefs = x,y,dxdt,dydt,dxdt2,dydt2
! srcvals  = x,y,dxdt,dydt,dxdt2,dydt2, rnx, rny
complex*16 zpars(5)

! output
complex*16 wnear(4*nquad)

! intermediate variables
integer ipars, ndd, ndz, ier, i, ipv, ndi
real*8 dpars
complex*16 omega, ep0, ep1, k0, k1, mu0, mu1
complex*16 zpars1(6), zpars2(6), zpars3(6), zpars4(6)

procedure(), pointer :: fker
external h2d_transmission_dir, h2d_transmission_neu

! initialize the approriate kernel function
omega = zpars(1)
ep0 = zpars(2)
mu0 = zpars(3)
ep1 = zpars(4)
mu1 = zpars(5)

k0 = zpars(1)*sqrt(ep0*zpars(3))
k1 = zpars(1)*sqrt(ep1*zpars(5))

zpars1(1) = k0
zpars2(1) = k0
zpars3(1) = k0
zpars4(1) = k0

zpars1(2) = k1
zpars2(2) = k1
zpars3(2) = k1
zpars4(2) = k1

zpars1(3) = ep0**2
zpars1(4) = -ep1**2
zpars1(5) = 0
zpars1(6) = 0

zpars2(3) = 0
zpars2(4) = 0
zpars2(5) = ep0
zpars2(6) = -ep1

zpars3(3) = ep0
zpars3(4) = -ep1
zpars3(5) = 0
zpars3(6) = 0

zpars4(3) = 0
zpars4(4) = 0
zpars4(5) = 1
zpars4(6) = -1

ndd = 0
ndi = 0
ndz = 6

ipv = 0

fker => h2d_transmission_dir

call zgetoncurvequad_ggq2d(nch, norders, ixys, &
iptype, npts, srccoefs, srcvals, adjs, &
eps, ipv, fker, ndd, dpars, ndz, zpars1, &
ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, wnear, ier)

call zgetoncurvequad_ggq2d(nch, norders, ixys, &
iptype, npts, srccoefs, srcvals, adjs, &
eps, ipv, fker, ndd, dpars, ndz, zpars2, &
ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, wnear(nquad + 1), ier)

fker => h2d_transmission_neu
call zgetoncurvequad_ggq2d(nch, norders, ixys, &
iptype, npts, srccoefs, srcvals, adjs, &
eps, ipv, fker, ndd, dpars, ndz, zpars3, &
ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, wnear(2*nquad + 1), ier)

call zgetoncurvequad_ggq2d(nch, norders, ixys, &
iptype, npts, srccoefs, srcvals, adjs, &
eps, ipv, fker, ndd, dpars, ndz, zpars4, &
ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, wnear(3*nquad + 1), ier)

return
end subroutine getoncurvequad_helm_comb_trans_2d

subroutine getnearquad_helm_comb_trans_2d(nch, norders, &
 ixys, iptype, npts, srccoefs, srcvals, ndtarg, ntarg, &
 targs, ich_id, ts_targ, eps, zpars, &
 iquadtype, nnz, row_ptr, col_ind, iquad, nquad, wnear)

! this subroutine generates the near quadrature
! for the representation
!
! u1 = ep1^2 S_{k1}[\lambda]+ep1*D_{k1}[\rho] (interior representation)
! u0 = ep0^2 S_{k0}[\lambda]+ep0*D_{k0}[\rho] (exterior representation)
!
! wnear must be of size 4*nquad as 4 different layer potentials are returned
!  * ep0^2 S_{k0} - ep1^2 D_{k1}
!  * ep0 D_{k0}   - ep1 D_{k1}
!  * ep0 S_{k0}'  - ep1 S_{k1}'
!  * D_{k0}'      - D_{k1}'
!
! zpars(5) = omega, ep0, mu0, ep1, mu1
!
! k0 = omega * sqrt(ep0*mu0)
! k1 = omega * sqrt(ep1*mu1)
! --------------------------------
! Input arguments:
! --------------------------------
!    - nch: integer
!     number of chunks
!    - norders: integer(nch)
!        order of discretization on each patch
!    - ixys: integer(nch+1)
!        starting location of data on patch i
!    - iptype: integer(nch)
!        type of patch
!        iptype = 1 -> chunk discretized with Gauss Legendre nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (6,npts)
!        basis coefficients of x,y,dxdt,dydt,dxdt2,dydt2
!        if(iptype.eq.1) then basis = legendre polynomials
!    - srcvals: real *8 (8,npts)
!        x,y,dxdt,dydt,dxdt2,dydt2,rnx,rny at the discretization nodes
!    - ndtarg: integer
!        leading dimension of target array
!    - ntarg: integer
!        number of targets
!    - targs: real *8 (ndtarg, ntarg)
!        target information, the first two components must be
!        xy coordinate
!    - ich_id: integer(ntarg)
!        id of patch of target i, id = -1 if target is off-surface
!    - ts_targ: real *8(ntarg)
!        local t coordinate of chunk if on-surface
!    - eps: real *8
!        precision requested
!    - zpars: complex *16 (5)
!        kernel parameters, as in
!        zpars(1:5) = omega, ep0, mu0, ep1, mu1
!    - iquadtype - integer
!        quadrature type
!        iquadtype = 1, using ggq for on curve quadrature.
!    - nnz: integer
!        number of source patch-> target interactions in the near field
!    - row_ptr: integer(ntarg+1)
!        row_ptr(i) is the pointer to col_ind array where list of
!        relevant source patches for target i start
!    - col_ind: integer (nnz)
!        list of source patches relevant for all targets, sorted
!        by the target number
!    - iquad: integer(nnz+1)
!        location in wnear array where quadrature for col_ind(i)
!        starts
!    - nquad: integer
!        number of entries in wnear
!
!  Output
!  -------------------------------
!    - wnear: complex *16(4*nquad)
!        the desired near field quadrature

implicit none

! input
integer nch, norders(nch), npts, nquad, iquadtype, nnz, row_ptr(npts + 1), col_ind(nnz), iquad(nnz + 1)
integer ndtarg, ntarg
integer ixys(nch + 1), iptype(nch), adjs(2, nch), ich_id(ntarg)
real*8 srccoefs(6, npts), srcvals(8, npts), targs(ndtarg, ntarg), ts_targ(ntarg), eps
! srccoefs = x,y,dxdt,dydt,dxdt2,dydt2
! srcvals  = x,y,dxdt,dydt,dxdt2,dydt2, rnx, rny
complex*16 zpars(5)

! output
complex*16 wnear(4*nquad)

! intermediate variables

integer ipars
integer ndd, ndz, ndi
real*8 dpars

integer i, j, ipv

complex*16 omega, ep0, ep1, k0, k1, mu0, mu1
complex*16 zpars1(6), zpars2(6), zpars3(6), zpars4(6)

procedure(), pointer :: fker
external h2d_transmission_dir, h2d_transmission_neu

! initialize the approriate kernel function
omega = zpars(1)
ep0 = zpars(2)
mu0 = zpars(3)
ep1 = zpars(4)
mu1 = zpars(5)

k0 = zpars(1)*sqrt(ep0*zpars(3))
k1 = zpars(1)*sqrt(ep1*zpars(5))

zpars1(1) = k0
zpars2(1) = k0
zpars3(1) = k0
zpars4(1) = k0

zpars1(2) = k1
zpars2(2) = k1
zpars3(2) = k1
zpars4(2) = k1

zpars1(3) = ep0**2
zpars1(4) = -ep1**2
zpars1(5) = 0
zpars1(6) = 0

zpars2(3) = 0
zpars2(4) = 0
zpars2(5) = ep0
zpars2(6) = -ep1

zpars3(3) = ep0
zpars3(4) = -ep1
zpars3(5) = 0
zpars3(6) = 0

zpars4(3) = 0
zpars4(4) = 0
zpars4(5) = 1
zpars4(6) = -1

ndz = 6
ndd = 0
ndi = 0
ipv = 0

fker => h2d_transmission_dir

call zgetnearquad_adap_guru2d(nch, norders, ixys, &
iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
ich_id, ts_targ, eps, ipv, fker, ndd, dpars, ndz, zpars1, &
ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, wnear)

call zgetnearquad_adap_guru2d(nch, norders, ixys, &
iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
ich_id, ts_targ, eps, ipv, fker, ndd, dpars, ndz, zpars2, &
ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, wnear(nquad + 1))

fker => h2d_transmission_neu

call zgetnearquad_adap_guru2d(nch, norders, ixys, &
iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
ich_id, ts_targ, eps, ipv, fker, ndd, dpars, ndz, zpars3, &
ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, wnear(2*nquad + 1))

call zgetnearquad_adap_guru2d(nch, norders, ixys, &
iptype, npts, srccoefs, srcvals, ndtarg, ntarg, targs, &
ich_id, ts_targ, eps, ipv, fker, ndd, dpars, ndz, zpars4, &
ndi, ipars, nnz, row_ptr, col_ind, iquad, nquad, wnear(3*nquad + 1))

return
end subroutine getnearquad_helm_comb_trans_2d

subroutine lpcomp_helm_comb_trans_add_sub_2d(nch, norders, ixys, &
    iptype, npts, srccoefs, srcvals, eps, zpars, &
    nnz, row_ptr, col_ind, iquad, nquad, wnear, sigma, novers, &
    nptso, ixyso, srcover, whtsover, pot)
!
!f2py intent(in) nch,norders,ixys,iptype,npts,srccoefs,srcvals
!f2py intent(in) ndtarg,ntarg,targs,eps,zpars
!f2py intent(in) sigma,nnz,row_ptr,col_ind,iquad,nquad,wnear,novers
!f2py intent(in) nptso,ixyso,srcover,whtsover
!f2py intent(out) pot
!
!
! ------------------------------
!   This subroutine evaluates the layer potential on surface for the representation
!
!    u1 = ep1^2 S_{k1}[\lambda]+ep1*D_{k1}[\rho] (interior representation)
!    u0 = ep0^2 S_{k0}[\lambda]+ep0*D_{k0}[\rho] (exterior representation)
!
!   and returns quantities related to u0-u1 , and 1/ep0 u0' - 1/ep1 u1'
!   on surface
!
!   On imposing the boundary condition, we get the following
!   sets of operators
!
!    u0-u1 = (ep0+ep1)/2 \rho + (ep0 D_{k0} - \ep1 D_{k1})[\rho] +
!            (ep0^2 S_{k1} - ep1^2 S_{k0})[\lambda]
!
!    1/ep0 u0' - 1/ep1 u1' = -(ep0 + ep1)/2 \lambda +
!            (D_{k0}'-D_{k1}')[\rho] + (ep0 S_{k0}' - ep1 S_{k1}')[\lambda]
!
!   where the near field is precomputed and stored
!   in the row sparse compressed format.
!
!   Note: For targets on the boundary, this routine only computes
!   the principal value part, the identity term corresponding to the jump
!   in the layer potential is not included in the layer potential.
!
!
!   Input arguments:
!
!     - nch: integer
!         number of patches
!     - norders: integer(nch)
!         order of discretization on each patch
!     - ixys: integer(nch+1)
!         ixys(i) denotes the starting location in srccoefs,
!         and srcvals array where information for patch i begins
!     - iptype: integer(nch)
!         type of patch
!     - npts: integer
!         total number of discretization points on the boundary
!     - srccoefs: double precision (6,npts)
!         koornwinder expansion coefficients of x, $\partial_{u} x$,
!         and $\partial_{v} x$.
!     - srcvals: double precision (8,npts)
!         x, $\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
!         discretization nodes
!     - eps: double precision
!         precision requested
!     - zpars: double complex (5)
!         the kernel parameter
!         zpars = (omega, ep0, mu0, ep1, mu1)
!     - nnz: integer
!         number of source patch-> target interactions in the near field
!     - row_ptr: integer(ntarg+1)
!         row_ptr(i) is the pointer to col_ind array where list of
!         relevant source patches for target i start
!     - col_ind: integer (nnz)
!         list of source patches relevant for all targets, sorted
!         by the target number
!     - iquad: integer(nnz+1)
!         location in wnear array where quadrature for col_ind(i)
!         starts
!     - nquad: integer
!         number of entries in wnear
!     - wnear: complex *16(4*nquad)
!         the desired near field quadrature
!     - sigma: double complex(2*npts)
!         density for layer potential
!         the first half is for lambda, the second half is for rho.
!     - novers: integer(nch)
!         order of discretization for oversampled sources and
!         density
!     - nptso: integer
!         number of oversampled points
!     - ixyso: integer(nch+1)
!         ixyso(i) denotes the starting location in srcover,
!         whtsover array, and other functions sampled at oversampled
!         nodes where information for patch i begins
!     - srcover: double precision (8,npts)
!         x, $\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
!         oversampled nodes
!     - whtsover: double precision (npts)
!         quadrature weights for integrating smooth functions sampled
!         at the oversampled nodes
!
!   Output arguments
!     - pot: double complex(2*ntarg)
!         layer potential evaluated at the target points
!
! -----------------------------------
!
!
implicit none
integer nch, npts
integer ndtarg, ntarg
integer norders(nch), ixys(nch + 1)
integer ixyso(nch + 1), iptype(nch)
real*8 srccoefs(6, npts), srcvals(8, npts), eps
complex*16 zpars(5)
complex*16 zk0, zk1
integer nnz, row_ptr(npts + 1), col_ind(nnz), nquad
integer iquad(nnz + 1)
complex*16 sigma(2*npts)
complex*16 wnear(4*nquad)

integer, intent(in) :: novers(nch)
integer, intent(in) :: nptso
real*8, intent(in) :: srcover(8, nptso), whtsover(nptso)
complex*16, intent(out) :: pot(2*npts)

integer norder, npols, nover, npolso
complex*16, allocatable :: potsort(:)

real*8, allocatable :: sources(:, :), srctmp(:, :)

!❗this might cause some problems?
complex*16, allocatable :: charges0(:), dipstr0(:), sigmaover(:)
complex*16, allocatable :: charges1(:), dipstr1(:)
real*8, allocatable :: dipvec(:, :)
integer ns, nt
complex*16 alpha, beta
integer ifcharge, ifdipole
integer ifpgh, ifpghtarg
complex*16 tmp(10), val, E(4)

integer i, j, jpatch, jquadstart, jstart
complex*16 zdotu, pottmp, gradtmp(2)
complex*16, allocatable :: pot_aux(:), grad_aux(:, :)
complex*16 ep0, ep1, ep0sq, ep1sq, ep0inv, ep1inv
real*8 radexp, epsfmm

! ❓
integer ipars
real*8 dpars, timeinfo(10), t1, t2, omp_get_wtime

real*8, allocatable :: radsrc(:)
real*8, allocatable :: srctmp2(:, :)
complex*16, allocatable :: ctmp0(:), dtmp0(:)
complex*16, allocatable :: ctmp1(:), dtmp1(:)

real*8, allocatable :: dipvec2(:, :)

real*8 thresh, ra
real*8 rr, rmin
real*8 over4pi
integer nss, ii, l, npover, ier
complex*16 ima, ztmp

integer nd, ntarg0, nmax

real*8 ttot, done, pi

integer iper
real*8 xmin, xmax, ymin, ymax, zmin, zmax, sizey, sizez, boxsize

integer ifaddsub

integer ntj

data ima/(0.0d0, 1.0d0)/
data over4pi/0.07957747154594767d0/

parameter(nd=1, ntarg0=1)

ns = nptso
ntarg = npts
done = 1
pi = atan(done)*4

! compute the interior and exterior wave number
zk0 = zpars(1)*sqrt(zpars(2)*zpars(3))
zk1 = zpars(1)*sqrt(zpars(4)*zpars(5))
ep0 = zpars(2)
ep1 = zpars(4)
ep0sq = ep0*ep0
ep1sq = ep1*ep1
ep0inv = 1.0d0/ep0
ep1inv = 1.0d0/ep1

ifpgh = 0
ifpghtarg = 2

allocate (sources(2, ns), srctmp(2, npts))
allocate (charges0(ns), dipstr0(ns), dipvec(2, ns))
allocate (charges1(ns), dipstr1(ns))
allocate (sigmaover(2*ns))
allocate (pot_aux(npts), grad_aux(2, npts))

!
!     estimate max number of sources in near field of
!     any target
!
nmax = 0
call get_near_corr_max2d(npts, row_ptr, nnz, col_ind, nch, ixyso, nmax)
allocate (srctmp2(2, nmax), ctmp0(nmax), dtmp0(nmax))
allocate (ctmp1(nmax), dtmp1(nmax))

!
!        oversample density
!

call oversample_fun_curv2d(2, nch, norders, ixys, iptype, &
npts, sigma, novers, ixyso, ns, sigmaover)
call oversample_fun_curv2d(2, nch, norders, ixys, iptype, &
npts, sigma(npts + 1), novers, ixyso, ns, sigmaover(npts + 1))

ra = 0

! ❗ the following few lines might be deleted
! !
! !        set relevatn parameters for the fmm
! !
!    alpha = zpars(2)
!    beta = zpars(3)

! extract source and target info

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
do i = 1, ns

sources(1, i) = srcover(1, i)
sources(2, i) = srcover(2, i)

dipvec(1, i) = srcover(7, i)
dipvec(2, i) = srcover(8, i)

!❗ Maybe I should add over4pi in the following equations?

charges0(i) = sigmaover(ns + i)*whtsover(i)*ep0sq
charges1(i) = sigmaover(ns + i)*whtsover(i)*ep1sq

dipstr0(i) = sigmaover(i)*whtsover(i)*ep0
dipstr1(i) = sigmaover(i)*whtsover(i)*ep1

end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
do i = 1, ntarg

srctmp(1, i) = srcvals(1, i)
srctmp(2, i) = srcvals(2, i)
pot_aux(i) = 0
grad_aux(1, i) = 0
grad_aux(2, i) = 0
pot(i) = 0
pot(npts + i) = 0

end do
!$OMP END PARALLEL DO

!  Compute
!    ep0^2*S_{k0}[\lambda] + ep0*D_{k0}[\rho],
!    ep0^2*S_{k0}'[\lambda] + ep0*D_{k0}'[\rho]

ifcharge = 1
ifdipole = 1

! if(alpha.eq.0) ifcharge = 0
! if(beta.eq.0) ifdipole = 0

!        call the fmm
!
call cpu_time(t1)
!$  t1 = omp_get_wtime()

! ❗ the following line is for helm_comb_dir_2d, might be deleted
! call hfmm2d(nd,eps,zpars(1),ns,sources,ifcharge,charges,&
!    ifdipole,dipstr,dipvec,iper,ifpgh,tmp,tmp,tmp,ntarg,&
!    targvals,ifpghtarg,pot,tmp,tmp,ier)
! ❗the following line is from helm_comb_trans, might be deleted

! call hfmm3d_t_cd_g(eps,zk0,ns,sources,charges0,dipvec0,npts, &
! srctmp,pot_aux,grad_aux,ier)

call hfmm2d_t_cd_g(eps, zk0, ns, sources, charges0, &
dipstr0, dipvec, ntarg, srctmp, pot_aux, grad_aux, ier)

call cpu_time(t2)
!$  t2 = omp_get_wtime()

timeinfo(1) = t2 - t1

!
!   Add ep0^2*S_{k0}[\lambda]+ep0 D_{k0}[\rho] to pot and
!     1/ep0*(ep0^2*S_{k0}' + ep0 D_{k0}[\rho] to pot(npts+:)
!
!$OMP PARALLEL DO DEFAULT(SHARED)
do i = 1, npts
pot(i) = pot_aux(i)
pot(npts + i) = ep0inv*(grad_aux(1, i)*srcvals(7, i) + grad_aux(2, i)*srcvals(8, i))
end do
!$OMP END PARALLEL DO

!
!  Compute ep1^2*S_{k1}[\lambda] + ep1*D_{k1}[\rho],
!  ep1^2*S_{k1}'[\lambda] + ep1*D_{k1}'[\rho]
!
call hfmm2d_t_cd_g(eps, zk1, ns, sources, charges1, &
dipstr1, dipvec, npts, srctmp, pot_aux, grad_aux, ier)
!
!   subtract ep1^2*S_{k1}[\lambda]+ep1 D_{k1}[\rho] from pot and
!     1/ep1*(ep1^2*S_{k1}' + ep0 D_{k1}[\rho] from pot(npts+:)
!
!$OMP PARALLEL DO DEFAULT(SHARED)
do i = 1, npts
pot(i) = pot(i) - pot_aux(i)
pot(npts + i) = pot(npts + i) - ep1inv*(grad_aux(1, i)*srcvals(7, i) + grad_aux(2, i)*srcvals(8, i))
end do
!$OMP END PARALLEL DO

!
!         compute threshold for ignoring local computation
!
call get_fmm2d_thresh(2, ns, sources, 2, ntarg, srcvals, thresh)

!
!        add in precomputed quadrature

call cpu_time(t1)
!$  t1 = omp_get_wtime()

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart) &
!$OMP PRIVATE(jstart,npols,l)
do i = 1, npts
do j = row_ptr(i), row_ptr(i + 1) - 1
jpatch = col_ind(j)
npols = ixys(jpatch + 1) - ixys(jpatch)
jquadstart = iquad(j)
jstart = ixys(jpatch)
do l = 1, npols
pot(i) = pot(i) + &
wnear(jquadstart + l - 1)*sigma(jstart + l - 1 + npts)
pot(i) = pot(i) + &
wnear(nquad + jquadstart + l - 1)*sigma(jstart + l - 1)
pot(i + npts) = pot(i + npts) + &
wnear(2*nquad + jquadstart + l - 1)*sigma(jstart + l - 1 + npts)
pot(i + npts) = pot(i + npts) + &
wnear(3*nquad + jquadstart + l - 1)*sigma(jstart + l - 1)
end do
end do
end do
!$OMP END PARALLEL DO

! Remove near contribution of the FMM
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2) &
!$OMP PRIVATE(ctmp0,dtmp0,ctmp1,dtmp1,dipvec,l,jstart,nss,pottmp,gradtmp)
do i = 1, npts
nss = 0
do j = row_ptr(i), row_ptr(i + 1) - 1
jpatch = col_ind(j)
do l = ixyso(jpatch), ixyso(jpatch + 1) - 1
nss = nss + 1
srctmp2(1, nss) = srcover(1, l)
srctmp2(2, nss) = srcover(2, l)
srctmp2(3, nss) = srcover(3, l)

ctmp0(nss) = charges0(l)
ctmp1(nss) = charges1(l)

dtmp0(nss) = dipstr0(l)
dtmp1(nss) = dipstr1(l)
end do
end do

pottmp = 0
gradtmp(1) = 0
gradtmp(2) = 0

call h2d_directcdg(nd, zk0, srctmp2, nss, ctmp0, dtmp0, dipvec, &
srctmp(1, i), ntarg0, pottmp, gradtmp, thresh)

pot(i) = pot(i) - pottmp
pot(npts + i) = pot(npts + i) - &
(gradtmp(1)*srcvals(7, i) + gradtmp(2)*srcvals(8, i))*ep0inv

pottmp = 0
gradtmp(1) = 0
gradtmp(2) = 0

call h2d_directcdg(nd, zk1, srctmp2, nss, ctmp1, dtmp1, dipvec, &
srctmp(1, i), ntarg0, pottmp, gradtmp, thresh)
pot(i) = pot(i) + pottmp
pot(npts + i) = pot(npts + i) + &
(gradtmp(1)*srcvals(7, i) + gradtmp(2)*srcvals(8, i))*ep1inv
!
!  flip sign of pot(npts+i)
!
pot(npts + i) = -pot(npts + i)
end do
!$OMP END PARALLEL DO
call cpu_time(t2)

!$  t2 = omp_get_wtime()

timeinfo(2) = t2 - t1

!       call prin2('lpcomp timeinfo=*',timeinfo,2)
!       call prinf('nfmm *',ns,1)
!       call prinf('nquad *',nquad,1)
ttot = timeinfo(1) + timeinfo(2)
!c      call prin2('time in lpcomp=*',ttot,1)

return
end

! subroutine helm_comb_dir_solver_2d(nch,norders,ixys,
!    1    iptype,npts,srccoefs,srcvals,adjs,eps,zpars,numit,ifinout,
!    2    rhs,eps_gmres,niter,errs,rres,soln)

! subroutine helm_comb_trans_solver(npatches,norders,ixyzs, &
!    iptype,npts,srccoefs,srcvals,eps,zpars,numit, &
!    rhs,eps_gmres,niter,errs,rres,soln)

subroutine helm_comb_trans_solver_2d(nch, norders, ixys, &
iptype, npts, srccoefs, srcvals, adjs, eps, zpars, numit, &
rhs, eps_gmres, niter, errs, rres, soln)
!
!  Solve the Helmholtz transmission problem using the combined
!  field integral equation
!
!     u1 = ep1^2 S_{k1}[\lambda]+ep1*D_{k1}[\rho] (interior representation)
!     u0 = ep0^2 S_{k0}[\lambda]+ep0*D_{k0}[\rho] (exterior representation)
!
!  Input arguments:
!    - nch: integer
!        number of chunks
!    - norders: integer(nch)
!        order of discretization on each patch
!    - ixys: integer(nch+1)
!        starting location of data on patch i
!    - iptype: integer(nch)
!        type of patch
!        iptype = 1 -> chunk discretized with Gauss Legendre nodes
!    - npts: integer
!        total number of discretization points on the boundary
!    - srccoefs: real *8 (6,npts)
!        basis coefficients of x,y,dxdt,dydt,dxdt2,dydt2
!        if(iptype.eq.1) then basis = legendre polynomials
!    - srcvals: real *8 (8,npts)
!        x,y,dxdt,dydt,dxdt2,dydt2,rnx,rny at the discretization nodes
!    - eps: real *8
!        precision requested
!    - zpars: complex *16 (5)
!        kernel parameters
!        zpars = (omega, ep0, mu0, ep1, mu1)
!    - numit: integer
!        max number of gmres iterations
!    - rhs: complex *16(2*npts)
!        right hand side
!        rhs(1:npts) = f
!        rhs(npts+1:2*npts) = g
!    - eps_gmres: real *8
!        gmres tolerance requested
!
!  Output arguments:
!    - niter: integer
!        number of gmres iterations required for relative residual
!    - errs: real *8 (1:niter)
!        relative residual as a function of iteration number
!    - rres: real *8
!        relative residual for computed solution
!    - soln: complex *16(2*npts)
!        density which solves the dirichlet problem
!-----------------------------------
!

implicit none
integer nch, npts
integer norders(nch), ixys(nch + 1), iptype(nch)
real*8 srccoefs(6, npts), srcvals(8, npts)
integer adjs(2, nch)
integer numit, niter
real*8 eps, eps_gmres
real*8 errs(numit + 1), rres

complex*16 zpars(5)
complex*16 rhs(2*npts), soln(2*npts)

real*8, allocatable :: targs(:, :)
integer, allocatable :: ich_id(:)
real*8, allocatable :: ts_targ(:)

integer ndtarg, ntarg

integer nover, npolso, nptso
integer nnz, nquad
integer, allocatable :: row_ptr(:), col_ind(:), iquad(:)
complex*16, allocatable :: wnear(:)

real*8 rho
integer :: ier, ising, npolyfac
real*8, allocatable :: srcover(:, :), wover(:), rects(:, :, :)
integer, allocatable :: ixyso(:), novers(:)

real*8, allocatable :: cms(:, :), rads(:), rad_near(:)
integer i, j, jpatch, jquadstart, jstart

integer ipars(2)
real*8 dpars, timeinfo(10), t1, t2, omp_get_wtime
complex*16 zk0, zk1, zkuse

real*8 ttot, done, pi
real*8 rfac, rfac0, t0
integer iptype_avg, norder_avg
integer iquadtype, npts_over
integer n_var

! gmres variables
complex*16 zid, ztmp
real*8 rb, wnrm2
integer it, iind, it1, k, l
real*8 rmyerr
complex*16 temp
complex*16, allocatable :: vmat(:, :), hmat(:, :)
complex*16, allocatable :: cs(:), sn(:)
complex*16, allocatable :: svec(:), yvec(:), wtmp(:)

complex*16 ima
ima = (0.0d0, 1.0d0)
done = 1
pi = atan(done)*4

n_var = npts*2

allocate (vmat(n_var, numit + 1), hmat(numit, numit))
allocate (cs(numit), sn(numit))
allocate (wtmp(n_var), svec(numit + 1), yvec(numit + 1))

ndtarg = 8
ntarg = npts
allocate (targs(ndtarg, ntarg))

!$OMP PARALLEL DO DEFAULT(SHARED)
do i = 1, ntarg
targs(:, i) = srcvals(:, i)
! ipatch_id(i) = -1
! uvs_targ(1,i) = 0
! uvs_targ(2,i) = 0
end do
!$OMP END PARALLEL DO

!!!!! near field and oversampling
!!!!! ❕ there is no far field measurement!

rho = 2.0d0
npolyfac = 2

allocate (novers(nch), ixyso(nch + 1))

call ellipse_nearfield2d_getnovers(eps, rho, npolyfac, &
nch, norders, ising, novers, ier)

ixyso(1) = 1

do i = 1, nch
ixyso(i + 1) = ixyso(i) + novers(i)
end do

npts_over = ixyso(nch + 1) - 1

allocate (rects(2, 4, nch))
call ellipse_nearfield2d_definerects(nch, norders, &
ixys, iptype, npts, srccoefs, srcvals, rho, rects)

call findinrectangle_mem(nch, rects, npts, ndtarg, targs, nnz, ier)

allocate (row_ptr(npts + 1), col_ind(nnz))

call findinrectangle(nch, rects, npts, ndtarg, targs, row_ptr, nnz, col_ind, ier)

allocate (iquad(nnz + 1))

call get_iquad_rsc2d(nch, ixys, npts, nnz, row_ptr, col_ind, iquad)

!!!!! oversample geometry and get over sampled weights for FMM
allocate (srcover(8, npts_over), wover(npts_over))

call oversample_geom2d(nch, norders, ixys, iptype, npts, &
srccoefs, srcvals, novers, ixyso, npts_over, srcover)

call get_qwts2d(nch, novers, ixyso, iptype, npts_over, srcover, wover)

!!!!! compute the near quadrature correction

nquad = iquad(nnz + 1) - 1
allocate (wnear(4*nquad))

!$OMP PARALLEL DO DEFAULT(SHARED)
do i = 1, 4*nquad
wnear(i) = 0
end do
!$OMP END PARALLEL DO

iquadtype = 1
call cpu_time(t0)
!$  t1 = omp_get_wtime()
call getoncurvequad_helm_comb_trans_2d(nch, norders, ixys, iptype, &
  npts, srccoefs, srcvals, adjs, eps, zpars, iquadtype, nnz, &
  row_ptr, col_ind, iquad, nquad, wnear)

call cpu_time(t1)

call prin2('time generating near quad *', t1 - t0, 1)
print *, "done generating near quadrature, now starting gmres"

!
!
!     start gmres code here
!
!     NOTE: matrix equation should be of the form (z*I + K)x = y
!       the identity scaling (z) is defined via zid below,
!       and K represents the action of the principal value
!       part of the matvec
!

zid = (zpars(2) + zpars(4))/2.0d0

niter = 0

! compute the norm of rhs and initialize v

rb = 0

do i = 1, numit
cs(i) = 0
sn(i) = 0
end do

!$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rb)
do i = 1, n_var
rb = rb + abs(rhs(i))**2
end do
!$OMP END PARALLEL DO
rb = sqrt(rb)

!$OMP PARALLEL DO DEFAULT(SHARED)
do i = 1, npts
vmat(i, 1) = rhs(i)/rb
vmat(i + npts, 1) = -rhs(i + npts)/rb
end do
!$OMP END PARALLEL DO
svec(1) = rb

do it = 1, numit
it1 = it + 1

call lpcomp_helm_comb_trans_add_sub_2d(nch, norders, ixys, &
      iptype, npts, srccoefs, srcvals, eps, zpars, nnz, row_ptr, &
      col_ind, iquad, nquad, wnear, vmat(1, it), novers, &
      npts_over, ixyso, srcover, wover, wtmp)
do k = 1, it
ztmp = 0
!$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:ztmp)
do j = 1, n_var
ztmp = ztmp + wtmp(j)*conjg(vmat(j, k))
end do
!$OMP END PARALLEL DO
hmat(k, it) = ztmp

!$OMP PARALLEL DO DEFAULT(SHARED)
do j = 1, n_var
wtmp(j) = wtmp(j) - hmat(k, it)*vmat(j, k)
end do
!$OMP END PARALLEL DO
end do

hmat(it, it) = hmat(it, it) + zid
wnrm2 = 0
!$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:wnrm2)
do j = 1, n_var
wnrm2 = wnrm2 + abs(wtmp(j))**2
end do
!$OMP END PARALLEL DO
wnrm2 = sqrt(wnrm2)

!$OMP PARALLEL DO DEFAULT(SHARED)
do j = 1, n_var
vmat(j, it1) = wtmp(j)/wnrm2
end do
!$OMP END PARALLEL DO

do k = 1, it - 1
temp = cs(k)*hmat(k, it) + conjg(sn(k))*hmat(k + 1, it)
hmat(k + 1, it) = -sn(k)*hmat(k, it) + cs(k)*hmat(k + 1, it)
hmat(k, it) = temp
end do

ztmp = wnrm2

call zrotmat_gmres2d(hmat(it, it), ztmp, cs(it), sn(it))

hmat(it, it) = cs(it)*hmat(it, it) + conjg(sn(it))*wnrm2
svec(it1) = -sn(it)*svec(it)
svec(it) = cs(it)*svec(it)
rmyerr = abs(svec(it1))/rb
errs(it) = rmyerr
print *, "iter=", it, errs(it)

if (rmyerr .le. eps_gmres .or. it .eq. numit) then

!
!            solve the linear system corresponding to
!            upper triangular part of hmat to obtain yvec
!
!            y = triu(H(1:it,1:it))\s(1:it);
!
do j = 1, it
iind = it - j + 1
yvec(iind) = svec(iind)
do l = iind + 1, it
yvec(iind) = yvec(iind) - hmat(iind, l)*yvec(l)
end do
yvec(iind) = yvec(iind)/hmat(iind, iind)
end do

!
!          estimate x
!
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
do j = 1, n_var
soln(j) = 0
do i = 1, it
soln(j) = soln(j) + yvec(i)*vmat(j, i)
end do
end do
!$OMP END PARALLEL DO

rres = 0
!$OMP PARALLEL DO DEFAULT(SHARED)
do i = 1, n_var
wtmp(i) = 0
end do
!$OMP END PARALLEL DO
!
!        NOTE:
!        replace this routine by appropriate layer potential
!        evaluation routine
!
call lpcomp_helm_comb_trans_add_sub_2d(nch, norders, ixys, &
      iptype, npts, srccoefs, srcvals, eps, zpars, nnz, row_ptr, &
      col_ind, iquad, nquad, wnear, vmat(1, it), novers, &
      npts_over, ixyso, srcover, wover, wtmp)


!$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rres)
do i = 1, npts
rres = rres + abs(zid*soln(i) + wtmp(i) - rhs(i))**2
rres = rres + abs(zid*soln(npts + i) + &
wtmp(npts + i) + rhs(npts + i))**2
end do
!$OMP END PARALLEL DO
rres = sqrt(rres)/rb
niter = it
return

end if
end do

end
