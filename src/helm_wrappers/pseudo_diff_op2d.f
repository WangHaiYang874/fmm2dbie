      subroutine pd_helm2d(zk, nch, norders, ixys, 
     1     npts, srcvals, adjs, dmat, smat)
cf2py intent(in) zk, nch, norders, ixys, npts, srcvals, adjs
cf2py intent(out) dmat 
cf2py intent(out) smat

    !   This function returns the double layer and single layer potential
    !   matrices evaluate on the boundary of a 2D domain.

      implicit none
    ! input parameters
      integer *8, intent(in) :: nch, npts, norders(nch)
      integer *8, intent(in) :: ixys(nch+1), adjs(2,nch)
      real    *8, intent(in) :: srcvals(8,npts)
      complex*16, intent(in) :: zk
      
      ! output parameters
      complex*16, allocatable, intent(out):: dmat(:,:), smat(:,:)


    ! intermediate parameters
      integer *8 :: ndi, ndd, ndz, iptype(nch), opdims(2)
      integer *8 :: ifrobust, ier, ising, iquad
      real *8 , allocatable :: dpars(:)
      complex *16 , allocatable :: zpars(:)
      integer*8, allocatable :: ipars(:)
      procedure (), pointer :: fker
      external h2d_slp, h2d_dlp
      integer *8 :: i, itype, k
      real    *8 :: srccoefs(6,npts) 
      real    *8 :: ds(npts), alpha, beta
      real    *8, allocatable :: xs(:), ys(:), u(:,:), v(:,:), ws(:)

      alpha = 1.0d0
      beta  = 0.0d0
      itype = 2
      k = norders(1)

      allocate(xs(k), ys(k), u(k,k))
      allocate(v(k,k), ws(k))

      call legeexps(itype,k, xs, ys, u, v, ws)

      ds = sqrt((srcvals(3,:)**2 + srcvals(4,:)**2))
      ! computing srccoefs 
        
      do i = 1, nch
        call dgemm('n','t',6,k,k, alpha, 
     1    srcvals(1,ixys(i)),8,u,k,beta,srccoefs(1,ixys(i)),6)
      enddo
      
      do i = 1, nch
        iptype(i) = 1
      enddo
      ising = 0
      iquad = 0
      ifrobust = 0
      opdims(1) = 1
      opdims(2) = 1
      ndz = 1
      ndd = 0
      ndi = 0
      allocate(zpars(ndz), dpars(ndd), ipars(ndi))
      allocate(dmat(npts,npts), smat(npts,npts))
      zpars(1) = zk
      fker => h2d_slp
      call zchunk_matbuild_ggq(nch,norders,ixys,iptype,
     1       npts,srccoefs,srcvals,adjs,fker,
     2       ndd,dpars,ndz,zpars,ndi,ipars,opdims,
     3       ising, iquad, ifrobust, smat, ier)


      fker => h2d_dlp
      call zchunk_matbuild_ggq(nch,norders,ixys,iptype,
     1       npts,srccoefs,srcvals,adjs,fker,
     2       ndd,dpars,ndz,zpars,ndi,ipars,opdims,
     3       ising, iquad, ifrobust, dmat, ier)

      end subroutine pd_helm2d