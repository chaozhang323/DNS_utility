module MTSTimeData
  implicit none
  real, dimension(:), pointer :: tstime
  logical, private :: IsTimeRead = .false.
  integer  ::  ntpoint
  character(200) :: flowpath  ! flowpath of the timeseries output by DNS
contains
    integer function GetPointRead()
      GetPointRead = ntpoint
    end function GetPointRead

    real(8) function Getdts()
      if(IsTimeRead) then
         Getdts = tstime(2) - tstime(1)
      else
         print *, 'TSTIME needs to be read first ... STOP'
         stop
      endif
    end function Getdts

    subroutine ReadTSTime(filepath,ntime,ns)
      character(*),intent(in) :: filepath
      integer, intent(in) :: ntime
      integer, intent(in), optional :: ns   ! # of data skipped
      integer :: nt, i, j, n, nskip
      real(8) :: tmp
      character(200) :: fn

      ntpoint = ntime
      flowpath = filepath
      print *, 'ntpoint=',ntpoint
      allocate(tstime(ntpoint))
      nskip = 0
      if (present(ns)) nskip=ns
      fn = trim(flowpath)//'series_time.dat'
      open(21, file=fn, form='unformatted', status='old')
      do nt=1,nskip
        read(21,end=101) tmp
      end do
      do nt = 1, ntpoint
        read(21,end=101) tstime(nt)
      end do
101   close(21)
      print*,'Actaully read time point = ', ntpoint
      print*,'tstart=',tstime(1),',tend=',tstime(ntpoint)
      do nt = 2, ntpoint
         if(tstime(nt)-tstime(nt-1).le.0.d0) then
            print *, 'nt-1=',nt-1, 'tstime(nt-1)=',tstime(nt-1)
            print *, 'nt  =',nt  , 'tstime(nt)  =',tstime(nt)
            print *, 'T+dt should be great than T. Time series corrupted. STOP'
            stop
         endif
      enddo
      IsTimeRead = .true.
    end subroutine ReadTSTime

    subroutine ReadTsData(fn,nvarin,nx,ny,nz,timeseries,iasc,ns)
      implicit none
      character(*), intent(in) :: fn
      integer, intent(in) :: nvarin, nx, ny, nz
      real(8), intent(out) :: timeseries(nvarin,nx,ny,nz,ntpoint)
      integer, intent(in), optional :: iasc   ! data type, whether ascii or not 
      integer, intent(in), optional :: ns   ! # of data skipped
      integer :: nt, i, j, k, n, nskip, iascii
      real(8) :: tmp
      nskip = 0
      if (present(iasc)) then
         iascii = iasc
      else
         iascii = 0
      endif
      if (present(ns)) then
          nskip=ns
      else
          nskip = 0
      endif
      if(iascii.eq.1) then
         open(21, file=fn, status='old')
         do nt=1,nskip
            read(21,*) ((((tmp,n=1,nvarin),i=1,nx),j=1,ny),k=1,nz)
         end do
         do nt = 1, ntpoint
            read(21,*) ((((timeseries(n,i,j,k,nt),n=1,nvarin),i=1,nx),j=1,ny),k=1,nz)
         end do
         close(21)
      else
         open(21, file=fn, form='unformatted', status='old')
         do nt=1,nskip
            read(21) ((((tmp,n=1,nvarin),i=1,nx),j=1,ny),k=1,nz)
         end do

!     print *, 'readig file: ', trim(fn)
!     print *, 'nvarin = ', nvarin
!     print *, 'nx = ', nx
!     print *, 'ny = ', ny
!     print *, 'nz = ', nz

         do nt = 1, ntpoint
            read(21) ((((timeseries(n,i,j,k,nt),n=1,nvarin),i=1,nx),j=1,ny),k=1,nz)
         end do
         close(21)

!         print *, 'p2 = ', timeseries(4,1,1:ny,1,1)

      endif
    end subroutine ReadTsData

end module MTSTimeData

module MTSKplane
  use MTSTimeData
  implicit none
  type node_info
      logical :: IsNodeIncluded = .false.
      integer :: myid, nv, nx, ny, nz
      integer, dimension(:), allocatable :: iloc, jloc, kloc
  end type node_info
  type(node_info), dimension(:), allocatable, private :: node_kp
  integer, private :: iminloc, imaxloc
  integer, private :: nxpoint, nypoint
  integer, private :: ioutputmin, ioutputmax, joutputmin, joutputmax, koutputmin, koutputmax   ! Index ranges of timeseries output by DNS code
  integer, private :: iskip, jskip
  integer, private :: ibe_ave, iend_ave, jbe_ave, jend_ave
  integer, private :: nvar, nxpt, nypt, nodes_k
  integer, private :: iascii, nskip
  parameter(nvar=26)
  character(3) :: varname_Kplane(nvar)
  parameter(varname_Kplane = (/'u','v','w','p','t','ux','uy','uz','vx','vy','vz','wx','wy','wz','uxx','uxy','uxz', &
                               'vxy','vyy','vyz','wxz','wyz','wzz','uk','vk','wk'/))
  
  contains

    subroutine GetDNSInfo_kplane(nodes, iasc, ns, dxs, dys, ibe_k, iend_k, iskip_k, jbe_k, jend_k, jskip_k)
      integer, intent(in) :: nodes
      integer, intent(in) :: iasc, ns
      real(8), intent(out) :: dxs, dys
      integer, intent(in) :: ibe_k, iend_k, jbe_k, jend_k
      integer, intent(out) :: iskip_k, jskip_k
      integer :: n, i, j, ibe, jbe, tmp
      integer :: imintmp, imaxtmp, jmintmp, jmaxtmp

      nskip = ns
      iascii = iasc
      nodes_k = nodes
      if(allocated(node_kp)) deallocate(node_kp)
      allocate(node_kp(nodes_k))
      open(20,file=trim(flowpath)//'kplane_index.dat')
      rewind(20)
      read(20,*)
      
      nxpt = 0; nypt = 0
      ioutputmin = 100000; ioutputmax = 0
      joutputmin = 100000; joutputmax = 0
      do n = 1, nodes_k
!         read(20,*) node_kp(n)%myid,node_kp(n)%nv,node_kp(n)%nx,node_kp(n)%ny, &
!                    ibe,tmp,iskip, &
!                    jbe,tmp,jskip, dxs, dys
         read(20,*) node_kp(n)%myid,node_kp(n)%nx,node_kp(n)%ny, &
                    ibe,tmp,iskip, &
                    jbe,tmp,jskip, dxs, dys
         node_kp(n)%nv = 23
         allocate(node_kp(n)%iloc(node_kp(n)%nx))
         allocate(node_kp(n)%jloc(node_kp(n)%ny))
         if(nxpt.lt.node_kp(n)%nx) nxpt = node_kp(n)%nx
         if(nypt.lt.node_kp(n)%ny) nypt = node_kp(n)%ny
         do i = 1, node_kp(n)%nx
            node_kp(n)%iloc(i) = ibe + (i-1)*iskip
         enddo
         do j = 1, node_kp(n)%ny
            node_kp(n)%jloc(j) = jbe + (j-1)*jskip
         enddo
         imintmp = node_kp(n)%iloc(1)
         imaxtmp = node_kp(n)%iloc(node_kp(n)%nx)
         jmintmp = node_kp(n)%jloc(1)
         jmaxtmp = node_kp(n)%jloc(node_kp(n)%ny)
         if(ioutputmin.gt.imintmp) ioutputmin = imintmp
         if(joutputmin.gt.jmintmp) joutputmin = jmintmp
         if(ioutputmax.lt.imaxtmp) ioutputmax = imaxtmp
         if(joutputmax.lt.jmaxtmp) joutputmax = jmaxtmp
       enddo
       if(any(node_kp%nv.ne.nvar)) then
           print *, 'nvar_input, nvar_read NOT match for Kplanes'
       endif
       node_kp%nz = 1
       close(20)

        ibe_ave = ibe_k
        iend_ave = iend_k
        jbe_ave = jbe_k
        jend_ave = jend_k

       if(ibe_ave .lt.ioutputmin.or.iend_ave.gt.ioutputmax.or.jbe_ave .lt.joutputmin.or.jend_ave.gt.joutputmax)  then
          print *, 'ibe_inp,iend_inp, jbe_inp, jend_inp out of range. STOP'
          stop
       endif

      print *, 'DNS Output timeseries spatial index range (K-plane)'
      print *, '   ioutputmin=', ioutputmin, 'ioutputmax=',ioutputmax
      print *, '   joutputmin=', joutputmin, 'joutputmax=',joutputmax

      print *, 'Physical index ranges for reading TS (K-plane)'
      print *, '   ibe =', ibe_ave, 'iend =', iend_ave
      print *, '   jbe =', jbe_ave, 'jend =', jend_ave

        iskip_k = iskip
        jskip_k = jskip

    end subroutine GetDNSInfo_kplane

    subroutine InitTSKplane(nxp, nyp, ibe_k, iend_k)
    ! ns: number of data points to be skipped
    ! iasc: wether the data is iascii (iasc = 1) or binary/unformatted (iasc = 0)
      integer, intent(in) :: ibe_k, iend_k
      integer, intent(out) :: nxp, nyp
      integer :: n, i, j, ibe, jbe, tmp
 
        ibe_ave = ibe_k
        iend_ave = iend_k

      print *, 'Physical index ranges for reading TS (K-plane)'
      print *, '   ibe =', ibe_ave, 'iend =', iend_ave
      print *, '   jbe =', jbe_ave, 'jend =', jend_ave

      do n = 1, nodes_k
         if(any(node_kp(n)%iloc.ge.ibe_ave.and.node_kp(n)%iloc.le.iend_ave).and. &
            any(node_kp(n)%jloc.ge.jbe_ave.and.node_kp(n)%jloc.le.jend_ave)) then
            node_kp(n)%IsNodeIncluded = .true.
         endif
      enddo

      if(.not.any(node_kp%IsNodeIncluded)) then
          print *, 'None of the nodes include points required. STOP'
          stop
      endif

! Dimension range for tsdata 

      nxpoint = (iend_ave-ibe_ave)/iskip + 1
      nypoint = (jend_ave-jbe_ave)/jskip + 1
      nxp = nxpoint
      nyp = nypoint

    end subroutine InitTSKplane

    subroutine ReadTSKplane(kindex, nvar_output, varindex_output, varout, zloc, dzdk)
     implicit none
     integer, intent(in) :: kindex, nvar_output
     integer, intent(in) :: varindex_output(nvar_output)
     real(8), intent(out) :: varout(nvar_output,nxpoint,nypoint,ntpoint)
     real(8), intent(out) :: zloc, dzdk
     integer :: n, i, j, ii, jj, nn, k
     character(4) :: fnum4, fnum5
     character(300) :: fileinput
     real(8) :: tsdata_pe(nvar_output,nxpt,nypt,1,ntpoint)
     integer :: istat, k1

     write(unit=fnum4,fmt='(I04.4)') kindex
!     open(100,file=trim(flowpath)//'kplane_zloc.dat',status='old')
!     rewind(100)
!     read(100,*)
!     do ! k1 = 1, 10 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        read(100,*,iostat=istat) k, zloc, dzdk
!        if(k.eq.kindex.or.istat.lt.0) exit  ! istat < 0: endfile condition
!     enddo
!     if(k.ne.kindex) then
!        print *, 'kindex should match one of those in the file kplane_zloc.dat'
!        stop
!     endif
     do n = 1, nodes_k
        if(node_kp(n)%IsNodeIncluded) then
            write(unit=fnum5,fmt='(I04.4)') node_kp(n)%myid
            fileinput = trim(flowpath)//'series_kplane'//fnum4//'_nid'//fnum5
            call ReadTsData(trim(fileinput),node_kp(n)%nv,node_kp(n)%nx,node_kp(n)%ny,1,&
                           tsdata_pe(1:node_kp(n)%nv,1:node_kp(n)%nx,1:node_kp(n)%ny,:,:), iascii, nskip)
            do j = 1, node_kp(n)%ny
            do i = 1, node_kp(n)%nx
                 if(node_kp(n)%iloc(i).ge.ibe_ave.and.node_kp(n)%iloc(i).le.iend_ave.and. &
                    node_kp(n)%jloc(j).ge.jbe_ave.and.node_kp(n)%jloc(j).le.jend_ave) then
                    ii = (node_kp(n)%iloc(i)-ibe_ave)/iskip+1
                    jj = (node_kp(n)%jloc(j)-jbe_ave)/jskip+1

                    do nn = 1, nvar_output                 !!!!!!
                       varout(nn,ii,jj,:) = tsdata_pe(varindex_output(nn),i,j,1,:) 
                    enddo
                 endif
             enddo
             enddo
        endif
     enddo
  end subroutine ReadTSKplane

end module MTSKplane

module MTSIplane
  use MTSTimeData
  implicit none
  type node_info
      logical :: IsNodeIncluded = .false.
      integer :: myid, nv, nx, ny, nz
      integer, dimension(:), allocatable :: iloc, jloc, kloc
  end type node_info
  type(node_info), dimension(:), allocatable, private :: node_ip
  integer, private :: nodes_i
  integer, private :: nvar, nypt, nzpt
  integer, private :: nypoint, nzpoint
  integer, private :: joutputmin, joutputmax, koutputmin, koutputmax   ! Index ranges of timeseries output by DNS code
  integer, private :: jskip, kskip
  integer, private :: jbe_ave, jend_ave, kbe_ave, kend_ave
  integer, private :: iascii, nskip
  parameter(nvar=26)
  character(3) :: varname_Iplane(nvar)
  parameter(varname_Iplane = (/'u','v','w','p','t','ux','uy','uz','vx','vy','vz','wx','wy','wz','uxx','uxy','uxz', &
                               'vxy','vyy','vyz','wxz','wyz','wzz','ui','vi','wi'/))
contains
   
    subroutine GetDNSInfo_iplane(iloc, nodes, iasc, ns, dxs, dys, jbe_i, jend_i, jskip_i, kbe_i, kend_i, kskip_i)
      integer, intent(in) :: iloc, nodes
      integer, intent(in) :: iasc, ns
      integer, intent(in) :: jbe_i, jend_i, kbe_i, kend_i
      integer, intent(out) :: jskip_i, kskip_i
      real(8), intent(out) :: dxs, dys
      integer :: n, i, j, k, jbe, kbe, tmp
      integer :: jmintmp, jmaxtmp
      character(4) :: fnum
      character(200) :: filename

      iascii = iasc
      nskip = ns
      nodes_i = nodes
      if(allocated(node_ip)) deallocate(node_ip)
      allocate(node_ip(nodes))
      write(unit=fnum,fmt='(I04.4)') iloc
      filename = trim(flowpath)//'iplane_index_i'//fnum//'.dat'
      open(20,file=trim(filename))
      rewind(20)
      read(20,*)
      nypt = 0; nzpt = 0
      joutputmin = 100000; joutputmax = 0
      do n = 1, nodes
!         read(20,*) node_ip(n)%myid,node_ip(n)%nv,node_ip(n)%ny,node_ip(n)%nz, &
!                    jbe,tmp,jskip, &
!                    kbe,tmp,kskip, dxs, dys
         read(20,*) node_ip(n)%myid,node_ip(n)%ny,node_ip(n)%nz, &
                    jbe,tmp,jskip, &
                    kbe,tmp,kskip, dys ! dxs
         node_ip(n)%nv = 23
         if(allocated(node_ip(n)%jloc)) deallocate(node_ip(n)%jloc)
         if(allocated(node_ip(n)%kloc)) deallocate(node_ip(n)%kloc)
         allocate(node_ip(n)%jloc(node_ip(n)%ny))
         allocate(node_ip(n)%kloc(node_ip(n)%nz))
         if(nypt.lt.node_ip(n)%ny) nypt = node_ip(n)%ny
         if(nzpt.lt.node_ip(n)%nz) nzpt = node_ip(n)%nz
         do j = 1, node_ip(n)%ny
            node_ip(n)%jloc(j) = jbe + (j-1)*jskip
         enddo
         do k = 1, node_ip(n)%nz
            node_ip(n)%kloc(k) = kbe + (k-1)*kskip
         enddo
         jmintmp = node_ip(n)%jloc(1)
         jmaxtmp = node_ip(n)%jloc(node_ip(n)%ny)
         if(joutputmin.gt.jmintmp) joutputmin = jmintmp
         if(joutputmax.lt.jmaxtmp) joutputmax = jmaxtmp
      enddo
      koutputmin = node_ip(1)%kloc(1)
      koutputmax = node_ip(1)%kloc(node_ip(1)%nz)
       if(any(node_ip%nv.ne.nvar)) then
           print *, 'nvar_input, nvar_read NOT match for Iplanes'
       endif
      node_ip%nx = 1
      close(20)


      jbe_ave = jbe_i
      jend_ave = jend_i
      kbe_ave = kbe_i
      kend_ave = kend_i

       if(jbe_ave .lt.joutputmin.or.jend_ave.gt.joutputmax.or.kbe_ave.lt.koutputmin.or.kend_ave.gt.koutputmax)  then
          print *, 'jbe_inp,jend_inp, kbe_inp, kend_inp out of range. STOP'
          stop
       endif
 
      print *, 'DNS Output time series index range (i plane)'
      print *, '   joutputmin=', joutputmin, 'joutputmax=',joutputmax
      print *, '   koutputmin=', koutputmin, 'koutputmax=',koutputmax

      print *, 'Physical index ranges for reading TS (i plane)'
      print *, '   jbe =', jbe_ave, 'jend =', jend_ave
      print *, '   kbe =', kbe_ave, 'kend =', kend_ave

      jskip_i = jskip
      kskip_i = kskip 

    end subroutine GetDNSInfo_iplane


    subroutine InitTSIplane(nyp, nzp, jbe_i, jend_i)
    ! ns: number of data points to be skipped
    ! iasc: wether the data is iascii (iasc = 1) or binary/unformatted (iasc = 0)
      integer, intent(in) :: jbe_i, jend_i
      integer, intent(out) :: nyp, nzp
      integer :: n, i, j, k, jbe, kbe

      jbe_ave = jbe_i
      jend_ave = jend_i

      print *, 'Physical index ranges for reading TS (i plane)'
      print *, '   jbe =', jbe_ave, 'jend =', jend_ave
      print *, '   kbe =', kbe_ave, 'kend =', kend_ave


      nypoint = (jend_ave-jbe_ave)/jskip + 1
      nzpoint = (kend_ave-kbe_ave)/kskip + 1
      nyp = nypoint
      nzp = nzpoint      

      do n = 1, nodes_i
         if(any(node_ip(n)%jloc.ge.jbe_ave.and.node_ip(n)%jloc.le.jend_ave)) then
            node_ip(n)%IsNodeIncluded = .true.
         endif
      enddo

      if(.not.any(node_ip%IsNodeIncluded)) then
          print *, 'None of the nodes include points required. STOP'
          stop
      endif

    end subroutine InitTSIplane

    subroutine ReadTSIplane(iindex, nvar_output, varindex_output, varout)
       implicit none
       integer, intent(in) :: iindex, nvar_output
       integer, intent(in) :: varindex_output(nvar_output)
       real(8), intent(out) :: varout(nvar_output,nypoint,nzpoint,ntpoint)               
       integer :: n, k, j, jj, kk, nn
       real(8) :: tsdata_pe(nvar_output,1,nypt,nzpt,ntpoint) !!!!!!!!!!!
       character(4) :: fnum4, fnum5
       character(300) :: fileinput

       write(unit=fnum4,fmt='(I04.4)') iindex
       do n = 1, nodes_i                     !!!!!!!!!!!!!!!!!!!!!! node_i
          if(node_ip(n)%IsNodeIncluded) then
             write(unit=fnum5,fmt='(I04.4)') node_ip(n)%myid
             fileinput = trim(flowpath)//'series_iplane'//fnum4//'_nid'//fnum5
             call ReadTsData(trim(fileinput),node_ip(n)%nv,1,node_ip(n)%ny,node_ip(n)%nz,&
                              tsdata_pe(:,:,1:node_ip(n)%ny,1:node_ip(n)%nz,:),iascii,nskip)

!   if(n.eq.1) print *, 'p = ', tsdata_pe(4,1,1:5,1,1)

              do k = 1, node_ip(n)%nz
              do j = 1, node_ip(n)%ny
                 if(node_ip(n)%jloc(j).ge.jbe_ave.and.node_ip(n)%jloc(j).le.jend_ave.and. &
                    node_ip(n)%kloc(k).ge.kbe_ave.and.node_ip(n)%kloc(k).le.kend_ave) then
                    jj = (node_ip(n)%jloc(j)-jbe_ave)/jskip+1
                    kk = (node_ip(n)%kloc(k)-kbe_ave)/kskip+1
                    do nn = 1, nvar_output                  !!!!!!!!
                       varout(nn,jj,kk,:) = tsdata_pe(varindex_output(nn),1,j,k,:) 
                    enddo
                 endif
              enddo  ! end j loop
              enddo  ! end k loop
          endif
       enddo

   end subroutine ReadTSIplane
end module MTSIplane

module MTSJplane
  use MTSTimeData
  implicit none
  type node_info
      logical :: IsNodeIncluded = .false.
      integer :: myid, nv, nx, ny, nz
      integer, dimension(:), allocatable :: iloc, jloc, kloc
  end type node_info
  type(node_info), dimension(:), allocatable, private :: node_jp
  integer, private :: nxpoint, nzpoint
  integer, private :: ioutputmin, ioutputmax, joutputmin, joutputmax, koutputmin, koutputmax   ! Index ranges of timeseries output by DNS code
  integer, private :: iskip, kskip
  integer, private :: ibe_ave, iend_ave, kbe_ave, kend_ave
  integer, private :: nvar, nxpt, nzpt, nodes_j
  integer, private :: iascii, nskip
  parameter(nvar=26)
  character(3) :: varname_Jplane(nvar)
  parameter(varname_Jplane = (/'u','v','w','p','t','ux','uy','uz','vx','vy','vz','wx','wy','wz','uxx','uxy','uxz', &
                               'vxy','vyy','vyz','wxz','wyz','wzz','uj','vj','wj'/))
  contains

    subroutine GetDNSInfo_jplane(jloc, nodes, iasc, ns, dxs, dys, ibe_j, iend_j, iskip_j, kbe_j, kend_j, kskip_j)
      integer, intent(in) :: jloc, nodes
      integer, intent(in) :: iasc, ns
      integer, intent(in) :: ibe_j, iend_j, kbe_j, kend_j
      integer, intent(out) :: iskip_j, kskip_j
      real(8), intent(out) :: dxs, dys
      integer :: n, i, k, ibe, kbe, tmp
      integer :: imintmp, imaxtmp
      integer :: temp(1)
      character(4) :: fnum
      character(200) :: filename

      iascii = iasc
      nskip = nskip
      nodes_j = nodes
      if(allocated(node_jp)) deallocate(node_jp)
      allocate(node_jp(nodes_j))
      write(unit=fnum,fmt='(I04.4)') jloc
      filename = trim(flowpath)//'jplane_index_j'//fnum//'.dat'  
      open(20,file=trim(filename))
      rewind(20)
      read(20,*)
      nxpt = 0; nzpt = 0
      ioutputmin = 100000; ioutputmax = 0
      do n = 1, nodes_j
!         read(20,*) node_jp(n)%myid,node_jp(n)%nv,node_jp(n)%nx,node_jp(n)%nz, &
!                    ibe,tmp,iskip, &
!                    kbe,tmp,kskip, dxs, dys
         read(20,*) node_jp(n)%myid,node_jp(n)%nx,node_jp(n)%nz, &
                    ibe,tmp,iskip, &
                    kbe,tmp,kskip, dxs, dys
         node_jp(n)%nv = 23
         allocate(node_jp(n)%iloc(node_jp(n)%nx))
         allocate(node_jp(n)%kloc(node_jp(n)%nz))
         if(nxpt.lt.node_jp(n)%nx) nxpt = node_jp(n)%nx
         if(nzpt.lt.node_jp(n)%nz) nzpt = node_jp(n)%nz
         do i = 1, node_jp(n)%nx
            node_jp(n)%iloc(i) = ibe + (i-1)*iskip
         enddo
         do k = 1, node_jp(n)%nz
            node_jp(n)%kloc(k) = kbe + (k-1)*kskip
         enddo
         imintmp = node_jp(n)%iloc(1)
         imaxtmp = node_jp(n)%iloc(node_jp(n)%nx)
         if(ioutputmin.gt.imintmp) ioutputmin = imintmp
         if(ioutputmax.lt.imaxtmp) ioutputmax = imaxtmp
       enddo
       koutputmin = node_jp(1)%kloc(1)
       koutputmax = node_jp(1)%kloc(node_jp(1)%nz)
       if(any(node_Jp%nv.ne.nvar)) then
           print *, 'nvar_input, nvar_read NOT match for Jplanes'
       endif
       node_jp%ny = 1
       close(20)

       ibe_ave = ibe_j
       iend_ave = iend_j
       kbe_ave = kbe_j
       kend_ave = kend_j

      if(ibe_ave .lt.ioutputmin.or.iend_ave.gt.ioutputmax.or.kbe_ave.lt.koutputmin.or.kend_ave.gt.koutputmax)  then
          print *, 'ibe_inp,iend_inp, kbe_inp, kend_inp out of range. STOP'
          stop
       endif 

       print *, 'DNS Output timeseries spatial index range (J Plane)'
       print *, '   ioutputmin=', ioutputmin, 'ioutputmax=',ioutputmax
       print *, '   koutputmin=', koutputmin, 'koutputmax=',koutputmax

       print *, 'Physical index ranges for reading TS (J-plane)'
       print *, '   ibe =', ibe_ave, 'iend =', iend_ave
       print *, '   kbe =', kbe_ave, 'kend =', kend_ave

       iskip_j = iskip
       kskip_j = kskip

    end subroutine GetDNSInfo_jplane

    subroutine InitTSJplane(nxp, nzp, ibe_j, iend_j)
    ! ns: number of data points to be skipped
    ! iasc: wether the data is iascii (iasc = 1) or binary/unformatted (iasc = 0)
      integer, intent(in) :: ibe_j, iend_j
      integer, intent(out) :: nxp, nzp
      integer :: n, i, k, ibe, kbe, tmp
     
       ibe_ave = ibe_j
       iend_ave = iend_j

       print *, 'Physical index ranges for reading TS (J-plane)'
       print *, '   ibe =', ibe_ave, 'iend =', iend_ave
       print *, '   kbe =', kbe_ave, 'kend =', kend_ave


      ! begin, end, and number of spatial points for outputting

      nxpoint = (iend_ave-ibe_ave)/iskip + 1
      nzpoint = (kend_ave-kbe_ave)/kskip + 1
      ! Dimension range for tsdata 
      nxp = nxpoint
      nzp = nzpoint          

      do n = 1, nodes_j
         if(any(node_jp(n)%iloc.ge.ibe_ave.and.node_jp(n)%iloc.le.iend_ave)) then
            node_jp(n)%IsNodeIncluded = .true.
         endif
      enddo

      if(.not.any(node_jp%IsNodeIncluded)) then
          print *, 'None of the nodes include points required. STOP'
          stop
      endif

    end subroutine InitTSJplane

    subroutine ReadTSJplane(jindex, nvar_output, varindex_output, varout)
     implicit none
     integer, intent(in) :: jindex, nvar_output
     integer, intent(in) :: varindex_output(nvar_output)
     real(8), intent(out) :: varout(nvar_output,nxpoint,nzpoint,ntpoint)
     integer :: n, i, k, ii, kk, nn
     character(4) :: fnum4, fnum5
     character(300) :: fileinput
     real(8) :: tsdata_pe(nvar,nxpt,1,nzpt,ntpoint)

     write(unit=fnum4,fmt='(I04.4)') jindex
     do n = 1, nodes_j
        if(node_jp(n)%IsNodeIncluded) then
            write(unit=fnum5,fmt='(I04.4)') node_jp(n)%myid
            fileinput = trim(flowpath)//'series_jplane'//fnum4//'_nid'//fnum5
            call ReadTsData(trim(fileinput),node_jp(n)%nv,node_jp(n)%nx,1,node_jp(n)%nz,&
                           tsdata_pe(1:node_jp(n)%nv,1:node_jp(n)%nx,1,1:node_jp(n)%nz,:), iascii, nskip)
            do k = 1, node_jp(n)%nz
            do i = 1, node_jp(n)%nx
                 if(node_jp(n)%iloc(i).ge.ibe_ave.and.node_jp(n)%iloc(i).le.iend_ave.and.&
                    node_jp(n)%kloc(k).ge.kbe_ave.and.node_jp(n)%kloc(k).le.kend_ave) then
                    ii = (node_jp(n)%iloc(i)-ibe_ave)/iskip+1
                    kk = (node_jp(n)%kloc(k)-kbe_ave)/kskip+1
                    do nn = 1, nvar_output                  !!!!!!
                       varout(nn,ii,kk,:) = tsdata_pe(varindex_output(nn),i,1,k,:) 
                    enddo
                 endif
             enddo
             enddo
        endif
     enddo
  end subroutine ReadTSJplane

end module MTSJplane


