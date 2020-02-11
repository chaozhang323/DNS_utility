!> 

program tsplanes
  use decomp_2d
  use decomp_2d_io
  use modVolume
  use modRWPHDF5
  implicit none

  character(400) :: fname
  character(400) :: fname_grid
  character(200) :: datapath
  character(8) :: fnum
  character(4) :: fnum1
  integer :: idim, jdim, kdim
  integer :: iil, iih, stride
  integer :: ivolume, iplane_xz, iplane_yz, iplane_xy
  type(index_output) :: idx_vol, idx_xyp, idx_yzp, idx_xzp
  integer, dimension(:),  allocatable :: planes_xy, planes_yz, planes_xz
  real(8), dimension(:,:,:), allocatable :: x, y, z
  real(8), dimension(:,:,:), allocatable :: xh, yh, zh
  real(8), dimension(:,:,:,:), allocatable :: buffer_grid
  real(8), dimension(:,:,:), allocatable :: buffer_grid_tmp
  real(8), dimension(:,:,:,:), allocatable :: buffer_gridh
  real(8), dimension(:,:,:,:,:), allocatable :: buffer_flow
  real(8), dimension(:,:,:), allocatable :: buffer_flow_tmp
  real(8), dimension(:,:,:,:,:), allocatable :: buffer_flowh
  real(8), dimension(:,:,:), allocatable :: grad_p, grad_rho

  integer :: n,i,j,k,m, nn
  integer :: hdferr, rank, ierror
  integer :: ilen, jlen, klen, ibe, jbe, kbe
  integer :: ist_tsjplane, iend_tsjplane, isp_tsjplane, kst_tsjplane, kend_tsjplane, ksp_tsjplane, nt, nt_skip
  integer :: jst_tsiplane, jend_tsiplane, jsp_tsiplane, kst_tsiplane, kend_tsiplane, ksp_tsiplane

  real(8), dimension(:,:,:,:), allocatable :: drdj, drdi  !!!!
  real(8), dimension(:,:,:), allocatable :: drdj_tmp, drdi_tmp
  real(8), dimension(:,:,:,:), allocatable :: drdjh, drdih

  integer :: p_row, p_col !!!!!!!

  call MPI_INIT(ierror)

  call InitHDF5()
  call Input()

  ilen = xsize(1)
  jlen = xsize(2)
  klen = xsize(3)
  ibe = xstart(1)
  jbe = xstart(2)
  kbe = xstart(3)

!print *, 'nrank = ', nrank
 !       print *, 'xsize(1) = ', xsize(1)
 !       print *, 'xsize(2) = ', xsize(2)
 !       print *, 'xsize(3) = ', xsize(3)
 !       print *, 'xstart(1) = ', xstart(1)
!        print *, 'xstart(2) = ', xstart(2)
!        print *, 'xstart(3) = ', xstart(3)

  fname_grid = trim(datapath)//'timeseries_GridMetrics.h5'
  if(nrank.eq.0) print *, 'reading grid: ', trim(fname_grid)

  if(iplane_yz.eq.1) then
    call ReadtsGridMetrics_P(fname_grid,"/iplane",xsize(1),xsize(2),xsize(3),xstart(1),xstart(2),xstart(3),buffer_grid)

    ! update halo-cell
    do i=1, 12
      call update_halo(buffer_grid(:,:,:,i),buffer_grid_tmp,2)
!    if(nrank.eq.0) print *, 'buffer_grid_tmp(1,2,1) = ', buffer_grid_tmp(1,2,1)

      if(i.eq.1) allocate(buffer_gridh(size(buffer_grid_tmp,1),size(buffer_grid_tmp,2),size(buffer_grid_tmp,3),12))
      buffer_gridh(1:size(buffer_grid_tmp,1),1:size(buffer_grid_tmp,2),1:size(buffer_grid_tmp,3),i) = &
                buffer_grid_tmp(1:size(buffer_grid_tmp,1),1:size(buffer_grid_tmp,2),1:size(buffer_grid_tmp,3))
      deallocate(buffer_grid_tmp)
    enddo ! end i loop

!    if(nrank.eq.0) print *, 'buffer_gridh(1,2,1,1) = ', buffer_gridh(1,2,1,1)

    ! for debug
!    print *, 'size(buffer_grid,1) = ', size(buffer_grid,1)
!    print *, 'size(buffer_grid,2) = ', size(buffer_grid,2)
!    print *, 'size(buffer_grid,3) = ', size(buffer_grid,3)
    if(nrank.eq.0) then
      print *, 'size(buffer_gridh,1) = ', size(buffer_gridh,1)
      print *, 'size(buffer_gridh,2) = ', size(buffer_gridh,2)
      print *, 'size(buffer_gridh,3) = ', size(buffer_gridh,3)
    endif

    do i=1, xsize(1)
      do j=1, xsize(2)
        do k=1, xsize(3)
          y(k,i,j) = buffer_grid(i,j,k,2)
          z(k,i,j) = buffer_grid(i,j,k,3)
        enddo
      enddo
    enddo

    ! for debug
    if(nrank.eq.0) then
!        print *, 'xsize(1) = ', xsize(1)
!        print *, 'xsize(2) = ', xsize(2)
!        print *, 'xsize(3) = ', xsize(3)
!        print *, 'xstart(1) = ', xstart(1)
!        print *, 'xstart(2) = ', xstart(2)
!        print *, 'xstart(3) = ', xstart(3)
!        print *, ' ', buffer_grid(1,1,1,1)
!        print *, 'z(1,1,1) = ', z(1,1,1)
!        print *, 'z(1,1,2) = ', z(1,1,2)
        print *, 'buffer_grid(1,1,1,3) = ', buffer_grid(1,1,1,3)
        print *, 'buffer_grid(1,2,1,3) = ', buffer_grid(1,2,1,3)
        print *, 'buffer_gridh(1,1,1,3) = ', buffer_gridh(1,1,1,3)
        print *, 'buffer_gridh(1,2,1,3) = ', buffer_gridh(1,2,1,3)

    endif

!    ilen = 1
!    jlen = (idx_yzp%jend - idx_yzp%jst)/idx_yzp%jsp + 1
!    klen = (idx_yzp%kend - idx_yzp%kst)/idx_yzp%ksp + 1

!    idx_yzp%ist=1
!    idx_yzp%iend=1
!    idx_yzp%isp=1

    do n=iil,iih,stride
      write(unit=fnum,fmt='(I08.8)')n
      fname = trim(datapath)//'timeseries_'//fnum//'.h5'

      if(nrank.eq.0) print *, 'reading file: ', trim(fname)

!!      call Readtsflow(fname,"/iplane",nt,jdim,kdim,1,buffer_flow(:,:,:,:,1:10),2)
      call Readtsflow_P(fname,"/iplane",nt,xsize(1),xsize(2),xsize(3),xstart(1),xstart(2),xstart(3),buffer_flow(:,:,:,:,1:10),2)


      if(nrank.eq.0) then
        print *, 'buffer_flow(1,1,2,1,1) = ', buffer_flow(1,1,2,1,1)
      endif

      do nn=1, nt, nt_skip
       do j=1, xsize(1) !jdim
          do k=1, xsize(2) !kdim
            drdi(nn,j,k,1) = 1/287.0*(1.0/buffer_flow(nn,j,k,1,5)*buffer_flow(nn,j,k,1,9) - &
                            buffer_flow(nn,j,k,1,4)*buffer_flow(nn,j,k,1,10)/(buffer_flow(nn,j,k,1,5)**2)) ! drhodi
            buffer_flow(nn,j,k,1,11) = buffer_flow(nn,j,k,1,4)/dble(287.0*buffer_flow(nn,j,k,1,5)) ! rho
          enddo
        enddo
      enddo

      do nn=1, nt, nt_skip
        call update_halo(drdi(nn,:,:,:), drdi_tmp, 2)
        if(nn.eq.1) allocate(drdih(nt,size(drdi_tmp,1), size(drdi_tmp,2),size(drdi_tmp,3)) )
        drdih(nn,1:size(drdi_tmp,1),1:size(drdi_tmp,2),1:size(drdi_tmp,3)) = &
             drdi_tmp(1:size(drdi_tmp,1),1:size(drdi_tmp,2),1:size(drdi_tmp,3))
        deallocate(drdi_tmp)
      enddo

      do nn=1, nt, nt_skip
        do i=1, 11
          call update_halo(buffer_flow(nn,:,:,:,i), buffer_flow_tmp,2)
          if(nn.eq.1.and.i.eq.1) allocate(buffer_flowh(nt,size(buffer_flow_tmp,1),size(buffer_flow_tmp,2),size(buffer_flow_tmp,3),11))
          buffer_flowh(nn,1:size(buffer_flow_tmp,1),1:size(buffer_flow_tmp,2),1:size(buffer_flow_tmp,3),i) = &
                       buffer_flow_tmp(1:size(buffer_flow_tmp,1),1:size(buffer_flow_tmp,2),1:size(buffer_flow_tmp,3))
          deallocate(buffer_flow_tmp)
        enddo ! end i loop
      enddo ! end nn loop

     if(nrank.eq.0) then
       print *, 'buffer_flowh(1,1,2,1,1)= ', buffer_flowh(1,1,2,1,1)
     endif

      do nn=1, nt, nt_skip

      ! write grad_rho
!        call CalGradient_tsiplane(idx_yzp,jdim,kdim,buffer_grid,buffer_flow(nn,:,:,1,11),drdi(nn,:,:,1),grad_rho,jst_tsiplane,kst_tsiplane)

     print *, 'nrank = ', nrank
        call CalGradient_tsiplane(idx_yzp,xsize(1),xsize(2)+4,buffer_gridh(:,:,1:1,1:12),buffer_flowh(nn,:,:,1,11),drdih(nn,:,:,1),grad_rho,jst_tsiplane,kst_tsiplane)

        write(unit=fnum1,fmt='(I04.4)') nn
        fname = 'tsiplane_'//fnum//'-'//fnum1//'.h5'
        if(nrank.eq.0) print *, 'writing file : ', trim(fname)

        if(nrank.eq.0) then
          print *, 'grad_rho(1,1,1) = ', grad_rho(1,1,1)
        endif

!        print *, 'nrank = ', nrank
!        print *, 'size(grad_rho,1) = ', size(grad_rho,1)
!        print *, 'size(grad_rho,2) = ', size(grad_rho,2)
!        print *, 'size(grad_rho,3) = ', size(grad_rho,3)
!        print *, 'xsize(1) = ', xsize(1)
!        print *, 'xsize(2) = ', xsize(2)
!        print *, 'xsize(3) = ', xsize(3)
!        print *, 'xstart(1) = ', xstart(1)
!        print *, 'xstart(2) = ', xstart(2)
!        print *, 'xstart(3) = ', xstart(3)

        call Writetsplanes_P(fname,idim,jdim,kdim,xstart(3),xstart(1),xstart(2),grad_rho(1:xsize(3),1:xsize(1),1:xsize(2)) )


      enddo ! end nn loop
    enddo ! end n loop
  endif ! end read tsiplane

  call FinalizeHDF5()
  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)


 contains
    subroutine Input()
      integer ::i

!      read(*,*)
!      read(*,'(A)') datapath
!      read(*,*)
!      read(*,*) iil, iih, stride
!      read(*,*)
!      read(*,*) p_row, p_col
!      read(*,*)
!      read(*,*) iplane_yz, idx_yzp%jst, idx_yzp%jend, idx_yzp%jsp, idx_yzp%kst, idx_yzp%kend, idx_yzp%ksp
!      read(*,*)
!      read(*,*) jst_tsiplane, jend_tsiplane, jsp_tsiplane, kst_tsiplane, kend_tsiplane, ksp_tsiplane, nt, nt_skip
!      read(*,*)
!      read(*,*)
!      read(*,*) iplane_xz, idx_xzp%ist, idx_xzp%iend, idx_xzp%isp, idx_xzp%kst, idx_xzp%kend, idx_xzp%ksp
!      read(*,*)
!      read(*,*) ist_tsjplane, iend_tsjplane, isp_tsjplane, kst_tsjplane, kend_tsjplane, ksp_tsjplane, nt, nt_skip


      datapath = "/share/duan/duanl/TestDNS/testPHDF5/M14_fromPleiades/TIMESERIES_1/"
      iil = 311000; iih = 311000; stride = 1000
      p_row = 4
      p_col = 1
      iplane_yz = 1
      idx_yzp%jst=1; idx_yzp%jend=460; idx_yzp%jsp=1; idx_yzp%kst=1; idx_yzp%kend=460; idx_yzp%ksp=1

      jst_tsiplane=1; jend_tsiplane=460; jsp_tsiplane=1; kst_tsiplane=1; kend_tsiplane=460; ksp_tsiplane=1; nt=100; nt_skip=100




      if(iplane_yz.eq.1) then
        idim = 1
        jdim = jend_tsiplane-jst_tsiplane+1
        kdim = kend_tsiplane-kst_tsiplane+1
!        allocate(buffer_flow(nt,jdim,kdim,idim,11))
!        allocate(drdi(nt,jdim,kdim,idim))
!        allocate(buffer_grid(jdim,kdim,idim,12))

!        call decomp_2d_init(idim,jdim,kdim,p_row,p_col)
        call decomp_2d_init(jdim,kdim,idim,p_row,p_col)
        ! debug
!        print *, 'nrank = ', nrank
!        print *, 'xsize(1) = ', xsize(1)
!        print *, 'xsize(2) = ', xsize(2)
!        print *, 'xsize(3) = ', xsize(3)

        allocate(buffer_flow(nt,xsize(1),xsize(2),xsize(3),11))
        allocate(drdi(nt,xsize(1),xsize(2),xsize(3)))
        allocate(buffer_grid(xsize(1),xsize(2),xsize(3),12))
        allocate(x(xsize(3),xsize(1),xsize(2)), &
                 y(xsize(3),xsize(1),xsize(2)), &
                 z(xsize(3),xsize(1),xsize(2)) )
      endif

!      if(iplane_xz.eq.1) then
!        idim = iend_tsjplane-ist_tsjplane+1
!        jdim = 1
!        kdim = kend_tsjplane-kst_tsjplane+1
!        allocate(buffer_flow(nt,idim,kdim,jdim,11))
!        allocate(drdj(nt,idim,kdim,jdim))
!        allocate(buffer_grid(idim,kdim,jdim,12))
!      endif


!      allocate(x(idim,jdim,kdim),y(idim,jdim,kdim),z(idim,jdim,kdim))


    end subroutine Input


end program tsplanes
