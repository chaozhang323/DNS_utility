

program tsjplane
  use MFileIO
  use modVolume
  use modReadHDF5
  use modWriteHDF5
  use modTecbin
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
  real, dimension(:,:,:), allocatable :: x, y, z
  real, dimension(:,:,:,:), allocatable :: buffer_grid
!  real, dimension(:,:,:,:), allocatable :: buffer_flow
  real, dimension(:,:,:,:,:), allocatable :: buffer_flow
  real, dimension(:,:,:), allocatable :: grad_p, grad_rho
  real, dimension(:,:,:), allocatable :: omx, omy, omz
  real, dimension(:,:,:), allocatable :: swirl
  real, dimension(:,:,:), allocatable :: div
  real, dimension(:,:,:), allocatable :: Q
  integer :: n,i,j,k,m, nn
  real(8) :: uinf, delta, utau, ztau
  integer(HSIZE_T), dimension(:), allocatable :: dimsf
  integer :: hdferr, rank
  integer :: ilen, jlen, klen
  integer :: ist_tsjplane, iend_tsjplane, isp_tsjplane, kst_tsjplane, kend_tsjplane, ksp_tsjplane, nt, nt_skip
  integer :: jst_tsiplane, jend_tsiplane, jsp_tsiplane, kst_tsiplane, kend_tsiplane, ksp_tsiplane

  real(8), dimension(:,:,:,:), allocatable :: drdj, drdi  !!!!

  call Input()
  call InitHDF5()

  fname_grid = trim(datapath)//'timeseries_GridMetrics.h5'
  print *, 'reading grid: ', trim(fname_grid)

  if(iplane_xz.eq.1) then
    call ReadtsGridMetrics(fname_grid,"/jplane",idim,kdim,jdim,buffer_grid)
    do i=1,idim
      do j=1, jdim
        do k=1, kdim
          x(i,j,k) = buffer_grid(i,k,1,1)
          z(i,j,k) = buffer_grid(i,k,1,3)
        enddo
      enddo
    enddo

    ilen = (idx_xzp%iend - idx_xzp%ist)/idx_xzp%isp + 1
    jlen = 1
    klen = (idx_xzp%kend - idx_xzp%kst)/idx_xzp%ksp + 1

    rank=3
    allocate(dimsf(rank))
    dimsf(1) = ilen
    dimsf(2) = jlen
    dimsf(3) = klen

    idx_xzp%jst=1
    idx_xzp%jend=1
    idx_xzp%jsp=1
    ! init Tecplot file *.plt
    call InitTec(nt/nt_skip,ilen,jlen,klen,3,"x z grad_rho", 2)

    do n=iil,iih,stride
      write(unit=fnum,fmt='(I08.8)')n
      fname = trim(datapath)//'timeseries_'//fnum//'.h5'
      print *, 'reading file: ', trim(fname)
!      call Readtsflow_jplane(fname,nt,idim,kdim,10,buffer_flow(:,:,:,:,1:10))
      call Readtsflow(fname,"/jplane",nt,idim,kdim,jdim,buffer_flow(:,:,:,:,1:10),3)

      do nn=1, nt, nt_skip
        do i=1, idim
          do k=1, kdim
            drdj(nn,i,k,1) = 1/287.0*(1.0/buffer_flow(nn,i,k,1,5)*buffer_flow(nn,i,k,1,9) - &
                            buffer_flow(nn,i,k,1,4)*buffer_flow(nn,i,k,1,10)/(buffer_flow(nn,i,k,1,5)**2)) ! drhodj
            buffer_flow(nn,i,k,1,11) = buffer_flow(nn,i,k,1,4)/dble(287.0*buffer_flow(nn,i,k,1,5)) ! rho
          enddo
        enddo
      enddo

      print *, 'min(rho) = ', minval(buffer_flow(:,:,:,:,11))

      m = 0
      do nn=1, nt, nt_skip
        write(unit=fnum1,fmt='(I04.4)') nn
        fname = 'tsjplane_'//fnum//'-'//fnum1//'.h5'
        print *, 'writing file: ', trim(fname)
      ! write x, z
  !    call Write3DHDF5(fname, rank, dimsf,"x",x,0)
  !    call Write3DHDF5(fname, rank, dimsf,"z",z,1)
      ! write grad_rho
        call CalGradient_tsjplane(idx_xzp,idim,kdim,buffer_grid,buffer_flow(nn,:,:,1,11),drdj(nn,:,:,1),grad_rho,ist_tsjplane,kst_tsjplane)
  !   call Write3DHDF5(fname,rank,dimsf,"grad_rho",grad_rho,1)
      ! write grad_p
  !    call CalGradient_tsjplane(idx_xzp,idim,kdim,buffer_grid,buffer_flow(nn,:,:,1,4),buffer_flow(nn,:,:,1,9),grad_p,ist_tsjplane,kst_tsjplane)
  !    call Write3DHDF5(fname,rank,dimsf,"grad_p",grad_p,1)

      ! write omx, omy, omz
  !    call CalVorticity_tsjplane(idx_xzp,idim,kdim,buffer_grid,buffer_flow(nn,:,:,1,1),buffer_flow(nn,:,:,1,2),buffer_flow(nn,:,:,1,3), &
  !                      buffer_flow(nn,:,:,1,6),buffer_flow(nn,:,:,1,7),buffer_flow(nn,:,:,1,8), &
 !                       omx,omy,omz,ist_tsjplane,kst_tsjplane)
 !     call Write3DHDF5(fname,rank,dimsf,"omx",omx,1)
 !     call Write3DHDF5(fname,rank,dimsf,"omy",omy,1)
 !     call Write3DHDF5(fname,rank,dimsf,"omz",omz,1)

      ! write div
 !     call CalDivergence_tsjplane(idx_xzp,idim,kdim,buffer_grid,buffer_flow(nn,:,:,1,1),buffer_flow(nn,:,:,1,2),buffer_flow(nn,:,:,1,3), &
 !                       buffer_flow(nn,:,:,1,6),buffer_flow(nn,:,:,1,7),buffer_flow(nn,:,:,1,8), &
 !                       div,ist_tsjplane,kst_tsjplane)
 !     call Write3DHDF5(fname,rank,dimsf,"div",div,1)

        fname = 'tsjplane_'//fnum//'.plt'
        m = m +1
        call WriteTec(fname,x((idx_xzp%ist-ist_tsjplane+1):(idx_xzp%iend-idx_xzp%ist+1):idx_xzp%isp,1:1:1,(idx_xzp%kst-kst_tsjplane+1):(idx_xzp%kend-idx_xzp%kst+1):idx_xzp%ksp), &
                            z((idx_xzp%ist-ist_tsjplane+1):(idx_xzp%iend-idx_xzp%ist+1):idx_xzp%isp,1:1:1,(idx_xzp%kst-kst_tsjplane+1):(idx_xzp%kend-idx_xzp%kst+1):idx_xzp%ksp), &
                            grad_rho,m)
      enddo ! end nn loop
    enddo ! end n loop
  endif ! end read tsjplane


  if(iplane_yz.eq.1) then
    call ReadtsGridMetrics(fname_grid,"/iplane",jdim,kdim,idim,buffer_grid)
!    call ReadtsGridMetrics(fname_grid,"/iplane",kdim,jdim,idim,buffer_grid)
    do i=1,idim
      do j=1, jdim
        do k=1, kdim
          y(i,j,k) = buffer_grid(j,k,i,2)
          z(i,j,k) = buffer_grid(j,k,i,3)
        enddo
      enddo
    enddo

    ilen = 1
    jlen = (idx_yzp%jend - idx_yzp%jst)/idx_yzp%jsp + 1
    klen = (idx_yzp%kend - idx_yzp%kst)/idx_yzp%ksp + 1

    rank=3
    allocate(dimsf(rank))
    dimsf(1) = ilen
    dimsf(2) = jlen
    dimsf(3) = klen

    idx_yzp%ist=1
    idx_yzp%iend=1
    idx_yzp%isp=1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! init Tecplot file *.plt
    call InitTec(nt/nt_skip,ilen,jlen,klen,3,"y z grad_rho", 2)

    do n=iil,iih,stride
      write(unit=fnum,fmt='(I08.8)')n
      fname = trim(datapath)//'timeseries_'//fnum//'.h5'
      print *, 'reading file: ', trim(fname)
!      call Readtsflow_jplane(fname,nt,idim,kdim,10,buffer_flow(:,:,:,:,1:10)) !!!!!!
      call Readtsflow(fname,"/iplane",nt,jdim,kdim,1,buffer_flow(:,:,:,:,1:10),2)

      do nn=1, nt, nt_skip
        do j=1, jdim
          do k=1, kdim
            drdi(nn,j,k,1) = 1/287.0*(1.0/buffer_flow(nn,j,k,1,5)*buffer_flow(nn,j,k,1,9) - &
                            buffer_flow(nn,j,k,1,4)*buffer_flow(nn,j,k,1,10)/(buffer_flow(nn,j,k,1,5)**2)) ! drhodi
            buffer_flow(nn,j,k,1,11) = buffer_flow(nn,j,k,1,4)/dble(287.0*buffer_flow(nn,j,k,1,5)) ! rho
          enddo
        enddo
      enddo

      print *, 'min(rho) = ', minval(buffer_flow(:,:,:,:,11))

      m = 0
      do nn=1, nt, nt_skip
        write(unit=fnum1,fmt='(I04.4)') nn
        fname = 'tsiplane_'//fnum//'-'//fnum1//'.h5'
        print *, 'writing file: ', trim(fname)
      ! write x, z
  !    call Write3DHDF5(fname, rank, dimsf,"x",x,0)
  !    call Write3DHDF5(fname, rank, dimsf,"z",z,1)
      ! write grad_rho
        call CalGradient_tsiplane(idx_yzp,jdim,kdim,buffer_grid,buffer_flow(nn,:,:,1,11),drdi(nn,:,:,1),grad_rho,jst_tsiplane,kst_tsiplane)

     print *, ' grad_rho(1,1,1) = ', grad_rho(1,1,1)
     print *, ' grad_rho(1,2,1) = ', grad_rho(1,2,1)

        fname = 'tsiplane_'//fnum//'.plt'
        m = m +1
        call WriteTec(fname,y(1:1:1,(idx_yzp%jst-jst_tsiplane+1):(idx_yzp%jend-idx_yzp%jst+1):idx_yzp%jsp,(idx_yzp%kst-kst_tsiplane+1):(idx_yzp%kend-idx_yzp%kst+1):idx_yzp%ksp), &
                            z(1:1:1,(idx_yzp%jst-jst_tsiplane+1):(idx_yzp%jend-idx_yzp%jst+1):idx_yzp%jsp,(idx_yzp%kst-kst_tsiplane+1):(idx_yzp%kend-idx_yzp%kst+1):idx_yzp%ksp), &
                            grad_rho,m)
      enddo ! end nn loop
    enddo ! end n loop
  endif ! end read tsiplane






  call FinalizeHDF5()

 contains
    subroutine Input()
      integer ::i

      read(*,*)
      read(*,'(A)') datapath
      read(*,*)
      read(*,*) uinf, delta, utau, ztau
      read(*,*)
      read(*,*) iil, iih, stride
      read(*,*)
      read(*,*) iplane_yz, idx_yzp%jst, idx_yzp%jend, idx_yzp%jsp, idx_yzp%kst, idx_yzp%kend, idx_yzp%ksp
      read(*,*)
      read(*,*) jst_tsiplane, jend_tsiplane, jsp_tsiplane, kst_tsiplane, kend_tsiplane, ksp_tsiplane, nt, nt_skip
      read(*,*)
      read(*,*) iplane_xz, idx_xzp%ist, idx_xzp%iend, idx_xzp%isp, idx_xzp%kst, idx_xzp%kend, idx_xzp%ksp
      read(*,*)
      read(*,*) ist_tsjplane, iend_tsjplane, isp_tsjplane, kst_tsjplane, kend_tsjplane, ksp_tsjplane, nt, nt_skip

      if(iplane_yz.eq.1) then
        idim = 1
        jdim = jend_tsiplane-jst_tsiplane+1
        kdim = kend_tsiplane-kst_tsiplane+1
        allocate(buffer_flow(nt,jdim,kdim,idim,11))
        allocate(drdi(nt,jdim,kdim,idim))
        allocate(buffer_grid(jdim,kdim,idim,12))
      endif

      if(iplane_xz.eq.1) then
        idim = iend_tsjplane-ist_tsjplane+1
        jdim = 1
        kdim = kend_tsjplane-kst_tsjplane+1
        allocate(buffer_flow(nt,idim,kdim,jdim,11))
        allocate(drdj(nt,idim,kdim,jdim))
        allocate(buffer_grid(idim,kdim,jdim,12))
      endif


      allocate(x(idim,jdim,kdim),y(idim,jdim,kdim),z(idim,jdim,kdim))


    end subroutine Input


end program tsjplane
