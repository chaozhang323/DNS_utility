

program tsjplane
  use MFileIO
  use modVolume
  use modReadHDF5
  use modWriteHDF5
  implicit none
  character(400) :: fname
  character(400) :: fname_grid
  character(200) :: datapath
  character(8) :: fnum
  character(4) :: fnum1
  integer :: imax, jmax, kmax
  integer :: iil, iih, stride
  integer :: ivolume, iplane_xz, iplane_yz, iplane_xy
  type(index_output) :: idx_vol, idx_xyp, idx_yzp, idx_xzp
  integer, dimension(:),  allocatable :: planes_xy, planes_yz, planes_xz
  real, dimension(:,:,:), allocatable :: x, y, z
  real, dimension(:,:,:,:), allocatable :: buffer_grid
  real, dimension(:,:,:,:), allocatable :: buffer_flow
  real, dimension(:,:,:,:,:), allocatable :: buffer_flow_tmp
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

  real(8), dimension(:,:,:,:), allocatable :: drdj

  call Input()
  call InitHDF5()

  fname_grid = trim(datapath)//'timeseries_GridMetrics.h5'
  print *, 'reading grid: ', trim(fname_grid)
  call ReadtsGrid_jplane(fname_grid,imax,kmax,12,buffer_grid)

  do i=1,imax
    do j=1, jmax
      do k=1, kmax
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



  do n=iil,iih,stride
    write(unit=fnum,fmt='(I08.8)')n
    fname = trim(datapath)//'timeseries_'//fnum//'.h5'
    print *, 'reading file: ', trim(fname)
    call Readtsflow_jplane(fname,nt,iend_tsjplane-ist_tsjplane+1,kend_tsjplane-kst_tsjplane+1,10,buffer_flow_tmp(:,:,:,:,1:10))

    do nn=1, nt, nt_skip
      do i=1, imax
        do k=1, kmax
          drdj(nn,i,k,1) = 1/287.0*(1.0/buffer_flow_tmp(nn,i,k,1,5)*buffer_flow_tmp(nn,i,k,1,9) - &
                          buffer_flow_tmp(nn,i,k,1,4)*buffer_flow_tmp(nn,i,k,1,10)/(buffer_flow_tmp(nn,i,k,1,5)**2)) ! drhodj
          buffer_flow_tmp(nn,i,k,1,11) = buffer_flow_tmp(nn,i,k,1,4)/dble(287.0*buffer_flow_tmp(nn,i,k,1,5)) ! rho
        enddo
      enddo
    enddo

    print *, 'min(rho) = ', minval(buffer_flow_tmp(:,:,:,:,11))

    do nn=1, nt, nt_skip
      write(unit=fnum1,fmt='(I04.4)') nn
      fname = 'tsjplane_'//fnum//'-'//fnum1//'.h5'
      print *, 'writing file: ', trim(fname)
      ! write x
      call Write3DHDF5(fname, rank, dimsf,"x",x,0)
      call Write3DHDF5(fname, rank, dimsf,"z",z,1)
      ! write grad_rho
      call CalGradient_tsjplane(idx_xzp,imax,kmax,buffer_grid,buffer_flow_tmp(nn,:,:,1,11),drdj(nn,:,:,1),grad_rho,ist_tsjplane,kst_tsjplane)
      call Write3DHDF5(fname,rank,dimsf,"grad_rho",grad_rho,1)


 !  print *, 'buffer_flow_tmp(nn,1,1,1,9) = ', buffer_flow_tmp(nn,1,1,1,9)

      ! write grad_p
      call CalGradient_tsjplane(idx_xzp,imax,kmax,buffer_grid,buffer_flow_tmp(nn,:,:,1,4),buffer_flow_tmp(nn,:,:,1,9),grad_p,ist_tsjplane,kst_tsjplane)
      call Write3DHDF5(fname,rank,dimsf,"grad_p",grad_p,1)

      ! write omx, omy, omz
      call CalVorticity_tsjplane(idx_xzp,imax,kmax,buffer_grid,buffer_flow_tmp(nn,:,:,1,1),buffer_flow_tmp(nn,:,:,1,2),buffer_flow_tmp(nn,:,:,1,3), &
                        buffer_flow_tmp(nn,:,:,1,6),buffer_flow_tmp(nn,:,:,1,7),buffer_flow_tmp(nn,:,:,1,8), &
                        omx,omy,omz,ist_tsjplane,kst_tsjplane)
      call Write3DHDF5(fname,rank,dimsf,"omx",omx,1)
      call Write3DHDF5(fname,rank,dimsf,"omy",omy,1)
      call Write3DHDF5(fname,rank,dimsf,"omz",omz,1)

      ! write div
      call CalDivergence_tsjplane(idx_xzp,imax,kmax,buffer_grid,buffer_flow_tmp(nn,:,:,1,1),buffer_flow_tmp(nn,:,:,1,2),buffer_flow_tmp(nn,:,:,1,3), &
                        buffer_flow_tmp(nn,:,:,1,6),buffer_flow_tmp(nn,:,:,1,7),buffer_flow_tmp(nn,:,:,1,8), &
                        div,ist_tsjplane,kst_tsjplane)
      call Write3DHDF5(fname,rank,dimsf,"div",div,1)


    enddo ! end nn loop

  enddo ! end n loop

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
      read(*,*) iplane_xz, idx_xzp%ist, idx_xzp%iend, idx_xzp%isp, idx_xzp%kst, idx_xzp%kend, idx_xzp%ksp
      read(*,*)
      read(*,*) ist_tsjplane, iend_tsjplane, isp_tsjplane, kst_tsjplane, kend_tsjplane, ksp_tsjplane, nt, nt_skip

      imax = iend_tsjplane-ist_tsjplane+1
      jmax = 1
      kmax = kend_tsjplane-kst_tsjplane+1

      allocate(buffer_grid(imax,kmax,1,12))
      allocate(x(imax,1,kmax),y(imax,1,kmax),z(imax,1,kmax))
      allocate(buffer_flow_tmp(nt,imax,kmax,1,11))
      allocate(drdj(nt,imax,kmax,1))

    end subroutine Input


end program tsjplane
