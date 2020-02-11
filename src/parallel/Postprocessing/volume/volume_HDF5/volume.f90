program volume
  use decomp_2d
  use MFileIO
  use modVolume
  use modPRWHDF5

  implicit none

  character(400) :: fname, finfo
  character(200) :: datapath
  character(20) :: surffix
  character(8) :: fnum 
  character(4) :: fnum1
  integer :: imax, jmax, kmax
  integer :: iil, iih, stride
  integer :: ivolume, iplane_xz, iplane_yz, iplane_xy
  type(index_output) :: idx_vol, idx_xyp, idx_yzp, idx_xzp
  integer, dimension(:),  allocatable :: planes_xy, planes_yz, planes_xz
  integer, parameter :: ns=1
  real, dimension(:,:,:), allocatable :: x, y, z
  real, dimension(:,:,:), allocatable :: xx, yy, zz, xxh, yyh, zzh
  real, dimension(:,:,:,:), allocatable :: buffer_grid
  real, dimension(:,:,:,:), allocatable :: buffer_flow
  real, dimension(:,:,:), allocatable :: buffer_flowh     !!!!!!!
  real, dimension(:,:,:), allocatable :: buffer_flowh1,  buffer_flowh2, buffer_flowh3
  real, dimension(:,:,:), allocatable :: grad_p, grad_rho
  real, dimension(:,:,:), allocatable :: omx, omy, omz
  real, dimension(:,:,:), allocatable :: swirl
  real, dimension(:,:,:), allocatable :: div
  real, dimension(:,:,:), allocatable :: Q
  real, dimension(:,:,:), allocatable :: dudn, dtdn, ds
  integer :: n, i, j, k, m, nn
  real(8) :: uinf, delta, utau, ztau
!  integer, dimension(:), allocatable :: varindex_vol_output, varindex_xyp_output, varindex_yzp_output, varindex_xzp_output
  integer, dimension(17) :: varindex_vol_output, varindex_xyp_output, varindex_yzp_output, varindex_xzp_output
!  character(10), dimension(:), allocatable :: dname_vol_output, dname_xyp_output, dname_yzp_output, dname_xzp_output
  character(10), dimension(17) :: dname_vol_output, dname_xyp_output, dname_yzp_output, dname_xzp_output

  integer :: nvar, nvar_vol_output, nvar_xyp_output, nvar_yzp_output, nvar_xzp_output
  parameter(nvar=17)
  character(10) :: dname(nvar)
  parameter(dname = (/'x','y','z','u','v','w','p','t','rho','grad_p','grad_rho','omx','omy','omz','div','swil','Q'/))
  character(200) :: varname_output
  integer, dimension(:), allocatable :: dims, dims_total
  integer :: rank
  integer :: ilen, jlen, klen
  real(8), dimension(:,:,:,:), allocatable :: buffer_tec

  integer :: p_row, p_col
  integer :: ierror

  ! initialize MPI
  call MPI_INIT(ierror)
  call InitHDF5()

  !! dimension information
!  imax = 400; jmax = 300; kmax = 130
!  p_row = 2;  p_col = 1

  call Input()

!! *****************************************
  !! test output
   if(nrank.eq.0) then
     print *, 'nrank = ', nrank
     print *, 'xsize(1) = ', xsize(1)
     print *, 'xsize(2) = ', xsize(2)
     print *, 'xsize(3) = ', xsize(3)
     print *, 'xstart(1) = ', xstart(1)
     print *, 'xstart(2) = ', xstart(2)
     print *, 'xstart(3) = ', xstart(3)
     print *, '***************************'
   endif
!! *****************************************

  if(nrank.eq.0) then
    print *, 'reading grid file'
    print *, '***************************'
  endif
  ! read grid file
  fname = trim(datapath)//'grid.h5'
  call ReadHDF5grid_P(fname,xsize(2),xsize(3),xsize(1),xstart(2),xstart(3),xstart(1),xx,yy,zz)

  if(nrank.eq.0) then
    print *, 'finish reading grid file'
    print *, '***************************'
  endif

  ! read flowdata files
  do n=iil, iih, stride

    write(unit=fnum,fmt='(I08.8)') n
    fname = trim(datapath)//'flowdata_'//fnum//'.h5'
    if(nrank.eq.0)  print *, 'Reading file: ', trim(fname)

    call ReadHDF5sol_P(fname,xsize(2),xsize(3),xsize(1),xstart(2),xstart(3),xstart(1),buffer_flow(:,:,:,1:5))
    buffer_flow(:,:,:,6) = buffer_flow(:,:,:,4)/(287.0*buffer_flow(:,:,:,5))

    ! check info
   ! if(nrank.eq.0) then
   !   print *, 'nrank =  ', nrank
   !   print *, 'xx(1,2,1) = ', xx(1,2,1)
   !   print *, 'xx(1,2,63) = ', xx(1,2,63)
   !   print *, 'xx(1,2,65) = ', xx(1,2,65)
   !   print *, 'yy(1,1,2) = ', yy(1,1,2)
   !   print *, 'zz(2,1,1) = ', zz(2,1,1)
   ! endif

    if (ivolume.gt.0) then

      rank = 3
      allocate( dims(rank), dims_total(rank) )
      dims(1) = xsize(1)
      dims(2) = xsize(2)
      dims(3) = xsize(3)
      dims_total(1) = kmax
      dims_total(2) = imax
      dims_total(3) = jmax

      ! check info
      ! print *, 'dims = ', dims
      ! print *, 'dims_total = ', dims_total

      do nn=1, nvar_vol_output
        dname_vol_output(nn) = dname(varindex_vol_output(nn))
      enddo
      if(nrank.eq.0) print *, 'dname_output = ' , (trim(dname_vol_output(nn))//',', nn=1, nvar_vol_output)

      fname = 'volume_'//fnum//'.h5'

      ! call Write3dHDF5_svariable_P(fname,dims,dims_total,xstart(1),xstart(2),xstart(3),"x",xx,0)
      ! call Write3dHDF5_svariable_P(fname,dims,dims_total,xstart(1),xstart(2),xstart(3),"y",yy,1)
      ! call Write3dHDF5_svariable_P(fname,dims,dims_total,xstart(1),xstart(2),xstart(3),"z",zz,1)

      call update_halo(xx,xxh,2)
      call update_halo(yy,yyh,2)
      call update_halo(zz,zzh,2)

      do nn=1, nvar_vol_output
        ! write grad_p
        if(varindex_vol_output(nn).eq.10) then
        ! call update_halo( buffer_flow(:,:,:,4),buffer_flowh,2 )
        ! call CalGradient(idx_vol,xsize(2)+4,xsize(3)+4,xsize(1),xxh,yyh,zzh,buffer_flowh(1:xsize(1),-1:xsize(2)+2,-1:xsize(3)+2),grad_p)
        ! call Write3dHDF5_svariable_P(fname,dims,dims_total,xstart(1),xstart(2),xstart(3),"grad_p",grad_p,0) !!!
        ! deallocate(buffer_flowh)
        endif
        ! write grad_rho
        if(varindex_vol_output(nn).eq.11) then
         ! call update_halo( buffer_flow(:,:,:,6),buffer_flowh,2 )
          ! check the size for idim, jdim, kdim
         ! call CalGradient(idx_vol,xsize(2)+4,xsize(3)+4,xsize(1),xxh,yyh,zzh,buffer_flowh(1:xsize(1),-1:xsize(2)+2,-1:xsize(3)+2),grad_rho)
         ! call Write3dHDF5_svariable_P(fname,dims,dims_total,xstart(1),xstart(2),xstart(3),"grad_rho",grad_rho,0) !!!
         ! deallocate(buffer_flowh)
        endif

        call update_halo( buffer_flow(:,:,:,1),buffer_flowh1,2 )
        call update_halo( buffer_flow(:,:,:,2),buffer_flowh2,2 )
        call update_halo( buffer_flow(:,:,:,3),buffer_flowh3,2 )

        ! write omx, omy, omz
        if(varindex_vol_output(nn).eq.12) then

         ! call CalVorticity(idx_vol,xsize(2)+4,xsize(3)+4,xsize(1),xxh,yyh,zzh,buffer_flowh1(1:xsize(1),-1:xsize(2)+2,-1:xsize(3)+2), &
         !                   buffer_flowh2(1:xsize(1),-1:xsize(2)+2,-1:xsize(3)+2),buffer_flowh3(1:xsize(1),-1:xsize(2)+2,-1:xsize(3)+2),omx,omy,omz)
         ! call Write3dHDF5_svariable_P(fname,dims,dims_total,xstart(1),xstart(2),xstart(3),"omx",omx,0)
         ! call Write3dHDF5_svariable_P(fname,dims,dims_total,xstart(1),xstart(2),xstart(3),"omy",omy,0)
         ! call Write3dHDF5_svariable_P(fname,dims,dims_total,xstart(1),xstart(2),xstart(3),"omz",omz,0)


        endif

        ! write div
        if(varindex_vol_output(nn).eq.15) then
         ! call CalDivergence(idx_vol,xsize(2)+4,xsize(3)+4,xsize(1),xxh,yyh,zzh,buffer_flowh1(1:xsize(1),-1:xsize(2)+2,-1:xsize(3)+2), &
         !                    buffer_flowh2(1:xsize(1),-1:xsize(2)+2,-1:xsize(3)+2),buffer_flowh3(1:xsize(1),-1:xsize(2)+2,-1:xsize(3)+2),div)
         ! call Write3dHDF5_svariable_P(fname,dims,dims_total,xstart(1),xstart(2),xstart(3),"div",div,0)
        endif

        ! write swirl
        if(varindex_vol_output(nn).eq.15) then
          call CalSwirl(idx_vol,xsize(2)+4,xsize(3)+4,xsize(1),xxh,yyh,zzh,buffer_flowh1(1:xsize(1),-1:xsize(2)+2,-1:xsize(3)+2), &
                            buffer_flowh2(1:xsize(1),-1:xsize(2)+2,-1:xsize(3)+2),buffer_flowh3(1:xsize(1),-1:xsize(2)+2,-1:xsize(3)+2),swirl)
          call Write3dHDF5_svariable_P(fname,dims,dims_total,xstart(1),xstart(2),xstart(3),"swirl",swirl,0)

        endif


        deallocate(buffer_flowh1)
        deallocate(buffer_flowh2)
        deallocate(buffer_flowh3)

      enddo ! end nn loop

    endif ! end ivolume.gt.0

  enddo ! end n loop


  call FinalizeHDF5()
  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)

  contains
    subroutine Input()
      integer :: i

     call MPI_COMM_RANK(MPI_COMM_WORLD, nrank, ierror)

     if(nrank.eq.0) then
       read(*,*)
       read(*,*) imax, jmax, kmax, p_row, p_col
       print *, 'imax = ', imax, 'jmax = ', jmax, 'kmax = ', kmax
       print *,  'node_1 = ', p_row, 'node_2 = ', p_col
     endif
     call MPI_Bcast(imax,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
     call MPI_Bcast(jmax,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
     call MPI_Bcast(kmax,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
     call MPI_Bcast(p_row, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
     call MPI_Bcast(p_col, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

     ! init
     call decomp_2d_init(kmax,imax,jmax,p_row,p_col)
     call MPI_Barrier(MPI_COMM_WORLD,ierror)

      if(nrank.eq.0) then
        read(*,*)
        read(*,'(A)')datapath  !path where data files are stored
        read(*,*)
        read(*,*)uinf,delta,utau,ztau
        read(*,*)
        read(*,*)iil,iih,stride
        read(*,*)
        read(*,*)ivolume, idx_vol%ist, idx_vol%iend, idx_vol%isp, idx_vol%jst, idx_vol%jend, idx_vol%jsp, idx_vol%kst, idx_vol%kend, idx_vol%ksp
        read(*,*)
        read(*,*) nvar_vol_output
        ! allocate(varindex_vol_output(nvar_vol_output), dname_vol_output(nvar_vol_output))
        read(*,*) (varindex_vol_output(i), i=1,nvar_vol_output)
        ! xy planes
        read(*,*)
        read(*,*)iplane_xy,idx_xyp%ist, idx_xyp%iend, idx_xyp%isp, idx_xyp%jst, idx_xyp%jend, idx_xyp%jsp


        if (iplane_xy.gt.0) then
          allocate(planes_xy(1:iplane_xy))
          read(*,*)(planes_xy(i),i=1,iplane_xy)
        else
          read(*,*)
        end if
        read(*,*)
        read(*,*) nvar_xyp_output
        ! allocate(varindex_xyp_output(nvar_xyp_output), dname_xyp_output(nvar_xyp_output))
        read(*,*) (varindex_xyp_output(i), i=1,nvar_xyp_output)

        ! yz planes
        read(*,*)
        read(*,*)iplane_yz, idx_yzp%jst, idx_yzp%jend, idx_yzp%jsp, idx_yzp%kst, idx_yzp%kend, idx_yzp%ksp
        if (iplane_yz.gt.0) then
          allocate(planes_yz(1:iplane_yz))
          read(*,*)(planes_yz(i),i=1,iplane_yz)
        else
          read(*,*)
        end if
        read(*,*)
        read(*,*) nvar_yzp_output
        ! allocate(varindex_yzp_output(nvar_yzp_output), dname_yzp_output(nvar_yzp_output))
        read(*,*) (varindex_yzp_output(i), i=1,nvar_yzp_output)
        ! xz planes
        read(*,*)
        read(*,*)iplane_xz, idx_xzp%ist, idx_xzp%iend, idx_xzp%isp, idx_xzp%kst, idx_xzp%kend, idx_xzp%ksp
        if (iplane_xz.gt.0) then
          allocate(planes_xz(1:iplane_xz))
          read(*,*)(planes_xz(i),i=1,iplane_xz)
        else
          read(*,*)
        end if
        read(*,*)
        read(*,*) nvar_xzp_output
        ! allocate(varindex_xzp_output(nvar_xzp_output), dname_xzp_output(nvar_xzp_output))
        read(*,*) (varindex_xzp_output(i), i=1,nvar_xzp_output)

        ! call ReadPlot3DGridPlane(trim(datapath)//'gridp3d.grd',imax,jmax,kmax,x,y,z)

!! ****************************************
!        ! check
!        print *, 'imax = ', imax
!        print *, 'jmax = ', jmax
!        print *, 'kmax = ', kmax
!        print *, 'x(2,63,1) = ', x(2,63,1)
!        print *, 'x(2,65,1) = ', x(2,65,1)
!! ****************************************

      endif ! nrank=0

      call MPI_Bcast(datapath, 200, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(uinf,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(delta,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(utau,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(ztau,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(iil,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(iih,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(stride, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(ivolume,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(idx_vol%ist,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(idx_vol%iend, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(idx_vol%isp,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(idx_vol%jst,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(idx_vol%jend, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(idx_vol%jsp,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(idx_vol%kst,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(idx_vol%kend, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(idx_vol%ksp,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(nvar_vol_output,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
      call MPI_Bcast(varindex_vol_output,  nvar_vol_output, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)


      if(nrank.eq.0) then
        print *, 'iil = ', iil
        print *, 'iih = ', iih
        print *, 'stride = ', stride
      endif


      ! dimension order k,i,j
      allocate( buffer_flow(xsize(1),xsize(2),xsize(3),6) )
      allocate( xx(xsize(1),xsize(2),xsize(3)), yy(xsize(1),xsize(2),xsize(3)), zz(xsize(1),xsize(2),xsize(3)) )

      return
    end subroutine Input


end program volume    
