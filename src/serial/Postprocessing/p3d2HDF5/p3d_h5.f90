!>
!> \file p3d_h5.f90
!> @brief
!> @details

PROGRAM p3d_h5
  use modRWHDF5
  use MFileIO
  implicit none

  integer ::   i, j, k, kk, k2, kp, n, nn, imax, jmax, kmax, ires, icurrent
  integer, dimension(:), allocatable :: isize,jsize
  real(8), parameter :: Univ = 8314.3D0
  real(8), dimension(:,:,:), allocatable :: X, Y, Z
  real(8), dimension(:,:,:,:), allocatable :: vars, buffer_grid, buffer_flow, buffer_flow_tmp
  real(8), dimension(:,:,:), allocatable :: vars2d
  real(8) ::time, rbar
  integer :: iConvert_grid, iConvert_flow, iConvert_inlet
  integer :: iFormat_grid, iFormat_flow, iFormat_inlet
  character(300) :: gridfile_rd, gridfile_wt, flowpath_rd, flowpath_wt
  integer :: file_be, file_end, file_skip, num_file
  character(400) :: fname
  character(8) :: fnum
  character(4) :: fnum4_1, fnum4_2
  type(tp_rdwt_hdf5) :: grd, fsol, fsolRescale, grd1V, fsol1V
  integer :: irhoinlet  ! 'inlet.sol' varialbe list (=1: u, v, w, rho, T)
  real(8) :: Mw
  integer :: ireadRescalemean, istencilsize, kioplane, nblks, kdim, idim, nvar, inode, jnode
  real(8) :: theta_r
  integer :: nsample_rescale
  character(5) :: gridname(3) = (/'x','y','z'/)
  character(5) :: flowdataname(5) = (/'u','v','w','p','T'/)

  call InitHDF5()
  call Initial()
  rbar = Univ/Mw
  print *, 'Gas constant =', rbar

  if(iConvert_grid.eq.1) then ! Convert grid file
    call InitGridHDF5(grd)
!    grd%dimsf = (/kmax,imax,jmax/)
!    allocate(buffer_grid(kmax,imax,jmax,3))

    if(iFormat_grid.eq.0) then ! Convert plot3d grid file to HDF5
      grd%dimsf = (/kmax,imax,jmax/)
      allocate(buffer_grid(kmax,imax,jmax,3))
      print *, 'reading file ... ', trim(gridfile_rd)
      allocate(X(imax,jmax,kmax),Y(imax,jmax,kmax),Z(imax,jmax,kmax) )
      call ReadPlot3DGridPlane(trim(gridfile_rd),imax,jmax,kmax,X,Y,Z)
      do i=1, imax
        do j=1, jmax
          do k=1, kmax
            buffer_grid(k,i,j,1) = X(i,j,k)
            buffer_grid(k,i,j,2) = Y(i,j,k)
            buffer_grid(k,i,j,3) = Z(i,j,k)
          enddo
        enddo
      enddo
      grd%fname = trim(gridfile_wt)
      print *, 'writing file ... ', trim(grd%fname)
      call WriteHDF5_3D(grd, buffer_grid)
      deallocate(X,Y,Z)
    elseif(iFormat_grid.eq.1) then ! transfer HDF5 grid file to plot3d
      grd%dimsf = (/kmax,imax,jmax/)
      allocate(buffer_grid(kmax,imax,jmax,3))
      grd%fname = trim(gridfile_rd)
      print *, 'reading file ... ', trim(grd%fname)
      call ReadHDF5_3D(grd,buffer_grid)
      allocate(X(imax,jmax,kmax),Y(imax,jmax,kmax),Z(imax,jmax,kmax) )
      do i=1, imax
        do j=1, jmax
          do k=1, kmax
            X(i,j,k) = buffer_grid(k,i,j,1)
            Y(i,j,k) = buffer_grid(k,i,j,2)
            Z(i,j,k) = buffer_grid(k,i,j,3)
          enddo
        enddo
      enddo
      print *, 'writing file ... ', trim(gridfile_wt)
      call WritePlot3DGridPlane(trim(gridfile_wt),imax,jmax,kmax,X,Y,Z)
      deallocate(X,Y,Z)
    elseif(iFormat_grid.eq.2) then ! separate single HDF5 file into multiple files
      grd%dimsf = (/kmax,imax,jmax/)
      allocate(buffer_grid(kmax,imax,jmax,3))
      grd%fname = trim(gridfile_rd)
      print *, 'reading file ... ', trim(grd%fname)
      call InitHDF5_1V(grd1V)
      grd1V%dimsf = (/kmax,imax,jmax/)
      do n=1, 3
        call ReadHDF5_3D_1V(grd,buffer_grid(:,:,:,n),n)
        grd1V%fname = 'grid_'//trim(gridname(n))//'.h5'
        grd1V%dname = trim(gridname(n))
        print *, 'writing file ... ', trim(grd1V%fname)
        call WriteHDF5_3D(grd1V, buffer_grid(:,:,:,n:n))
      enddo

    elseif(iFormat_grid.eq.3) then
      do kk=1, kmax, kioplane
        k2 = min(kk+kioplane-1,kmax)
        kp = k2-kk+1
        grd%dimsf = (/kp,imax,jmax/)
        grd%offset = (/kk-1,0,0/)
        if(kk.eq.1) allocate(buffer_grid(kp,imax,jmax,3))
        grd%fname = trim(gridfile_rd)
        if(kk.eq.1) print *, 'reading file ... ', trim(grd%fname)
        call ReadHDF5_3D(grd,buffer_grid(1:kp,1:imax,1:jmax,1:3))

        grd%offset = 0
        write(unit=fnum4_1,fmt='(I04.4)') kk
        write(unit=fnum4_2,fmt='(I04.4)') k2
        grd%fname = 'grid_kplane'//fnum4_1//'_'//fnum4_2//'.h5'
        print *, 'writing file ... ', trim(grd%fname)
        call WriteHDF5_3D(grd, buffer_grid(1:kp,1:imax,1:jmax,1:3))
      enddo ! end kk loop

    else
      print *, 'unknown grid format. Stop!'
      stop
    endif
    deallocate(buffer_grid)
  endif

  if(iConvert_flow.eq.1) then ! Convert flow file
    call InitFlowHDF5(fsol)
    num_file = (file_end - file_be)/file_skip + 1

    if(iFormat_flow.eq.0) then ! convert plot3d solution file to HDF5
      fsol%sname = 'time'
      fsol%dimsf = (/kmax,imax,jmax/)
      allocate(vars(imax,jmax,kmax,6))
      allocate(buffer_flow(kmax,imax,jmax,5))
      do n=1, num_file
        write(unit=fnum,fmt='(I08.8)') (file_be + (n-1)*file_skip)
        fname = trim(flowpath_rd)//'flowdata_'//fnum//'.sol'
        print *, 'reading file ', trim(fname)
        call ReadPlot3DSolPlane(trim(fname),imax,jmax,kmax,6,vars,time)  ! Read plot3d sol file

        do i=1, imax
          do j=1, jmax
            do k=1, kmax
              buffer_flow(k,i,j,1) = vars(i,j,k,1)
              buffer_flow(k,i,j,2) = vars(i,j,k,2)
              buffer_flow(k,i,j,3) = vars(i,j,k,3)
              buffer_flow(k,i,j,4) = vars(i,j,k,4)
              buffer_flow(k,i,j,5) = vars(i,j,k,5)
            end do
          end do
        end do

        fsol%fname = trim(flowpath_wt)//'flowdata_'//fnum//'.h5'
        print *, 'writing file ', trim(fsol%fname)
        call WriteHDF5_3D(fsol,buffer_flow)
        call WriteHDF5_scalar(fsol,time)
      enddo ! end n loop

    elseif(iFormat_flow.eq.1) then  ! convert HDF5 solution file to plot3d
      fsol%sname = 'time'
      fsol%dimsf = (/kmax,imax,jmax/)
      allocate(vars(imax,jmax,kmax,6))
      allocate(buffer_flow(kmax,imax,jmax,5))
      do n=1, num_file
        write(unit=fnum,fmt='(I08.8)') (file_be + (n-1)*file_skip)
        fsol%fname = trim(flowpath_rd)//'flowdata_'//fnum//'.h5'
        print *, 'reading file ', trim(fsol%fname)
        call ReadHDF5_3D(fsol,buffer_flow)
        call ReadHDF5_scalar(fsol,time)

        do i=1, imax
          do j=1, jmax
            do k=1, kmax
              vars(i,j,k,1) = buffer_flow(k,i,j,1)
              vars(i,j,k,2) = buffer_flow(k,i,j,2)
              vars(i,j,k,3) = buffer_flow(k,i,j,3)
              vars(i,j,k,4) = buffer_flow(k,i,j,4)
              vars(i,j,k,5) = buffer_flow(k,i,j,5)
            end do
          end do
        end do
        vars(:,:,:,6) = vars(:,:,:,4)/(rbar*vars(:,:,:,5))

        fname = trim(flowpath_wt)//'flowdata_'//fnum//'.sol'
        print *, 'writing file ', trim(fname)
        call WritePlot3DSolPlane(trim(fname),imax,jmax,kmax,6,vars,time)
      enddo ! end n loop

    elseif(iFormat_flow.eq.2) then
      fsol%sname = 'time'
      fsol%dimsf = (/kmax,imax,jmax/)
      allocate(buffer_flow(kmax,imax,jmax,5))
      call InitHDF5_1V(fsol1V)
      fsol1V%dimsf = (/kmax,imax,jmax/)
      do n=1, num_file
        write(unit=fnum,fmt='(I08.8)') (file_be + (n-1)*file_skip)
        fsol%fname = trim(flowpath_rd)//'flowdata_'//fnum//'.h5'
        print *, 'reading file ', trim(fsol%fname)
        do nn=1, 5
          call ReadHDF5_3D_1V(fsol,buffer_flow(:,:,:,nn),nn)
          if(nn.eq.1) call ReadHDF5_scalar(fsol,time)
          fsol1V%fname = trim(flowpath_wt)//'flowdata_'//trim(flowdataname(nn))//'_'//fnum//'.h5'
          fsol1V%dname = trim(flowdataname(nn))
          print *, 'writing file ... ', trim(fsol1V%fname)
          call WriteHDF5_3D(fsol1V, buffer_flow(:,:,:,nn:nn))
          call WriteHDF5_scalar(fsol1V,time)
        enddo
      enddo ! end n loop

    elseif(iFormat_flow.eq.3) then
      do n=1, num_file
        write(unit=fnum,fmt='(I08.8)') (file_be + (n-1)*file_skip)
        do kk=1, kmax, kioplane
          k2 = min(kk+kioplane-1,kmax)
          kp = k2-kk+1
          fsol%dimsf = (/kp,imax,jmax/)
          fsol%offset = (/kk-1,0,0/)
          if(kk.eq.1.and.n.eq.1) allocate(buffer_flow(kp,imax,jmax,5))
          fsol%fname = trim(flowpath_rd)//'flowdata_'//fnum//'.h5'
          if(kk.eq.1) print *, 'reading file ', trim(fsol%fname)
          call ReadHDF5_3D(fsol,buffer_flow(1:kp,1:imax,1:jmax,1:5))
          if(kk.eq.1) call ReadHDF5_scalar(fsol,time)

          fsol%offset = 0
          write(unit=fnum4_1,fmt='(I04.4)') kk
          write(unit=fnum4_2,fmt='(I04.4)') k2
          fsol%fname = 'flowdata_'//fnum//'_kplane'//fnum4_1//'_'//fnum4_2//'.h5'
          print *, 'writing file ... ', trim(fsol%fname)
          call WriteHDF5_3D(fsol, buffer_flow(1:kp,1:imax,1:jmax,1:5))
          call WriteHDF5_scalar(fsol,time)
        enddo ! end kk loop
      enddo ! end n loop
    elseif(iFormat_flow.eq.4) then      ! separate single HDF5 flowdata into multiple files using inode
      allocate(isize(0:inode-1))
      i=imax/inode
      ires = imax - i*inode
      do n=0, inode-1
        isize(n) = i
        if(n.lt.ires) isize(n) = i+1
      enddo

      do n=1, num_file
        write(unit=fnum,fmt='(I08.8)') (file_be + (n-1)*file_skip)
        icurrent = 0
        do nn=0, inode-1
          icurrent = icurrent + isize(nn)
          fsol%dimsf = (/kmax,isize(nn),jmax/)
          fsol%offset = (/0,icurrent-isize(nn),0/)
          if(nn.eq.0.and.n.eq.1) allocate(buffer_flow(kmax,isize(nn),jmax,5))
          fsol%fname = trim(flowpath_rd)//'flowdata_'//fnum//'.h5'
          if(nn.eq.0) print *, 'reading file ', trim(fsol%fname)
          call ReadHDF5_3D(fsol,buffer_flow(1:kmax,1:isize(nn),1:jmax,1:5))
          if(nn.eq.0) call ReadHDF5_scalar(fsol,time)

          fsol%offset = 0
          write(unit=fnum4_1,fmt='(I04.4)') icurrent-isize(nn)+1
          write(unit=fnum4_2,fmt='(I04.4)') icurrent
          fsol%fname = 'flowdata_i'//fnum4_1//'_'//fnum4_2//'_'//fnum//'.h5'
          print *, 'writing file ... ', trim(fsol%fname)
          call WriteHDF5_3D(fsol, buffer_flow(1:kmax,1:isize(nn),1:jmax,1:5))
          call WriteHDF5_scalar(fsol,time)
        enddo ! end nn loop
      enddo ! end n loop
      deallocate(isize)
    elseif(iFormat_flow.eq.5) then      ! separate single HDF5 flowdata into multiple files using jnode
      allocate(jsize(0:jnode-1))
      i=jmax/jnode
      ires = jmax - i*jnode
      do n=0, jnode-1
        jsize(n) = i
        if(n.lt.ires) jsize(n) = i+1
      enddo

      do n=1, num_file
        write(unit=fnum,fmt='(I08.8)') (file_be + (n-1)*file_skip)
        icurrent = 0
        do nn=0, jnode-1
          icurrent = icurrent + jsize(nn)
          fsol%dimsf = (/kmax,imax,jsize(nn)/)
          fsol%offset = (/0,0,icurrent-jsize(nn)/)
          if(nn.eq.0.and.n.eq.1) allocate(buffer_flow(kmax,imax,jsize(nn),5))
          fsol%fname = trim(flowpath_rd)//'flowdata_'//fnum//'.h5'
          if(nn.eq.0) print *, 'reading file ', trim(fsol%fname)
          call ReadHDF5_3D(fsol,buffer_flow(1:kmax,1:imax,1:jsize(nn),1:5))
          if(nn.eq.0) call ReadHDF5_scalar(fsol,time)

          fsol%offset = 0
          write(unit=fnum4_1,fmt='(I04.4)') icurrent-jsize(nn)+1
          write(unit=fnum4_2,fmt='(I04.4)') icurrent
          fsol%fname = 'flowdata_j'//fnum4_1//'_'//fnum4_2//'_'//fnum//'.h5'
          print *, 'writing file ... ', trim(fsol%fname)
          call WriteHDF5_3D(fsol, buffer_flow(1:kmax,1:imax,1:jsize(nn),1:5))
          call WriteHDF5_scalar(fsol,time)
        enddo ! end nn loop
      enddo ! end n loop
      deallocate(jsize)
    else
      print *, 'unknown flow format. Stop!'
      stop
    endif
    deallocate(buffer_flow)

  endif  ! end if iConvert_flow

  if(iConvert_inlet.gt.0) then
    call InitFlowHDF5(fsol)
!    fsol%dimsf = (/kmax,1,jmax/)
    allocate(vars2d(kmax,jmax,5))
    allocate(buffer_flow(kmax,1,jmax,5))

    if(iFormat_inlet.eq.0) then     ! convert plot3d file to HDF5
      
      select case(irhoinlet)
      case(0) ! Take 1st i-plane from flowdata file as inlet
        write(unit=fnum,fmt='(I08.8)') file_be
        fsol%dimsf = (/kmax,imax,jmax/)
        fsol%fname = trim(flowpath_wt)//'flowdata_'//fnum//'.h5'
        print *, 'reading file for generating inlet: ', trim(fsol%fname)
        allocate(buffer_flow_tmp(kmax,imax,jmax,5))
        call ReadHDF5_3D(fsol,buffer_flow_tmp)      
!        forall(k=1:kmax, j=1:jmax, n = 1:5)
!            buffer_flow(k,1,j,n) = buffer_flow_tmp(k,1,j,n)
!        end forall   
        do i = 1, 1
           buffer_flow(:,i,:,:) = buffer_flow_tmp(:,i,:,:)
        enddo  
        deallocate(buffer_flow_tmp)
      case(1) ! Read Plot2D 'inlet.sol' with u, v, w, p, T       
         fname = trim(flowpath_rd)//'inlet.sol'
         call ReadPlot2DSolution(trim(fname),kmax,jmax,5,vars2d)
         do i=1, 1
            do j=1, jmax
               do k=1, kmax
                  buffer_flow(k,i,j,1) = vars2d(k,j,1)
                  buffer_flow(k,i,j,2) = vars2d(k,j,2)
                  buffer_flow(k,i,j,3) = vars2d(k,j,3)
                  buffer_flow(k,i,j,4) = vars2d(k,j,4)
                  buffer_flow(k,i,j,5) = vars2d(k,j,5)
               enddo
            enddo
         enddo
      case(2)  ! Read Plot2D 'inlet.sol' with u, v, w, rho, T
         fname = trim(flowpath_rd)//'inlet.sol'
         call ReadPlot2DSolution(trim(fname),kmax,jmax,5,vars2d)
         do i=1, 1
            do j=1, jmax
               do k=1, kmax
                  buffer_flow(k,i,j,1) = vars2d(k,j,1)
                  buffer_flow(k,i,j,2) = vars2d(k,j,2)
                  buffer_flow(k,i,j,3) = vars2d(k,j,3)
                  buffer_flow(k,i,j,4) = vars2d(k,j,4)*rbar*vars2d(k,i,5)  ! p = rho*R*T
                  buffer_flow(k,i,j,5) = vars2d(k,j,5)
               enddo
            enddo
         enddo
      case default
         print *, 'Inlet conversion fail: Unknown irhoinlet type ... STOP!!!'   
         stop
      end select     
      fsol%dimsf = (/kmax,1,jmax/)
      fsol%fname = trim(flowpath_wt)//'inlet.h5'
      print *, 'writing file ', trim(fsol%fname)
      call WriteHDF5_3D(fsol,buffer_flow)
    elseif(iFormat_inlet.eq.1) then  ! convert 'inlet.h5' to 'inlet.sol'
      print *, 'iFormat_inlet type NOT implemented'
      stop          
    elseif(iFormat_inlet.eq.2) then     ! convert plot3d to plot3d (p/t -- rho/t)
        ! read in 'inlet.sol'
        fname = trim(flowpath_rd)//'inlet.sol'
        call ReadPlot2DSolution(trim(fname),kmax,jmax,5,vars2d)
        if(irhoinlet.eq.1) then     !   irhoinlet=1 convert 'u,v,w,p,t' to 'u,v,w,rho,t'
            do i=1, 1
            do j=1, jmax
               do k=1, kmax
                  buffer_flow(k,i,j,1) = vars2d(k,j,1)  ! u
                  buffer_flow(k,i,j,2) = vars2d(k,j,2)  ! v
                  buffer_flow(k,i,j,3) = vars2d(k,j,3)  ! w
                  buffer_flow(k,i,j,4) = vars2d(k,j,4)/rbar/vars2d(k,i,5)  ! rho = p/r/t
                  buffer_flow(k,i,j,5) = vars2d(k,j,5)  ! t
               enddo
            enddo
            enddo
            fname = trim(flowpath_wt)//'inlet.sol_rt'
            call WritePlot2DSolution(trim(fname),kmax,imax,5,buffer_flow(1:kmax,1,1:jmax,1:5))
        elseif(irhoinlet.eq.2) then !   irhoinlet=2 convert 'u,v,w,rho,t' to 'u,v,w,p,t'
            do i=1, 1
            do j=1, jmax
               do k=1, kmax
                  buffer_flow(k,i,j,1) = vars2d(k,j,1)  ! u
                  buffer_flow(k,i,j,2) = vars2d(k,j,2)  ! v
                  buffer_flow(k,i,j,3) = vars2d(k,j,3)  ! w
                  buffer_flow(k,i,j,4) = vars2d(k,j,4)*rbar*vars2d(k,i,5)  ! p = rho*R*T
                  buffer_flow(k,i,j,5) = vars2d(k,j,5)  ! t
               enddo
            enddo
            enddo
            fname = trim(flowpath_wt)//'inlet.sol_pt'
            call WritePlot2DSolution(trim(fname),kmax,jmax,5,buffer_flow(1:kmax,1,1:jmax,1:5))
        else
            print*,'Unknown PLOT3D format of inlet.sol, choose irhoinlet=1/2...'
            stop
        endif
    endif  ! end if iFormat_inlet
    deallocate(vars2d, buffer_flow)
  endif  ! end if iConvert_inlet

  if(ireadRescalemean.eq.1) then
    call InitRescalemean(fsolRescale)
    allocate(buffer_flow(kmax,istencilsize,1,4))
    fsolRescale%dimsf = (/kmax,istencilsize,1/)

    num_file = (file_end - file_be)/file_skip + 1
    do n=1, num_file
      write(unit=fnum,fmt='(I08.8)') (file_be + (n-1)*file_skip)
      fname = trim(flowpath_rd)//'rescalemean_'//fnum//'.sol'
      print *, 'reading file ', trim(fname)
      open(21,file=trim(fname),form='unformatted',status='old')
        read(21) nblks
        if(nblks.ne.1) then
            print*,'number of blocks must be 1 in file: ', trim(fname)
            stop
        end if
        read(21) kdim, idim, nvar
        if(kdim.ne.kmax.or.idim.ne.istencilsize.or.nvar.ne.4) then
          print*,'dimension insistent in rescale mean file'
          stop
        end if
        read(21)((buffer_flow(k,i,1,1),   k=1,kmax),i=1,istencilsize),&
                ((buffer_flow(k,i,1,2),   k=1,kmax),i=1,istencilsize),&
                ((buffer_flow(k,i,1,3),   k=1,kmax),i=1,istencilsize),&
                ((buffer_flow(k,i,1,4),   k=1,kmax),i=1,istencilsize)
        read(21)nsample_rescale, theta_r
      close(21)

      fsolRescale%fname = trim(flowpath_wt)//'rescalemean_'//fnum//'.h5'
      print *, 'Writing file: ', trim(fsolRescale%fname)
      call WriteHDF5_3D(fsolRescale, buffer_flow)
      fsolRescale%sname = 'nsample_rescale'
      call WriteHDF5_scalar(fsolRescale,real(nsample_rescale))
      fsolRescale%sname = 'theta_r'
      call WriteHDF5_scalar(fsolRescale,real(theta_r))
    enddo
  endif ! end iradRescalemean.eq.1

  call FinalizeHDF5()

contains

  subroutine Initial()
     implicit none
     read(*,*)
     read(*,*) imax, jmax, kmax, kioplane, inode, jnode, Mw
     read(*,*)
     read(*,*) iConvert_grid, iFormat_grid
     read(*,*)
     read(*,*) iConvert_flow, iFormat_flow
     read(*,*)
     read(*,*) iConvert_inlet, iFormat_inlet
     read(*,*)
     read(*,*) ireadRescalemean, irhoinlet, istencilsize
     read(*,*)
     read(*,*)
     read(*,*)
     read(*,'(a)') gridfile_rd
     read(*,*)
     read(*,'(a)') gridfile_wt
     read(*,*)
     read(*,*)
     read(*,*)
     read(*,'(a)') flowpath_rd
     read(*,*)
     read(*,'(a)') flowpath_wt
     read(*,*)
     read(*,*)
     read(*,*)
     read(*,*) file_be, file_end, file_skip

     if(iConvert_grid.eq.1) then
       print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
       if(iFormat_grid.eq.0) then
         print *, 'Convert grid plot3d file to HDF5 ... '
       elseif(iFormat_grid.eq.1.or.iFormat_flow.eq.1) then
         print *, 'Convert grid HDF5 file to plot3d ... '
       elseif(iFormat_grid.eq.2.or.iFormat_flow.eq.2) then
         print *, 'separate single HDF5 grid and flowdata into multiple files, each has only one variable ... '
       elseif(iFormat_grid.eq.3.or.iFormat_flow.eq.3) then
         print *, 'separate single HDF5 grid and flowdata into multiple files using kioplane size'
       elseif(iFormat_grid.eq.4.or.iFormat_flow.eq.4) then
         print *, 'separate single HDF5 grid and flowdata into multiple files using inode'
       elseif(iFormat_flow.eq.5) then
         print *, 'separate single HDF5 flowdata into multiple files using jnode'
       else
         print *, 'unknown grid format. Stop!'
         stop
       endif
       print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
     endif

!     if(iConvert_grid.eq.1.or.iConvert_flow.eq.1) then
!       print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
!       if(iFormat_grid.eq.0.or.iFormat_flow.eq.0) then
!         print *, 'Convert grid plot3d file to HDF5 ... '
!       elseif(iFormat_grid.eq.1.or.iFormat_flow.eq.1) then
!         print *, 'Convert grid HDF5 file to plot3d ... '
!       elseif(iFormat_grid.eq.2.or.iFormat_flow.eq.2) then
!         print *, 'separate single HDF5 grid and flowdata into multiple files, each has only one variable ... '
!       elseif(iFormat_grid.eq.3.or.iFormat_flow.eq.3) then
!         print *, 'separate single HDF5 grid and flowdata into multiple files using kioplane size'
!       elseif(iFormat_grid.eq.4.or.iFormat_flow.eq.4) then
!         print *, 'separate single HDF5 grid and flowdata into multiple files using inode'
!       elseif(iFormat_flow.eq.5) then
!         print *, 'separate single HDF5 flowdata into multiple files using jnode'
!       else
!         print *, 'unknown grid format. Stop!'
!         stop
!       endif
!       print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
!     endif

     print *, 'File dimension: imax = ', imax, 'jmax = ', jmax, 'kmax = ', kmax

  end subroutine Initial

end program p3d_h5


