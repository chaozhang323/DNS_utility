!>
!> Create single volume Virtual Dataset
!>
!>

program VDS
  use modRWHDF5
  use modVDSHDF5
  implicit none

  integer ::   i, j, k, kk, k2, kp, n, nn, imax, jmax, kmax, ires, jres, icurrent, jcurrent, kcurrent
  integer, dimension(:), allocatable :: isize, jsize
  real(8), parameter :: Univ = 8314.3D0
  integer :: iConvert_grid, iConvert_flow
  integer :: iFormat_grid, iFormat_flow
  character(300) :: filepath_rd, filepath_wt
  integer :: file_be, file_end, file_skip, num_file, num_links
  character(400) :: fname, VDSname
  character(8) :: fnum
  character(4) :: fnum4_1, fnum4_2
  type(tp_rdwt_hdf5) :: grd, fsol
  type(tp_rdwt_hdf5), dimension(:), allocatable :: FILES
  integer :: kioplane, inode, jnode
  real(8) :: time

  call InitHDF5()
  call Initial()
  
  ! Convert grid file
  if(iConvert_grid.gt.0) then
    call InitGridHDF5(grd)
    grd%dimsf=(/kmax,imax,jmax/)

    ! create VDS file using kioplane
    if(iFormat_grid.eq.1) then
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !finding amount of files
      k=kmax/kioplane
      ires=kmax-k*kioplane
      num_links = k
      if(ires.gt.0) num_links=k+1
      allocate(FILES(num_links))

      grd%fname= trim(filepath_wt)//'grid.h5'
      print *, 'writing file: ', trim(grd%fname)
      kcurrent=0
      do kk=1, kmax, kioplane !defining info of each file to be linked
        kcurrent=kcurrent+1
        call InitGridHDF5(FILES(kcurrent))

        k2 = min(kk+kioplane-1,kmax)
        kp = k2-kk+1

        FILES(kcurrent)%dimsf = (/kp,imax,jmax/)
        FILES(kcurrent)%offset = (/kk-1,0,0/)

        !generating expected file name
        write(unit=fnum4_1,fmt='(I04.4)') kk
        write(unit=fnum4_2,fmt='(I04.4)') k2
        FILES(kcurrent)%fname = trim(filepath_rd)//'grid_kplane'//fnum4_1//'_'//fnum4_2//'.h5'
      enddo!end kk loop
      call CreateHDF5_VDS_3D(grd,num_links,FILES)
      deallocate(FILES)
    elseif(iFormat_grid.eq.2) then
      num_links=inode
      allocate(FILES(num_links))
      allocate(isize(0:inode-1))

      !defining size of each file
      i=imax/inode
      ires = imax - i*inode
      do n=0, inode-1
        isize(n) = i
        if(n.lt.ires) isize(n) = i+1
      enddo
      grd%fname= trim(filepath_wt)//'grid.h5'
      print *, 'writing file: ', trim(grd%fname)
      icurrent = 0
      do nn=0,inode-1 !loop each file to generate information for VDS
        call InitGridHDF5(FILES(nn+1))
        FILES(nn+1)%offset=(/0,icurrent,0/)
        icurrent = icurrent + isize(nn)
        FILES(nn+1)%dimsf = (/kmax,isize(nn),jmax/)
        write(unit=fnum4_1,fmt='(I04.4)') icurrent-isize(nn)+1
        write(unit=fnum4_2,fmt='(I04.4)') icurrent
        FILES(nn+1)%fname = trim(filepath_rd)//'grid_i'//fnum4_1//'_'//fnum4_2//'.h5'
      enddo!end nn loop defining each file data
      call CreateHDF5_VDS_3D(grd,num_links,FILES)
      deallocate(isize)
      deallocate(FILES)

    elseif(iFormat_grid.eq.3) then
      num_links=jnode
      allocate(FILES(num_links))
      allocate(jsize(0:jnode-1))

      !defining size of each file
      j=jmax/jnode
      jres = jmax - j*jnode
      do n=0, jnode-1
        jsize(n) = j
        if(n.lt.jres) jsize(n) = j+1
      enddo
      grd%fname= trim(filepath_wt)//'grid.h5'
      print *, 'writing file: ', trim(grd%fname)
      jcurrent = 0
      do nn=0,jnode-1 !loop each file to generate information for VDS
        call InitGridHDF5(FILES(nn+1))
        FILES(nn+1)%offset=(/0,0,jcurrent/)
        jcurrent = jcurrent + jsize(nn)
        FILES(nn+1)%dimsf = (/kmax,imax,jsize(nn)/)
        write(unit=fnum4_1,fmt='(I04.4)') jcurrent-jsize(nn)+1
        write(unit=fnum4_2,fmt='(I04.4)') jcurrent
        FILES(nn+1)%fname = trim(filepath_rd)//'grid_j'//fnum4_1//'_'//fnum4_2//'.h5'
      enddo!end nn loop defining each file data
      call CreateHDF5_VDS_3D(grd,num_links,FILES)
      deallocate(jsize)
      deallocate(FILES)
    else
      print *, 'unknown iFormat_grid ... STOP !'
      stop
    endif


  endif

  if(iConvert_flow.gt.0) then ! Convert flow file
    call InitFlowHDF5(fsol)
    fsol%dimsf=(/kmax,imax,jmax/)
    num_file = (file_end - file_be)/file_skip + 1
    print *, 'Total number of files: ', num_file

    ! create VDS file using kioplane
    if(iFormat_flow.eq.1) then
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !finding amount of files
      k=kmax/kioplane
      ires=kmax-k*kioplane
      num_links = k
      if(ires.gt.0) num_links=k+1
      allocate(FILES(num_links))

      do n=1, num_file
        write(unit=fnum,fmt='(I08.8)') (file_be + (n-1)*file_skip)
        fsol%fname= trim(filepath_wt)//'flowdata_'//fnum//'.h5'
        print *, 'writing file: ', trim(fsol%fname)
        kcurrent=0
        do kk=1, kmax, kioplane !defining info of each file to be linked
          kcurrent=kcurrent+1    
          call InitFlowHDF5(FILES(kcurrent))
           
          k2 = min(kk+kioplane-1,kmax)
          kp = k2-kk+1
          
          FILES(kcurrent)%dimsf = (/kp,imax,jmax/)
          FILES(kcurrent)%offset = (/kk-1,0,0/)
          
          !generating expected file name
          write(unit=fnum4_1,fmt='(I04.4)') kk
          write(unit=fnum4_2,fmt='(I04.4)') k2
          FILES(kcurrent)%fname = trim(filepath_rd)//'flowdata_'//fnum//'_kplane'//fnum4_1//'_'//fnum4_2//'.h5'
        enddo!end kk loop
        call CreateHDF5_VDS_3D(fsol,num_links,FILES)
        call ReadHDF5_scalar(FILES(1),time)
        call WriteHDF5_scalar(fsol,time)
      enddo ! end n loop
      deallocate(FILES)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! create VDS file using i-node
    elseif(iFormat_flow.eq.2) then !inode
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      num_links=inode
      allocate(FILES(num_links))
      allocate(isize(0:inode-1))
      
      !defining size of each file
      i=imax/inode
      ires = imax - i*inode
      do n=0, inode-1
        isize(n) = i
        if(n.lt.ires) isize(n) = i+1
      enddo

      do n=1, num_file !loop each file group
        write(unit=fnum,fmt='(I08.8)') (file_be + (n-1)*file_skip)
        fsol%fname= trim(filepath_wt)//'flowdata_'//fnum//'.h5'
        print *, 'writing file: ', trim(fsol%fname)
        icurrent = 0
        do nn=0,inode-1 !loop each file to generate information for VDS
          call InitFlowHDF5(FILES(nn+1))
          FILES(nn+1)%offset=(/0,icurrent,0/)
          icurrent = icurrent + isize(nn)
          FILES(nn+1)%dimsf = (/kmax,isize(nn),jmax/)
          write(unit=fnum4_1,fmt='(I04.4)') icurrent-isize(nn)+1
          write(unit=fnum4_2,fmt='(I04.4)') icurrent
          FILES(nn+1)%fname = trim(filepath_rd)//'flowdata_i'//fnum4_1//'_'//fnum4_2//'_'//fnum//'.h5'
        enddo!end nn loop defining each file data
        call CreateHDF5_VDS_3D(fsol,num_links,FILES)
        call ReadHDF5_scalar(FILES(1),time)
        call WriteHDF5_scalar(fsol,time)
      enddo ! end n loop
      
      deallocate(isize)        
      deallocate(FILES)
    
    elseif(iFormat_flow.eq.3) then !jnode
      num_links=jnode
      allocate(FILES(num_links))
      allocate(jsize(0:jnode-1))
      
      !defining size of each file
      j=jmax/jnode
      jres = jmax - j*jnode
      do n=0, jnode-1
        jsize(n) = j
        if(n.lt.jres) jsize(n) = j+1
      enddo

      do n=1, num_file
        write(unit=fnum,fmt='(I08.8)') (file_be + (n-1)*file_skip)
        fsol%fname= trim(filepath_wt)//'flowdata_'//fnum//'.h5'
        print *, 'writing file: ', trim(fsol%fname)
        jcurrent = 0
        do nn=0,jnode-1 !loop each file to generate information for VDS
          call InitFlowHDF5(FILES(nn+1))
          FILES(nn+1)%offset=(/0,0,jcurrent/)
          jcurrent = jcurrent + jsize(nn)
          FILES(nn+1)%dimsf = (/kmax,imax,jsize(nn)/)
          write(unit=fnum4_1,fmt='(I04.4)') jcurrent-jsize(nn)+1
          write(unit=fnum4_2,fmt='(I04.4)') jcurrent
          FILES(nn+1)%fname = trim(filepath_rd)//'flowdata_j'//fnum4_1//'_'//fnum4_2//'_'//fnum//'.h5'
        enddo!end nn loop defining each file data
        call CreateHDF5_VDS_3D(fsol,num_links,FILES)
        call ReadHDF5_scalar(FILES(1),time)
        call WriteHDF5_scalar(fsol,time)
      enddo ! end n loop
      
      deallocate(jsize)        
      deallocate(FILES)

    else
      print *, 'unknown iFormat_flow ... STOP !'
      stop
    endif!

  endif  ! end if iConvert_flow

 
  call FinalizeHDF5()

contains

  subroutine Initial()
     implicit none
     read(*,*)
     read(*,*) imax, jmax, kmax, kioplane, inode, jnode
     read(*,*)
     read(*,*) iConvert_grid, iFormat_grid
     read(*,*)
     read(*,*) iConvert_flow, iFormat_flow
     read(*,*)
     read(*,*)
     read(*,*)
     read(*,'(a)') filepath_rd
     read(*,*)
     read(*,'(a)') filepath_wt
     read(*,*)
     read(*,*)
     read(*,*)
     read(*,*) file_be, file_end, file_skip

     if(iConvert_grid.eq.1.or.iConvert_flow.eq.1) then
       print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
       if(iFormat_grid.eq.1.or.iFormat_flow.eq.1) then
         print *, 'Make VDS file from multiple HDF5 files using kioplane'
       elseif(iFormat_grid.eq.2.or.iFormat_flow.eq.2) then
         print *, 'Make VDS file from multiple HDF5 files using inode'
       elseif(iFormat_grid.eq.3.or.iFormat_flow.eq.3) then
         print *, 'Make VDS file from multiple HDF5 files using jnode'
       else
         print *, 'unknown grid format. Stop!'
         stop
       endif
       print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
     endif

     print *, 'File dimension for the Virtual dataset'
     print *, 'imax = ', imax, 'jmax = ', jmax, 'kmax = ', kmax
     print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

  end subroutine Initial
  
end program VDS




