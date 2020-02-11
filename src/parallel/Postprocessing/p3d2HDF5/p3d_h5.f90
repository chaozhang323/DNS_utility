 !>
!> \file p3d_h5.f90
!> @brief
!> @details

PROGRAM p3d_h5
  use decomp_2d
  use MPRWHDF5
  use MFileIO
  implicit none

  integer ::   i, j, k, kk, k2, kp, n, nn, imax, jmax, kmax, ires, icurrent
  integer, dimension(:), allocatable :: isize
  real(8), parameter :: Univ = 8314.3D0
  real(8), dimension(:,:,:), allocatable :: X, Y, Z
  real(8), dimension(:,:,:,:), allocatable :: vars, buffer_grid, buffer_flow, buffer_flow_tmp
  real(8), dimension(:,:,:), allocatable :: vars2d, buffer
  real(8) ::time, rbar
  integer :: iConvert_grid, iConvert_flow, iConvert_inlet
  integer :: iFormat_grid, iFormat_flow, iFormat_inlet
  character(300) :: gridfile_rd, gridfile_wt, flowpath_rd, flowpath_wt
  integer :: file_be, file_end, file_skip, num_file
  character(400) :: fname
  character(8) :: fnum
  character(4) :: fnum4_1, fnum4_2
  type(tp_rdwt_hdf5) :: grd, fsol, fsolRescale, grd1V, fsol1V,fsol_new,fsol_new2
  integer :: irhoinlet  ! 'inlet.sol' varialbe list (=1: u, v, w, rho, T)
  real(8) :: Mw
  integer :: ireadRescalemean, istencilsize, kioplane, nblks, kdim, idim, nvar, inode
  real(8) :: theta_r
  integer :: nsample_rescale
  character(5) :: gridname(3) = (/'x','y','z'/)
  character(5) :: flowdataname(5) = (/'u','v','w','p','T'/)
  !add for multi proc
  type(DECOMP_INFO) :: decomp
  integer :: ierr, numprocs, myid ,p_row,p_col
  logical :: isGrdHDF5initiated = .False., isFlowHDF5initiated = .False.
  
  
  !testting garbo
  !open(unit=1,file="data")

  ! initialize MPI
  call MPI_INIT(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  
  if(myid.eq.0) print *, 'Starting reading parameters'
  call Initial()
  if(myid.eq.0) print *, 'Finished reading parameters'
  call InitHDF5()
  
  if(myid.eq.0) then
    rbar = Univ/Mw
    print *, 'Gas constant =', rbar
  endif
  

  
  if(iConvert_flow.eq.1) then ! Convert flow file
    call InitFlowHDF5(fsol)
    call InitFlowHDF5(fsol_new)
    
    
    
    if(iFormat_flow.eq.4) then

    call decomp_2d_init(kmax,imax,jmax,p_row,p_col)
    call decomp_info_init(kmax,imax,jmax,decomp)
    
    !if(myid.eq.3) print *, 'xst(1)=', decomp%xen(1),  'xst(2)=', decomp%xen(2),  'xst(3)=', decomp%xen(3)
  


      do nn=1, num_file !would stay
      
        write(unit=fnum,fmt='(I08.8)') (file_be + (nn-1)*file_skip)
        

        fsol%dimsf = (/kmax,decomp%xsz(2),jmax/) !not used currently
        fsol_new%dimsf = (/kmax,decomp%xsz(2),jmax/)
        
        fsol%offset = (/0,imax-decomp%xen(2)+1,0/) !not used currently
        if(nn.eq.1) allocate(buffer_flow(kmax,decomp%xsz(2),jmax,5)) !not used currently 
        if(nn.eq.1) allocate(buffer(kmax,decomp%xsz(2),jmax)) 
        fsol%fname = trim(flowpath_rd)//'flowdata_'//fnum//'.h5' 
        if(myid.eq.0) print *, 'reading file ', trim(fsol%fname) 
        
       !reading single varibles  
        do n=1, fsol%dnum

             if(myid.eq.0) print *, 'Reading variable: ', trim(fsol%dname(n))
           call Read3dHDF5_svariable_P(fsol%fname,fsol%dname(n),kmax,decomp%xsz(2),jmax, &
                                                                decomp%xst(1),imax-decomp%xen(2)+1,decomp%xst(3), &
                                                                buffer)
           if(n.eq.1)call ReadHDF5_scalar(fsol,time)
             


           
           
                                                                
           if(n.eq.1) then             !set .eq.0 to .eq.1 for file spliting                                        
            write(unit=fnum4_1,fmt='(I04.4)') imax-decomp%xen(2)+1!decomp%xst(2)
            write(unit=fnum4_2,fmt='(I04.4)') imax-decomp%xst(2)+1!decomp%xen(2)
            fsol_new%fname = 'flowdata_i'//fnum4_1//'_'//fnum4_2//'_'//fnum//'.h5'
            print *, "making file: ",fsol_new%fname                                                   
           endif                                                   
                                             
                                                                

           call Write3dHDF5_svariable_S(fsol_new,buffer,n)  !change imax to decomp%xsz(2) and decomp%xst to 1 for multi file
           if(n.eq.1) call WriteHDF5_scalar(fsol_new,time)
        enddo ! end n loop    
        
        
        !endif !myid testing  
        !close(1)
            
        !call ReadHDF5_3D(fsol,buffer_flow)!(1:kmax,1:decomp%xsz(2),1:jmax),1:5)) !decompxsz(1)

            
        
        !if(myid.eq.0) call ReadHDF5_scalar(fsol,time) !myid

        !fsol%offset = 0
        !write(unit=fnum4_1,fmt='(I04.4)') decomp%xst(2)
        !write(unit=fnum4_2,fmt='(I04.4)') decomp%xen(2)
        !fsol%fname = 'flowdata_i'//fnum4_1//'_'//fnum4_2//'_'//fnum//'.h5'
        !print *, 'writing file ... ', trim(fsol%fname)
        !call WriteHDF5_3D_S(fsol, buffer_flow)!(1:kmax,1:decomp%xsz(2),1:jmax),1:5))
        !if(myid.eq.0) call WriteHDF5_scalar(fsol,time)
        
      enddo ! end nn loop
    endif

    deallocate(buffer_flow)
    deallocate(buffer)
    call decomp_2d_finalize
    call decomp_info_finalize(decomp)
    
  endif  ! end if iConvert_flow

 
  call FinalizeHDF5()
  call MPI_FINALIZE(ierr)

contains

  subroutine Initial()
     implicit none
     
     if(myid.eq.0) then
     
       read(*,*)
       read(*,*) imax, jmax, kmax, kioplane, inode, Mw
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

       if(iConvert_grid.eq.1.or.iConvert_flow.eq.1) then
         print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
         if(iFormat_grid.eq.0.or.iFormat_flow.eq.0) then
           print *, 'Convert grid plot3d file to HDF5 ... '
         elseif(iFormat_grid.eq.1.or.iFormat_flow.eq.1) then
           print *, 'Convert grid HDF5 file to plot3d ... '
         elseif(iFormat_grid.eq.2.or.iFormat_flow.eq.2) then
           print *, 'separate single HDF5 grid and flowdata into multiple files, each has only one variable ... '
         elseif(iFormat_grid.eq.3.or.iFormat_flow.eq.3) then
           print *, 'separate single HDF5 grid and flowdata into multiple files using kioplane size'
         elseif(iFormat_grid.eq.4.or.iFormat_flow.eq.4) then
           print *, 'separate single HDF5 grid and flowdata into multiple files using inode'
         else
           print *, 'unknown grid format. Stop!'
           stop
         endif
         print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
       endif

       print *, 'File dimension: imax = ', imax, 'jmax = ', jmax, 'kmax = ', kmax
       
       p_row= inode !this could be a function of numprocs
       p_col= 1
     
       num_file = (file_end - file_be)/file_skip + 1 !slightly out of place but makes more sense here
       
       
       

       




     endif !end myid=0
     
     call MPI_Bcast(imax,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(jmax,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(kmax,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(kioplane,          1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(inode,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(MW,                1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(iConvert_grid,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(iFormat_grid,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(iConvert_flow,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(iFormat_flow,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(iConvert_inlet,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(iFormat_inlet,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr) 
     call MPI_Bcast(ireadRescalemean,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(irhoinlet,         1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(istencilsize,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(gridfile_rd,       300, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(gridfile_wt,       300, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(flowpath_rd,       300, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(flowpath_wt,       300, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(file_be,           1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(file_end,          1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(file_skip,         1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(num_file,          1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(p_row,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(p_col,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)




  end subroutine Initial

end program p3d_h5


