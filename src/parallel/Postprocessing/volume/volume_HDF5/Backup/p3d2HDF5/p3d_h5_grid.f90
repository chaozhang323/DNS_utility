!>
!> \file p3d_to_h5_grid.f90
!> @brief
!> Read in plot3d file and Write out HDF5 file.
!> @details
!> Write out HDF5 grid file.
!> The basic subroutines for Write are:
!> CALL h5open_f(error) 
!> CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
!> CALL h5screate_simple_f(rank, dims, dspace_id, error)
!> CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
!> CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, X, dims, error) 
!> CALL h5sclose_f(dspace_id, error) 
!> CALL h5dclose_f(dset_id, error)
!> CALL h5fclose_f(file_id, error) 
!> CALL h5close_f(error) 


PROGRAM p3d_h5_grid

  USE HDF5 ! This module contains all necessary modules

  IMPLICIT NONE

  CHARACTER(LEN=8), PARAMETER :: filename = "grid.h5" ! File name
  CHARACTER(LEN=5), PARAMETER :: dsetname1 = "x"     ! Dataset name
  CHARACTER(LEN=5), PARAMETER :: dsetname2 = "y"     ! Dataset name
  CHARACTER(LEN=5), PARAMETER :: dsetname3 = "z"     ! Dataset name

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier


  INTEGER(HSIZE_T), DIMENSION(3) :: dims           ! Dataset dimensions
  INTEGER     ::   rank = 3                        ! Dataset rank

  INTEGER     ::   error ! Error flag
  INTEGER     ::   i, j, k, imax, jmax, kmax

  REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: X, Y, Z
  REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE :: vars

  REAL(8) :: rtime_begin, rtime_end, dtime_begin, dtime_end
  CHARACTER(200) :: fn1

  CALL Initial()

  CALL CPU_TIME(rtime_begin) 
  CALL ReadPlot3DGridPlane(fn1,imax,jmax,kmax,X,Y,Z)  ! Read plot3d grid file.
  CALL CPU_TIME(rtime_end)    

  Do i=1, imax
    Do j=1, jmax
      Do k=1, kmax
        vars(k,i,j,1) = x(i,j,k)
        vars(k,i,j,2) = y(i,j,k)
        vars(k,i,j,3) = z(i,j,k)
     
      End Do
    End Do
  End Do

  dims(1) = kmax
  dims(2) = imax
  dims(3) = jmax

  ! Write HDF5 grid file.      

  CALL CPU_TIME(dtime_begin)

  CALL h5open_f(error)                                                               ! Initialize FORTRAN interface.
  CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)                        ! Create a new file (filename2) using default properties.

  CALL h5screate_simple_f(rank, dims, dspace_id, error)                              ! Create the dataspace.
  CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error) ! Create the dataset (dsetname1) with default properties.
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vars(:,:,:,1), dims, error)                        ! Write dataset
  CALL h5sclose_f(dspace_id, error)                                                  ! Terminate access to the data space.
  CALL h5dclose_f(dset_id, error)                                                    ! End access to the dataset and release resources used by it.

  CALL h5screate_simple_f(rank, dims, dspace_id, error)
  CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vars(:,:,:,2), dims, error)
  CALL h5sclose_f(dspace_id, error)
  CALL h5dclose_f(dset_id, error)

  CALL h5screate_simple_f(rank, dims, dspace_id, error)
  CALL h5dcreate_f(file_id, dsetname3, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vars(:,:,:,3), dims, error)
  CALL h5sclose_f(dspace_id, error)
  CALL h5dclose_f(dset_id, error)

  CALL h5fclose_f(file_id, error)                                                    ! Close the file.
  CALL h5close_f(error)                                                              ! Close FORTRAN interface.

  CALL CPU_TIME(dtime_end)



  print *, 'Time of operation for reading p3d grid file was ', rtime_end - rtime_begin, ' seconds '
  print *, 'Time of operatoin for writing H5 grid file was ', dtime_end - dtime_begin, ' seconds'


contains

     Subroutine ReadPlot3DGridPlane(fn,idim,jdim,kdim,x,y,z)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim
       real(8), intent(out), dimension(idim,jdim,kdim) :: x, y, z
       integer :: id,jd,kd,nb
       integer :: i, j, k 

       open(unit=11,file=fn,form='unformatted',status='old')
       read(11)nb
       if (nb.ne.1) then
         print*,'number of blocks must be 1 in file: ', trim(fn)
         stop
       end if
       read(11)id,jd,kd
       if (id.ne.idim.or.jd.ne.jdim.or.kd.ne.kdim) then
         print*,'size inconsistent with input in file: ', trim(fn)
         print*,'input idim jdim kdim', idim, jdim, kdim
         print*,'file idim jdim kdim', id, jd, kd
         print*,'can not proceed due to smaller file size'
         stop
       end if
         do k=1,kdim
           read(11)((x(i,j,k),i=1,idim),j=1,jdim)&
                 , ((y(i,j,k),i=1,idim),j=1,jdim)&
                 , ((z(i,j,k),i=1,idim),j=1,jdim)
         end do
       close(11)

  
     End Subroutine ReadPlot3DGridPlane


  
  Subroutine Initial()
       implicit none
       Read(*,*)
       Read(*,*) imax
       Read(*,*)
       Read(*,*) jmax
       Read(*,*)
       Read(*,*) kmax
       Read(*,*)
       Read(*,'(A)') fn1

      print *, 'imax = ', imax
      print *, 'jmax = ', jmax
      print *, 'kmax = ', kmax


       Allocate(X(imax,jmax,kmax), Y(imax,jmax,kmax), Z(imax,jmax,kmax))
       Allocate(vars(kmax,imax,jmax,3))
  End Subroutine Initial




END PROGRAM p3d_h5_grid



