!>
!> \file p3d_to_h5_sol.f90
!> @brief
!> Read in plot3d file and Write out HDF5 file.
!> @details
!> Write out HDF5 sol file.
!> The basic subroutines for Write are:
!> CALL h5open_f(error) 
!> CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
!> CALL h5screate_simple_f(rank, dims, dspace_id, error)
!> CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
!> CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, X, dims, error) 
!> CALL h5sclose_f(dspace_id, error) 
!> CALL h5dclose_f(dset_id, error)
!> CALL h5fclose_f(file_id, error) 
!> CALL h5close_f(error) 


PROGRAM p3d_h5_sol

  USE HDF5 ! This module contains all necessary modules

  IMPLICIT NONE

  CHARACTER(LEN=100) :: filename            ! File name
  CHARACTER(LEN=5)   :: dsetname1 = "u"     ! Dataset name
  CHARACTER(LEN=5)   :: dsetname2 = "v"     ! Dataset name
  CHARACTER(LEN=5)   :: dsetname3 = "w"     ! Dataset name
  CHARACTER(LEN=5)   :: dsetname4 = "p"     ! Dataset name
  CHARACTER(LEN=5)   :: dsetname5 = "T"     ! Dataset name
  CHARACTER(LEN=5)   :: dsetname6 = "time"   ! Dataset name

  INTEGER(HID_T) :: file_id       ! File identifier
  INTEGER(HID_T) :: dset_id       ! Dataset identifier
  INTEGER(HID_T) :: dspace_id     ! Dataspace identifier


  INTEGER(HSIZE_T), DIMENSION(3) :: dims           ! Dataset dimensions
  INTEGER     ::   rank = 3                        ! Dataset rank

  INTEGER     ::   error ! Error flag
  INTEGER     ::   i, j, k, n, imax, jmax, kmax

  REAL(8), DIMENSION(:,:,:,:), ALLOCATABLE :: vars, vars2
  REAL(8) :: rtime_begin, rtime_end, dtime_begin, dtime_end, time
  CHARACTER(200) :: fn1
  CHARACTER(8) :: fnum
  INTEGER :: file_be, file_end, file_skip, num_file
  CHARACTER(200) :: filepath !!!!!!

  CALL Initial()

  num_file = (file_end - file_be)/file_skip + 1
  
  Do n=1, num_file
    write(unit=fnum,fmt='(I08.8)') (file_be + (n-1)*file_skip)
    fn1 = trim(filepath)//'flowdata_'//fnum//'.sol'
    CALL ReadPlot3DSolPlane(fn1,imax,jmax,kmax,6,vars,time)  ! Read plot3d sol file
    dims(1) = kmax
    dims(2) = imax
    dims(3) = jmax

    Do i=1, imax
      Do j=1, jmax
        Do k=1, kmax
          vars2(k,i,j,1) = vars(i,j,k,1)
          vars2(k,i,j,2) = vars(i,j,k,2)
          vars2(k,i,j,3) = vars(i,j,k,3)
          vars2(k,i,j,4) = vars(i,j,k,4)
          vars2(k,i,j,5) = vars(i,j,k,5)
        End Do
      End Do
    End Do 

    filename = trim(filepath)//'flowdata_'//fnum//'.h5'

    ! Write HDF sol file

    CALL CPU_TIME(dtime_begin)

    CALL h5open_f(error)                                                                    ! Initialize FORTRAN interface.
    CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)                             ! Create a new file (filename2) using default properties.

    ! Write variable vars(1:imax,1:jmax,1:kmax,1) to the file
    CALL h5screate_simple_f(rank, dims, dspace_id, error)                                   ! Create the dataspace. 
    CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)      ! Create the dataset (dsetname1) with default properties.
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vars2(:,:,:,1), dims, error)  ! Write dataset
    CALL h5sclose_f(dspace_id, error)                                                       ! Terminate access to the data space.
    CALL h5dclose_f(dset_id, error)                                                         ! End access to the dataset and release resources used by it.

    ! Write variable vars(1:imax,1:jmax,1:kmax,2) to the file
    CALL h5screate_simple_f(rank, dims, dspace_id, error)
    CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vars2(:,:,:,2), dims, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5dclose_f(dset_id, error)

    ! Write variable vars(1:imax,1:jmax,1:kmax,3) to the file
    CALL h5screate_simple_f(rank, dims, dspace_id, error)
    CALL h5dcreate_f(file_id, dsetname3, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vars2(:,:,:,3), dims, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5dclose_f(dset_id, error)

    ! Write variable vars(1:imax,1:jmax,1:kmax,4) to the file
    CALL h5screate_simple_f(rank, dims, dspace_id, error)
    CALL h5dcreate_f(file_id, dsetname4, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vars2(:,:,:,4), dims, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5dclose_f(dset_id, error)

    ! Write variable vars(1:imax,1:jmax,1:kmax,5) to the file
    CALL h5screate_simple_f(rank, dims, dspace_id, error)
    CALL h5dcreate_f(file_id, dsetname5, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vars2(:,:,:,5), dims, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5dclose_f(dset_id, error)

  ! Write variable vars(1:imax,1:jmax,1:kmax,6) to the file
!  CALL h5screate_simple_f(rank, dims, dspace_id, error)
!  CALL h5dcreate_f(file_id, dsetname6, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
!  CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, vars(:,:,:,6), dims, error)
!  CALL h5sclose_f(dspace_id, error)
!  CALL h5dclose_f(dset_id, error)


      CALL h5screate_f(H5S_SCALAR_F, dspace_id, error)
      CALL h5dcreate_f(file_id, dsetname6, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
      CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, time, dims, error) 
      CALL h5sclose_f(dspace_id, error)
      CALL h5dclose_f(dset_id, error)



    CALL h5fclose_f(file_id, error)                                                         ! Close the file.
    CALL h5close_f(error)                                                                   ! Close FORTRAN interface.

    print *, 'finish file: ',  'flowdata_'//fnum//'.h5'  

  enddo ! end n loop




contains

     subroutine ReadPlot3DSolPlane(fn,idim,jdim,kdim,nvar,vars,time)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim,nvar
       real(8), intent(out), dimension(idim,jdim,kdim,nvar) :: vars
       real(8), intent(out), optional :: time
       real(8) :: rh
       integer :: id,jd,kd, nv,nb
       integer :: i, j, k, n

       open(unit=11,file=fn,form='unformatted',status='old')
       read(11)nb
       if (nb.ne.1) then
         print*,'number of blocks must be 1 in file: ', trim(fn)
         stop
       end if
       read(11)id,jd,kd,nv
       if (id.ne.idim.or.jd.ne.jdim.or.kd.ne.kdim.or.nv.ne.nvar) then
         print*,'size inconsistent with input in file: ', trim(fn)
         print*,'input idim jdim kdim nvar', idim, jdim, kdim, nvar
         print*,'file idim jdim kdim nvar', id, jd, kd, nv
         print*,'can not proceed due to smaller file size'
         stop
       end if
       do k=1,kdim
         read(11)(((vars(i,j,k,n),i=1,idim),j=1,jdim),n=1,nvar)
       end do
       if (present(time)) then
         Read(11)time
       end if
       close(11)

     end subroutine ReadPlot3DSolPlane


 
  Subroutine Initial()
       implicit none
       Read(*,*)
       Read(*,*) imax, jmax, kmax
       Read(*,*)
       Read(*,'(a)') filepath
      read(*,*)
      read(*,*) file_be, file_end, file_skip


      print *, 'imax = ', imax
      print *, 'jmax = ', jmax
      print *, 'kmax = ', kmax

       Allocate(vars(imax,jmax,kmax,6),vars2(kmax,imax,jmax,5))

  End Subroutine Initial


END PROGRAM p3d_h5_sol



