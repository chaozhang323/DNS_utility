
!> parallel version to init grid/flowdata file

   program InitGridFlow_HDF5
     use decomp_2d
     use MPRWHDF5
     implicit none

     integer :: initGrid, initFlow, i, j, k
     integer :: imaxo, jmaxo, kmaxo, imax, jmax, kmax
     integer :: ixWavelengthConst, iyWavelengthConst
     real(8) :: time

     real(8), dimension(:), allocatable :: dim1input_index, dim1grid_index
     real(8), dimension(:), allocatable :: dim2input_index, dim2grid_index
     real(8), dimension(:), allocatable :: dim3input_index, dim3grid_index
     real(8), dimension(:,:,:), allocatable :: buffer, buffer_new

     type(tp_rdwt_hdf5) :: grd, fsol, grd_new1, grd_new2, fsol_new1, fsol_new2
     logical :: isGrdHDF5initiated = .False., isFlowHDF5initiated = .False.

     integer :: ierr
     integer :: myid, numprocs
     integer :: p_row, p_col
     integer :: errcode
     type(DECOMP_INFO) :: decompo, decomp

     ! initialize MPI
     call MPI_INIT(ierr)
     call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
     call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
     if(myid.eq.0) print *, 'Started reading parameters'
     call Input()
     if(myid.eq.0) print *, 'Finished reading parameters'
     call InitHDF5()

     call ReadGridFlow()

     call FinalizeHDF5()
     call MPI_FINALIZE(ierr)

contains

     subroutine Input()
       implicit none

       if(myid.eq.0) then
         read(*,*)
         read(*,*) initGrid, initFlow
         read(*,*)
         read(*,*) imaxo, jmaxo, kmaxo
         read(*,*)
         read(*,*) imax, jmax, kmax
         read(*,*)
         read(*,*) ixWavelengthConst, iyWavelengthConst
         read(*,*)
         read(*,'(a)') grd%fname
         read(*,*)
         read(*,'(a)') fsol%fname
         p_row = numprocs
         p_col = 1

       endif ! end myid.eq.0

       call MPI_Bcast(initGrid,           1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(initFlow,           1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(imaxo,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(jmaxo,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(kmaxo,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(imax,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(jmax,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(kmax,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(ixWavelengthConst,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(iyWavelengthConst,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(p_row,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(p_col,              1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(grd%fname,  300, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
       call MPI_Bcast(fsol%fname, 300, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

       if(myid.eq.0) then
         print *, 'imaxo = ', imaxo, 'jmaxo = ', jmaxo, 'kmaxo = ', kmaxo
         print *, 'imax  = ', imax,  'jmax  = ', jmax,  'kmax  = ', kmax
         print *, 'p_row = ', p_row, 'p_col = ', p_col
       endif

       if(initGrid.gt.0) then
         grd%dnum = 3
         allocate(grd%dname(grd%dnum))
         grd%dname = (/'x','y','z'/)
         isGrdHDF5initiated = .true.

         grd_new1%fname  = 'grid_new1.h5'
         grd_new2%fname  = 'grid_new2.h5'

       endif
       if(initFlow.gt.0) then
         fsol%dnum = 5
         allocate(fsol%dname(fsol%dnum))
         fsol%dname = (/'u','v','w','p','T'/)
         fsol%sname = 'time'
         isFlowHDF5initiated = .true.

         fsol_new1%fname = 'flowdata_new1.h5'
         fsol_new1%sname = 'time'
         fsol_new2%fname = 'flowdata_new2.h5'
         fsol_new2%sname = 'time'
       endif

     end subroutine Input

     subroutine GenIndex(dim1o,dim2o,dim3o,dim1,dim2,dim3)
       implicit none
       integer, intent(in) :: dim1o,dim2o,dim3o,dim1,dim2,dim3

       if(allocated(dim2input_index)) deallocate(dim2input_index, dim2grid_index)
       allocate(dim2input_index(dim2o), dim2grid_index(dim2))
       dim2input_index = 1.d0; dim2grid_index = 1.d0
       if(ixWavelengthConst.eq.0) then
          forall(i=1:dim2o, dim2o.gt.1)
             dim2input_index(i) = dble(i-1)/dble(dim2o-1)
          end forall
          forall(i=1:dim2,dim2.gt.1)
             dim2grid_index(i) = dble(i-1)/dble(dim2-1)
          end forall
        else
          forall(i=1:dim2o, dim2o.gt.1)
             dim2input_index(i) = dble(i-1)/dble(dim2o)
          end forall
          forall(i=1:dim2,dim2.gt.1)
             dim2grid_index(i) = dble(i-1)/dble(dim2)
          end forall
        endif

        if(allocated(dim3input_index)) deallocate(dim3input_index, dim3grid_index)
        allocate(dim3input_index(dim3o), dim3grid_index(dim3))
        dim3input_index = 1.d0; dim3grid_index = 1.d0
        if(iyWavelengthConst.eq.0) then
           forall(j=1:dim3o,dim3o.gt.1)
               dim3input_index(j) = dble(j-1)/dble(dim3o-1)
            end forall
            forall(j=1:dim3,dim3.gt.1)
               dim3grid_index(j) = dble(j-1)/dble(dim3-1)
             end forall
         else
           forall(j=1:dim3o,dim3o.gt.1)
               dim3input_index(j) = dble(j-1)/dble(dim3o)
           end forall
           forall(j=1:dim3,dim3.gt.1)
               dim3grid_index(j) = dble(j-1)/dble(dim3)
           end forall
         endif

        if(allocated(dim1input_index)) deallocate(dim1input_index, dim1grid_index)
        allocate(dim1input_index(dim1o), dim1grid_index(dim1))
        dim1input_index = 1.d0; dim1grid_index = 1.d0
        forall(k=1:dim1o)
           dim1input_index(k) = dble(k-1)/dble(dim1o-1)
        end forall
        forall(k=1:dim1)
           dim1grid_index(k) = dble(k-1)/dble(dim1-1)
        end forall

     end subroutine GenIndex

     subroutine ReadGridFlow()
       implicit none
       integer :: n

       !! init decomp
       call decomp_2d_init(kmaxo,imaxo,jmaxo,p_row,p_col)
       call decomp_info_init(kmaxo,imaxo,jmaxo,decompo)
       call decomp_2d_finalize
       call decomp_2d_init(kmaxo,imax,jmax,p_row,p_col)
       call decomp_info_init(kmaxo,imax,jmax,decomp)

       call GenIndex(decompo%zsz(1),imaxo,jmaxo,decompo%zsz(1),imax,jmax)

       allocate(buffer(decompo%zsz(1),decompo%zsz(2),decompo%zsz(3)))
       allocate(buffer_new(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)))
       if(initGrid.gt.0) then
         if(myid.eq.0) print *, 'Reading file: ', trim(grd%fname)
         do n=1, grd%dnum
           if(myid.eq.0) print *, 'Reading variable: ', trim(grd%dname(n))
           call Read3dHDF5_svariable_P(grd%fname,grd%dname(n),decompo%zsz(1),decompo%zsz(2),decompo%zsz(3), &
                                                              decompo%zst(1),decompo%zst(2),decompo%zst(3), &
                                                              buffer)
           if(myid.eq.0) print *, 'Interpolating ...'
           call InitVar(buffer,buffer_new)
           call Write3dHDF5_svariable_P(trim(grd_new1%fname),trim(grd%dname(n)),kmaxo,imax,jmax, &
                                                      decomp%zst(1),decomp%zst(2),decomp%zst(3), &
                                                      buffer_new,n-1)
          enddo ! end n loop
       endif ! initGrid.gt.0

       if(initFlow.gt.0) then
         if(myid.eq.0) print *, 'Reading file: ', trim(fsol%fname)
         do n=1, fsol%dnum
           if(myid.eq.0) print *, 'Reading variable: ', trim(fsol%dname(n))
           call Read3dHDF5_svariable_P(fsol%fname,fsol%dname(n),decompo%zsz(1),decompo%zsz(2),decompo%zsz(3), &
                                                                decompo%zst(1),decompo%zst(2),decompo%zst(3), &
                                                                buffer)
           call InitVar(buffer,buffer_new)
           if(myid.eq.0) print *, 'Interpolating ... '
           call Write3dHDF5_svariable_P(trim(fsol_new1%fname),trim(fsol%dname(n)),kmaxo,imax,jmax, &
                                                          decomp%zst(1),decomp%zst(2),decomp%zst(3), &
                                                          buffer_new,n-1)
         enddo ! end n loop
       endif ! initFlow.gt.0

       deallocate(buffer,buffer_new)

       if(myid.eq.0) then
         call ReadHDF5_scalar(fsol,time)
         call WriteHDF5_scalar(fsol_new1,time)
       endif

       call decomp_info_finalize(decompo)
       call decomp_info_finalize(decomp)
       call decomp_2d_finalize

       if(kmaxo.ne.kmax) then ! kmaxo.ne.kmax
         call decomp_2d_init(kmaxo,imax,jmax,p_row,p_col)
         call decomp_info_init(kmaxo,imax,jmax,decompo)
         call decomp_2d_finalize
         call decomp_2d_init(kmax,imax,jmax,p_row,p_col)
         call decomp_info_init(kmax,imax,jmax,decomp)

         call GenIndex(kmaxo,decompo%xsz(2),jmax,kmax,decompo%xsz(2),jmax)
         allocate(buffer(decompo%xsz(1),decompo%xsz(2),decompo%xsz(3)))
         allocate(buffer_new(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)))
         if(initGrid.gt.0) then
           if(myid.eq.0) print *, 'Reading file: ', trim(grd_new1%fname)
           do n=1, grd%dnum
             if(myid.eq.0) print *, 'Reading variable: ', trim(grd%dname(n))
             call Read3dHDF5_svariable_P(trim(grd_new1%fname),grd%dname(n),decompo%xsz(1),decompo%xsz(2),decompo%xsz(3), &
                                                                decompo%xst(1),decompo%xst(2),decompo%xst(3), &
                                                                buffer)
             if(myid.eq.0) print *, 'Interpolating ...'
             call InitVar(buffer,buffer_new)
             call Write3dHDF5_svariable_P(trim(grd_new2%fname),trim(grd%dname(n)),kmax,imax,jmax, &
                                                        decomp%xst(1),decomp%xst(2),decomp%xst(3), &
                                                        buffer_new,n-1)
            enddo ! end n loop
         endif ! initGrid.gt.0

         if(initFlow.gt.0) then
           if(myid.eq.0) print *, 'Reading file: ', trim(fsol_new1%fname)
           do n=1, fsol%dnum
             if(myid.eq.0) print *, 'Reading variable: ', trim(fsol%dname(n))
             call Read3dHDF5_svariable_P(trim(fsol_new1%fname),fsol%dname(n),decompo%xsz(1),decompo%xsz(2),decompo%xsz(3), &
                                                                               decompo%xst(1),decompo%xst(2),decompo%xst(3), &
                                                                               buffer)
             call InitVar(buffer,buffer_new)
             if(myid.eq.0) print *, 'Interpolating ... '
             call Write3dHDF5_svariable_P(trim(fsol_new2%fname),trim(fsol%dname(n)),kmax,imax,jmax, &
                                                            decomp%xst(1),decomp%xst(2),decomp%xst(3), &
                                                            buffer_new,n-1)
           enddo ! end n loop
         endif ! initFlow.gt.0

         deallocate(buffer,buffer_new)
         if(myid.eq.0) then
           call ReadHDF5_scalar(fsol,time)
           call WriteHDF5_scalar(fsol_new2,time)
         endif

         call decomp_info_finalize(decompo)
         call decomp_info_finalize(decomp)
         call decomp_2d_finalize

       endif ! end kmaxo.eq.kmax

     end subroutine ReadGridFlow

     subroutine InitVar(varinput,varoutput)
       implicit none
       real(8), intent(in) :: varinput(:,:,:)
       real(8), intent(out) :: varoutput(:,:,:)

       call Interp3D(dim2input_index,dim3input_index,dim1input_index,varinput, &
                     dim2grid_index,dim3grid_index,dim1grid_index,varoutput)

     end subroutine InitVar

     subroutine Interp3D(xindex_old,yindex_old,zindex_old,profilevars_old, &
                         xindex_new,yindex_new,zindex_new,profilevars_new)
       implicit none
       real(8),dimension(:), intent(in) :: xindex_old, yindex_old, zindex_old
       real(8),dimension(:), intent(in) :: xindex_new, yindex_new, zindex_new
       real(8),dimension(:,:,:), intent(in) :: profilevars_old
       real(8),dimension(:,:,:), intent(out) :: profilevars_new
       integer :: imax_old, jmax_old, kmax_old
       integer :: imax_new, jmax_new, kmax_new
       integer :: i, j, k, n
       real(8),dimension(:,:,:), allocatable :: profilevars_new1
       real(8),dimension(:,:,:), allocatable :: profilevars_new2

       ! Check the dimensions of input arrays
       imax_old = size(xindex_old)
       jmax_old = size(yindex_old)
       kmax_old = size(zindex_old)
       imax_new = size(xindex_new)
       jmax_new = size(yindex_new)
       kmax_new = size(zindex_new)

       allocate( profilevars_new1(kmax_old,imax_new,jmax_old) )
       if(imax_old.eq.1) then
         forall(i=1:imax_new,j=1:jmax_old,k=1:kmax_old)
           profilevars_new1(k,i,j) = profilevars_old(k,1,j)
         end forall
       elseif(imax_new.ne.imax_old) then
         do k = 1, kmax_old
           do j = 1, jmax_old
             call SplineInterp(xindex_old, profilevars_old(k,:,j), imax_old, &
                               xindex_new, profilevars_new1(k,:,j),imax_new,2)
           enddo
         enddo
        else
         profilevars_new1 = profilevars_old
       endif

       allocate( profilevars_new2(kmax_old,imax_new,jmax_new) )
       if(jmax_old.eq.1) then
         forall(i=1:imax_new,j=1:jmax_new,k=1:kmax_old)
           profilevars_new2(k,i,j) = profilevars_new1(k,i,1)
         end forall
       elseif(jmax_old.ne.jmax_new) then
         do k = 1, kmax_old
           do i = 1, imax_new
             call SplineInterp(yindex_old,profilevars_new1(k,i,:), jmax_old, &
                               yindex_new,profilevars_new2(k,i,:), jmax_new,2)
           enddo
         enddo
       else
         profilevars_new2 = profilevars_new1
       endif
       deallocate( profilevars_new1)

       if(kmax_old.eq.1) then
         forall(i=1:imax_new,j=1:jmax_new,k=1:kmax_new)
           profilevars_new(k,i,j) = profilevars_new2(1,i,j)
         end forall
       elseif(kmax_old.ne.kmax_new) then
         do j = 1, jmax_new
           do i = 1, imax_new
             call SplineInterp(zindex_old, profilevars_new2(:,i,j),kmax_old, &
                               zindex_new, profilevars_new(:,i,j), kmax_new,1)
           enddo
         enddo
       else
         profilevars_new = profilevars_new2
       endif
       deallocate( profilevars_new2 )

     end subroutine Interp3D


   end program InitGridFlow_HDF5

