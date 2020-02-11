
!> calculate statistical variables in kplane, jplane and iplane

    program CalAngle
      use MTSHDF5
      use MTSkplane
      use MTSjplane
      use MTSiplane
      use MTSAngle
      implicit none

      integer :: nvar_output
      integer, dimension(:), allocatable :: varindex_output
      integer :: ntpoint_total, npath, ireadVOL
      integer :: ireadgrid, ixwindowl, ixwindowr
      real(8) :: dt_sample
      real(8), dimension(:,:,:), allocatable :: bufferout, bufferout_Int, bufferout_AN

      integer :: ireadip, jbe_i, jend_i, kbe_i, kend_i, num_iplane
      integer, dimension(:), allocatable :: iplane
      real(8), dimension(:,:,:,:,:), allocatable :: buffer_iplane
      real(8), dimension(:,:,:,:), allocatable :: buffer_igrid
      integer :: iplane_be, iplane_end

      type(tp_DNSIndex) :: DNSIndex_i
      type(fprop), dimension(:), allocatable :: fileprop
      integer :: i, j, k, n, nsample, num_file
      character(400) :: filename
      character(8) :: fnum
      character(4) :: fnum_1, fnum_2
      integer :: iCalAngle

      call Input()
      call InitHDF5()

     if(ireadip.eq.1.and.num_iplane.gt.0) then
       ixwindowl = 0; ixwindowr = 0
       select case(ireadVOL)
          case(0)
            call InitTSiplane(fileprop(1)%ntpoint, npath, DNSIndex_i, jbe_i, jend_i, kbe_i, kend_i, &
                              nvar_output, varindex_output, fileprop, ireadgrid)
            nsample = (fileprop(1)%ntpoint)*(nypoint_i)

            allocate(buffer_iplane(fileprop(1)%ntpoint,nypoint_i,nzpoint_i,1,nvar_output))
            allocate(buffer_igrid(nypoint_i,nzpoint_i,1,12) )


              call InitTS_AN_iplane(nsample,nzpoint_i)

              allocate(bufferout_AN(1,nzpoint_i,nv_Angle+1))

              do i=1, num_iplane

                call ReadTSiplane(iplane(i),buffer_iplane,buffer_igrid)

                call CalTS_AN_iplane(nzpoint_i,buffer_iplane(:,:,:,1,:),bufferout_AN(1,:,1:nv_Angle))

                bufferout_AN(1,:,nv_Angle+1) = buffer_igrid(1,:,1,3)

                filename = 'timeseries_iplane-Angle.h5'
                if(i.eq.1) print *, 'Writing file: ', trim(filename)
                print *, 'Writing the iplane: ', iplane(i)
                call WriteTSAngle_iplane(trim(filename),iplane(i),nzpoint_i,bufferout_AN(:,1:nzpoint_i,1:nv_Angle+1))

              enddo ! end i loop
              deallocate(bufferout_AN)


          case(1)

!            call InitVOLTSiplane(fileprop(1)%ntpoint, npath, DNSIndex_i, jbe_i, jend_i, kbe_i, kend_i, &
!                      ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop, ireadgrid)
!            nsample = (fileprop(1)%ntpoint-4)*(nypoint_i-4)
!
!            allocate(buffer_iplane(fileprop(1)%ntpoint,nypoint_i,nzpoint_i,1,nvar_output))
!            allocate(buffer_igrid(nypoint_i,nzpoint_i,1,12))
!
!            if(iCalAcousticSource.eq.0) then
!              call InitTSStat_iplane(nsample,nzpoint_i)
!
!              allocate(bufferout(1,nzpoint_i,nv2d_ipStat+1))
!              allocate(bufferout_Int(1,nzpoint_i,nv2d_ipInt+1))
!
!              num_file = (fileprop(1)%file_end-fileprop(1)%file_be)/fileprop(1)%file_skip + 1
!              do n=1, num_file
!                do i=iplane_be, iplane_end
!                  write(unit=fnum,fmt='(I08.8)') fileprop(1)%file_be + (n-1)*fileprop(1)%file_skip
!                  filename = trim(fileprop(1)%filepath)//'timeseries_'//fnum//'.h5'
!                  if(i-iplane_be.eq.0) print *, 'Reading file: ', trim(filename)
!                  call ReadVOLTSiplane(filename, i, buffer_iplane, buffer_igrid)
!
!                  bufferout_Int(:,:,nv2d_ipInt+1) = buffer_igrid(:,:,1,3)
!                  bufferout(1,:,nv2d_ipStat+1) = buffer_igrid(1,:,1,3)      !!!!!!!!!!!!
!                  call CalTSDeriv_iplane(fileprop(1)%ntpoint, nypoint_i, nzpoint_i,nvar_output, buffer_iplane(:,:,:,1,:),buffer_igrid(:,:,1,4:12),dt_sample)
!                  call CalTSStat_ipInt(nzpoint_i, bufferout_Int(:,1,nv2d_ipInt+1),bufferout_Int(1,1:nzpoint_i,1:nv2d_ipInt),buffer_igrid(1,:,1,4:12))
!                  call CalTSStat_ipStat(nzpoint_i, bufferout)
!
!                  filename = 'timeseries_iplane-ave'//fnum//'.h5'
!                  if(i.eq.1) print *, 'Writing file: ', trim(filename)
!                  print *, 'Writing the iplane: ', iplane(i)
!                  call WriteTSStat_ipInt(trim(filename),iplane(i),nzpoint_i,bufferout_Int)
!                  call WriteTSStat_ipStat(trim(filename),iplane(i),nzpoint_i,bufferout(1,:,:))
!                enddo ! end i loop
!              enddo ! end n loop
!              deallocate(bufferout,bufferout_Int)
!
!            else ! iCalAcousticSource.eq.1
!              call InitTS_AS_iplane(nsample,nzpoint_i)
!              allocate(src_Philip_iplane(fileprop(1)%ntpoint,nypoint_i,nzpoint_i))
!              allocate(bufferout_ASource(1,nzpoint_i,nv_ASource+2))
!
!              num_file = (fileprop(1)%file_end-fileprop(1)%file_be)/fileprop(1)%file_skip + 1
!              do n=1, num_file
!                do i=1, num_iplane
!                  write(unit=fnum,fmt='(I08.8)') fileprop(1)%file_be + (n-1)*fileprop(1)%file_skip
!                  filename = trim(fileprop(1)%filepath)//'timeseries_'//fnum//'.h5'
!                  if(i.eq.1) print *, 'Reading file: ', trim(filename)
!                  call ReadVOLTSiplane(filename, i, buffer_iplane, buffer_igrid)
!
!
!                  call CalTS_AS_Philip_iplane(fileprop(1)%ntpoint,nypoint_i,nzpoint_i,nvar_output, &
!                                              buffer_iplane(:,1:nypoint_i,1:nzpoint_i,1,1:nvar_output), &
!                                              buffer_igrid(1:nypoint_i,1:nzpoint_i,1,4:12), &
!                                              src_Philip_iplane(:,1:nypoint_i,1:nzpoint_i))
!                  call CalTSDeriv_AS_iplane(fileprop(1)%ntpoint,nypoint_i,nzpoint_i,nvar_output, &
!                                            buffer_iplane(:,1:nypoint_i,1:nzpoint_i,1,1:nvar_output), &
!                                            buffer_igrid(1:nypoint_i,1:nzpoint_i,1,4:12), dt_sample, &
!                                            src_Philip_iplane(:,1:nypoint_i,1:nzpoint_i))
!                  call CalTS_AS_iplane(nzpoint_i,bufferout_ASource(1,:,1:nv_ASource))
!                  bufferout_ASource(1,:,nv_ASource+1) = buffer_igrid(1,:,1,2)
!                  bufferout_ASource(1,:,nv_ASource+2) = buffer_igrid(1,:,1,3)
!
!                  filename = 'timeseries_iplane-ASource'//fnum//'.h5'
!                  if(i.eq.1) print *, 'Writing file: ', trim(filename)
!                  print *, 'Writing the iplane: ', iplane(i)
!                  call WriteTSASource_iplane(trim(filename),iplane(i),nzpoint_i,bufferout_ASource(:,1:nzpoint_i,1:nv_ASource+2))
!                enddo ! end i loop
!              enddo ! end n loop
!              deallocate(src_Philip_iplane,bufferout_ASource)
!            endif ! end iCalAcousticSource.eq.0

          case default
            print *, 'ireadVOL parameter out of range'
            stop
        end select

        deallocate(buffer_iplane,buffer_igrid)
      endif ! end ireadjp.eq.1


  contains
      subroutine Input()
        implicit none
        integer :: i, n

        npath = 1
        read(*,*)
        read(*,*) dt_sample, uinf, rhoinf, Tinf
        read(*,*)
        read(*,*) ireadVOL
        allocate(fileprop(npath))
        do n=1, npath
          read(*,*)
          read(*,*) fileprop(n)%ntpoint, fileprop(n)%nskip
          read(*,*)
          read(*,*) fileprop(n)%file_be, fileprop(n)%file_end, fileprop(n)%file_skip
          read(*,*)
          read(*,'(a)') fileprop(n)%filepath
        enddo
        read(*,*)
        read(*,*)
        read(*,*) ireadip
        read(*,*)
        read(*,*) jbe_i, jend_i, kbe_i, kend_i
        read(*,*)
        read(*,*) DNSIndex_i%jbe, DNSIndex_i%jend, DNSIndex_i%jskip, &
                  DNSIndex_i%kbe, DNSIndex_i%kend, DNSIndex_i%kskip, DNSIndex_i%jbuffer
        read(*,*)
        read(*,*) num_iplane, iplane_be, iplane_end
        allocate(iplane(num_iplane))
        read(*,*) (iplane(i), i=1, num_iplane)

        ireadgrid = 1
        nvar_output = 3
        allocate(varindex_output(nvar_output))

        varindex_output(1) = 1 ! u
        varindex_output(2) = 3 ! w
        varindex_output(3) = 4 ! p

      end subroutine Input

    end program CalAngle
