
!> calculate statistical variables in kplane, jplane and iplane

    program Calstat
      use MTSHDF5
      use MTSkplane
      use MTSjplane
      use MTSiplane
      use MTSStat_kplane
      use MTSStat_jplane
      use MTSStat_iplane
      use MTSAcousticSource
      implicit none

      real(8), parameter :: R=8314.3D0
      integer :: nvar_output
      integer, dimension(:), allocatable :: varindex_output
      integer :: ntpoint_total, npath
      integer :: ireadgrid, ixwindowl, ixwindowr
      real(8) :: dt_sample
      real(8), dimension(:,:,:), allocatable :: bufferout, bufferout_Int, bufferout_ASource

      integer :: ireadkp, ibe_k, iend_k, jbe_k, jend_k, num_kplane
      integer, dimension(:), allocatable :: kplane
      real(8), dimension(:,:,:,:,:), allocatable :: buffer_kplane
      real(8), dimension(:,:,:,:), allocatable :: buffer_kgrid
      integer :: kplane_be, kplane_end

      integer :: ireadjp, ibe_j, iend_j, kbe_j, kend_j, num_jplane
      integer, dimension(:), allocatable :: jplane
      real(8), dimension(:,:,:,:,:), allocatable :: buffer_jplane
      real(8), dimension(:,:,:,:), allocatable :: buffer_jgrid
      integer :: jplane_be, jplane_end

      integer :: ireadip, jbe_i, jend_i, kbe_i, kend_i, num_iplane
      integer, dimension(:), allocatable :: iplane
      real(8), dimension(:,:,:,:,:), allocatable :: buffer_iplane
      real(8), dimension(:,:,:,:), allocatable :: buffer_igrid
      integer :: iplane_be, iplane_end

      type(tp_DNSIndex) :: DNSIndex_k, DNSIndex_j, DNSIndex_i
      type(fprop), dimension(:), allocatable :: fileprop
      integer :: i, j, k, n, nsample, num_file
      character(400) :: filename
      character(8) :: fnum
      character(4) :: fnum_1, fnum_2
      real(8), dimension(:,:,:), allocatable :: src_Philip_kplane, src_Philip_iplane, src_Philip_jplane
      integer :: iCalAcousticSource
      real(8) :: fl, fu, rbar_rd, Rm

      call Input()
      call InitHDF5()

     if(ireadip.eq.1.and.num_iplane.gt.0) then
       ixwindowl = 0; ixwindowr = 0
       select case(iFileType)
          case(0)
            call InitTSiplane(fileprop(1)%ntpoint, npath, DNSIndex_i, jbe_i, jend_i, kbe_i, kend_i, &
                              nvar_output, varindex_output, fileprop, ireadgrid)
            nsample = (fileprop(1)%ntpoint-4)*(nypoint_i-4)

            allocate(buffer_iplane(fileprop(1)%ntpoint,nypoint_i,nzpoint_i,1,nvar_output))
            allocate(buffer_igrid(nypoint_i,nzpoint_i,1,12) )

            if(iCalAcousticSource.eq.0) then
!              print *, 'iCalAcousticSource should be set to 1. STOP'
!              stop
              call InitTSStat_iplane(nsample,nzpoint_i,rbar_rd)
              allocate(bufferout(1,nzpoint_i,nv2d_ipStat+1))
              allocate(bufferout_Int(1,nzpoint_i,nv2d_ipInt+1))

              do i=1, num_iplane
                call ReadTSiplane(iplane(i),buffer_iplane,buffer_igrid)
                bufferout_Int(:,:,nv2d_ipInt+1) = buffer_igrid(:,:,1,3)
                bufferout(1,:,nv2d_ipStat+1) = buffer_igrid(1,:,1,3)      !!!!!!!!!!!!

                !!! FFT filter
               ! do n=1, nvar_output
               ! do j=1, nypoint_i
               !   do k=1, nzpoint_i
               !     call FFTFilter(fileprop(1)%ntpoint,buffer_iplane(:,j,k,1,n),dt_sample,fl,fu,buffer_iplane(:,j,k,1,n))
               !   enddo
               ! enddo
               ! enddo
                !! end FFT filter

                call CalTSDeriv_iplane(fileprop(1)%ntpoint, nypoint_i, nzpoint_i,nvar_output, buffer_iplane(:,:,:,1,:),buffer_igrid(:,:,1,4:12),dt_sample)
                call CalTSStat_ipInt(nzpoint_i, bufferout_Int(:,1,nv2d_ipInt+1),bufferout_Int(1,1:nzpoint_i,1:nv2d_ipInt),buffer_igrid(1,:,1,4:12))
                call CalTSStat_ipStat(nzpoint_i, bufferout)

                filename = 'timeseries_iplane-ave.h5'
                if(i.eq.1) print *, 'Writing file: ', trim(filename)
                print *, 'Writing the iplane: ', iplane(i)
                call WriteTSStat_ipInt(trim(filename),iplane(i),nzpoint_i,bufferout_Int)
                call WriteTSStat_ipStat(trim(filename),iplane(i),nzpoint_i,bufferout(1,:,:))

              enddo ! end i loop

            else ! iCalAcousticSource.eq.1
              call InitTS_AS_iplane(nsample,nzpoint_i)
              allocate(src_Philip_iplane(fileprop(1)%ntpoint,nypoint_i,nzpoint_i))
              allocate(bufferout_ASource(1,nzpoint_i,nv_ASource+2))

              do i=1, num_iplane

                call ReadTSiplane(iplane(i),buffer_iplane,buffer_igrid)

                call CalTS_AS_Philip_iplane(fileprop(1)%ntpoint, nypoint_i, nzpoint_i, nvar_output, &
                                            buffer_iplane(:,1:nypoint_i,1:nzpoint_i,1,1:nvar_output), &
                                            buffer_igrid(1:nypoint_i, 1:nzpoint_i, 1, 4:12), &
                                            src_Philip_iplane(:,1:nypoint_i, 1:nzpoint_i) )
                call CalTSDeriv_AS_iplane(fileprop(1)%ntpoint, nypoint_i, nzpoint_i, nvar_output, &
                                          buffer_iplane(:,1:nypoint_i, 1:nzpoint_i, 1, 1:nvar_output), &
                                          buffer_igrid(1:nypoint_i, 1:nzpoint_i, 1, 4:12), dt_sample, &
                                          src_Philip_iplane(:,1:nypoint_i, 1:nzpoint_i) )
                call CalTS_AS_iplane(nzpoint_i,bufferout_ASource(1,:,1:nv_ASource))
                bufferout_ASource(1,:,nv_ASource+1) = buffer_igrid(1,:,1,2)
                bufferout_ASource(1,:,nv_ASource+2) = buffer_igrid(1,:,1,3)

                filename = 'timeseries_iplane-ASource.h5'
                if(i.eq.1) print *, 'Writing file: ', trim(filename)
                print *, 'Writing the iplane: ', iplane(i)
                call WriteTSASource_iplane(trim(filename),iplane(i),nzpoint_i,bufferout_ASource(:,1:nzpoint_i,1:nv_ASource+2))

              enddo ! end i loop
              deallocate(src_Philip_iplane,bufferout_ASource)
            endif !

          case(1,2)

            call InitVOLTSiplane(fileprop(1)%ntpoint, npath, DNSIndex_i, jbe_i, jend_i, kbe_i, kend_i, &
                      ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop, ireadgrid)
            nsample = (fileprop(1)%ntpoint-4)*(nypoint_i-4)

            allocate(buffer_iplane(fileprop(1)%ntpoint,nypoint_i,nzpoint_i,1,nvar_output))
            allocate(buffer_igrid(nypoint_i,nzpoint_i,1,12))

            if(iCalAcousticSource.eq.0) then
              call InitTSStat_iplane(nsample,nzpoint_i,rbar_rd)

              allocate(bufferout(1,nzpoint_i,nv2d_ipStat+1))
              allocate(bufferout_Int(1,nzpoint_i,nv2d_ipInt+1))

              num_file = (fileprop(1)%file_end-fileprop(1)%file_be)/fileprop(1)%file_skip + 1
              do n=1, num_file
                do i=iplane_be, iplane_end
                  write(unit=fnum,fmt='(I08.8)') fileprop(1)%file_be + (n-1)*fileprop(1)%file_skip
                  filename = trim(fileprop(1)%filepath)//'timeseries_'//fnum//'.h5'
                  if(i-iplane_be.eq.0) print *, 'Reading file: ', trim(filename)
                  call ReadVOLTSiplane(filename, i, buffer_iplane, buffer_igrid)

                  bufferout_Int(:,:,nv2d_ipInt+1) = buffer_igrid(:,:,1,3)
                  bufferout(1,:,nv2d_ipStat+1) = buffer_igrid(1,:,1,3)      !!!!!!!!!!!!
                  call CalTSDeriv_iplane(fileprop(1)%ntpoint, nypoint_i, nzpoint_i,nvar_output, buffer_iplane(:,:,:,1,:),buffer_igrid(:,:,1,4:12),dt_sample)
                  call CalTSStat_ipInt(nzpoint_i, bufferout_Int(:,1,nv2d_ipInt+1),bufferout_Int(1,1:nzpoint_i,1:nv2d_ipInt),buffer_igrid(1,:,1,4:12))
                  call CalTSStat_ipStat(nzpoint_i, bufferout)

                  filename = 'timeseries_iplane-ave'//fnum//'.h5'
                  if(i.eq.1) print *, 'Writing file: ', trim(filename)
                  print *, 'Writing the iplane: ', iplane(i)
                  call WriteTSStat_ipInt(trim(filename),iplane(i),nzpoint_i,bufferout_Int)
                  call WriteTSStat_ipStat(trim(filename),iplane(i),nzpoint_i,bufferout(1,:,:))
                enddo ! end i loop
              enddo ! end n loop
              deallocate(bufferout,bufferout_Int)
!        allocate(buffer_iplane(fileprop(1)%ntpoint,nypoint_i,nzpoint_i,1,nvar_output))
!        allocate(buffer_igrid(nypoint_i,nzpoint_i,1,12))
            else ! iCalAcousticSource.eq.1
              call InitTS_AS_iplane(nsample,nzpoint_i)
              allocate(src_Philip_iplane(fileprop(1)%ntpoint,nypoint_i,nzpoint_i))
              allocate(bufferout_ASource(1,nzpoint_i,nv_ASource+2))

              num_file = (fileprop(1)%file_end-fileprop(1)%file_be)/fileprop(1)%file_skip + 1
              do n=1, num_file
                do i=1, num_iplane
                  write(unit=fnum,fmt='(I08.8)') fileprop(1)%file_be + (n-1)*fileprop(1)%file_skip
                  filename = trim(fileprop(1)%filepath)//'timeseries_'//fnum//'.h5'
                  if(i.eq.1) print *, 'Reading file: ', trim(filename)
                  call ReadVOLTSiplane(filename, i, buffer_iplane, buffer_igrid)

           ! test
           ! print *, 'T = ', buffer_iplane(1,1,1,1,5)
           ! print *, '', buffer_igrid(1,1,1,4)
           ! print *, '', buffer_igrid(1,1,1,5)
           ! print *, '', buffer_igrid(1,1,1,6)
           ! print *, '', buffer_igrid(1,1,1,7)
           ! print *, '', buffer_igrid(1,1,1,8)
           ! print *, '', buffer_igrid(1,1,1,9)
           ! print *, '', buffer_igrid(1,1,1,10)
           ! print *, '', buffer_igrid(1,1,1,11)
           ! print *, '', buffer_igrid(1,1,1,12)


                  call CalTS_AS_Philip_iplane(fileprop(1)%ntpoint,nypoint_i,nzpoint_i,nvar_output, &
                                              buffer_iplane(:,1:nypoint_i,1:nzpoint_i,1,1:nvar_output), &
                                              buffer_igrid(1:nypoint_i,1:nzpoint_i,1,4:12), &
                                              src_Philip_iplane(:,1:nypoint_i,1:nzpoint_i))
                  call CalTSDeriv_AS_iplane(fileprop(1)%ntpoint,nypoint_i,nzpoint_i,nvar_output, &
                                            buffer_iplane(:,1:nypoint_i,1:nzpoint_i,1,1:nvar_output), &
                                            buffer_igrid(1:nypoint_i,1:nzpoint_i,1,4:12), dt_sample, &
                                            src_Philip_iplane(:,1:nypoint_i,1:nzpoint_i))
                  call CalTS_AS_iplane(nzpoint_i,bufferout_ASource(1,:,1:nv_ASource))
                  bufferout_ASource(1,:,nv_ASource+1) = buffer_igrid(1,:,1,2)
                  bufferout_ASource(1,:,nv_ASource+2) = buffer_igrid(1,:,1,3)

                  filename = 'timeseries_iplane-ASource'//fnum//'.h5'
                  if(i.eq.1) print *, 'Writing file: ', trim(filename)
                  print *, 'Writing the iplane: ', iplane(i)
                  call WriteTSASource_iplane(trim(filename),iplane(i),nzpoint_i,bufferout_ASource(:,1:nzpoint_i,1:nv_ASource+2))
                enddo ! end i loop
              enddo ! end n loop
              deallocate(src_Philip_iplane,bufferout_ASource)
            endif ! end iCalAcousticSource.eq.0

          case default
            print *, 'iFileType parameter out of range'
            stop
        end select

        deallocate(buffer_iplane,buffer_igrid)
      endif ! end ireadjp.eq.1


      if(ireadjp.eq.1.and.num_jplane.gt.0) then
        ixwindowl = 2; ixwindowr = 2
        select case(iFileType)
          case(0)
            print *, 'iFileType should be set to 1. '
          case(1,2)
            call InitVOLTSJplane(fileprop(1)%ntpoint, npath, DNSIndex_j, ibe_j, iend_j, kbe_j, kend_j, &
                      ixwindowl, ixwindowr, 0, 0, nvar_output, varindex_output, fileprop, ireadgrid)
            nsample = (fileprop(1)%ntpoint-4)
            allocate(buffer_jplane(fileprop(1)%ntpoint,imints_j:imaxts_j,nzpoint_j,1,nvar_output))
            allocate(buffer_jgrid(imints_j:imaxts_j,nzpoint_j,1,12))

            if(iCalAcousticSource.eq.0) then
              call InitTSStat_jplane(nsample, nzpoint_j)

              allocate(bufferout(nzpoint_j,nxpoint_j,nv2d_jpStat+2))
              allocate(bufferout_Int(nzpoint_j,nxpoint_j,nv2d_jpInt+2))

              num_file = (fileprop(1)%file_end-fileprop(1)%file_be)/fileprop(1)%file_skip + 1

              do n=1, num_file
                do j=jplane_be, jplane_end
                  write(unit=fnum,fmt='(I08.8)') fileprop(1)%file_be + (n-1)*fileprop(1)%file_skip
                  filename = trim(fileprop(1)%filepath)//'timeseries_'//fnum//'.h5'
                  if(j-jplane_be.eq.0) print *, 'Reading file: ', trim(filename)

                  call ReadVOLTSJplane(filename, j, buffer_jplane, buffer_jgrid)

                  do i=1, nxpoint_j
                    bufferout(:,i,nv2d_jpStat+1)  = buffer_jgrid(i,:,1,1)
                    bufferout(:,i,nv2d_jpStat+2)  = buffer_jgrid(i,:,1,3)
                    bufferout_Int(:,i,nv2d_jpInt+1) = buffer_jgrid(i,:,1,1)
                    bufferout_Int(:,i,nv2d_jpInt+2) = buffer_jgrid(i,:,1,3)
                    call CalTSDeriv_jplane(fileprop(1)%ntpoint, nzpoint_j,nvar_output, buffer_jplane(:,i-2:i+2,:,1,:), buffer_jgrid(i,:,1,4:12), dt_sample)

                    call CalTSStat_jpInt(nzpoint_j, bufferout_Int(:,i,nv2d_jpInt+2),bufferout_Int(1:nzpoint_j,i,1:nv2d_jpInt),buffer_jgrid(i,:,1,4:12))
                    call CalTSStat_jpStat(nzpoint_j, bufferout(:,i,:))
                  enddo ! end i loop
                  filename = 'timeseries_jplane-ave'//fnum//'.h5'
                  if(j.eq.1) print *, 'Writing file: ', trim(filename)
                  print *, 'Writing jplane: ', jplane(j)
                  call WriteTSStat_jpInt(trim(filename),jplane(j),nxpoint_j,nzpoint_j,bufferout_Int)
                  call WriteTSSTat_jpStat(trim(filename), jplane(j),nxpoint_j,nzpoint_j,bufferout)
                enddo ! end j loop
              enddo ! end n loop
              deallocate(bufferout,bufferout_Int)
            else ! iCalAcousticSource.eq.1
              call InitTS_AS_jplane(nsample,nzpoint_j)
              allocate(src_Philip_jplane(fileprop(1)%ntpoint,1:nxpoint_j,nzpoint_j))
              allocate(bufferout_ASource(nzpoint_j,3:nxpoint_j-2,nv_ASource+2))         !!!!!!!!!

              num_file = (fileprop(1)%file_end-fileprop(1)%file_be)/fileprop(1)%file_skip + 1

              do n=1, num_file
                do j=1, num_jplane
                  write(unit=fnum,fmt='(I08.8)') fileprop(1)%file_be + (n-1)*fileprop(1)%file_skip
                  filename = trim(fileprop(1)%filepath)//'timeseries_'//fnum//'.h5'
                  if(j.eq.1) print *, 'Reading file: ', trim(filename)

                  call ReadVOLTSJplane(filename, j, buffer_jplane, buffer_jgrid)
                  call CalTS_AS_Philip_jplane(fileprop(1)%ntpoint, nxpoint_j+4, nzpoint_j,nvar_output, &
                                buffer_jplane(:,imints_j:imaxts_j,1:nzpoint_j,1,1:nvar_output), buffer_jgrid(imints_j:imaxts_j,1:nzpoint_j,1,4:12), &
                                src_Philip_jplane(:,1:nxpoint_j,1:nzpoint_j))

                  do i=3, nxpoint_j-2
                    bufferout_ASource(:,i,nv_ASource+1)  = buffer_jgrid(i,:,1,1)
                    bufferout_ASource(:,i,nv_ASource+2)  = buffer_jgrid(i,:,1,3)
                    call CalTSDeriv_AS_jplane(fileprop(1)%ntpoint, nzpoint_j, nvar_output, buffer_jplane(:,i-2:i+2,1:nzpoint_j,1,:), &
                          buffer_jgrid(i,1:nzpoint_j,1,4:12), dt_sample, src_Philip_jplane(:,i-2:i+2,1:nzpoint_j) )
                    call CalTS_AS_jplane(nzpoint_j,bufferout_ASource(:,i,1:nv_ASource))
                  enddo ! end i loop
                  filename = 'timeseries_jplane-ASource'//fnum//'.h5'
                  if(j.eq.1) print *, 'Writing file: ', trim(filename)
                  print *, 'Writing jplane: ', jplane(j)
                  call WriteTSASource_jplane(trim(filename), jplane(j),nxpoint_j-4,nzpoint_j,bufferout_ASource(1:nzpoint_j,3:nxpoint_j-2,1:nv_ASource+2))
                enddo ! end j loop
              enddo ! end n loop
              deallocate(src_Philip_jplane, bufferout_ASource)
            endif ! end iCalAcousticSource.eq.0

          case default
            print *, 'iFileType parameter out of range'
            stop
        end select

        deallocate(buffer_jplane,buffer_jgrid)
      endif ! end ireadjp.eq.1

      if(ireadkp.eq.1.and.num_kplane.gt.0) then

!        ixwindowl = 2; ixwindowr = 2
        ixwindowl = 2; ixwindowr = 2

        select case(iFileType)
          case(0) ! read data from the converted kplane data
           call InitTSkplane(fileprop(1)%ntpoint,npath, DNSIndex_k, ibe_k, iend_k, jbe_k, jend_k, &
                  ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop, ireadgrid)
            nsample = (nypoint_k-4)*(fileprop(1)%ntpoint-4)

            allocate(buffer_kplane(fileprop(1)%ntpoint, nypoint_k, imints_k:imaxts_k, 1, nvar_output))
            allocate(buffer_kgrid(nypoint_k, imints_k:imaxts_k, 1, 12))

            if(iCalAcousticSource.eq.0) then
              call InitTSStat(nsample)
              allocate(bufferout_Int(nxpoint_k,1,nv2d_kpInt+1))
              allocate(bufferout(nxpoint_k,1,nv2d_kpStat+2))

              write(unit=fnum_1,fmt='(I04.4)') ibe_k
              write(unit=fnum_2,fmt='(I04.4)') iend_k

              do k=kplane_be, kplane_end
                call ReadTSkplane(kplane(k),buffer_kplane,buffer_kgrid,k-1)

                do i=1, nxpoint_k
                  call CalTSDeriv_kplane(fileprop(1)%ntpoint, nypoint_k, nvar_output, buffer_kplane(:,:,i-2:i+2,1,:), buffer_kgrid(:,i,1,4:12), dt_sample)
                  if(kplane(k).eq.1) then
                    bufferout_Int(i,:,nv2d_kpInt+1) = buffer_kgrid(1,i,1,1) ! x

                    call CalTSStat_kpInt(bufferout_Int(i,1,1:nv2d_kpInt),buffer_kgrid(1,i,1,4:12))

                  endif
                  bufferout(i,1,nv2d_kpStat+1) = buffer_kgrid(1,i,1,1) ! x
                  bufferout(i,1,nv2d_kpStat+2) = buffer_kgrid(1,i,1,3) ! z
                  call CalTSStat_kpStat(bufferout(i,1,1:nv2d_kpStat))
                enddo ! end i loop
                filename = 'timeseries_kplane'//fnum_1//'-'//fnum_2//'.h5'
                if(k-kplane_be.eq.0) print *, 'Writing file: ', trim(filename)
                print *, 'Writing kplane: ', kplane(k)
                if(kplane(k).eq.1) then
                  call WriteTSStat_kpInt(trim(filename), kplane(k), nxpoint_k, bufferout_Int(1:nxpoint_k,1,1:nv2d_kpInt+1))
                endif
                call WriteTSStat_kpStat(trim(filename), kplane(k), nxpoint_k, bufferout(1:nxpoint_k,1,1:nv2d_kpStat+2))

              enddo ! end k loop

            else ! iCalAcousticSource.eq.1
              print *, 'iCalAcousticSource should be set to 0. STOP '
              stop
            endif ! end iCalAcousticSource.eq.0

!            call InitTSStat(nsample)
!!           allocate( buffer_kplane(ntpoint, nypoint_k, imints_k:imaxts_k, 1, nvar_output) )
!            allocate(buffer_kplane(ntpoint_total, nypoint_k, imints_k:imaxts_k, 1, nvar_output))
!            allocate(buffer_kgrid(nypoint_k, imints_k:imaxts_k, 1, 12))
!            allocate(bufferout(nxpoint_k,num_kplane,nv2dp2))
!            ! Reading k plane data
!            do k=1, num_kplane
!              call ReadTSkplane(kplane(k),buffer_kplane,buffer_kgrid)
!
!              do i = 1, nxpoint_k
!
!!                 call CalTSDeriv_kplane(ntpoint_total, nypoint_k, nvar_output, buffer_kplane(:,:,i-2:i+2,1,:), buffer_kgrid(:,i,1,4:12), dt_sample)
!                 call CalTSStat(bufferout(i,k,1:nv2d))
!                 bufferout(i,k,nv2d+1) = buffer_kgrid(1,i,1,1) ! x
!                 bufferout(i,k,nv2d+2) = buffer_kgrid(1,i,1,3) ! z
!              enddo
!            enddo
!            call  WriteTSStat_kplane(nxpoint_k, num_kplane, bufferout)
!
          case(1,2) ! directly read the data from the timeseries volume

            call InitVOLTSkplane(ntpoint_total, npath, DNSIndex_k, ibe_k, iend_k, jbe_k, jend_k, &
                     ixwindowl, ixwindowr, nvar_output, varindex_output, fileprop, ireadgrid)
            nsample = (nypoint_k-4)*(fileprop(1)%ntpoint-4)
            allocate(buffer_kplane(fileprop(1)%ntpoint, nypoint_k, imints_k:imaxts_k, 1, nvar_output))
            allocate(buffer_kgrid(nypoint_k, imints_k:imaxts_k, 1, 12))

            if(iCalAcousticSource.eq.0) then  ! iCalAcousticSource.eq.0
              call InitTSStat(nsample)
              allocate(bufferout_Int(nxpoint_k,1,nv2d_kpInt+1))
              allocate(bufferout(nxpoint_k,1,nv2d_kpStat+2))

              num_file = (fileprop(1)%file_end-fileprop(1)%file_be)/fileprop(1)%file_skip + 1
              do n=1, num_file
                do k=kplane_be, kplane_end
                  write(unit=fnum,fmt='(I08.8)') fileprop(1)%file_be + (n-1)*fileprop(1)%file_skip
                  filename = trim(fileprop(1)%filepath)//'timeseries_'//fnum//'.h5'
                  if(k-kplane_be.eq.0) print *, 'Reading file: ', trim(filename)
                  call ReadVOLTSkplane_perFile(filename,k,buffer_kplane,buffer_kgrid)

                  do i = 1, nxpoint_k
                    call CalTSDeriv_kplane(fileprop(1)%ntpoint, nypoint_k, nvar_output, buffer_kplane(:,:,i-2:i+2,1,:),buffer_kgrid(:,i,1,4:12),dt_sample)
                    if(kplane(k).eq.1) then
                      bufferout_Int(i,:,nv2d_kpInt+1) = buffer_kgrid(1,i,1,1) ! x
                      call CalTSStat_kpInt(bufferout_Int(i,1,1:nv2d_kpInt),buffer_kgrid(1,i,1,4:12))
                    endif
                    bufferout(i,1,nv2d_kpStat+1) = buffer_kgrid(1,i,1,1) ! x
                    bufferout(i,1,nv2d_kpStat+2) = buffer_kgrid(1,i,1,3) ! z
                    call CalTSStat_kpStat(bufferout(i,1,1:nv2d_kpStat))
                  enddo ! end i loop

                  filename = 'timeseries_kplane-ave'//fnum//'.h5'
                  if(k-kplane_be.eq.0) print *, 'Writing file: ', trim(filename)
                  print *, 'Writing kplane: ', kplane(k)
                  if(kplane(k).eq.1) then
                    call WriteTSStat_kpInt(trim(filename), kplane(k), nxpoint_k, bufferout_Int(1:nxpoint_k,1,1:nv2d_kpInt+1))
                  endif
                  call WriteTSStat_kpStat(trim(filename), kplane(k), nxpoint_k, bufferout(1:nxpoint_k,1,1:nv2d_kpStat+2))
                enddo ! end k loop
              enddo ! end n loop

              deallocate(bufferout,bufferout_Int)

            else  ! iCalAcousticSource.eq.1

!             allocate(buffer_kplane(fileprop(1)%ntpoint, nypoint_k, imints_k:imaxts_k, 1, nvar_output))
              call InitTS_AS_kplane(nsample)
              allocate(src_Philip_kplane(fileprop(1)%ntpoint,nypoint_k,1:nxpoint_k))
              allocate(bufferout_ASource(3:nxpoint_k-2,1,nv_Asource+2) )

              num_file = (fileprop(1)%file_end-fileprop(1)%file_be)/fileprop(1)%file_skip + 1

              do n=1, num_file
                do k=1, num_kplane
                  write(unit=fnum,fmt='(I08.8)') fileprop(1)%file_be + (n-1)*fileprop(1)%file_skip
                  filename = trim(fileprop(1)%filepath)//'timeseries_'//fnum//'.h5'
                  if(k.eq.1) print *, 'Reading file: ', trim(filename)
                  call ReadVOLTSkplane_perFile(filename,k,buffer_kplane,buffer_kgrid)
                  call CalTS_AS_Philip_kplane(fileprop(1)%ntpoint, nxpoint_k+4, nypoint_k, nvar_output, &
                                              buffer_kplane(:,1:nypoint_k,imints_k:imaxts_k,1,:), buffer_kgrid(1:nypoint_k,imints_k:imaxts_k,1,4:12), &
                                              src_Philip_kplane(:,1:nypoint_k,1:nxpoint_k))

                  do i = 3, nxpoint_k-2 !
                    bufferout_ASource(i,1,nv_ASource+1) = buffer_kgrid(1,i,1,1) ! x
                    bufferout_ASource(i,1,nv_ASource+2) = buffer_kgrid(1,i,1,3) ! z
                    call CalTSDeriv_AS_kplane(fileprop(1)%ntpoint, nypoint_k, nvar_output, buffer_kplane(:,:,i-2:i+2,1,:), &
                                              buffer_kgrid(:,i,1,4:12),dt_sample,src_Philip_kplane(:,:,i-2:i+2))
                    call CalTS_AS_kplane(bufferout_ASource(i,1,1:nv_ASource))
                  enddo ! end i loop

                  filename = 'timeseries_kplane-ASource'//fnum//'.h5'
                  if(k.eq.1) print *, 'Writing file: ', trim(filename)
                  print *, 'Writing kplane: ', kplane(k)
                  call WriteTSASource_kplane( trim(filename), kplane(k), nxpoint_k-4, bufferout_ASource(3:nxpoint_k-2,1,1:nv_ASource+2) )
                enddo ! end k loop
              enddo ! end n loop
              deallocate(src_Philip_kplane,bufferout_ASource)
            endif ! end iCalAcousticSource.eq.0

          case default
            print *, 'iFileType parameter out of range'
            stop
        end select
        deallocate(buffer_kplane,buffer_kgrid)
      endif

  contains
      subroutine Input()
        implicit none
        integer :: i, n

        npath = 1
        read(*,*)
        read(*,*) dt_sample, uinf, rhoinf, fl, fu
        read(*,*)
        read(*,*) Rm
        read(*,*)
        read(*,*) iFileType, iCalAcousticSource
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
        read(*,*) ireadkp
        read(*,*)
        read(*,*) ibe_k, iend_k, jbe_k, jend_k  !, dzdk_kp
        read(*,*)
        read(*,*) DNSIndex_k%ibe, DNSIndex_k%iend, DNSIndex_k%iskip, &
                  DNSIndex_k%jbe, DNSIndex_k%jend, DNSIndex_k%jskip, DNSIndex_k%ibuffer
        read(*,*)
        read(*,*) num_kplane, kplane_be, kplane_end
        allocate(kplane(num_kplane))
        read(*,*) (kplane(i), i=1, num_kplane)
        read(*,*)
        read(*,*)
        read(*,*) ireadjp
        read(*,*)
        read(*,*) ibe_j, iend_j, kbe_j, kend_j
        read(*,*)
        read(*,*) DNSIndex_j%ibe, DNSIndex_j%iend, DNSIndex_j%iskip, &
                  DNSIndex_j%kbe, DNSIndex_j%kend, DNSIndex_j%kskip
        read(*,*)
        read(*,*) num_jplane, jplane_be, jplane_end
        allocate(jplane(num_jplane))
        read(*,*) (jplane(i), i=1, num_jplane)
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
!        ntpoint_total = 0
!        if(ireadVOL.eq.0) then
!          do n=1, npath
!            fileprop(n)%ntpoint = fileprop(n)%ntpoint - fileprop(n)%nskip
!            ntpoint_total = ntpoint_total + fileprop(n)%ntpoint
!          enddo
!        elseif(ireadVOL.eq.1) then
!          do n=1, npath
!            ntpoint_total = ntpoint_total + ( (fileprop(n)%file_end-fileprop(n)%file_be)/fileprop(n)%file_skip + 1 )*fileprop(n)%ntpoint
!          enddo
!        endif

        ireadgrid = 1
        nvar_output = 10
        allocate(varindex_output(nvar_output))
        do i = 1, nvar_output
          varindex_output(i) = i
        enddo

        if(ireadkp.eq.1) then
          if(ibe_k.lt.DNSIndex_k%ibe+2.or.iend_k.gt.DNSIndex_k%iend-2) then
            print *, 'Streamwise index out of bound. iwindowl or iwindowr need to be at least 2. STOP'
            print *, 'ibe_k =', ibe_k, 'iend_k =', iend_k
            print *, 'ibe_DNS =', DNSIndex_k%ibe, 'iend_DNS =', DNSIndex_k%iend
            STOP
          endif
        endif
        if(ireadjp.eq.1) then
          if(ibe_j.lt.DNSIndex_j%ibe+2.or.iend_j.gt.DNSIndex_j%iend-2) then !!!!! iwin!!!!!!!!!!!!!!!!!!!!!!!
            print *, 'Streamwise index out of bound. iwindowl or iwindowr need to be at least 2. STOP'
            print *, 'ibe_j =', ibe_j, 'iend_j =', iend_j
            print *, 'ibe_DNS =', DNSIndex_j%ibe, 'iend_DNS =', DNSIndex_j%iend
            STOP
          endif
        endif
!        if(ireadip.eq.1) then
!          if(jbe_i.lt.DNSIndex_i%ibe+0.or.iend_i.gt.DNSIndex_i%iend-0) then !!!!! iwin!!!!!!!!!!!!!!!!!!!!!!!
!            print *, 'Streamwise index out of bound. iwindowl or iwindowr need to be at least 2. STOP'
!            print *, 'ibe_i =', ibe_i, 'iend_i =', iend_i
!            print *, 'ibe_DNS =', DNSIndex_i%ibe, 'iend_DNS =', DNSIndex_i%iend
!            STOP
!          endif
!        endif
         rbar_rd = R/Rm
         print *, 'rbar = ', rbar_rd
      end subroutine Input

    end program Calstat
