program probe_ts
use modRWHDF5
implicit none
! parameters
real(8), parameter :: R = 8314.3d0
!integer,parameter :: nvars_basic = 6
!character(10) :: dname(nvars_basic)
!parameter(dname = (/'u', 'v', 'w', 'p', 'T', 'rho'/))
integer,parameter :: nvars_basic = 5
character(10) :: dname(nvars_basic)
parameter(dname = (/'u', 'v', 'w', 'p', 'T'/))

integer :: kmax_rd
real(8) ::  rbar, Mw
character(10), dimension(:), allocatable :: dname_vol_output
character(500) :: varname_output

! input parameters
integer :: iprobe, iindx, iflow, file_st, file_ed, file_sp, nvars, npt
integer :: ist_jp,ied_jp,kst_jp,ked_jp,nt_jp,num_jp,indx_jp
integer :: jst_ip,jed_ip,kst_ip,ked_ip,nt_ip,num_ip,indx_ip
integer, dimension(:), allocatable :: varindex_vol_output
integer, dimension(:,:), allocatable :: points
character(500) :: fn_path, fn_grd


call InitHDF5()
call InitInput()
call DoSaveProbe()

contains
    subroutine InitInput()
        integer :: i
        read(*,*)
        read(*,'(a)') fn_path
        read(*,*)
        read(*,'(a)') fn_grd
        read(*,*)
        read(*,*) iprobe, iindx, iflow
        read(*,*)
        read(*,*) file_st, file_ed, file_sp
        read(*,*)
        read(*,*) nvars
        allocate(varindex_vol_output(nvars), dname_vol_output(nvars))
        read(*,*) (varindex_vol_output(i), i=1,nvars)
        read(*,*)   ! jplane
        read(*,*)
        read(*,*) ist_jp,ied_jp,kst_jp,ked_jp,num_jp,indx_jp,nt_jp
        read(*,*)   ! iplane
        read(*,*)
        read(*,*) jst_ip,jed_ip,kst_ip,ked_ip,num_ip,indx_ip,nt_ip
        read(*,*)
        read(*,*) npt
        allocate(points(3,npt))
        print*,'Number of saving points = ',npt
        print*,'The probing points are (i,j,k):'
        varname_output = 'time'
        do i=1, npt
          read(*,*) points(:,i)    ! global index probed in the whole flowfield
          print*,points(:,i)
        enddo
        kmax_rd = maxval(points(3,:))
        print*,'kmax_rd: ',kmax_rd
        ! finish reading input

        if(nvars.lt.1 .or. nvars.gt.nvars_basic) then
            print*,'Max # of variables SHOULD be: ',nvars_basic
!            print*,'1, 2, 3, 4, 5, 6'
!            print*,'u, v, w, p, T, rho'
            print*,'1, 2, 3, 4, 5'
            print*,'u, v, w, p, T'
            stop
        endif
        do i = 1,nvars
          dname_vol_output(i) = dname(varindex_vol_output(i))
          varname_output = trim(adjustl(varname_output))//' '//trim(adjustl(dname_vol_output(i)))
        enddo
        print*,'The probing variables are: ', trim(varname_output)
        ! flow property
        Mw = 28.97  ! air
        if(iFlow .eq. 2) Mw = 28.01     ! nitrogen
        rbar = R/Mw
    end subroutine InitInput

    subroutine DoSaveProbe()
        integer :: n,nf,nv
        integer :: idim_jp,kdim_jp ! jplane
        integer :: jdim_ip,kdim_ip ! iplane
        integer,dimension(:),allocatable :: idx_out,jdx_out,kdx_out
        real(8),dimension(:),allocatable :: buffer_time
        real(8),dimension(:,:),allocatable :: buffer_flow
        real(8),dimension(:,:,:,:),allocatable :: buffer_tmp
        real(8),dimension(:,:,:,:,:),allocatable :: buffer_rd
        character(4) :: fnum1, fnum2, fnum3
        character(8) :: fnum08
        character(500) :: fname_rd
        character(500),dimension(:),allocatable :: fn_ts
        type(tp_rdwt_hdf5) :: ftimeseries,fgrd_ts,fsol_ts

        select case(iprobe)
        case(1) ! read in flowdata_xxxxxxxx.h5
            print*,'Not implemented! STOP!'
            stop
        case(2) ! read timeseries_xxxxxxxx.h5/jplane
            ! check index range
            if( maxval(points(1,:)).gt.ied_jp .or. minval(points(1,:)).lt.ist_jp &
            .or.maxval(points(3,:)).gt.ked_jp .or. minval(points(3,:)).lt.kst_jp) then
                print*,'Probe point index out of range. STOP!'
                stop
            endif
            idim_jp = ied_jp - ist_jp + 1
            kdim_jp = ked_jp - kst_jp + 1
            ! initialize grid for timeseries data
            call InitGridHDF5(fgrd_ts)
            fgrd_ts%fname = trim(fn_path)//'timeseries_GridMetrics.h5'
            fgrd_ts%gname = "/jplane"
            fgrd_ts%dimsf(1) = idim_jp
            fgrd_ts%dimsf(2) = kdim_jp
            fgrd_ts%dimsf(3) = num_jp

            fgrd_ts%dimsm(1) = idim_jp
!            fgrd_ts%dimsm(2) = kdim_jp
            fgrd_ts%dimsm(2) = kmax_rd
            fgrd_ts%dimsm(3) = 1

            fgrd_ts%block = fgrd_ts%dimsm

            fgrd_ts%offset = 0
            fgrd_ts%offset(3) = indx_jp - 1
            ! read grid file
            allocate(buffer_tmp(fgrd_ts%dimsm(1),fgrd_ts%dimsm(2),fgrd_ts%dimsm(3),fgrd_ts%dnum))    ! x,y,z
            call ReadtsGridMetrics_Subset(fgrd_ts,buffer_tmp)

            allocate(fn_ts(npt),idx_out(npt),kdx_out(npt))
            do n = 1,npt    ! create probe files
                write(unit=fnum1,fmt='(I04.4)')points(1,n)
                write(unit=fnum2,fmt='(I04.4)')points(2,n)
                write(unit=fnum3,fmt='(I04.4)')points(3,n)
                fn_ts(n) = 'series_i'//fnum1//'_j'//fnum2//'_k'//fnum3//'.dat'
                open(20+n,file = trim(fn_ts(n)),status = 'unknown')
                write(20+n,*) 'Variables = '//trim(varname_output)
                idx_out(n) = points(1,n) - ist_jp + 1
                kdx_out(n) = points(3,n) - kst_jp + 1
                write(20+n,'(A,3E10.3)') '# x,y,z (m): ',buffer_tmp(idx_out(n),kdx_out(n),1,1),buffer_tmp(idx_out(n),kdx_out(n),1,2),&
                                         buffer_tmp(idx_out(n),kdx_out(n),1,3)
                close(20+n)
            enddo
            deallocate(buffer_tmp)

            ! initialize the time-series probe flow data
            call InitFlowHDF5_tsprobe(fsol_ts,varindex_vol_output,nvars)
            fsol_ts%gname = "/jplane"
            fsol_ts%dimsf(1) = nt_jp
            fsol_ts%dimsf(2) = idim_jp
            fsol_ts%dimsf(3) = kdim_jp
            fsol_ts%dimsf(4) = num_jp
            fsol_ts%dimsm = fsol_ts%dimsf
            fsol_ts%dimsm(3) = kmax_rd
            fsol_ts%dimsm(4) = 1
            fsol_ts%block = fsol_ts%dimsm
            fsol_ts%offset = 0
            fsol_ts%offset(4) = indx_jp - 1 ! choose which plane to read

            call InitFlow_1D(ftimeseries)
            ftimeseries%dimsf = (/nt_jp/) !!
            allocate(buffer_time(nt_jp))
            allocate(buffer_rd(fsol_ts%dimsm(1),fsol_ts%dimsm(2),fsol_ts%dimsm(3),fsol_ts%dimsm(4),fsol_ts%dnum))!,buffer_flow(nt_jp,nvars))
            do nf = file_st,file_ed,file_sp  ! loop within available files
                write(unit=fnum08,fmt='(I08.8)') nf
                fname_rd = trim(fn_path)//'timeseries_'//fnum08//'.h5'
                print*,'Reading timeseries file: ',trim(fname_rd)
                ! read flow data
                fsol_ts%fname = trim(fname_rd)
                call Readtsflow_Subset(fsol_ts,buffer_rd)

                ftimeseries%fname = trim(fname_rd)
                call ReadHDF5_1D(ftimeseries,buffer_time)
                do n = 1,npt
                    call WriteDat_timeseries(fn_ts(n),nt_jp,buffer_time,buffer_rd(:,idx_out(n),kdx_out(n),1,:))
                enddo   ! enddo n = 1,npt
            enddo   ! end do nf = file_st,file_ed,file_sp
            deallocate(idx_out,kdx_out)
            if(allocated(buffer_rd)) deallocate(buffer_rd)
            if(allocated(buffer_time)) deallocate(buffer_time)
        case(3) ! read timeseries_xxxxxxxx.h5/iplane
            ! check index range
            if( maxval(points(3,:)).gt.ked_ip .or. minval(points(3,:)).lt.kst_ip &
            .or.maxval(points(2,:)).gt.jed_ip .or. minval(points(2,:)).lt.jst_ip) then
                print*,'Probe point index out of range. STOP!'
                stop
            endif
            jdim_ip = jed_ip - jst_ip + 1
            kdim_ip = ked_ip - kst_ip + 1
            ! initialize grid for timeseries data
            call InitGridHDF5(fgrd_ts)
            fgrd_ts%fname = trim(fn_path)//'timeseries_GridMetrics.h5'
            fgrd_ts%gname = "/iplane"
            fgrd_ts%dimsf(1) = jdim_ip
            fgrd_ts%dimsf(2) = kdim_ip
            fgrd_ts%dimsf(3) = num_ip

            fgrd_ts%dimsm(1) = jdim_ip
!            fgrd_ts%dimsm(2) = kdim_ip
            fgrd_ts%dimsm(2) = kmax_rd
            fgrd_ts%dimsm(3) = 1

            fgrd_ts%block = fgrd_ts%dimsm

            fgrd_ts%offset = 0
            fgrd_ts%offset(3) = indx_ip - 1
            ! read grid file
            allocate(buffer_tmp(fgrd_ts%dimsm(1),fgrd_ts%dimsm(2),fgrd_ts%dimsm(3),fgrd_ts%dnum))    ! x,y,z
            call ReadtsGridMetrics_Subset(fgrd_ts,buffer_tmp)

            allocate(fn_ts(npt),jdx_out(npt),kdx_out(npt))
            do n = 1,npt    ! create probe files
                write(unit=fnum1,fmt='(I04.4)') points(1,n)
                write(unit=fnum2,fmt='(I04.4)') points(2,n)
                write(unit=fnum3,fmt='(I04.4)') points(3,n)
                fn_ts(n) = 'series_i'//fnum1//'_j'//fnum2//'_k'//fnum3//'.dat'
                open(20+n,file = trim(fn_ts(n)),status = 'unknown')
                write(20+n,*) 'Variables = '//trim(varname_output)
                jdx_out(n) = points(2,n) - jst_ip + 1
                kdx_out(n) = points(3,n) - kst_ip + 1
                write(20+n,'(A,3E10.3)') '# x,y,z (m): ',buffer_tmp(jdx_out(n),kdx_out(n),1,1),buffer_tmp(jdx_out(n),kdx_out(n),1,2),&
                                         buffer_tmp(jdx_out(n),kdx_out(n),1,3)
                close(20+n)
            enddo
            deallocate(buffer_tmp)

            ! initialize the time-series probe flow data
            call InitFlowHDF5_tsprobe(fsol_ts,varindex_vol_output,nvars)
            fsol_ts%gname = "/iplane"
            fsol_ts%dimsf(1) = nt_ip
            fsol_ts%dimsf(2) = jdim_ip
            fsol_ts%dimsf(3) = kdim_ip
            fsol_ts%dimsf(4) = num_ip
            fsol_ts%dimsm = fsol_ts%dimsf
            fsol_ts%dimsm(3) = kmax_rd
            fsol_ts%dimsm(4) = 1
            fsol_ts%block = fsol_ts%dimsm
            fsol_ts%offset = 0
            fsol_ts%offset(4) = indx_ip - 1 ! choose which plane to read

            ! initialize the time variable
            call InitFlow_1D(ftimeseries)
            ftimeseries%dimsf = (/nt_ip/) !!
            allocate(buffer_time(fsol_ts%dimsm(1)))
            allocate(buffer_rd(fsol_ts%dimsm(1),fsol_ts%dimsm(2),fsol_ts%dimsm(3),fsol_ts%dimsm(4),fsol_ts%dnum))
            do nf = file_st,file_ed,file_sp  ! loop within available files
                write(unit=fnum08,fmt='(I08.8)') nf
                fname_rd = trim(fn_path)//'timeseries_'//fnum08//'.h5'
                print*,'Reading timeseries file: ',trim(fname_rd)
                ! read flow data
                fsol_ts%fname = trim(fname_rd)
                call Readtsflow_Subset(fsol_ts,buffer_rd)

                ! read time series
                ftimeseries%fname = trim(fname_rd)
                call ReadHDF5_1D(ftimeseries,buffer_time)
                do n = 1,npt
                    call WriteDat_timeseries(fn_ts(n),nt_ip,buffer_time,buffer_rd(:,jdx_out(n),kdx_out(n),1,:))
                enddo   ! enddo n = 1,npt
            enddo   ! end do nf = file_st,file_ed,file_sp
            deallocate(jdx_out,kdx_out)
            if(allocated(buffer_rd)) deallocate(buffer_rd)
            if(allocated(buffer_time)) deallocate(buffer_time)
        case default
            print*,'Unknown case! STOP!'
            stop
        end select
    end subroutine DoSaveProbe

    subroutine WriteDat_timeseries(fn,itmax,buf_time,buf_flow)
        character(*),intent(in) :: fn
        integer,intent(in) :: itmax
        real(8),dimension(itmax),intent(in) :: buf_time
        real(8),dimension(:,:),intent(in) :: buf_flow
        integer :: nvar_buf,n,i
        character(100) :: format_str

        if(size(buf_flow,dim = 1).ne.itmax) then
            print*,'Dimension inconsistent. STOP!'
            stop
        endif
        nvar_buf = size(buf_flow,dim = 2)
        write(format_str,'(A2,I1,A6,A2)') '(',nvar_buf+1,'E23.15',')'
        open(20,file=trim(fn),position='append')
        do i = 1,itmax
!            write(20,'(E23.15)') buf_time(i),(buf_flow(i,n),n=1,nvar_buf)
            write(20,format_str) buf_time(i),(buf_flow(i,n),n=1,nvar_buf)
!            write(20,*) buf_time(i),(buf_flow(i,n),n=1,nvar_buf)
        enddo
        close(20)
    end subroutine WriteDat_timeseries

end program probe_ts
