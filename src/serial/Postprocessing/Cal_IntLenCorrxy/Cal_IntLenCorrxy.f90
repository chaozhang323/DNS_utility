program Cal_IntLenCorrxy
use modRWHDF5
implicit none
real(8), parameter :: rbar = 287.d0, gamma = 1.4
real(8) :: pi
! input parameter
integer :: nzones,nheader,imax,jmax,nvars,nvars_all,&
           ilenx,ileny,ithres,iref_grd, icylin
integer,dimension(:),allocatable :: indx,kref_grd
real(8),dimension(:),allocatable :: thres
character(500) :: filename,gridname

call InitHDF5()
call Input()
call Cal_IntLen_by_Corrxy()

contains

subroutine Input()
    implicit none
    integer, parameter :: nid=5
    integer :: n

    read(*,*)
    read(*,*) ilenx, ileny
    read(*,*)
    read(*,*) ithres
    allocate(thres(ithres))
    read(*,*) (thres(n),n=1,ithres)
    read(*,*)
    read(*,'(a)') filename
    read(*,*)
    read(*,*) imax, jmax, iref_grd, icylin
    read(*,*)
    read(*,*) nzones, nvars_all, nheader
    allocate(kref_grd(nzones))
    read(*,*) (kref_grd(n),n=1,nzones)
    read(*,*)
    read(*,*) nvars
    allocate(indx(nvars))
    read(*,*) (indx(n),n=1,nvars)
    read(*,*)
    read(*,'(a)') gridname
    if (nid.ne.5) close(nid)

    pi = 4.0d0*atan(1.d0)
end subroutine Input

subroutine Cal_IntLen_by_Corrxy()
    integer :: n,nh,i,j, iref,jref
    real(8) :: dx,dy
    real(8),dimension(:),allocatable :: vars_tmp
    real(8),dimension(:,:),allocatable :: corr_coef,lenx_scale,leny_scale
    real(8),dimension(:,:,:),allocatable :: vars_sz,z2d
    integer,dimension(:,:),allocatable :: ibe_x,ied_x,ibe_y,ied_y
    character(2) :: fthres
    character(500) :: file_wt
    type(tp_rdwt_hdf5) :: grd

    ! read z2d
    grd%fname = trim(gridname)
    grd%gname='/'
    call DetectHDF5(grd)
    print*,'Reading grid file: ',trim(grd%fname)
    print*,'imax, kmax: ',grd%dimsf(2), grd%dimsf(1)
    grd%rank = 2
    grd%dnum = 1
    if(allocated(grd%dname)) deallocate(grd%dname)
    allocate(grd%dname(grd%dnum))
    grd%dname(1)='z'
    grd%IsHSInitialized = .true.
    allocate(z2d(grd%dimsf(1),grd%dimsf(2),1))
    call ReadHDF5_2D(grd,z2d)

    allocate(vars_tmp(nvars_all))
    allocate(corr_coef(imax,jmax),vars_sz(imax,jmax,nvars))
    allocate(lenx_scale(size(thres),nzones),ibe_x(size(thres),nzones),ied_x(size(thres),nzones))
    allocate(leny_scale(size(thres),nzones),ibe_y(size(thres),nzones),ied_y(size(thres),nzones))

    iref = imax/2+mod(imax,2)
    jref = jmax/2+mod(jmax,2)

    ! initial nheader line: nheader
    ! read the data
    open(33,file=filename,status='old')
    do n=1,nzones
        print*,'Zone #: ',kref_grd(n)
        do nh = 1,nheader
            read(33,*)
        enddo
        do j=1,jmax
        do i=1,imax
            read(33,*) vars_tmp
            do nh=1,nvars
                vars_sz(i,j,nh) = vars_tmp(indx(nh))
            enddo
        enddo
        enddo
        ! calculate correlation coefficient
        do j=1,jmax
        do i=1,imax
            corr_coef(i,j) = vars_sz(i,j,3)/sqrt(abs(vars_sz(iref,jref,4)*vars_sz(i,jref,5))+1.d-30)
        enddo
        enddo
        dx = vars_sz(2,1,1)-vars_sz(1,1,1)
        dy = vars_sz(1,2,2)-vars_sz(1,1,2)
        ! for cylinder, dy=dtheta*r
        if(icylin.eq.1) dy = dy*z2d(kref_grd(n),iref_grd,1)

        ! streamwise integral length scale
        if(ilenx.eq.1) then
            ! calculate integral length scale for each threshold
            call Cal_IntegralLength(corr_coef(:,jref),thres,lenx_scale(:,n),ibe_x(:,n),ied_x(:,n))
            lenx_scale(:,n) = lenx_scale(:,n)*dx
        endif   ! endif ilenx=1

        ! spanwise integral length scale
        if(ileny.eq.1) then
            ! calculate integral length scale for each threshold
            call Cal_IntegralLength(corr_coef(iref,:),thres,leny_scale(:,n),ibe_y(:,n),ied_y(:,n))
            leny_scale(:,n) = leny_scale(:,n)*dy
        endif   ! endif ileny=1
        nheader = 1
    enddo   ! end n=1,nzones
    close(33)
    deallocate(vars_tmp,vars_sz)

    ! write to file
    if(ilenx.eq.1) then
        do n=1,size(thres)
            write(fthres,fmt='(I02.2)') int(thres(n)*100)    !!!
            file_wt = 'IntLength_x_th0p'//fthres//'.dat'
            print*,'Writing ',trim(file_wt)
            open(33,file=file_wt,status='unknown')
            write(33,*) 'Variables=z,length,ibe,iend'
            write(33,'(A,f)') '# Threshold: ',thres(n)
            do nh=1,nzones
                write(33,*) z2d(kref_grd(nh),iref_grd,1),lenx_scale(n,nh),ibe_x(n,nh),ied_x(n,nh)
            enddo
            close(33)
        enddo   ! n=1,size(thres)
    endif   ! endif ilenx=1

    ! write to file
    if(ileny.eq.1) then
        do n=1,size(thres)
            write(fthres,fmt='(I02.2)') int(thres(n)*100)    !!!
            file_wt = 'IntLength_y_th0p'//fthres//'.dat'
            print*,'Writing ',trim(file_wt)
            open(33,file=file_wt,status='unknown')
            write(33,*) 'Variables=z,length,ibe,iend'
            write(33,'(A,f)') '# Threshold: ',thres(n)
            do nh=1,nzones
                write(33,*) z2d(kref_grd(nh),iref_grd,1),leny_scale(n,nh),ibe_y(n,nh),ied_y(n,nh)
            enddo
            close(33)
        enddo   ! n=1,size(thres)
    endif   ! endif ileny=1

end subroutine Cal_IntLen_by_Corrxy

! calculate integral length scale with threshold
subroutine Cal_IntegralLength(coef,thres_in,len_out,ibeg,iend)
    real(8),dimension(:),intent(in) :: coef,thres_in
    real(8),dimension(:),intent(out) :: len_out
    integer,dimension(:),intent(out) :: ibeg,iend
    integer :: i,n,idim,iref,ithres_in

    idim = size(coef)
    iref = idim/2+mod(idim,2)

    ithres_in = size(thres_in)

    if(size(len_out).ne.ithres_in .or. size(ibeg).ne.ithres_in .or. size(iend).ne.ithres_in) then
        print*,'Cal_IntegralLength: array size not consistent.'
        stop
    endif

    do n=1,ithres_in
        ! find the integral limite
        do i=iref,idim
            if(thres_in(n).ge.coef(i) .and. thres_in(n).le.coef(i-1)) then
                iend(n) = i
                ibeg(n) = 2*iref-i
                exit
            endif
        enddo
        ! do integration
        len_out(n) = 0.d0
        do i=ibeg(n),iend(n)-1
            len_out(n) = len_out(n) + 0.5*(coef(i)+coef(i+1))
        enddo
    enddo   ! end do n=1,ithres_in
end subroutine Cal_IntegralLength

end program Cal_IntLenCorrxy
