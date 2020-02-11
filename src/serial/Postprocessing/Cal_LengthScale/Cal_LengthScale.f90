!> using the correlation files to calculate the corresponding length scales.

program Cal_LengthScale
  use modRWHDF5
  implicit none

  integer :: nvarcorrxy = 2, nvarcorrxz = 5
  integer :: i, j, k, ilen, jlen, klen,  num_zone
  integer :: icorrxy, icorrxz, ireadGrid
  integer :: index_i, index_j
  integer, dimension(:), allocatable :: index_be, index_end, index_k
  real(8), dimension(:), allocatable :: th_dim
  real(8), dimension(:), allocatable :: xdiff, ydiff
  real(8), dimension(:,:,:,:), allocatable :: buffer
  real(8), dimension(:,:,:), allocatable :: buffer_grid2d
  real(8) :: dlen
  character(400) :: gridname, fname
  real(8), dimension(:), allocatable :: zgrid

  type tp_fileInfo
    character(400) :: filename
    integer :: imax, jmax, kmax
    integer :: TotZone
    integer :: variable_index
    real(8) :: th
  end type tp_fileInfo

  type(tp_rdwt_hdf5) :: grd
  type(tp_fileInfo) :: fcorrxy, fcorrxz

  call InitHDF5()
  call Input()

  !! read data from corrxy.dat
  if(icorrxy.eq.1) then
    open(7, file=trim(fcorrxy%filename), status='unknown')
      read(7,*)
      do k=1, num_zone
        read(7,'("Zone T=kloc_",I04.4, ",I=",I4, ",J=",I4, ", AUXDATA pave=""",E16.9, """, AUXDATA uave=""",E16.9, """, AUXDATA vave=""",E16.9, """, AUXDATA wave=""",E16.9, """, AUXDATA prms=""",E16.9, """, AUXDATA urms=""",E16.9, """, AUXDATA vrms=""",E16.9, """, AUXDATA wrms=""",E16.9,"""")') &
              index_k(k)
        do j=1, jlen
          do i=1, ilen
            read(7,*) xdiff(i), ydiff(j), buffer(i,j,k,1:nvarcorrxy)
          enddo
        enddo
      enddo
    close(7)

    !! find the location for the maximum value in each wall-normal plane
    call findIndex(xdiff,fcorrxy%imax,index_i)
    call findIndex(ydiff,fcorrxy%jmax,index_j)
    print *, 'The location for the maximum value: ', 'index_i = ', index_i, 'index_j = ', index_j

    !! dimensionalize the threshold value using the maximum value in each wall-normal plane
    th_dim(1:num_zone) = fcorrxy%th*buffer(index_i,index_j,1:num_zone,fcorrxy%variable_index)

    !! read GridMetric fro AveAcoustic_GridMetric.h5
    if(ireadGrid.eq.1) then
      call Read2dGrid(gridname)
    else
      call Read1dGrid(gridname)
    endif

    !! find the threshold location in streamwise direction
    call FindThLoc(ilen,buffer(1:ilen,index_j,1:num_zone,fcorrxy%variable_index),fcorrxy%variable_index)
    dlen = xdiff(2) - xdiff(1)
    fname = 'streamwizeLength.dat'
    call WriteFile(fname,fcorrxy%variable_index,dlen)

    !! find the threshold location in spanwise direction
    call FindThLoc(jlen,buffer(index_i,1:jlen,1:num_zone,fcorrxy%variable_index),fcorrxy%variable_index)
    dlen = ydiff(2) - ydiff(1)
    fname = 'spanwizeLength.dat'
    call WriteFile(fname,fcorrxy%variable_index,dlen)

  endif

  call FinalizeHDF5()

contains

  subroutine Read1dGrid(fn)
    implicit none
    character(*), intent(in) :: fn

    open(7,file=trim(fn),status='unknown')
      read(7,*)
      do k=1, num_zone
        read(7,*) zgrid(k)
      enddo
    close(7)

  end subroutine Read1dGrid

  subroutine Read2dGrid(fn)
    implicit none
    character(*), intent(in) :: fn

    grd%fname = trim(fn)
    grd%gname = '/'
    call DetectHDF5(grd)
    allocate(buffer_grid2d(grd%dimsf(1),grd%dimsf(2),2))
    grd%dnum = 2
    grd%dname(1) = 'x'
    grd%dname(2) = 'z'
    print *, 'Reading grid file: ', trim(grd%fname)
    call ReadHDF5_2D(grd, buffer_grid2d)

    do k=1, num_zone
      zgrid(k) = buffer_grid2d(index_k(k),1,2)
    enddo

  end subroutine Read2dGrid

  subroutine WriteFile(fn,index_in,dlen)
    implicit none
    character(*), intent(in) :: fn
    integer, intent(in) :: index_in
    real(8), intent(in) :: dlen
    integer :: kbe

    kbe = 1
    if(index_in.eq.2) kbe=2

    print *, 'Writing file: ', trim(fn)
    open(7,file=trim(fn),status='unknown')
      rewind(7)
      write(7,'(a)') 'variables=z, length'
      do k=kbe, num_zone
        write(7,*) zgrid(k), (index_end(k)-index_be(k)+1)*(dlen)
      enddo
    close(7)

!    print *, 'Writing file: ', trim('spanwizeLength.dat')
!    open(7,file='spanwizeLength.dat',status='unknown')
!      rewind(7)
!      write(7,'(a)') 'variables=z, length'
!      do k=kbe, num_zone
!        write(7,*) zgrid(k), (index_end(k)-index_be(k)+1)*(ydiff(2)-ydiff(1))
!      enddo
!    close(7)

  end subroutine WriteFile

  subroutine FindThLoc(dim_in,varin, var_index)
    implicit none
    integer, intent(in) :: dim_in, var_index
    real(8), dimension(dim_in,num_zone) :: varin
    integer :: kbe

    index_be  = 0
    index_end = 0
    kbe = 1

    if(var_index.eq.2) then
      kbe = 2
    endif
    do k=kbe, num_zone
      do i=2, dim_in
        if(th_dim(k).ge.varin(i-1,k).and.th_dim(k).lt.varin(i,k)) then
          index_be(k) = i
        elseif(th_dim(k).le.varin(i-1,k).and.th_dim(k).gt.varin(i,k)) then
          index_end(k) = i
        endif
      enddo ! end i loop
    enddo ! end k loop

    ! check results
    do k=kbe, num_zone
      if(index_end(k).eq.0.or.index_be(k).eq.0.or.index_end(k).le.index_be(k)) then
        print *, 'error in finding th location'
        stop
      endif
    enddo

  end subroutine FindThLoc

!  subroutine FindThLoc1()
!    implicit none
!    integer :: kbe
!
!    index_be  = 0
!    index_end = 0
!!    if() then
!!      kbe = 2
!!    endif
!
!!    do k=1, num_zone
!    do k=2, num_zone
!      do j=2, fcorrxy%jmax
!        if(th_dim(k).ge.buffer(index_i,j-1,k,fcorrxy%variable_index).and.&
!           th_dim(k).lt.buffer(index_i,j,k,fcorrxy%variable_index)) then
!          index_be(k) = j
!
!          print *, 'k = ', k
!          print *, 'index_be(k) = ', index_be(k)
!          print *, buffer(index_i,j-1,k,fcorrxy%variable_index), th_dim(k), buffer(index_i,j,k,fcorrxy%variable_index)
!
!        elseif(th_dim(k).le.buffer(index_i,j-1,k,fcorrxy%variable_index).and.&
!               th_dim(k).gt.buffer(index_i,j,k,fcorrxy%variable_index)) then
!          index_end(k) = j
!
!          print *, 'index_end(k)', index_end(k)
!          print *, buffer(index_i,j-1,k,fcorrxy%variable_index), th_dim(k), buffer(index_i,j,k,fcorrxy%variable_index)
!          print *, '#######################'
!        endif
!      enddo ! end j loop
!    enddo ! end k loop
!
!    ! check results
!    do k=2, num_zone
!!      print *, 'k = ', k, 'index_be(k) = ', index_be(k), 'index_end(k)', index_end(k)
!      if(index_end(k).eq.0.or.index_be(k).eq.0.or.index_end(k).le.index_be(k)) then
!        print *, 'error in finding th location'
!        stop
!      endif
!    enddo
!
!  end subroutine FindThLoc1

  subroutine findIndex(var,length,ii)
    implicit none
    integer, intent(in) :: length
    real(8), dimension(length), intent(in) :: var
    integer, intent(out) :: ii

    do i=1, length
      if(var(i).eq.0) then
        ii = i
        exit
      elseif(i.eq.length) then
        print *, 'error in finding index for var(i) = 0 '
        stop
      endif
    enddo

  end subroutine findIndex

  subroutine Input()
    implicit none

!    read(*,*)
!    read(*,*) delta, ztau
    read(*,*)
    read(*,*) icorrxy, icorrxz, ireadGrid
    read(*,*)
    read(*,'(a)') gridname
    read(*,*)
    read(*,*)
    read(*,'(a)') fcorrxy%filename
    read(*,*)
    read(*,*) fcorrxy%imax, fcorrxy%jmax, fcorrxy%TotZone
    read(*,*)
    read(*,*) fcorrxy%variable_index
    read(*,*)
    read(*,*) fcorrxy%th
    read(*,*)
    read(*,*)
    read(*,'(a)') fcorrxz%filename
    read(*,*)
    read(*,*) fcorrxz%jmax, fcorrxz%kmax, fcorrxz%TotZone
    read(*,*)
    read(*,*) fcorrxz%variable_index
    read(*,*)
    read(*,*) fcorrxz%th


    if(icorrxy.eq.1) then
      num_zone = fcorrxy%TotZone
      ilen = fcorrxy%imax
      jlen = fcorrxy%jmax
      allocate(th_dim(num_zone),index_be(num_zone), index_end(num_zone))
      allocate(xdiff(ilen), ydiff(jlen), buffer(ilen,jlen,num_zone,nvarcorrxy))
      allocate(index_k(num_zone))
      allocate(zgrid(num_zone))
    endif
!    if(icorrxz.eq.1) then
!      num_zone = fcorrxz%zone_end - fcorrxz%zone_be + 1
!      if(num_zone.gt.fcorrxz%TotZone) then
!        print *, 'Input zone_be or zone_end # error ! '
!        stop
!      endif
!
!
!    endif

  end subroutine Input




end program Cal_LengthScale





