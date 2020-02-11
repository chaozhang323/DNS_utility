!> using the streamwise-spanwise correlation (corrxy) files to calculate the corresponding length scales.
!> using the streanwise-wall-normal correlation (corrxz) files to calculate the angle

program Cal_corrxz_angle
  use modRWHDF5
  implicit none

  integer :: i, j, k, ilen, jlen, klen, ntrans, num_zone
  integer :: icorrxy, icorrxz, icohLength, ireadGrid
  integer :: index_i, index_j
  integer, dimension(:), allocatable :: index_be, index_end, index_k
  real(8), dimension(:), allocatable :: xdiff, ydiff, z, ff
  real(8), dimension(:,:,:,:), allocatable :: buffer
  real(8), dimension(:,:,:), allocatable :: buffer_grid2d
  real(8) :: dlen, pi, tmp
  character(400) :: gridname, fname
  real(8), dimension(:), allocatable :: zgrid
  real(8), dimension(:), allocatable :: IntLen

  real(8), dimension(:), allocatable   :: WORK
  integer :: M, N, NRHS, INFO, LWORK
  integer :: LDA, LDB, LWMAX

  type tp_fileInfo
    character(1000) :: filename
    integer :: nvar, dim1, dim2, dim3
    integer :: TotZone, nheaders
    integer :: variable_index
    real(8) :: th, iIntegralLen
  end type tp_fileInfo

  real(8), dimension(:,:,:), allocatable :: coh

  type(tp_rdwt_hdf5) :: grd
  type(tp_fileInfo) :: fcorrxy, fcorrxz, fcoh

  integer :: kbe = 3, flimit = 300, index_total = 10000  !!!!
  real(8) :: Lz
  real(8), dimension(:,:), allocatable :: coh_new
  real(8), dimension(:), allocatable :: xdiff_new

  call InitHDF5()
  call Input()

  pi = 4.d0*atan(1.)

  !! read coherence data
  if(icohLength.eq.1) then
    num_zone = fcoh%TotZone
    ntrans  = fcoh%dim1
    ilen    = fcoh%dim2
    jlen    = fcoh%dim3
    allocate(index_be(ntrans),index_end(ntrans))
    allocate(ff(ntrans),xdiff(ilen),ydiff(jlen),buffer(ntrans,ilen,jlen,1))
    allocate(coh(ntrans,ilen,jlen))

    print *, 'Reading file: ', trim(fcoh%filename)
    open(7,file=trim(fcoh%filename),status='old')
      do i=1, fcoh%nheaders
        read(7,*)
      enddo
      do i=1, ilen
        do j=1, jlen
          do k=1, ntrans
            read(7,*) ff(k), xdiff(i), coh(k,i,j)
          enddo
        enddo
      enddo ! end k loop
    close(7)

    print *, 'x = ', xdiff(1:10)

    !! find the location for the maximum value in each frequency
    call findIndex(xdiff,ilen,index_i)
    print *, 'The location for the maximum value: ', 'index = ', index_i

    print *, 'coh(5,1,1) = ', coh(5,1,1)
    print *, 'coh(5,2,1) = ', coh(5,2,1)
    print *, 'coh(3,index_i:index_i+5,1) = ', coh(3,index_i:index_i+5,1)

    allocate(IntLen(flimit))
    allocate(coh_new(flimit,index_total))
    allocate(xdiff_new(index_total))

    xdiff_new = 0.d0
    coh_new = 0.d0
    dlen = xdiff(2) - xdiff(1)

    !print *, '', ilen-index_i+1
    print *, 'dlen = ', dlen

    do i=1, index_total
      xdiff_new(i) = (i-1)*dlen
    enddo

    do k=kbe, flimit
      Lz = Integral(k,ilen)
      !Lz = Integral(k,index_i+10)       !!!!! change
      print *, 'Lz = ', Lz
      do i=1, index_total
        coh_new(k,i) = exp( -(xdiff_new(i)/Lz))
      enddo ! end i loop
      IntLen(k) = Integral_new(k,index_total)
    enddo ! end k loop
    !print *, 'xdiff_new(657) = ', xdiff_new(657)

    fname = 'IntLen.dat'
    call WriteFile_Int(trim(fname),IntLen)

    print *, 'Writing file: coh_new.dat '
    open(7,file='coh_new.dat',status='unknown')
      write(7,'(a)')  'variables=f,x,coh'
      write(7,'("Zone T=coh",",I=",I4, ",J=",I8 )') flimit-kbe+1, index_total
      do i=1, index_total
        do k=kbe, flimit
          write(7,*) ff(k), xdiff_new(i), coh_new(k,i)
        enddo
      enddo
    close(7)

  endif

  call FinalizeHDF5()

contains

  function Integral(dim1,dim2)
    integer, intent(in) :: dim1, dim2
    real(8) :: Integral
    real(8) :: tmp
    integer :: i

    tmp = 0.d0
    do i=index_i+1, dim2-1
      tmp = tmp + ( coh(dim1,i,1) + coh(dim1,i-1,1) )*dlen/2.d0
    enddo
    Integral = tmp

  end function Integral

  function Integral_new(dim1,dim2)
    integer, intent(in) :: dim1, dim2
    real(8) :: Integral_new
    real(8) :: tmp
    integer :: i

    tmp = 0.d0
    do i=2, dim2
      tmp = tmp + ( coh_new(dim1,i) + coh_new(dim1,i-1) )*dlen/2.d0
    enddo

    Integral_new = tmp

  end function Integral_new

  subroutine WriteFile_Int(fn,IntLen)
    character(*), intent(in) :: fn
    real(8), dimension(flimit), intent(in) :: IntLen

    print *, 'Writing file: ', trim(fn)
    open(7,file=trim(fn),status='unknown')
      write(7,'(a)') 'variables=f, length'
      do k=kbe, flimit
        write(7,*) ff(k), IntLen(k)
      enddo
    close(7)
  end subroutine WriteFile_Int

  subroutine findIndex(var,length,ii)
    implicit none
    integer, intent(in) :: length
    real(8), dimension(length), intent(in) :: var
    integer, intent(out) :: ii

    do i=1, length
      if(var(i).ge.0) then
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

    read(*,*)
    read(*,*) icorrxy, icorrxz, icohLength, ireadGrid
    read(*,*)
    read(*,'(a)') gridname
    read(*,*)
    read(*,*)
    read(*,'(a)') fcorrxy%filename
    read(*,*)
    read(*,*) fcorrxy%dim1, fcorrxy%dim2, fcorrxy%nvar, fcorrxy%TotZone, fcorrxy%nheaders
    read(*,*)
    read(*,*) fcorrxy%variable_index
    read(*,*)
    read(*,*) fcorrxy%th, fcorrxy%iIntegralLen
    read(*,*)
    read(*,*)
    read(*,'(a)') fcorrxz%filename
    read(*,*)
    read(*,*) fcorrxz%dim1, fcorrxz%dim2, fcorrxz%nvar, fcorrxz%TotZone, fcorrxz%nheaders
    read(*,*)
    read(*,*) fcorrxz%variable_index
    read(*,*)
    read(*,*) fcorrxz%th
    read(*,*)
    read(*,*)
    read(*,'(a)') fcoh%filename
    read(*,*)
    read(*,*) fcoh%dim1, fcoh%dim2, fcoh%dim3, fcoh%nvar, fcoh%TotZone, fcoh%nheaders
    read(*,*)
    read(*,*) fcoh%variable_index
    read(*,*)
    read(*,*) fcoh%th

  end subroutine Input




end program Cal_corrxz_angle





