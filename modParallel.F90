!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  module that initializes MPI parallel computing                    !!
!!  public variables: myid, nid_left, nid_right, nid_bottom, nid_top  !!
!!                    iaux, iauxall, itype_all, jtype_all, COMM_2D    !!
!!                    COMM_jline, COMM_iline, myid_jline, myid_iline  !!
!!  public subroutines: InitMPI, InitParallel                         !!
!!  Description:                                                      !!
!!    This module contains parallel informations that are needed      !!
!!    for communication between processes. It contains subroutines    !!
!!    to decompose the domain and to generate derived data types.     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module MParallel
  !use MPI
  use MFlowData
  implicit none
  include 'mpif.h'
  integer :: myid, numprocs, nid_left,nid_right,nid_bottom,nid_top
  integer :: iaux(10) !, iauxall(10,0:nodes-1)
  integer, dimension(:,:), allocatable :: iauxall
  integer :: itype_single,itype_all,jtype_single,jtype_all
  integer :: COMM_2D, COMM_jline,myid_jline,numprocs_jline&
           , COMM_iline, myid_iline, numprocs_iline
  ! temp variables
  integer :: status(MPI_STATUS_SIZE)
  integer :: req(4),status_array(MPI_STATUS_SIZE,4),ierr
  contains
    subroutine InitMPI(t0)
      real(8), intent(out) :: t0
      call MPI_Init(ierr)
      t0 = MPI_Wtime()
      
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr )

    end subroutine InitMPI

    subroutine InitParallel()
      call InitCoordinates()
      call InitDerivedDataType()
      call MPI_BARRIER( MPI_COMM_WORLD, ierr )
      return
    end subroutine InitParallel

    subroutine InitCoordinates()
      logical periods(2)
      integer dims(2),coords(2)
      data periods/2*.true./
      integer isize(0:inode),jsize(0:jnode)
      integer ierr,ires,jres
      integer i,j,n
      
      allocate(iauxall(10,0:nodes-1))
      if (inletbc*ioutletbc.ne.0) periods(2) = .false.
      dims(1) = jnode
      dims(2) = inode
      call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periods,.false.,COMM_2D,ierr)
      call MPI_CART_GET(COMM_2D,2,dims,periods,coords,ierr)
      call MPI_CART_SHIFT(COMM_2D,1,1,nid_left,nid_right,ierr)
      call MPI_CART_SHIFT(COMM_2D,0,1,nid_bottom,nid_top,ierr)

      call MPI_COMM_SPLIT(COMM_2D,coords(2),1,comm_jline,ierr)
      call MPI_COMM_SIZE(comm_jline,numprocs_jline,ierr)
      call MPI_COMM_RANK(comm_jline,myid_jline,ierr)

      call MPI_COMM_SPLIT(COMM_2D,coords(1),1,comm_iline,ierr)
      call MPI_COMM_SIZE(comm_iline,numprocs_iline,ierr)
      call MPI_COMM_RANK(comm_iline,myid_iline,ierr)

      iaux(1)=coords(1)
      iaux(2)=coords(2)
      iaux(3)=nid_left
      iaux(4)=nid_right
      iaux(5)=nid_bottom
      iaux(6)=nid_top
      ! compute start and end points of each PE
      ! x direction
      i=imax/dims(2)
      ires=imax-i*dims(2)
      do n=0,dims(2)-1
        isize(n)=i
        if(n.lt.ires) isize(n)=i+1
      end do
      iend = 0
      do i = 0,coords(2)
        iend = iend + isize(i)
      end do
      ibeg = iend - isize(coords(2)) + 1
      ! y direction
      j=jmax/dims(1)
      jres=jmax-j*dims(1)
      do n=0,dims(1)-1
        jsize(n)=j
        if(n.lt.jres) jsize(n)=j+1
      end do
      jend = 0
      do j = 0,coords(1)
        jend = jend + jsize(j)
      end do
      jbeg = jend - jsize(coords(1)) + 1

      iaux(7)  = ibeg
      iaux(8)  = iend
      iaux(9)  = jbeg
      iaux(10) = jend
      call MPI_Allgather(iaux,10,MPI_INTEGER,iauxall,10,MPI_INTEGER,COMM_2D,ierr)
      ! Compute the loop indexes
      ilen=iend-ibeg+1
      jlen=jend-jbeg+1
      iend_bc=ilen+istencilsize
      jend_bc=jlen+jstencilsize
      return
    end subroutine InitCoordinates

    subroutine InitDerivedDataType()
      integer :: icount, iblock,istride
      integer :: iblock_itype(nvars), iblock_jtype(nvars)
      integer(MPI_ADDRESS_KIND) :: indicesmpi(nvars), idist(0:3)
      integer :: indices(nvars)
      !integer :: indices(nvars), idist(0:3)

      call GetArrayAddress(indicesmpi, idist)
      indices = indicesmpi
      ! jtype_single (single variable)
      icount  = jstencilsize
      !iblock  = ilen*kmax
      iblock  = ilen*(kmax + 2*kstencilsize )
      istride = idist(1)-idist(0)

      call mpi_Type_Hvector(icount,iblock,istride,MPI_DOUBLE_PRECISION,jtype_single,ierr)
      call mpi_Type_commit(jtype_single, ierr)
      ! jtype_all (mutiple variables)
      iblock_jtype=1
      call mpi_Type_Hindexed(nvars,iblock_jtype,indices,jtype_single,jtype_all,ierr)
      call mpi_Type_commit(jtype_all, ierr)

      ! itype_single
      icount  = jdimlen_bc
      !iblock  = istencilsize*kmax
      iblock  = istencilsize*(kmax + 2*kstencilsize )
      istride = idist(1)-idist(0)
      call mpi_Type_Hvector(icount,iblock,istride,MPI_DOUBLE_PRECISION,itype_single,ierr)
      call mpi_Type_commit(itype_single, ierr)
      ! itype_all
      iblock_itype=1
      call mpi_Type_Hindexed(nvars,iblock_itype,indices,itype_single,itype_all,ierr)
      call mpi_Type_commit(itype_all, ierr)
    end subroutine InitDerivedDataType

!*********************************************************************
!****         compute offset arrays for derived datatypes         ****
!*********************************************************************
    subroutine GetArrayAddress(indx,idis)
      integer(MPI_ADDRESS_KIND), intent(out) :: indx(nvars), idis(0:3)
      !integer, intent(out) :: indx(nvars), idis(0:3)
      !integer :: in, ifoo,n
      integer(MPI_ADDRESS_KIND) :: in, ifoo
      integer :: n
      ! get offsets between variables to be sent
      call MPI_Get_ADDRESS(u(1,1,1),     indx(1), ierr)
      call MPI_Get_ADDRESS(v(1,1,1),     indx(2), ierr)
      call MPI_Get_ADDRESS(w(1,1,1),     indx(3), ierr)
      call MPI_Get_ADDRESS(p(1,1,1),     indx(4), ierr)
      call MPI_Get_ADDRESS(t(1,1,1),     indx(5), ierr)
      call MPI_Get_ADDRESS(rho(1,1,1),   indx(6), ierr)
      call MPI_Get_ADDRESS(e(1,1,1),     indx(7), ierr)
      call MPI_Get_ADDRESS(mu(1,1,1),    indx(8), ierr)
      call MPI_Get_ADDRESS(kappa(1,1,1), indx(9), ierr)
      call MPI_Get_ADDRESS(a   (1,1,1),  indx(10),ierr)

      in=indx(1)
      ifoo=indx(1)
      do n=1,nvars
        indx(n)=indx(n)-in
        if(myid.eq.0)WRITE(*,*)n,indx(n),indx(n)-ifoo
        ifoo=indx(n)
      end do
      ! get offsets between adjacent elements of arrays
      call MPI_Get_ADDRESS(u(1,1,1), idis(0),ierr)
      call MPI_Get_ADDRESS(u(1,1,2), idis(1),ierr)
      call MPI_Get_ADDRESS(u(1,2,1), idis(2),ierr)
      call MPI_Get_ADDRESS(u(2,1,1), idis(3),ierr)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      return
    end subroutine GetArrayAddress
end module MParallel
