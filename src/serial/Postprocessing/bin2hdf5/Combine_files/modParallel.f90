

    module modParallel
      use decomp_2d
      use decomp_2d_io
      use decomp_2d_fft
      implicit none
      include 'mpif.h'

      integer, parameter :: imax=200, jmax=200, kmax=130
      integer, parameter :: inode=40, jnode=1
      integer, parameter :: nodes=inode*jnode
      integer, parameter :: ibotbc=0, itopbc=1
      integer, parameter :: inletbc=0, ioutletbc=0, istreamwisegeo=0, ispanwisegeo=0
      integer, parameter :: kioplane=4
      logical, dimension(3) :: periodic_bc
      real(8), parameter :: c12i = 1.d0/12.d0

      integer, parameter :: imax_pe=imax/inode+1, jmax_pe=jmax/jnode+1
      integer, parameter :: p_row=inode, p_col=jnode
      integer, parameter :: istencilsize=2, jstencilsize=2, istencil_grid=istencilsize+2, &
                                                            jstencil_grid=jstencilsize+2

      integer, parameter :: idimbeg_var=1-istencilsize, idimend_var=imax_pe+istencilsize &
                          , jdimbeg_var=1-jstencilsize, jdimend_var=jmax_pe+jstencilsize
      integer, parameter :: idimbeg_grid=1-istencil_grid, idimend_grid=imax_pe+istencil_grid &
                          , jdimbeg_grid=1-jstencil_grid, jdimend_grid=jmax_pe+jstencil_grid

      integer :: COMM_2D, COMM_jline, myid_jline, numprocs_jline, COMM_iline, myid_iline, numprocs_iline
      integer :: nid_left, nid_right, nid_bottom, nid_top
      integer :: myid, ierr, errcode
      integer :: status(MPI_STATUS_SIZE)
      integer :: iaux(4), iauxall(4,0:nodes-1)


      type pencil
        integer :: ibeg, iend, jbeg, jend, kbeg, kend
        integer :: ilen, jlen, klen
      end type pencil
      type(pencil) :: xpen, ypen, zpen

contains

      subroutine Init2decomp()
        implicit none
        logical :: periods(2)
        integer :: dims(2), coords(2)
        data periods/2*.true./

        if(inletbc*ioutletbc.ne.0) periods(2) = .false.
        dims(1) = jnode
        dims(2) = inode
        call MPI_INIT(ierr)

        periodic_bc(1) = .false.
        periodic_bc(2) = .true.
        periodic_bc(3) = .true.

        ! init 2decomp with periodic_bc
        call decomp_2d_init(kmax,imax,jmax,p_row,p_col,periodic_bc)

        call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)

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


        xpen%ibeg = xstart(2); xpen%iend = xend(2)
        xpen%jbeg = xstart(3); xpen%jend = xend(3)
        xpen%kbeg = xstart(1); xpen%kend = xend(1)
        xpen%ilen = xsize(2)
        xpen%jlen = xsize(3)
        xpen%klen = xsize(1)

        ypen%ibeg = ystart(1); ypen%iend = yend(1)
        ypen%jbeg = ystart(3); ypen%jend = yend(3)
        ypen%kbeg = ystart(2); ypen%kend = yend(2)
        ypen%ilen = ysize(1)
        ypen%jlen = ysize(3)
        ypen%klen = ysize(2)

        zpen%ibeg = zstart(1); zpen%iend = zend(1)
        zpen%jbeg = zstart(2); zpen%jend = zend(2)
        zpen%kbeg = zstart(3); zpen%kend = zend(3)
        zpen%ilen = zsize(1)
        zpen%jlen = zsize(2)
        zpen%klen = zsize(3)

        iaux(1) = xpen%ibeg; iaux(2) = xpen%iend
        iaux(3) = xpen%jbeg; iaux(4) = xpen%jend
        call MPI_Allgather(iaux,4,MPI_INTEGER,iauxall,4,MPI_INTEGER,COMM_2D,ierr)

      end subroutine Init2decomp

      subroutine Finalize2decomp()
        implicit none

        call decomp_2d_finalize
        call MPI_FINALIZE(ierr)

      end subroutine Finalize2decomp

    end module modParallel
