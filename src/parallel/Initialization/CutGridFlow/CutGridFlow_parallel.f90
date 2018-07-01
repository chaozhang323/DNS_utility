program CutGridFlow_parallel
use decomp_2d
use MPRWHDF5
implicit none

! input parameters
integer :: iCutGrid, iCutFlow, idigfilter, isponge
integer :: imaxo,jmaxo,kmaxo, ibeg, iend, iskip, jbeg, jend, jskip, kbeg, kend, kskip
character(300) :: flowfilename_rd, gridfilename_rd, fn_digfilter_rd, fn_sponge_rd

integer :: niloc,njloc,nkloc,niloc_pe,njloc_pe
integer,dimension(:),allocatable :: iloc,jloc,kloc,iloc_pe,jloc_pe
type(tp_rdwt_hdf5) :: grd, fsol
logical :: IsWriteIncluded

! mpi related
integer :: myid, numprocs, ierr
integer :: p_row,p_col

! Initialize MPI
call MPI_Init(ierr)
! get the number of processes
call MPI_Comm_size(MPI_COMM_WORLD, numprocs, ierr)
! get the individual process ID
call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
!-------------------------------------------------

! initialize HDF5
call InitHDF5()
! Initialize input parameters
call Input()
! initialize 2decomp
call decomp_2d_init(imaxo,jmaxo,kmaxo,p_row,p_col)
if(myid.eq.0) then
    print*,'Z pencil: ',zsize(1),zsize(2),zsize(3)
    print*,' i range: ',zstart(1),zend(1)
    print*,' j range: ',zstart(2),zend(2)
    print*,' k range: ',zstart(3),zend(3)
endif

call InitCut()

if(iCutGrid.eq.1) call CutGrid()
if(iCutFlow.eq.1) call CutFlow()
if(idigfilter.eq.1) call GutDigFilter()
!if(isponge.ge.1) call CutSponge()

!-------------------------------------------------
! close mpi
call decomp_2d_finalize
call MPI_Finalize(ierr)

contains

    subroutine CutGrid()
        real(8),dimension(:,:,:),allocatable :: buffer_in, buffer_cut
        integer :: n,i,j,k,ii,jj
        integer :: nvars = 3
        type(tp_rdwt_hdf5) :: grd_wt

        call InitGridHDF5(grd)
        grd%fname = trim(gridfilename_rd)
        grd%dimsf = (/kmaxo,imaxo,jmaxo/)
        grd%stride = 1
        grd%count =  1
        grd%dimsm(1) = kmaxo
        grd%block(1) = kmaxo
        grd%offset(1) = 0
        grd%dimsm(2) = zsize(1)
        grd%block(2) = zsize(1)
        grd%offset(2) = zstart(1)-1
        grd%dimsm(3) = zsize(2)
        grd%block(3) = zsize(2)
        grd%offset(3) = zstart(2)-1

        call InitGridHDF5(grd_wt)
        grd_wt%fname = 'grid.h5'
        grd_wt%dimsf = (/nkloc,niloc,njloc/)
        if(IsWriteIncluded) then
            grd_wt%offset(1) = 0
            grd_wt%offset(2) = (iloc_pe(1)-ibeg)/iskip
            grd_wt%offset(3) = (jloc_pe(1)-jbeg)/jskip
            grd_wt%block(1) = nkloc
            grd_wt%block(2) = niloc_pe
            grd_wt%block(3) = njloc_pe
            grd_wt%stride = 1
            grd_wt%count = 1
        else
            niloc_pe = 0; njloc_pe = 0
            grd_wt%offset = 0
            grd_wt%stride = 1
            grd_wt%block = 0
            grd_wt%count = 0
        endif
        grd_wt%dimsm = (/nkloc,niloc_pe,njloc_pe/)

        allocate(buffer_in(zsize(3),zsize(1),zsize(2)))
        allocate(buffer_cut(nkloc,niloc_pe,njloc_pe))
        do n = 1,nvars
            if(myid.eq.0) then
                if(n.eq.1) print*,'Reading grid: ',trim(grd%fname)
                if(n.eq.1) print*,'    x ...'
                if(n.eq.2) print*,'    y ...'
                if(n.eq.3) print*,'    z ...'
            endif
            call ReadHDF5_3D_1V(grd,buffer_in,n)
            if(myid.eq.0 .and. n.eq.3) print*,'Finish reading grid file.'
            if(myid.eq.0 .and. n.eq.3) print*,'Writting new grid file.'
            if(IsWriteIncluded) then
                do jj = 1,njloc_pe
                do ii = 1,niloc_pe
                    if(zstart(1).le.iloc_pe(ii).and.zend(1).ge.iloc_pe(ii) .and. &
                       zstart(2).le.jloc_pe(jj).and.zend(2).ge.jloc_pe(jj)) then
                        i = iloc_pe(ii)-zstart(1)+1
                        j = jloc_pe(jj)-zstart(2)+1
                        do k = 1,nkloc
                            buffer_cut(k,ii,jj) = buffer_in(kloc(k),i,j)
                        enddo
                    endif
                enddo
                enddo
            endif
            call WriteHDF5_3D_1V(grd_wt,buffer_cut,n)
        enddo   ! enddo n = 1,nvars
!        if(myid.eq.0) then
!            print*,'    x(1:5): ',buffer_in(10,1:5,1,1)
!            print*,'    y(1:5): ',buffer_in(10,1:5,1,2)
!            print*,'    z(1:5): ',buffer_in(10,1:5,1,3)
!        endif

        deallocate(buffer_in,buffer_cut)
    end subroutine CutGrid

    subroutine CutFlow()
        real(8),dimension(:,:,:),allocatable :: buffer_in, buffer_cut
        real(8) :: time
        integer :: n,i,j,k,ii,jj
        integer :: nvars = 5
        type(tp_rdwt_hdf5) :: fsol_wt

        call InitFlowHDF5(fsol)
        fsol%fname = trim(flowfilename_rd)
        fsol%dimsf = (/kmaxo,imaxo,jmaxo/)
        fsol%stride = 1
        fsol%count =  1
        fsol%dimsm(1) = kmaxo
        fsol%block(1) = kmaxo
        fsol%offset(1) = 0
        fsol%dimsm(2) = zsize(1)
        fsol%block(2) = zsize(1)
        fsol%offset(2) = zstart(1)-1
        fsol%dimsm(3) = zsize(2)
        fsol%block(3) = zsize(2)
        fsol%offset(3) = zstart(2)-1

        call InitFlowHDF5(fsol_wt)
        fsol_wt%fname = 'flowdata_00000000.h5'
        fsol_wt%dimsf = (/nkloc,niloc,njloc/)
        if(IsWriteIncluded) then
            fsol_wt%offset(1) = 0
            fsol_wt%offset(2) = (iloc_pe(1)-ibeg)/iskip
            fsol_wt%offset(3) = (jloc_pe(1)-jbeg)/jskip
            fsol_wt%block(1) = nkloc
            fsol_wt%block(2) = niloc_pe
            fsol_wt%block(3) = njloc_pe
            fsol_wt%stride = 1
            fsol_wt%count = 1
        else
            niloc_pe = 0; njloc_pe = 0
            fsol_wt%offset = 0
            fsol_wt%stride = 1
            fsol_wt%block = 0
            fsol_wt%count = 0
        endif
        fsol_wt%dimsm = (/nkloc,niloc_pe,njloc_pe/)

        allocate(buffer_in(zsize(3),zsize(1),zsize(2)))
        allocate(buffer_cut(nkloc,niloc_pe,njloc_pe))
        do n = 1,nvars
            if(myid.eq.0) then
                if(n.eq.1) print*,'Reading flow data: ',trim(fsol%fname)
                if(n.eq.1) print*,'    u ...'
                if(n.eq.2) print*,'    v ...'
                if(n.eq.3) print*,'    w ...'
                if(n.eq.4) print*,'    p ...'
                if(n.eq.5) print*,'    T ...'
            endif
            call ReadHDF5_3D_1V(fsol,buffer_in,n)
            if(myid.eq.0 .and. n.eq.nvars) print*,'Finish reading flow data file.'
            if(myid.eq.0 .and. n.eq.nvars) print*,'Writting new flow data file.'
            if(IsWriteIncluded) then
                do jj = 1,njloc_pe
                do ii = 1,niloc_pe
                    if(zstart(1).le.iloc_pe(ii).and.zend(1).ge.iloc_pe(ii) .and. &
                       zstart(2).le.jloc_pe(jj).and.zend(2).ge.jloc_pe(jj)) then
                        i = iloc_pe(ii)-zstart(1)+1
                        j = jloc_pe(jj)-zstart(2)+1
                        do k = 1,nkloc
                            buffer_cut(k,ii,jj) = buffer_in(kloc(k),i,j)
                        enddo
                    endif
                enddo
                enddo
            endif
            call WriteHDF5_3D_1V(fsol_wt,buffer_cut,n)
        enddo   ! enddo n = 1,nvars
        deallocate(buffer_in,buffer_cut)
        if(myid.eq.0) then
            fsol%sname = 'time'
            call ReadHDF5_scalar(fsol,time)
            print*,'Write time variable, time = ',time
            call WriteHDF5_scalar(fsol_wt,time)
        endif
    end subroutine CutFlow

    subroutine GutDigFilter()
        implicit none
        real(8),dimension(:,:,:),allocatable :: vars2d,vars_wt
        type(tp_rdwt_hdf5) :: fsol2d,fsol_wt

        if(myid.eq.0) then
            print*,'Cut flowdata to generate digital filter inflow ...'
            print*,'Read file: ',trim(fn_digfilter_rd)
            print*,'inflow i index: ',ibeg
            fsol2d%fname = trim(fn_digfilter_rd)
            fsol2d%gname = '/'
            fsol2d%rank = 2
            fsol2d%dnum = 12
            allocate(fsol2d%dname(fsol2d%dnum),fsol2d%dimsf(fsol2d%rank))
            fsol2d%dname(1) = 'u'
            fsol2d%dname(2) = 'v'
            fsol2d%dname(3) = 'w'
            fsol2d%dname(4) = 'p'
            fsol2d%dname(5) = 'T'
            fsol2d%dname(6) = 'rho'
            fsol2d%dname(7) = 'uu'
            fsol2d%dname(8) = 'vv'
            fsol2d%dname(9) = 'ww'
            fsol2d%dname(10) = 'uv'
            fsol2d%dname(11) = 'uw'
            fsol2d%dname(12) = 'vw'
            fsol2d%dimsf = (/kmaxo,imaxo/)
            fsol2d%IsHSInitialized = .true.
            allocate(vars2d(kmaxo,imaxo,fsol2d%dnum))
            call ReadHDF5_2D_serial(fsol2d,vars2d)
            ! write digital filter
            call InitHDF5_DigFilter(fsol_wt)
            fsol_wt%fname = 'DigFilter_Stat.h5'
            fsol_wt%gname = '/'
            fsol_wt%dimsf = (/kmaxo,1/)     ! needs to be fixed for partition in k direction
            allocate(vars_wt(kmaxo,1,fsol_wt%dnum))
            ! assign data
            vars_wt(:,1,1:3) = vars2d(:,ibeg,1:3)
            vars_wt(:,1,4) = vars2d(:,ibeg,6)
            vars_wt(:,1,5) = vars2d(:,ibeg,5)
            vars_wt(:,1,6:11) = vars2d(:,ibeg,7:12)
            ! write data
            print*,'writing file: DigFilter_Stat.h5'
            call WriteHDF5_2D_serial(fsol_wt,vars_wt)
        endif
    end subroutine GutDigFilter

    subroutine Input()
        implicit none
        integer, parameter :: nid = 5

        if(myid.eq.0) then
            print *, 'Started reading parameters'
            read(nid,*)
            read(nid,*) p_row, p_col, iCutGrid, iCutFlow
            read(nid,*)
            read(nid,*) imaxo, jmaxo, kmaxo
            read(nid,*)
            read(nid,*) ibeg, iend, iskip, jbeg, jend, jskip, kbeg, kend, kskip
            read(nid,*)
            read(nid,'(a)') gridfilename_rd
            read(nid,*)
            read(nid,'(a)') flowfilename_rd
            read(nid,*)
            read(nid,*) idigfilter, isponge
            read(nid,*)
            read(nid,'(a)') fn_digfilter_rd
            read(nid,*)
            read(nid,'(a)') fn_sponge_rd
            if (nid.ne.5) close(nid)
            if(p_row*p_col .ne. numprocs) then
                print*,'Number of processors inconsistent with input value. STOP!'
                print*,'Input: ',p_row*p_col
                print*,'# of processors: ',numprocs
                stop
            endif
            if(ibeg.lt.1 .or. iend.gt.imaxo .or. jbeg.lt.1 .or. jend.gt.jmaxo .or. kbeg.lt.1 .or. kend.gt.kmaxo) then
                print*,'Selected index range exceed original dimension range. STOP!'
                stop
            endif
            print *, 'Finished reading parameters'
        endif

        call MPI_Bcast(imaxo,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(jmaxo,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(kmaxo,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(p_row,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(p_col,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(iCutGrid,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(iCutFlow,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(ibeg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(iend,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(iskip,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(jbeg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(jend,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(jskip,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(kbeg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(kend,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(kskip,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(idigfilter,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(isponge,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(gridfilename_rd, 300, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(flowfilename_rd, 300, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(fn_digfilter_rd, 300, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(fn_sponge_rd, 300, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    end subroutine Input

    subroutine InitCut()
        integer :: i,j,k, ii,jj

        niloc = (iend-ibeg)/iskip+1
        njloc = (jend-jbeg)/jskip+1
        nkloc = (kend-kbeg)/kskip+1
        niloc_pe = zsize(1)/iskip+1     ! max dimension per node
        njloc_pe = zsize(2)/jskip+1     ! max dimension per node
        allocate(iloc(niloc),jloc(njloc),kloc(nkloc))
        allocate(iloc_pe(niloc_pe),jloc_pe(njloc_pe))
        iloc_pe = 0; jloc_pe = 0
        do i = 1,niloc
            iloc(i) = (i-1)*iskip+ibeg
        enddo
        do j = 1,njloc
            jloc(j) = (j-1)*jskip+jbeg
        enddo
        do k = 1,nkloc
            kloc(k) = (k-1)*kskip+kbeg
        enddo

        IsWriteIncluded = .false.
        do i = 1,niloc
        do j = 1,njloc
            if(zstart(1).le.iloc(i).and.zend(1).ge.iloc(i) .and. &
               zstart(2).le.jloc(j).and.zend(2).ge.jloc(j)) IsWriteIncluded = .true.
        enddo
        enddo

        if(IsWriteIncluded) then
            ii = 0; iloc_pe = 0
            do i = 1,niloc
                if(zstart(1).le.iloc(i).and.zend(1).ge.iloc(i)) then
                    ii = ii+1
                    iloc_pe(ii) = iloc(i)
                endif
            enddo
            niloc_pe = ii

            jj = 0; jloc_pe = 0
            do j = 1,njloc
                if(zstart(2).le.jloc(j).and.zend(2).ge.jloc(j)) then
                    jj = jj+1
                    jloc_pe(jj) = jloc(j)
                endif
            enddo
            njloc_pe = jj
        endif

!        print*,myid,zstart(1),zend(1),zstart(2),zend(2),iloc_pe(1),jloc_pe(1)

    end subroutine InitCut

end program CutGridFlow_parallel
