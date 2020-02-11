! program to compute coherence using timeseries data

program Calcoh
  use decomp_2d
  use MPRWHDF5
  use modSpect
  implicit none

  integer :: myid, numprocs, ierr

  integer :: nskip, kloc_ref, ixwindowl, ixwindowr, jywindowr
  integer :: nperiod, ioverlap, ntrans, nsection
  integer :: num_kplane, iwindow
  integer :: kplane_be, kplane_end, kplane_skip
  integer, dimension(:), allocatable :: kplane
  integer :: n, ibe_k, iend_k, jbe_k, jend_k
  type(fprop), dimension(:), allocatable :: fileprop
  character(400) :: fname
  character(4) :: fnum4, fnum4_1, fnum4_2, fnum4_3, fnum4_4
  integer, parameter :: nvar_output = 1
  integer :: varindex_output(nvar_output)
  character(3) :: varname_output(nvar_output)
  integer :: ii, i, j, k, kk, nn, m, jj, jjj
  real(8) :: dt_sample, dx_kp, dy_kp
  complex, dimension(:,:,:), allocatable :: spect12_tmp1, spect11_tmp1, spect22_tmp1
  complex, dimension(:,:,:), allocatable :: spect12, spect11, spect22
  integer, dimension(:,:), allocatable ::  numpt
  real(8) :: fperiod, corrx
  integer :: numptx
  real(8) :: zloc_ref, dzdk_ref, zloc_kp, dzdk_kp
  integer :: ntpoint_total
  real(8), dimension(:,:,:,:,:), allocatable :: buffer_kplane
  type(tp_DNSIndex) :: DNSIndex_k
  integer :: npath
  integer :: kmax_rd, imax_rd, jmax_rd, ishift, ist_rd, iend_rd, ilen_rd
  integer :: iave_st, iave_end
  logical, dimension(3) :: periodic_bc
  integer :: p_row, p_col
  type(DECOMP_INFO) :: decomp
  integer :: dim1_line, dim2_line, ntpoint_tmp
  integer :: icorrxt, iCalCoh

   ! initialize MPI
   call MPI_INIT(ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
   call InitHDF5()
   if(myid.eq.0) print *, 'Started reading parameters'
   call Input()
   if(myid.eq.0) print *, 'Finished reading parameters'

   call InitSpect1d(ntpoint_total, nsection, iwindow, ioverlap, dt_sample, ntrans, nperiod)
   fperiod = (ntrans-1)*dt_sample

   imax_rd = iend_k - ibe_k + 1
   jmax_rd = jend_k - jbe_k + 1
   periodic_bc = (/.true.,.false.,.false./)

   if(myid.eq.0) then
     print *, 'imax_rd', imax_rd, 'jmax_rd = ', jmax_rd
     print *, '# of data points per segments =', ntrans
     print *, '# of segments =', nsection
     print *, 'Sampling frequency (Herz) =', 1./dt_sample
     print *, 'Interval of sampling (sec) =', dt_sample
     print *, 'Length per segments (sec) =', (ntrans-1)*dt_sample
     if(ioverlap.eq.0) print *, 'Segments with no overlap'
     if(ioverlap.eq.1) print *, 'Segments with half overlap'
     select case(iwindow)
      case(0)  !flat top window
        print*,'calculating flat top window'
      case(1)  !Hanning
        print*,'calculating hanning window'
      case(2)  !modified Hanning with flat top
        print*,'calculating modified hanning window'
      case(3) ! Hamming
        print*,'calculating hamming window'
      case default
        print*,'unknow window type with ntype = ', ntype
        stop
      end select
   endif


   if(icorrxt.eq.1) then
     p_row = numprocs
     p_col = 1

     call decomp_2d_init(jmax_rd, jmax_rd, jmax_rd, p_row, p_col)
     call decomp_info_init(jmax_rd, jmax_rd, jmax_rd, decomp)
     call MPI_Cart_Sub(DECOMP_2D_COMM_CART_X,(/.true.,.false./),dim2_line,ierr)

     ! Indexes ranges for reading TIMESERIES data
     ishift   = ibe_k - DNSIndex_k%ibe + 1
     ist_rd   = ishift - ixwindowl
     iend_rd  = ishift + imax_rd + ixwindowr - 1
     ilen_rd  = iend_rd - ist_rd + 1
     iave_st  = ishift
     iave_end = ishift + imax_rd - 1

     if(myid.eq.0) print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
     if(myid.eq.0) print *, 'ishift = ', ishift, ' iave_st = ', iave_st, 'iave_end = ', iave_end, 'ist_rd = ', ist_rd, 'iend_rd = ', iend_rd

     allocate(spect12_tmp1(ntrans,-ixwindowl:ixwindowr,1) )
     allocate(spect12(ntrans,-ixwindowl:ixwindowr,1) )
     allocate(numpt(-ixwindowl:ixwindowr,1))
     spect12_tmp1 = (0.d0,0.d0); spect12 = (0.d0,0.d0); numpt = 0

     allocate( buffer_kplane(ntpoint_total,jmax_rd,ist_rd:iend_rd,1,1) )

     do k=1, num_kplane
       write(unit=fnum4_1,fmt='(I04.4)') kplane(k)
       write(unit=fnum4_2,fmt='(I04.4)') ibe_k
       write(unit=fnum4_3,fmt='(I04.4)') iend_k
       write(unit=fnum4_4,fmt='(I04.4)') nsection
       do kk=1, DNSIndex_k%nkplane
         if(kplane(k).eq.DNSIndex_k%kplane(kk)) then
           kplane_be = kk
         endif
       enddo ! end kk loop
       call ReadFiles(kplane_be)
       if(myid.eq.0) print *, 'wall-normal location k =', kplane(k)
       if(myid.eq.0) print *, 'streamwise average range i =', ibe_k, iend_k

       do j=decomp%xst(2), decomp%xen(2)
         do i=iave_st, iave_end
           do ii=-ixwindowl, ixwindowr
             spect12_tmp1(:,ii,1) = spect12_tmp1(:,ii,1) + crossspect1d(buffer_kplane(1:ntpoint_total,j,i+ii,1,1),buffer_kplane(1:ntpoint_total,j,i,1,1))
              numpt(ii,1) = numpt(ii,1) + 1
           enddo ! end ii loop
         enddo
       enddo ! end j loop

       do ii=-ixwindowl, ixwindowr
         spect12_tmp1(:,ii,1) = spect12_tmp1(:,ii,1)/dble(numpt(ii,1))
       enddo

       call MPI_Allreduce(spect12_tmp1,spect12,(ntrans)*(ixwindowl+ixwindowr+1),&
                           MPI_DOUBLE_COMPLEX,MPI_SUM,dim2_line,ierr)

       spect12 = spect12/dble(numprocs)

       if(myid.eq.0) then
         fname = '_k'//fnum4_1//'_i'//fnum4_2//'-'//fnum4_3//'_nsection'//fnum4_4
         print *, 'writing file: ', trim(fname)
         print *, 'ixwindowl,ixwindowr = ', ixwindowl,ixwindowr
         call CalCorrxt(ixwindowl,ixwindowr,dx_kp,spect12(1:ntrans,-ixwindowl:ixwindowr,1),fname)
!         open(unit=12,file=trim(fname),status='unknown')
!           rewind(12)
!           write(12,'(a)') 'variables = f,x,y,spect12_r,spect12_i,spect11,spect22'
!           write(12,*) 'zone T = k001, i =', ntrans/2, ' j =', (ixwindowl + ixwindowr + 1), ' k =', jywindowr+1
!           do jj=0, jywindowr
!             do ii = -ixwindowl,ixwindowr
!               do i = 1,ntrans/2
!                 write(12,*)  dble(i-1)/fperiod, dble(ii)*dx_kp, dble(jj)*dy_kp, real(spect12(i,ii,jj)),imag(spect12(i,ii,jj)), abs(spect11(i,ii,jj)),abs(spect22(i,ii,jj))
!               enddo
!             enddo
!           enddo
!         close(12)
       endif

     enddo ! end k loop
   elseif(iCalCoh.eq.1) then
     p_row = numprocs
     p_col = 1
     !call decomp_2d_init(ntrans, jmax_rd, imax_rd, p_row, p_col)
     !call decomp_info_init(ntrans, jmax_rd, imax_rd, decomp)
     call decomp_2d_init(jmax_rd, jmax_rd, jmax_rd, p_row, p_col)
     call decomp_info_init(jmax_rd, jmax_rd, jmax_rd, decomp)
     call MPI_Cart_Sub(DECOMP_2D_COMM_CART_X,(/.true.,.false./),dim2_line,ierr)

     ! Indexes ranges for reading TIMESERIES data
     ishift   = ibe_k - DNSIndex_k%ibe + 1
     ist_rd   = ishift - ixwindowl
     iend_rd  = ishift + imax_rd + ixwindowr - 1
     ilen_rd  = iend_rd - ist_rd + 1
     iave_st  = ishift
     iave_end = ishift + imax_rd - 1

     if(myid.eq.0) print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
     if(myid.eq.0) print *, 'ishift = ', ishift, ' iave_st = ', iave_st, 'iave_end = ', iave_end, 'ist_rd = ', ist_rd, 'iend_rd = ', iend_rd

     allocate(spect12_tmp1(ntrans,-ixwindowl:ixwindowr,0:jywindowr), spect11_tmp1(ntrans,-ixwindowl:ixwindowr,0:jywindowr), spect22_tmp1(ntrans,-ixwindowl:ixwindowr,0:jywindowr) )
     allocate(spect12(ntrans,-ixwindowl:ixwindowr,0:jywindowr),  spect11(ntrans,-ixwindowl:ixwindowr,0:jywindowr),  spect22(ntrans,-ixwindowl:ixwindowr,0:jywindowr) )
     allocate(numpt(-ixwindowl:ixwindowr,0:jywindowr))
     spect12_tmp1 = (0.d0,0.d0); spect11_tmp1 = (0.d0,0.d0); spect22_tmp1 = (0.d0,0.d0); numpt = 0

     allocate( buffer_kplane(ntpoint_total,jmax_rd,ist_rd:iend_rd,1,1) )

     do k=1, num_kplane
       write(unit=fnum4_1,fmt='(I04.4)') kplane(k)
       write(unit=fnum4_2,fmt='(I04.4)') ibe_k
       write(unit=fnum4_3,fmt='(I04.4)') iend_k
       write(unit=fnum4_4,fmt='(I04.4)') nsection
       do kk=1, DNSIndex_k%nkplane
         if(kplane(k).eq.DNSIndex_k%kplane(kk)) then
           kplane_be = kk
         endif
       enddo ! end kk loop
       call ReadFiles(kplane_be)
       if(myid.eq.0) print *, 'wall-normal location k =', kplane(k)
       if(myid.eq.0) print *, 'streamwise average range i =', ibe_k, iend_k

       do j=decomp%xst(2), decomp%xen(2)
       do i=iave_st, iave_end
         do jj=0, jywindowr
           jjj = modulo(j+jj,jmax_rd)
           if(jjj.eq.0) jjj = jmax_rd
           do ii=-ixwindowl, ixwindowr
             spect12_tmp1(:,ii,jj) = spect12_tmp1(:,ii,jj) + crossspect1d(buffer_kplane(1:ntpoint_total,j,i,1,1),buffer_kplane(1:ntpoint_total,jjj,i+ii,1,1))
             spect11_tmp1(:,ii,jj) = spect11_tmp1(:,ii,jj) +  autospect1d(buffer_kplane(1:ntpoint_total,j,i,1,1))
             spect22_tmp1(:,ii,jj) = spect22_tmp1(:,ii,jj) +  autospect1d(buffer_kplane(1:ntpoint_total,jjj,i+ii,1,1))
             numpt(ii,jj) = numpt(ii,jj) + 1
           enddo ! end ii loop
         enddo ! end jj loop
       enddo
       enddo ! end j loop

       do jj=0, jywindowr
         do ii=-ixwindowl, ixwindowr
           spect12_tmp1(:,ii,jj) = spect12_tmp1(:,ii,jj)/dble(numpt(ii,jj))
           spect11_tmp1(:,ii,jj) = spect11_tmp1(:,ii,jj)/dble(numpt(ii,jj))
           spect22_tmp1(:,ii,jj) = spect22_tmp1(:,ii,jj)/dble(numpt(ii,jj))
         enddo
       enddo

       call MPI_Allreduce(spect12_tmp1,spect12,(ntrans)*(ixwindowl+ixwindowr+1)*(jywindowr+1),&
                           MPI_DOUBLE_COMPLEX,MPI_SUM,dim2_line,ierr)
       call MPI_Allreduce(spect11_tmp1,spect11,(ntrans)*(ixwindowl+ixwindowr+1)*(jywindowr+1),&
                           MPI_DOUBLE_COMPLEX,MPI_SUM,dim2_line,ierr)
       call MPI_Allreduce(spect22_tmp1,spect22,(ntrans)*(ixwindowl+ixwindowr+1)*(jywindowr+1),&
                           MPI_DOUBLE_COMPLEX,MPI_SUM,dim2_line,ierr)

       spect12 = spect12/dble(numprocs)
       spect11 = spect11/dble(numprocs)
       spect22 = spect22/dble(numprocs)

       if(myid.eq.0) then
         fname = 'coherence_k'//fnum4_1//'_i'//fnum4_2//'-'//fnum4_3//'_nsection'//fnum4_4//'.dat'
         print *, 'writing file: ', trim(fname)
         open(unit=12,file=trim(fname),status='unknown')
           rewind(12)
           write(12,'(a)') 'variables = f,x,y,spect12_r,spect12_i,spect11,spect22'
           write(12,*) 'zone T = k001, i =', ntrans/2, ' j =', (ixwindowl + ixwindowr + 1), ' k =', jywindowr+1
           do jj=0, jywindowr
             do ii = -ixwindowl,ixwindowr
               do i = 1,ntrans/2
                 write(12,*)  dble(i-1)/fperiod, dble(ii)*dx_kp, dble(jj)*dy_kp, real(spect12(i,ii,jj)),imag(spect12(i,ii,jj)), abs(spect11(i,ii,jj)),abs(spect22(i,ii,jj))
               enddo
             enddo
           enddo
         close(12)
       endif

     enddo ! end k loop

   endif ! end icorrxt.eq.1

   call decomp_info_finalize(decomp)
   call decomp_2d_finalize
   call FinalizeHDF5()
   call MPI_FINALIZE(ierr)


  contains

     subroutine ReadFiles(kshift)
       integer, intent(in) :: kshift
       integer :: n, nn, num_file, ntpoint_tmp
       type(tp_rdwt_hdf5) :: TSkplane
       integer :: imax_file, jmax_file, kmax_file
       character(8) :: fnum8

       TSkplane%comm = MPI_COMM_WORLD
       TSkplane%info = MPI_INFO_NULL
       TSkplane%gname = '/kplane'
       TSkplane%dnum = 1
       TSkplane%rank = 4
       allocate(TSkplane%dname(TSkplane%dnum), TSkplane%dimsf(TSkplane%rank), TSkplane%dimsm(TSkplane%rank))
       allocate(TSkplane%count(TSkplane%rank), TSkplane%offset(TSkplane%rank), TSkplane%block(TSkplane%rank), TSkplane%stride(TSkplane%rank))
       TSkplane%dname(1) = 'p'
       TSkplane%IsHSInitialized = .true.
       TSkplane%count  = 1
       TSkplane%stride = 1

       imax_file = DNSIndex_k%iend - DNSIndex_k%ibe + 1
       jmax_file = DNSIndex_k%jend - DNSIndex_k%jbe + 1
       kmax_file = DNSIndex_k%nkplane
       TSkplane%dimsf = (/fileprop(1)%ntpoint,jmax_file,imax_file,kmax_file/)
       TSkplane%dimsm = (/fileprop(1)%ntpoint,jmax_rd,ilen_rd,1/)
       TSkplane%block = (/fileprop(1)%ntpoint,jmax_rd,ilen_rd,1/)
       TSkplane%offset(1) = 0
       TSkplane%offset(2) = jbe_k - 1
       TSkplane%offset(3) = ist_rd - 1
       TSkplane%offset(4) = kshift - 1

       num_file = ( fileprop(1)%file_end - fileprop(1)%file_be)/fileprop(1)%file_skip + 1
       if(myid.eq.0) print *, 'Total number of files = ', num_file
       ntpoint_tmp = 0
       if(myid.eq.0) print *, 'begin reading files  '
       do n=1, num_file
         write(unit=fnum8,fmt='(I08.8)') fileprop(1)%file_be + (n-1)*fileprop(1)%file_skip
         TSkplane%fname = trim(fileprop(1)%filepath)//'timeseries_kplane_'//fnum8//'.h5'
         if(myid.eq.0) print *, 'reading file: ', trim(TSkplane%fname)
         call ReadHDF5_4D(TSkplane,buffer_kplane(ntpoint_tmp+1:ntpoint_tmp+TSkplane%dimsf(1),1:jmax_rd,ist_rd:iend_rd,:,:) )
         ntpoint_tmp = ntpoint_tmp + TSkplane%dimsf(1)
       enddo
       if(myid.eq.0) print *, 'finish reading files'

     end subroutine ReadFiles

    subroutine Input()
      integer, parameter :: nid=11
      integer :: i, j, k, kk, n, errcode
      logical :: PlaneNameMatch = .false.

      npath = 1
      allocate(fileprop(npath))
      if(myid.eq.0) then
        read(*,*)
        read(*,*) icorrxt, iCalCoh
        read(*,*)
        do n=1, npath
          read(*,*)
          read(*,*) fileprop(n)%file_be, fileprop(n)%file_end, fileprop(n)%file_skip
          read(*,*)
          read(*,'(a)') fileprop(n)%filepath
        enddo
        read(*,*)
        read(*,*)
        read(*,*) nsection, ioverlap, iwindow
        read(*,*)
        read(*,*)
        read(*,*) ibe_k, iend_k, jbe_k, jend_k, ixwindowl, ixwindowr, jywindowr
        read(*,*)
        read(*,*) num_kplane
        allocate(kplane(num_kplane))
        read(*,*) (kplane(i), i = 1, num_kplane)
        call ReadDNS_index_kplane(fileprop(1)%ntpoint, dx_kp, dy_kp, dt_sample, DNSIndex_k,fileprop(1)%filepath)
        fileprop(1)%nskip = 0

        ! check kplane location index
        if(num_kplane.le.DNSIndex_k%nkplane) then
          do k=1, num_kplane
            do kk=1, DNSIndex_k%nkplane
              if(kplane(k).eq.DNSIndex_k%kplane(kk)) then
                PlaneNameMatch = .true.
              endif
            enddo ! end kk loop
            if(.not.PlaneNameMatch) then
              print *, '#############################################################'
              print *, 'The input kplane name does not match the name in DNSIndex.h5 '
              print *, 'Input kplane name ', kplane(k) ,' does not exist in file DNSIndex.h5 '
              print *, 'Stop ... '
              stop
            else
              PlaneNameMatch = .false.
            endif
          enddo ! end k loop
        else ! num_kplane.gt.DNSIndex_k%nkplane
          print *, 'The input num_kplane does not match the num_kplane in DNSIndex.h5. '
          print *, 'The input num_kplane = ', num_kplane, 'num_kplane in DNSIndex.h5 = ', DNSIndex_k%nkplane
          stop
        endif

        ! calculate the total number of time point
        ntpoint_total = 0
        do n=1, npath
          ntpoint_total = ntpoint_total + ( (fileprop(n)%file_end-fileprop(n)%file_be)/fileprop(n)%file_skip + 1 )*fileprop(n)%ntpoint
        enddo

        print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        print *, 'ixwindowl =', ixwindowl,  'ixwindowr =', ixwindowr, 'jywindowr =', jywindowr
        if(ioverlap.ne.0.and.ioverlap.ne.1) then
          print*,'ioverlap can ONLY be 0 or 1'
          stop
        end if
        print *, 'number of wall-normal locations =', num_kplane
        print *, 'Wall-normal plane indexes k =', kplane
        write(*,*)
        print *, 'ntpoint_total =', ntpoint_total, 'dt_sample = ', dt_sample
        print *, '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      endif ! end myid.eq.0

      call MPI_Bcast(icorrxt,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(iCalCoh,               1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if(icorrxt.eq.1.and.iCalCoh.eq.1) then
        if(myid.eq.0) then
          print *, 'icorrxt and iCalCoh cannot be 1 at the same time... Stop '
          print *, 'icorrxt = ', icorrxt, 'iCalCoh = ', iCalCoh
        endif
        stop
      endif

      call MPI_Bcast(fileprop(1)%filepath,   400, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(fileprop(1)%file_be,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(fileprop(1)%file_end,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(fileprop(1)%file_skip,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(fileprop(1)%ntpoint,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      call MPI_Bcast(nsection,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(ioverlap,      1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(iwindow,       1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(ntpoint_total, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(ixwindowl,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(ixwindowr,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(jywindowr,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(dt_sample,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      call MPI_Bcast(ibe_k,         1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(iend_k,        1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(jbe_k,         1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(jend_k,        1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(num_kplane,             1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if(myid.ne.0) allocate(kplane(num_kplane))
      call MPI_Bcast(kplane,        num_kplane, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      call MPI_Bcast(DNSIndex_k%ibe,         1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(DNSIndex_k%iend,        1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(DNSIndex_k%iskip,       1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(DNSIndex_k%jbe,         1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(DNSIndex_k%jend,        1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(DNSIndex_k%jskip,       1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(DNSIndex_k%nkplane,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if(myid.ne.0) allocate(DNSIndex_k%kplane(DNSIndex_k%nkplane))
      call MPI_Bcast(DNSIndex_k%kplane, DNSIndex_k%nkplane, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(dx_kp,                  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      call MPI_Bcast(dy_kp,                  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      if(mod((jend_k-jbe_k+1),numprocs).ne.0) then
        if(myid.eq.0) print *, 'number of processors/cores must be a fraction of jmax ... '
        errcode = 337
        call MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
      endif

    end subroutine input


end program Calcoh
