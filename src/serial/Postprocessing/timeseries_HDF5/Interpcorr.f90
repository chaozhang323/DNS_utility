Program Interpcorr
   implicit none

   type fDimInterp
      integer :: idimold
      integer :: idimnew
   end type fDimInterp

   type fprop
      character(200) :: fname
      integer :: nvars, nheaders, is_block_data, ndim, IsCalPhasespeed
      type(fDimInterp), dimension(:), allocatable :: DimInterp
   end type fprop

   integer :: nfiles
   type(fprop), dimension(:), allocatable :: fave
   integer :: ii

   call Input()
   print *, 'Number of Files =', nfiles
   do ii = 1, nfiles
     print *, ' File Number =', ii
     print *, ' File Name =', trim(fave(ii)%fname)
     call CalInterp(fave(ii))
   enddo


contains
   subroutine Input()
     implicit none
     integer, parameter :: nid=5
     integer :: n, i

     read(nid,*)
     read(nid,*) nfiles
     allocate(fave(nfiles))
     do n = 1, nfiles
       read(nid,*)
       read(nid,*)
       read(nid,'(a)') fave(n)%fname
       read(nid,*)
       read(nid,*) fave(n)%nvars, fave(n)%nheaders, fave(n)%is_block_data, fave(n)%ndim, fave(n)%IsCalPhasespeed
       if(fave(n)%ndim.ge.1) allocate(fave(n)%DimInterp(fave(n)%ndim))
       do i = 1, fave(n)%ndim
          read(nid,*)
          read(nid,*) fave(n)%DimInterp(i)%idimold, fave(n)%DimInterp(i)%idimnew
       enddo
     enddo
     if (nid.ne.5) close(nid)
   end subroutine Input

   subroutine CalInterp(fp)
     type(fprop), intent(in) :: fp
     integer :: nvars, ilength, nheaders, ndim, IsCalPhasespeed
     character(200) :: fname
     integer :: isblock
     real(8), dimension(:,:,:), allocatable :: var3d, vartmp3d, vartmp3d2
     real(8), dimension(:,:), allocatable :: varin, vartmp2d, var2d
     integer, dimension(:), allocatable :: idimold, idimnew
     real(8), dimension(:), allocatable :: indexold1, indexold2, indexnew1, indexnew2
     integer :: n, i, j, ii
     character(10000) :: buff10000

     nvars = fp%nvars
     ndim = fp%ndim
     allocate(idimold(ndim),idimnew(ndim))
     do i = 1, ndim
        idimold(i) = fp%DimInterp(i)%idimold
        idimnew(i) = fp%DimInterp(i)%idimnew
     enddo
     ilength = product(idimold)
     isblock = fp%is_block_data
     fname = fp%fname
     nheaders = fp%nheaders
     IsCalPhasespeed = fp%IsCalPhasespeed
    
     allocate(varin(nvars,ilength))
     open(111,file=trim(fname), status='old')
     open(112,file=trim(fname)//'_Interp')
     do i = 1, nheaders
        read(111,'(a)') buff10000
        if(i.eq.1) write(112,'(a)') trim(buff10000)
     enddo
     if(isblock.eq.1) then
         do i = 1, nvars
            read(111, *) varin(i,:)
         enddo
     else
         do i = 1, ilength
            read(111, *) varin(:,i)
         enddo
     endif
     close(111)
     print *, 'ilength =', ilength

     select case(ndim)
     case(1)

        allocate(vartmp2d(nvars,idimold(1)),var2d(nvars,idimnew(1)))
        allocate(indexold1(idimold(1)), indexnew1(idimnew(1)))

        vartmp2d = varin

        indexold1 = 1.d0; indexnew1 = 1.d0
        forall(i=1:idimold(1), idimold(1).gt.1)
           indexold1(i) = dble(i-1)/dble(idimold(1)-1)
        end forall
        forall(i=1:idimnew(1),idimnew(1).gt.1)
           indexnew1(i) = dble(i-1)/dble(idimnew(1)-1)
         end forall

         if(idimold(1).eq.1) then
            forall(n=1:nvars,i=1:idimnew(1))
                var2d(n,i) = vartmp2d(n,1)
            end forall
         elseif(idimnew(1).ne.idimold(1)) then
            do n = 1, nvars
               call SplineInterp(indexold1,vartmp2d(n,:),idimold(1), &
                                 indexnew1,var2d(n,:),idimnew(1),2)
            enddo
         else
            var2d = vartmp2d
         endif

         write(112,*) var2d
         deallocate(vartmp2d,var2d,indexold1,indexnew1)

      case(2)

         allocate(vartmp3d(nvars,idimold(1),idimold(2)))
         do ii = 1, ilength
            i = mod(ii-1,idimold(1))+1
            j = (ii-1)/idimold(1)+1
            vartmp3d(:,i,j) = varin(:,ii)
         enddo

         allocate(indexold1(idimold(1)), indexnew1(idimnew(1)))
         indexold1 = 1.d0; indexnew1 = 1.d0
         forall(i=1:idimold(1), idimold(1).gt.1)
            indexold1(i) = dble(i-1)/dble(idimold(1)-1)
         end forall
         forall(i=1:idimnew(1), idimnew(1).gt.1)
            indexnew1(i) = dble(i-1)/dble(idimnew(1)-1)
         end forall

         allocate(indexold2(idimold(2)), indexnew2(idimnew(2)))
         indexold2 = 1.d0; indexnew2 = 1.d0
         forall(i=1:idimold(2), idimold(2).gt.1)
            indexold2(i) = dble(i-1)/dble(idimold(2)-1)
         end forall
         forall(i=1:idimnew(2), idimnew(2).gt.1)
            indexnew2(i) = dble(i-1)/dble(idimnew(2)-1)
         end forall

         allocate(   vartmp3d2(nvars,idimnew(1),idimold(2)))
         if(idimold(1).eq.1) then
            forall(i=1:idimnew(1),j=1:idimold(2),n=1:nvars)
               vartmp3d2(n,i,j) = vartmp3d(n,1,j)
            end forall
         elseif(idimnew(1).ne.idimold(1)) then
            do n = 1, nvars
               do j = 1, idimold(2)
                  call SplineInterp(indexold1, vartmp3d(n,:,j),idimold(1), &
                                    indexnew1,vartmp3d2(n,:,j),idimnew(1), 2)
               enddo
            enddo
         else
            vartmp3d2 = vartmp3d
         endif
         deallocate(vartmp3d, indexold1, indexnew1)

         allocate(   var3d(nvars,idimnew(1),idimnew(2)))
         if(idimold(2).eq.1) then
            forall(i=1:idimnew(1),j=1:idimnew(2),n=1:nvars)
               var3d(n,i,j) = vartmp3d2(n,1,j)
            end forall
         elseif(idimold(2).ne.idimnew(2)) then
            do n = 1, nvars
               do i = 1, idimnew(1)
                  call SplineInterp(indexold2, vartmp3d2(n,i,:), idimold(2), &
                                    indexnew2,     var3d(n,i,:), idimnew(2), 2)  
               enddo
            enddo
          else
            var3d = vartmp3d2
          endif
          deallocate(vartmp3d2, indexold2, indexnew2)
         
          write(112,*) 'Zone T=Interp2D, I=', idimnew(1), 'J=', idimnew(2)
          do j = 1, idimnew(2)
             do i = 1, idimnew(1)
                write(112,*) (var3d(n,i,j), n=1,nvars)
             enddo
          enddo
          call CalCorr1d(var3d,fname)
          if(IsCalPhasespeed.eq.1) then
             call CalPhasespeed_fixdt(var3d,fname)
             call CalPhasespeed_fixdx(var3d,fname)
          endif
          deallocate(var3d)    
      case(3)
         print *, 'ndim=3 NOT implemented yet. STOP...'
         stop
      end select
      close(112)

   end subroutine CalInterp

   subroutine CalCorr1d(var2d,fn)
   ! Calculate 1d correlation
     real(8), dimension(:,:,:), intent(in) :: var2d
     character(*), intent(in) :: fn
     integer :: n, nmid1(1), nmid2(1)
     real(8) :: ds1, ds2, lint

     ds1 = var2d(1,2,1) - var2d(1,1,1)
     ds2 = var2d(2,1,2) - var2d(2,1,1)
     nmid1 = minloc(abs(var2d(1,:,1)))
     nmid2 = minloc(abs(var2d(2,1,:)))
     print *, 'nmid1 =', nmid1(1), 'nmid2 =', nmid2(1)
     print *, 'var2d =', var2d(:,nmid1(1),nmid2(1))
        
!     print *, 'maxloc1 =', maxloc(var2d(3,:,:),1), 'maxloc1 =', maxloc(var2d(3,:,:),2)

     open(21, file=trim(fn)//'_dim1.dat')
     write(21,'(a)')'variables=s,corr,lint'
     lint = 0.d0
     do n = nmid1(1), size(var2d,2)-1
        lint = lint + 0.5*(var2d(3,n,nmid2(1))+var2d(3,n+1,nmid2(1)))*ds1
        write(21,*) var2d(1,n,nmid2(1)), var2d(3,n,nmid2(1)), lint
     enddo
     close(21)

     open(22, file=trim(fn)//'_dim2.dat')
     write(22,'(a)')'variables=s,corr,lint'
     lint = 0.d0
     do n = nmid2(1), size(var2d,3)-1
        lint = lint + 0.5*(var2d(3,nmid1(1),n)+var2d(3,nmid1(1),n+1))*ds2
        write(22,*) var2d(2,nmid1(1),n), var2d(3,nmid1(1),n), lint
     enddo
     close(22)
   end subroutine CalCorr1d

   subroutine CalPhasespeed_fixdt(varxt,fn)
     ! Calculate convection velocity
     ! Ref: Bernardini & Pirozzoli, PhysFluids, 23, 085102 (2011)
     ! varxt(3,tdim, xdim)
      real(8), dimension(:,:,:), intent(in) :: varxt
      character(*), intent(in) :: fn
      integer :: n, nopt, iloc1(1), iloc2(1), nloc_lower(1), nloc_upper(1), nmid(1)
      real(8) :: t1, t2, r1, r2,  uc, tt, rr
      real(8) :: dts, dxs
      real(8) :: corr1, corr2, cor

      dts = varxt(1,2,1) - varxt(1,1,1)
      dxs = varxt(2,1,2) - varxt(2,1,1)
      nmid = minloc(abs(varxt(1,:,1)))

      nloc_lower = minloc(varxt(1,:,1),varxt(1,:,1).ge.dts) ! omega*delta/uinf ~ 100

      print *, 'Calculating convection velocity for fixed temporal separation...'
      open(21, file=trim(fn)//'_uconv_fixdt.dat')
      write(21,'(a)')'variables=dt,ucorr,xloc,cormax'
      do n = nloc_lower(1), 2*nmid(1) - 1           !!! change
         nopt = 2*nmid(1) - n
         iloc1 = maxloc(varxt(3,n,:))
         iloc2 = maxloc(varxt(3,nopt,:))
         t1 = varxt(1,n,iloc1(1));    r1 = varxt(2,n,iloc1(1));   corr1 = varxt(3,n,iloc1(1))
         t2 = varxt(1,nopt,iloc2(1)); r2 = varxt(2,nopt,iloc2(1)); corr2 = varxt(3,nopt,iloc2(1))
         tt = 0.5d0*(abs(t1)+abs(t2))
         rr = 0.5d0*(abs(r1)+abs(r2))
         cor = 0.5d0*(abs(corr1)+abs(corr2))
         write(21,'(4E20.12)') tt, rr/(tt+1.e-30), rr, cor

      enddo
      close(21)
   end subroutine CalPhasespeed_fixdt

   subroutine CalPhasespeed_fixdx(varxt,fn)
     ! Calculate convection velocity
     ! Ref: Bernardini & Pirozzoli, PhysFluids, 23, 085102 (2011)
     ! varxt(3,tdim, xdim)
      real(8), dimension(:,:,:), intent(in) :: varxt
      character(*), intent(in) :: fn
      integer :: n, nopt, iloc1(1), iloc2(1), nloc_lower(1), nloc_upper(1), nmid(1)
      real(8) :: t1, t2, r1, r2,  uc, tt, rr
      real(8) :: dts, dxs
      real(8) :: corr1, corr2, cor

      dts = varxt(1,2,1) - varxt(1,1,1)
      dxs = varxt(2,1,2) - varxt(2,1,1)
      nmid = minloc(abs(varxt(2,1,:)))

      nloc_lower = minloc(varxt(2,1,:),varxt(2,1,:).ge.dxs) ! 
      print *, 'Calculating convection velocity for fixed spatial separation...'
      open(21, file=trim(fn)//'_uconv_fixdx.dat')
      write(21,'(a)')'variables=dx,ucorr,tloc,cormax'
      do n = nloc_lower(1), size(varxt,3)
         nopt = 2*nmid(1) - n
         iloc1 = maxloc(varxt(3,:,n))
         iloc2 = maxloc(varxt(3,:,nopt))
         t1 = varxt(1,iloc1(1),n);    r1 = varxt(2,iloc1(1),n);    corr1 = varxt(3,iloc1(1),n)
         t2 = varxt(1,iloc2(1),nopt); r2 = varxt(2,iloc2(1),nopt); corr2 = varxt(3,iloc2(1),nopt)
         tt = 0.5d0*(abs(t1)+abs(t2))
         rr = 0.5d0*(abs(r1)+abs(r2))
         cor = 0.5d0*(abs(corr1)+abs(corr2))
         write(21,'(4E20.12)') rr, rr/(tt+1.e-30), tt, cor

      enddo
      close(21)
   end subroutine CalPhasespeed_fixdx


end program InterpCorr

