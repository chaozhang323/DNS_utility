

    module modTecbin
      use MTecplotIO
      implicit none

      integer, private :: TotZone
      type(Pointer3D), dimension(:), allocatable, private :: vars
      real(8), dimension(:,:,:,:), allocatable, target, private :: vartmp
      integer, dimension(:), allocatable, private :: varidx
      integer, private :: imax, jmax, kmax, numshare
      integer, private :: nvars
      integer, dimension(:), allocatable :: ShareVar
      real(8) :: time
      character(600) :: varoutname
      logical :: IsTecInitialized = .false.

   contains

      ! TZ: total zone in *.plt file
      ! nvarout: # of variables
      ! varname: variable name
      ! nshare: # of share variables
      subroutine InitTec(TZ,idim,jdim,kdim,nvarout,varname,nshare)
        integer, intent(in) :: TZ, idim, jdim, kdim, nvarout, nshare
        character(*), intent(in) :: varname
        integer :: i, n

        if(allocated(varidx)) deallocate(varidx)
        if(allocated(vartmp)) deallocate(vartmp)
        if(allocated(vars)) deallocate(vars)
        if(allocated(ShareVar)) deallocate(ShareVar)

        TotZone = TZ
        imax = idim
        jmax = jdim
        kmax = kdim
        nvars = nvarout
        numshare = nshare
        varoutname = trim(varname)

        allocate(varidx(nvars))
        do i=1, nvars
          varidx(i) = i
        enddo
        allocate(vartmp(imax,jmax,kmax,nvars))
        allocate(vars(nvars))
        allocate(ShareVar(nvars))
        ShareVar = 0
        do n=1, nvars
          vars(n)%pt => vartmp(1:imax,1:jmax,1:kmax,varidx(n))
          if(n.le.numshare) ShareVar(n) = 1
        enddo

        IsTecInitialized = .true.

      end subroutine InitTec

      ! write 3 variables with dim imax, jmax, kmax
      subroutine WriteTec(fn,buffer1,buffer2,buffer3,CurrentZone)
        character(*), intent(in) :: fn
        integer, intent(in) :: CurrentZone
        real(8), intent(in) :: buffer1(:,:,:),buffer2(:,:,:),buffer3(:,:,:)
        logical :: IsFileNew, IsFileClose
        integer :: i,j,k

        if(.not.IsTecInitialized) then
          print *, 'Tec should be Initialized'
          stop
        endif

        do k=1, kmax
          do j=1, jmax
            do i=1, imax
              vartmp(i,j,k,1) = buffer1(i,j,k)
              vartmp(i,j,k,2) = buffer2(i,j,k)
              vartmp(i,j,k,3) = buffer3(i,j,k)
            enddo
          enddo
        enddo

        time = time + 1
        if(CurrentZone.eq.1) then
          IsFileNew=.true.; IsFileClose=.false.
          call WriteTecBin(trim(fn),trim(varoutname),imax,jmax,kmax,nvars,time,vars,IsFileNew,IsFileClose)
        elseif(CurrentZone.lt.Totzone) then
          IsFileNew=.false.; IsFileClose=.false.
          call WriteTecBin(trim(fn),trim(varoutname),imax,jmax,kmax,nvars-numshare,time,vars(numshare+1:nvars),IsFileNew,IsFileClose,ShareVar)
        else
          IsFileNew=.false.; IsFileClose=.true.
          call WriteTecBin(trim(fn),trim(varoutname),imax,jmax,kmax,nvars-numshare,time,vars(numshare+1:nvars),IsFileNew,IsFileClose,ShareVar)
        endif

      end subroutine WriteTec

      ! write 1d grid with 1d AveAcousticsInt
      subroutine WriteTec1dGrid(fn,buffer1,buffer2)
        character(*), intent(in) :: fn
        real(8), intent(in) :: buffer1(:)
        real(8), intent(in) :: buffer2(:,:)
        logical :: IsFileNew, IsFileClose
        integer :: i,j,k,n

        if(.not.IsTecInitialized) then
          print *, 'Tec should be Initialized'
          stop
        endif

        do k=1, kmax
          do j=1, jmax
            do i=1, imax
              vartmp(i,j,k,1) = buffer1(i)
              vartmp(i,j,k,2:nvars) = buffer2(i,1:(nvars-1))
            enddo
          enddo
        enddo
        time = time + 1
        IsFileNew=.true.; IsFileClose=.true.
        call WriteTecBin(trim(fn),trim(varoutname),imax,jmax,kmax,nvars,time,vars,IsFileNew,IsFileClose)

      end subroutine WriteTec1dGrid

      ! write 2d grid with 2d AveAcousticsStat
      subroutine WriteTec2dGrid(fn,buffer1,buffer2,buffer3,buffer4,buffer5)
        character(*), intent(in) :: fn
        real(8), intent(in) :: buffer1(:,:,:)
        real(8), intent(in) :: buffer2(:,:,:)
        real(8), intent(in) :: buffer3(:,:,:), buffer4(:,:,:), buffer5(:,:,:)
        logical :: IsFileNew, IsFileClose
        integer :: i,j,k,n

        if(.not.IsTecInitialized) then
          print *, 'Tec should be Initialized'
          stop
        endif
!       print *, 'nvars = ', nvars
!       print *, 'size(buffer3, dim=1) = ', size(buffer3, dim=1)
!       print *, 'size(buffer3, dim=2) = ', size(buffer3, dim=2)
!       print *, 'size(buffer3, dim=3) = ', size(buffer3, dim=3)

        do k=1, kmax
          do j=1, jmax
            do i=1, imax
              vartmp(i,j,k,1) = buffer1(k,i,1)
              vartmp(i,j,k,2) = buffer1(k,i,2)
              vartmp(i,j,k,3:(nvars-3)) = buffer2(k,i,1:(nvars-5))
              vartmp(i,j,k,nvars-2) = buffer3(k,i,1)
              vartmp(i,j,k,nvars-1) = buffer4(k,i,1)
              vartmp(i,j,k,nvars)   = buffer5(k,i,1)
            enddo
          enddo
        enddo

      time = time + 1
        IsFileNew=.true.; IsFileClose=.true.
        call WriteTecBin(trim(fn),trim(varoutname),imax,jmax,kmax,nvars,time,vars,IsFileNew,IsFileClose)

      end subroutine WriteTec2dGrid



      ! write 3d flowdata
      subroutine WriteTec3D(fn,buffer)
        character(*), intent(in) :: fn
        real(8), intent(in) :: buffer(:,:,:,:)
        logical :: IsFileNew, IsFileClose
        integer :: i,j,k,n

        if(.not.IsTecInitialized) then
          print *, 'Tec should be Initialized'
          stop
        endif
!       print *, 'nvars = ', nvars
!       print *, 'size(buffer3, dim=1) = ', size(buffer3, dim=1)
!       print *, 'size(buffer3, dim=2) = ', size(buffer3, dim=2)
!       print *, 'size(buffer3, dim=3) = ', size(buffer3, dim=3)

        do k=1, kmax
          do j=1, jmax
            do i=1, imax
              vartmp(i,j,k,1:nvars) = buffer(i,j,k,1:nvars)
            enddo
          enddo
        enddo

      time = time + 1
        IsFileNew=.true.; IsFileClose=.true.
        call WriteTecBin(trim(fn),trim(varoutname),imax,jmax,kmax,nvars,time,vars,IsFileNew,IsFileClose)

      end subroutine WriteTec3D

      subroutine WriteTec12(fn,buffer1,buffer2,buffer3,buffer4,buffer5,buffer6,buffer7,buffer8,buffer9,buffer10, &
                               buffer11,buffer12)
        character(*), intent(in) :: fn
        real(8), intent(in) :: buffer1(:,:,:),buffer2(:,:,:),buffer3(:,:,:),buffer4(:,:,:),buffer5(:,:,:),buffer6(:,:,:)
        real(8), intent(in) :: buffer7(:,:,:),buffer8(:,:,:),buffer9(:,:,:),buffer10(:,:,:),buffer11(:,:,:),buffer12(:,:,:)
        logical :: IsFileNew, IsFileClose
        integer :: i,j,k


        if(.not.IsTecInitialized) then
          print *, 'Tec should be Initialized'
          stop
        endif

        do k=1, kmax
          do j=1, jmax
            do i=1, imax
              vartmp(i,j,k,1) = buffer1(i,j,k)
              vartmp(i,j,k,2) = buffer2(i,j,k)
              vartmp(i,j,k,3) = buffer3(i,j,k)
              vartmp(i,j,k,4) = buffer4(i,j,k)
              vartmp(i,j,k,5) = buffer5(i,j,k)
              vartmp(i,j,k,6) = buffer6(i,j,k)
              vartmp(i,j,k,7) = buffer7(i,j,k)
              vartmp(i,j,k,8) = buffer8(i,j,k)
              vartmp(i,j,k,9) = buffer9(i,j,k)
              vartmp(i,j,k,10) = buffer10(i,j,k)
              vartmp(i,j,k,11) = buffer11(i,j,k)
              vartmp(i,j,k,12) = buffer12(i,j,k)
            enddo
          enddo
        enddo

        time = time + 1
        IsFileNew=.true.; IsFileClose=.true.
        call WriteTecBin(trim(fn),trim(varoutname),imax,jmax,kmax,nvars,time,vars,IsFileNew,IsFileClose)
 

      end subroutine WriteTec12



    end module modTecbin
