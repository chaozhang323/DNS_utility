

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
      character(100) :: varoutname

   contains

      ! TZ: total zone in *.plt file
      ! nvarout: # of variables
      ! varname: variable name
      ! nshare: # of share variables
      subroutine InitTec(TZ,idim,jdim,kdim,nvarout,varname,nshare)
        integer, intent(in) :: TZ, idim, jdim, kdim, nvarout, nshare
        character(*), intent(in) :: varname
        integer :: i, n

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

      end subroutine InitTec

      subroutine WriteTec(fn,buffer1,buffer2,buffer3,CurrentZone)
        character(*), intent(in) :: fn
        integer, intent(in) :: CurrentZone
        real(8), intent(in) :: buffer1(:,:,:),buffer2(:,:,:),buffer3(:,:,:)
        logical :: IsFileNew, IsFileClose
        integer :: i,j,k

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













    end module modTecbin
