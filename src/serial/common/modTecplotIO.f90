!module contains subroutins to read/write Tecplot ASCII format files
!types:
!    Pointer1D, Pointer2D, Pointer3D
!subroutines:
!  ReadAsciiPoint1D:
!      Read 1 dimensional point-format Tecplot file
!  ReadAsciiPoint2D:
!      Read 2 dimensional point-format Tecplot file
!  WriteAsciiPoint1D:
!      Write 1 dimensional point-format Tecplot file
!  WriteAsciiPoint2D:
!      Write 2 dimensional point-format Tecplot file
!  WriteAsciiPoint3D:
!      Write 3 dimensional point-format Tecplot file
!  ReadAsciiBlock2D:
!      Read 2 dimensional block-format Tecplot file
!  WriteAsciiBlock1D:
!      Write 1 dimensional block-format Tecplot file
!  WriteAsciiBlock2D:
!      Write 2 dimensional block-format Tecplot file
!  WriteAsciiBlock3D:
!      Write 3 dimensional block-format Tecplot file
module MTecplotIO
  implicit none
  type Pointer1D
    real, pointer :: pt(:)
  end type

  type Pointer2D
    real, pointer :: pt(:,:)
  end type

  type Pointer3D
    real, pointer :: pt(:,:,:)
  end type

  interface
      INTEGER(4) FUNCTION tecini112 &
       (Title, &
        Variables, &
        FName, &
        ScratchDir, &
        FileType, &
        Debug, &
        VIsDouble)
        !MS$ATTRIBUTES STDCALL :: tecini112
        !MS$ATTRIBUTES REFERENCE :: Title,Variables,FName
        !MS$ATTRIBUTES REFERENCE :: ScratchDir,FileType,Debug,VIsDouble
        CHARACTER(LEN=*) Title
        CHARACTER(LEN=*) Variables
        CHARACTER(LEN=*) FName
        CHARACTER(LEN=*) ScratchDir
        INTEGER(4)       FileType
        INTEGER(4)       Debug
        INTEGER(4)       VIsDouble
      END FUNCTION tecini112
   
      INTEGER(4) FUNCTION teczne112 &
       (ZoneTitle, &
        ZoneType, &
        IMxOrNumPts, &
        JMxOrNumElements, &
        KMxOrNumFaces, &
        ICellMax, &
        JCellMax, &
        KCellMax, &
        SolutionTime, &
        StrandID, &
        ParentZone, &
        IsBlock, &
        NumFaceConnections, &
        FaceNeighborMode, &
        TotalNumFaceNodes, &
        NumConnectedBoundaryFaces, &
        TotalNumBoundaryConnections, &
        PassiveVarList, &
        ValueLocation, &
        ShareVarFromZone, &
        ShareConnectivityFromZone)
        !MS$ATTRIBUTES STDCALL :: teczne112
        !MS$ATTRIBUTES REFERENCE :: ZoneTitle,ZoneType,IMxOrNumPts
        !MS$ATTRIBUTES REFERENCE :: JMxOrNumElements,KMxOrNumFaces
        !MS$ATTRIBUTES REFERENCE :: ICellMax,JCellMax,KCellMax
        !MS$ATTRIBUTES REFERENCE :: SolutionTime,StrandID,ParentZone
        !MS$ATTRIBUTES REFERENCE :: IsBlock,PassiveVarList
        !MS$ATTRIBUTES REFERENCE :: NumFaceConnections,FaceNeighborMode
        !MS$ATTRIBUTES REFERENCE :: TotalNumFaceNodes
        !MS$ATTRIBUTES REFERENCE :: NumConnectedBoundaryFaces
        !MS$ATTRIBUTES REFERENCE :: TotalNumBoundaryConnections
        !MS$ATTRIBUTES REFERENCE :: ValueLocation,ShareVarFromZone
        !MS$ATTRIBUTES REFERENCE :: ShareConnectivityFromZone
        CHARACTER(LEN=*) ZoneTitle
        INTEGER(4)       ZoneType
        INTEGER(4)       IMxOrNumPts
        INTEGER(4)       JMxOrNumElements
        INTEGER(4)       KMx
        INTEGER(4)       ICellMax
        INTEGER(4)       JCellMax
        INTEGER(4)       KCellMax
        REAL(8)          SolutionTime
        INTEGER(4)       StrandID
        INTEGER(4)       ParentZone
        INTEGER(4)       IsBlock
        INTEGER(4)       NumFaceConnections
        INTEGER(4)       FaceNeighborMode
        INTEGER(4)       TotalNumFaceNodes
        INTEGER(4)       NumConnectedBoundaryFaces
        INTEGER(4)       TotalNumBoundaryConnections
        INTEGER(4)       PassiveVarList(*)
        INTEGER(4)       ValueLocation(*)
        INTEGER(4)       ShareVarFromZone(*)
        INTEGER(4)       ShareConnectivityFromZone
      END FUNCTION teczne112

      INTEGER(4) FUNCTION tecdat112 &
       (N, &
        FieldData, &
        IsDouble)
        !MS$ATTRIBUTES STDCALL :: tecdat112
        !MS$ATTRIBUTES REFERENCE :: N,FieldData,IsDouble
        INTEGER(4)  N
        REAL(8)     FieldData(*)
        INTEGER(4)  IsDouble
      END FUNCTION tecdat112

      INTEGER(4) FUNCTION tecend112()
        !MS$ATTRIBUTES STDCALL :: tecend112
      END FUNCTION tecend112

  end interface

  contains
    subroutine ReadAsciiPoint1D(fn,idim,nvar,vars,nheaderline)
      character(*), intent(in) :: fn
      integer, intent(in) :: nvar,idim
      real, intent(out), dimension(idim,nvar) :: vars
      integer, intent(in), optional :: nheaderline
      integer :: n, i
      open(11, file=fn, access='sequential', status='old')
      if (present(nheaderline)) then
        do n=1, nheaderline
          read(11,*)
        end do
      end if
      do i=1,idim
        read(11,*)(vars(i,n),n=1,nvar)
      end do
      close(11)
    end subroutine ReadAsciiPoint1D

    subroutine ReadAsciiPoint2D(fn,idim,jdim,nvar,vars,nheaderline)
      character(*), intent(in) :: fn
      integer, intent(in) :: nvar,idim,jdim
      type(Pointer2D), dimension(nvar), intent(out) :: vars
      integer, intent(in), optional :: nheaderline
      integer :: n, i, j
      open(11, file=fn, access='sequential', status='old')
      if (present(nheaderline)) then
        do n=1, nheaderline
          read(11,*)
        end do
      end if
      do j=1,jdim
        do i=1,idim
          read(11,*)(vars(n)%pt(i,j),n=1,nvar)
        end do
      end do
      close(11)
    end subroutine ReadAsciiPoint2D

!    subroutine ReadAsciiBlock2D(fn, idim, jdim, nvar, vars, nheaderline)
!      character(*), intent(in) :: fn
!      integer, intent(in) :: idim, jdim, nvar
!      type(Pointer2D), dimension(nvar), intent(out) :: vars
!      integer, intent(in), optional :: nheaderline
!      integer :: i, j, n
!      print *,'fn=',fn, 'nheaderline=',nheaderline
!      open(20, file=fn, access='sequential', status='old')
!      if (present(nheaderline)) then
!        do n=1, nheaderline
!          read(20,*)
!        end do
!      end if
!      do n = 1, nvar
!        read(20,*)((vars(n)%pt(i,j), i=1,idim), j=1,jdim)
!      end do
!!      close(20)
!      return
!    end subroutine ReadAsciiBlock2D
    subroutine ReadAsciiBlock2D(fn, idim, jdim, nvar, vars, nheaderline)
      character(*), intent(in) :: fn
      integer, intent(in) :: idim, jdim, nvar
      real(8), dimension(idim,jdim,nvar), intent(out) :: vars
      integer, intent(in), optional :: nheaderline
      integer :: i, j, n
      open(20, file=fn, access='sequential', status='old')
      if (present(nheaderline)) then
        do n=1, nheaderline
          read(20,*)
        end do
      end if
      do n = 1, nvar
        read(20,*)((vars(i,j,n), i=1,idim), j=1,jdim)
      end do
      close(20)
      return
    end subroutine ReadAsciiBlock2D

    ! output procedures
    subroutine WriteAsciiBlock1D(fn, varname, zonename, idim, nvar, vars, AddZone)
      character(len=*), intent(in) :: fn, varname, zonename
      integer, intent(in) :: nvar, idim
      type(Pointer1D), dimension(nvar), intent(in) :: vars
      logical, intent(in), optional :: AddZone
      integer :: i, n
      if (.not.present(AddZone)) then
        open(20,file=fn, access='sequential', status='unknown')
        write(20,1551)trim(varname)
        write(20,1552)trim(zonename), idim, 1, 1
      else
        open(20,file=fn, access='sequential', position='append')
        write(20,1552)trim(zonename), idim, 1, 1
      end if
      do n = 1, nvar
        write(20,*)(vars(n)%pt(i), i=1, idim)
      end do
      close(20)
1551  format('Variables=',A)
1552  format('ZONE T=',A ,',I=', I4, ',J=', I4, ',K=', I4, ',F=BLOCK')
    end subroutine WriteAsciiBlock1D

    subroutine WriteAsciiBlock2D(fn, varname, zonename, idim, jdim, nvar, vars, AddZone)
      character(len=*), intent(in) :: fn, varname, zonename
      integer, intent(in) :: nvar, idim, jdim
      type(Pointer2D), dimension(nvar), intent(in) :: vars
      logical, intent(in), optional :: AddZone
      integer :: i, j, n
      if (.not.present(AddZone)) then
        open(20,file=fn, access='sequential', status='unknown')
        write(20,1551)trim(varname)
        write(20,1552)trim(zonename), idim, jdim, 1
      else
        open(20,file=fn, access='sequential', position='append')
        write(20,1552)trim(zonename), idim, jdim, 1
      end if

      do n = 1, nvar
        write(20,*)((vars(n)%pt(i,j), i=1, idim), j=1, jdim)
      end do
      close(20)
1551  format('Variables=',A)
1552  format('ZONE T=',A ,',I=', I4, ',J=', I4, ',K=', I4, ',F=BLOCK')
    end subroutine WriteAsciiBlock2D

    subroutine WriteAsciiBlock3D(fn, varname, zonename, idim, jdim, kdim, nvar, vars, AddZone)
      character(len=*), intent(in) :: fn, varname, zonename
      integer, intent(in) :: nvar, idim, jdim, kdim
      type(Pointer3D), dimension(nvar), intent(in) :: vars
      logical, intent(in), optional :: AddZone
      integer :: i, j, k, n
      if (.not.present(AddZone)) then
        open(20,file=fn, access='sequential', status='unknown')
        write(20,1551)trim(varname)
        write(20,1552)trim(zonename), idim, jdim, kdim
      else
        open(20,file=fn, access='sequential', position='append')
        write(20,1552)trim(zonename), idim, jdim, kdim
      end if
      do n = 1, nvar
        write(20,*)(((vars(n)%pt(i,j,k), i=1, idim), j=1, jdim), k=1,kdim)
      end do
      close(20)
1551  format('Variables=',A)
1552  format('ZONE T=',A ,',I=', I4, ',J=', I4, ',K=', I4, ',F=BLOCK')
    end subroutine WriteAsciiBlock3D

    subroutine WriteAsciiPoint1D(fn, varname, zonename, idim, nvar, vars, AddZone)
      character(len=*), intent(in) :: fn, varname, zonename
      integer, intent(in) :: nvar, idim
      type(Pointer1D), dimension(nvar), intent(in) :: vars
      logical, intent(in), optional :: AddZone
      integer :: i, n
      if (.not.present(AddZone)) then
        open(20,file=fn, access='sequential', status='unknown')
        write(20,1561)trim(varname)
        write(20,1562)trim(zonename), idim
      else
        open(20,file=fn, access='sequential', position='append')
        write(20,1562)trim(zonename), idim
      end if

      do i = 1, idim
        write(20,1563)(vars(n)%pt(i), n=1, nVar)
      end do
      close(20)
1561  format('Variables=',A)
1562  format('ZONE T=',A ,',I=', I4)
1563  format(1X, 20E23.15)
!1563  format(1X, <nvar>E23.15)
    end subroutine WriteAsciiPoint1D
  
    subroutine WriteAsciiPoint2D(fn, varname, zonename, idim, jdim, nvar, vars, AddZone)
      character(len=*), intent(in) :: fn, varname, zonename
      integer, intent(in) :: nvar, idim, jdim
      type(Pointer2D), dimension(nvar), intent(in) :: vars
      logical, intent(in), optional :: AddZone
      integer :: i, j, n
      if (.not.present(AddZone)) then
        open(20,file=fn, access='sequential', status='unknown')
        write(20,1561)trim(varname)
        write(20,1562)trim(zonename), idim, jdim
      else
        open(20,file=fn, access='sequential', position='append')
        write(20,1562)trim(zonename), idim, jdim
      end if

      do j = 1, jdim
        do i = 1, idim
          write(20,1563)(vars(n)%pt(i,j), n=1, nvar)
        end do
      end do
      close(20)
1561  format('Variables=',A)
1562  format('ZONE T=',A ,',I=',I4,',J=',I4)
1563  format(1X, 20E23.15)
!1563  format(1X, <nvar>E23.15)
    end subroutine WriteAsciiPoint2D

    subroutine WriteAsciiPoint3D(fn, varname, zonename, idim, jdim, kdim, nvar, vars, AddZone)
      character(len=*), intent(in) :: fn, varname, zonename
      integer, intent(in) :: nvar, idim, jdim, kdim
      type(Pointer3D), dimension(nvar), intent(in) :: vars
      logical, intent(in), optional :: AddZone
      integer :: i, j, k, n
      if (.not.present(AddZone)) then
        open(20,file=fn, access='sequential', status='unknown')
        write(20,1561)trim(varname)
        write(20,1562)trim(zonename), idim, jdim, kdim
      else
        open(20,file=fn, access='sequential', position='append')
        write(20,1562)trim(zonename), idim, jdim, kdim
      end if
      do k = 1, kdim
        do j = 1, jdim
          do i = 1, idim
            write(20,1563)(vars(n)%pt(i,j,k), n=1, nvar)
          end do
        end do
      end do
      close(20)
1561  format('Variables=',A)
1562  format('ZONE T=',A ,',I=',I4,',J=',I4,',K=',I4)
1563  format(1X, 20E23.15)
    end subroutine WriteAsciiPoint3D

    subroutine WriteTecBin(fn, varname, idim, jdim, kdim, nvar, SolTime, vars, IsNewFile, IsFileClose,Sharevar)
      character(len=*), intent(in) :: fn, varname
      integer, intent(in) :: nvar, idim, jdim, kdim
      real(8), intent(in) :: SolTime
      type(Pointer3D), dimension(nvar), intent(in) :: vars
      logical, intent(in) :: IsNewFile, IsFileClose
      integer :: i_INI, i_ZNE, i_DAT, i_END, num_tot, n
      character(1) :: NullChr
      pointer (NullPtr, Null)
      integer :: Null(*)     
      integer, optional, dimension(nvar), intent(in) :: Sharevar

      i_INI = 0; i_ZNE = 0; i_DAT = 0; i_END = 0
      NullChr = char(0)
      NullPtr = 0
      if (IsNewFile) then
          i_INI = TECINI112(Title='SIMPLE DATASET'//NullChr, &
                           Variables=trim(varname)//NullChr, &
                           FName=trim(fn)//NullChr, &
                           ScratchDir='.'//NullChr, &
                           FileType=0, &
                           Debug=0, &
                           VISDouble=1)
      endif

      if(i_INI.eq.-1) then
         print *, 'ERROR: TECINI112. STOP'
         stop
      else
         if (.not.present(Sharevar)) then
            i_ZNE = TECZNE112(ZoneTitle='Simple Zone'//NullChr, &
                           ZoneType =0, &
                           IMxOrNumPts=idim, &
                           JMxOrNumElements=jdim, &
                           KMxOrNumFaces=kdim, &
                           ICellMax=0, &
                           JCellMax=0, &
                           KCellMax=0, &
                           SolutionTime=SolTime, &
                           StrandID=0, &
                           ParentZone=0, &
                           IsBlock=1,&
                           NumFaceConnections=0, &
                           FaceNeighborMode=0, &
                           TotalNumFaceNodes=0, &
                           NumConnectedBoundaryFaces=0, &
                           TotalNumBoundaryConnections=0, &
                           PassiveVarList=Null, &
                           ValueLocation=Null, &
                           ShareVarFromZone=Null, &
                           ShareConnectivityFromZone=0)
         else
            i_ZNE = TECZNE112(ZoneTitle='Simple Zone'//NullChr, &
                           ZoneType =0, &
                           IMxOrNumPts=idim, &
                           JMxOrNumElements=jdim, &
                           KMxOrNumFaces=kdim, &
                           ICellMax=0, &
                           JCellMax=0, &
                           KCellMax=0, &
                           SolutionTime=SolTime, &
                           StrandID=0, &
                           ParentZone=0, &
                           IsBlock=1,&
                           NumFaceConnections=0, &
                           FaceNeighborMode=0, &
                           TotalNumFaceNodes=0, &
                           NumConnectedBoundaryFaces=0, &
                           TotalNumBoundaryConnections=0, &
                           PassiveVarList=Null, &
                           ValueLocation=Null, &
                           ShareVarFromZone=Sharevar, &
                           ShareConnectivityFromZone=0)
        endif
      endif

      num_tot = idim*jdim*kdim
      if(i_ZNE.eq.-1) then
         print *, 'ERROR: TECZNE112. STOP'
         stop
      else
         do n = 1, nvar
            
            i_DAT = TECDAT112(N=num_tot, &
                              FieldData=vars(n)%pt, &
                              IsDouble=1)
            if(i_DAT.eq.-1) then
               print *, 'ERROR: TECDAT112. STOP'
               stop
            endif
         enddo
      endif

      if(IsFileClose) then
         i_END = TECEND112() 
         if(i_END.eq.-1) print *, 'ERROR: TECEND112. STOP'
      endif
    end subroutine WriteTecBin
end module MTecplotIO
