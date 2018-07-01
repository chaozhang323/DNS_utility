!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Module contains subroutines to read/write Plot3D unformatted files
!! All subroutines contained here do NOT necessarily need
!! to be in a module. Putting them in a module removes the
!! necessity for declaring interface when using them (otherwise
!! optional arguments do not work for some reason).
!!
!! Subroutines:
!! ReadPlot3DSolPlane: 
!!     Read Plot3D plane format function file with dimensions specified
!! ReadPlot3DSolPlaneShift: 
!!     Read Plot3D plane format function file with dimensions specified
!!     The 3rd dimension in output array is shifted to the 1st dimension
!! ReadPlot3DSolPlaneGen:
!!     Read Plot3D plane format function file without specifying the dimensions
!! WritePlot3DSolPlane: 
!!     Write Plot3D plane format function file
!! WritePlot3DSolPlaneShift: 
!!     Write Plot3D plane format function file
!!     The 3rd dimension in output array is shifted to the 1st dimension
!! ReadPlot3DGridPlane:
!!     Read Plot3D plane format grid file with dimensions specified
!! ReadPlot3DGridPlaneShift:
!!     Read Plot3D plane format grid file with dimensions specified
!!     The 3rd dimension in output array is shifted to the 1st dimension
!! ReadPlot3DGridPlaneGen:
!!     Read Plot3D plane format grid file without specifying the dimensions
!! WritePlot3DGridPlane: 
!!     Write Plot3D plane format grid file
!! WritePlot3DGridPlaneShift: 
!!     Write Plot3D plane format grid file
!!     The 3rd dimension in output array is shifted to the 1st dimension
!! ReadPlot3DSolWhole: 
!!     Read Plot3D whole format function file with dimensions specified
!! WritePlot3DSolWhole: 
!!     Write Plot3D whole format function file
!! ReadPlot3DGridWhole: 
!!     Read Plot3D whole format grid file with dimensions specified
!! WritePlot3DGridWhole: 
!!     Write Plot3D whole format grid file
!! ReadPlot2DSolution: 
!!     Read Plot2D format function file with dimensions specified
!! WritePlot2DSolution: 
!!     Write Plot2D format function file
!! WritePlot2DSolutionShift: 
!!     Write Plot2D format function file with dimensions specified
!!     The 2nd dimension of output array is shited to the 1st dimension
!! ReadPlot2DGrid:
!!     Read Plot2D format grid file with dimensions specified
!! WritePlot2DGrid:
!!     Write Plot2D format grid file
!! WritePlot2DGridshift:
!!     Write Plot2D format grid file with dimensions specified
!!     The 2nd dimension of output array is shited to the 1st dimension
!! ReadCartGrid:
!!     Read cartesian grid file used in the DNS code with dimension specified
!! WriteCartGrid:
!!     Write cartesian grid file used in the DNS code
!! WriteCartGridEx:
!!     Write cartesian grid file with grid derivatives used in the DNS code
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module MFileIO
  implicit none
  contains
     !Read Plot3D plane format function file with dimensions specified
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !kdim (input): size of the 3rd dimension
     !nvar (input): number of variables
     !vars (output): 4 dimensional array to store the data
     !               dimensions are idimxjdimxkdimxnvar
     !time (output, optional): if present, an extra double precision number
     !               will be read, which is the time parameter for DNS data files
     subroutine ReadPlot3DSolPlane(fn,idim,jdim,kdim,nvar,vars,time)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim,nvar
       real(8), intent(out), dimension(idim,jdim,kdim,nvar) :: vars
       real(8), intent(out), optional :: time
       real(8) :: rh
       integer :: id,jd,kd, nv,nb
       integer :: i, j, k, n
       open(unit=11,file=fn,form='unformatted',status='old')
       read(11)nb
       if (nb.ne.1) then
         print*,'number of blocks must be 1 in file: ', trim(fn)
         stop
       end if
       read(11)id,jd,kd,nv
       if (id.ne.idim.or.jd.ne.jdim.or.kd.ne.kdim.or.nv.ne.nvar) then
         print*,'size inconsistent with input in file: ', trim(fn)
         print*,'input idim jdim kdim nvar', idim, jdim, kdim, nvar
         print*,'file idim jdim kdim nvar', id, jd, kd, nv
         print*,'can not proceed due to smaller file size'
         stop
       end if
       do k=1,kdim
         read(11)(((vars(i,j,k,n),i=1,idim),j=1,jdim),n=1,nvar)
       end do
       if (present(time)) then
         Read(11)time
       end if
       close(11)
     end subroutine ReadPlot3DSolPlane

     !Read Plot3D plane format function file with dimensions specified
     !The 3rd dimension in output array is shifted to the 1st dimension
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !kdim (input): size of the 3rd dimension
     !nvar (input): number of variables
     !vars (output): 4 dimensional array to store the data
     !               dimensions are kdimxidimxjdimxnvar
     !time (output, optional): if present, an extra double precision number
     !               will be read, which is the time parameter for DNS data files
     !apd2 (output, optional): if present, two extra double precision numbers
     !               will be read, including time
     subroutine ReadPlot3DSolPlaneShift(fn,idim,jdim,kdim,nvar,vars,time,apd2)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim,nvar
       real(8), intent(out), dimension(kdim,idim,jdim,nvar) :: vars
       real(8), intent(out), optional :: time, apd2
       real(8) :: rh
       integer :: id,jd,kd, nv,nb
       integer :: i, j, k, n
       open(unit=11,file=fn,form='unformatted',status='old')
       read(11)nb
       if (nb.ne.1) then
         print*,'number of blocks must be 1 in file: ', trim(fn)
         stop
       end if
       read(11)id,jd,kd,nv
       if (id.ne.idim.or.jd.ne.jdim.or.kd.ne.kdim.or.nv.ne.nvar) then
         print*,'size inconsistent with input in file: ', trim(fn)
         print*,'input idim jdim kdim nvar', idim, jdim, kdim, nvar
         print*,'file idim jdim kdim nvar', id, jd, kd, nv
         print*,'can not proceed due to smaller file size'
         stop
       end if
       do k=1,kdim
         read(11)(((vars(k,i,j,n),i=1,idim),j=1,jdim),n=1,nvar)
       end do
       if (present(time)) then
         if (present(apd2)) then
            read(11)time,apd2
         else
           read(11)time
         end if
       end if
       close(11)
     end subroutine ReadPlot3DSolPlaneShift

     !Read Plot3D plane format function file without specifing dimensions
     !fn (input): file name
     !idim (output): size of the 1st dimension
     !jdim (output): size of the 2nd dimension
     !kdim (output): size of the 3rd dimension
     !nvar (output): number of variables
     !vars (output): 4 dimensional allocatable array to store the data
     !               dimensions are idimxjdimxkdimxnvar if ish=0
     !               dimensions are kdimxidimxjdimxnvar if ish=1
     !ish  (input): whether the dimesions in vars will be shifted or not
     !time (output, optional): if present, an extra double precision number
     !               will be read, which is the time parameter for DNS data files
     !apd2 (output, optional): if present, two extra double precision numbers
     !               will be read, including time
     subroutine ReadPlot3DSolPlaneGen(fn,idim,jdim,kdim,nvar,vars,ish,time,apd2)
       character(*), intent(in) :: fn
       real(8), intent(out), dimension(:,:,:,:), allocatable :: vars
       real(8), intent(out), optional :: time, apd2
       integer, intent(out) :: idim,jdim,kdim,nvar
       integer, intent(in), optional :: ish
       integer :: nblks, ishift
       integer :: i, j, k, n
       if (present(ish)) then
         ishift = ish
       else
         ishift = 1
       end if
       open(unit=11,file=fn,form='unformatted',status='old')
       read(11)nblks
       if (nblks.ne.1) then
         print*,'number of blocks must be 1 in file: ', trim(fn)
         stop
       end if
       read(11)idim,jdim,kdim,nvar
       print*,'idim,jdim,kdim,nvar:',idim,jdim,kdim,nvar
       if (allocated(vars)) deallocate(vars)
       if (ishift.eq.1) then
         allocate(vars(kdim,idim,jdim,nvar))
       else
         allocate(vars(idim,jdim,kdim,nvar))
       end if
       do k=1,kdim
         if (ishift.eq.1) then
           read(11)(((vars(k,i,j,n),i=1,idim),j=1,jdim),n=1,nvar)
         else
           read(11)(((vars(i,j,k,n),i=1,idim),j=1,jdim),n=1,nvar)
         end if
       end do
       if (present(time)) then
         if (present(apd2)) then
            read(11)time,apd2
         else
           read(11)time
         end if
       end if
       close(11)
     end subroutine ReadPlot3DSolPlaneGen

     !Write Plot3D plane format function file
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !kdim (input): size of the 3rd dimension
     !nvar (input): number of variables
     !vars (input): 4 dimensional array to store the data
     !               dimensions are idimxjdimxkdimxnvar
     !time (input, optional): if present, an extra double precision number
     !               will be written, which is the time parameter for DNS data files
     !apd2 (input, optional): if present, two extra double precision numbers
     !               will be written, including time
     subroutine WritePlot3DSolPlane(fn,idim,jdim,kdim,nvar,vars,time,apd2)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim,nvar
       real(8), intent(in), dimension(idim,jdim,kdim,nvar) :: vars
       real(8), intent(in), optional :: time, apd2
       integer :: i, j, k, n
       open(unit=11,file=fn,form='unformatted',status='unknown')
       write(11)1 !1 block
       write(11)idim,jdim,kdim,nvar
       do k=1,kdim
         write(11)(((vars(i,j,k,n),i=1,idim),j=1,jdim),n=1,nvar)
       end do
       if (present(time)) then
         if (present(apd2)) then
           write(11)time,apd2
         else
           write(11)time
         end if
       end if
       close(11)
     end subroutine WritePlot3DSolPlane

     !Write Plot3D plane format function file
     !The 3rd dimension of vars is shifted to the 1st dimension
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !kdim (input): size of the 3rd dimension
     !nvar (input): number of variables
     !vars (input): 4 dimensional array to store the data
     !               dimensions are kdimxidimxjdimxnvar
     !time (input, optional): if present, an extra double precision number
     !               will be written, which is the time parameter for DNS data files
     !apd2 (input, optional): if present, two extra double precision numbers
     !               will be written, including time
     subroutine WritePlot3DSolPlaneShift(fn,idim,jdim,kdim,nvar,vars,time,apd2)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim,nvar
       real(8), intent(in), dimension(kdim,idim,jdim,nvar) :: vars
       real(8), intent(in), optional :: time,apd2
       integer :: i, j, k, n
       open(unit=11,file=fn,form='unformatted',status='unknown')
       write(11)1 !1 block
       write(11)idim,jdim,kdim,nvar
       do k=1,kdim
         write(11)(((vars(k,i,j,n),i=1,idim),j=1,jdim),n=1,nvar)
       end do
       if (present(time)) then
         if (present(apd2)) then
           write(11)time, apd2
         else
           write(11)time
         end if
       end if
       close(11)
     end subroutine WritePlot3DSolPlaneShift

     !Read Plot3D plane format grid file with dimensions specified
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !kdim (input): size of the 3rd dimension
     !x,y,z(output): 3 dimensional arrays to store the grid
     !               dimensions are idimxjdimxkdim
     subroutine ReadPlot3DGridPlane(fn,idim,jdim,kdim,x,y,z)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim
       real(8), intent(out), dimension(idim,jdim,kdim) :: x, y, z
       integer :: id,jd,kd,nb
       integer :: i, j, k
       open(unit=11,file=fn,form='unformatted',status='old')
       read(11)nb
       if (nb.ne.1) then
         print*,'number of blocks must be 1 in file: ', trim(fn)
         stop
       end if
       read(11)id,jd,kd
       if (id.ne.idim.or.jd.ne.jdim.or.kd.ne.kdim) then
         print*,'size inconsistent with input in file: ', trim(fn)
         print*,'input idim jdim kdim', idim, jdim, kdim
         print*,'file idim jdim kdim', id, jd, kd
         print*,'can not proceed due to smaller file size'
         stop
       end if
         do k=1,kdim
           read(11)((x(i,j,k),i=1,idim),j=1,jdim)&
                 , ((y(i,j,k),i=1,idim),j=1,jdim)&
                 , ((z(i,j,k),i=1,idim),j=1,jdim)
         end do
       close(11)
     end subroutine ReadPlot3DGridPlane

     !Read Plot3D plane format grid file with dimensions specified
     !The 3rd dimension of x,y,z are shifted to the 1st dimension
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !kdim (input): size of the 3rd dimension
     !x,y,z(output): 3 dimensional arrays to store the grid
     !               dimensions are kdimxidimxjdim
     subroutine ReadPlot3DGridPlaneShift(fn,idim,jdim,kdim,x,y,z)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim
       real(8), intent(out), dimension(kdim,idim,jdim) :: x, y, z
       integer :: id,jd,kd,nb
       integer :: i, j, k
       open(unit=11,file=fn,form='unformatted',status='old')
       read(11)nb
       if (nb.ne.1) then
         print*,'number of blocks must be 1 in file: ', trim(fn)
         stop
       end if
       read(11)id,jd,kd
       if (id.ne.idim.or.jd.ne.jdim.or.kd.ne.kdim) then
         print*,'size inconsistent with input in file: ', trim(fn)
         print*,'input idim jdim kdim', idim, jdim, kdim
         print*,'file idim jdim kdim', id, jd, kd
         print*,'can not proceed due to smaller file size'
         stop
       end if
         do k=1,kdim
           read(11)((x(k,i,j),i=1,idim),j=1,jdim)&
                 , ((y(k,i,j),i=1,idim),j=1,jdim)&
                 , ((z(k,i,j),i=1,idim),j=1,jdim)
         end do
       close(11)
     end subroutine ReadPlot3DGridPlaneShift

     !Read Plot3D plane format grid file without specifying the dimensions
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !kdim (input): size of the 3rd dimension
     !x,y,z(output): 3 dimensional allocatable arrays to store the grid
     !               dimensions are idimxjdimxkdim if ish=0
     !               dimensions are kdimxidimxjdim if ish=1
     !ish (input): whether the 3rd dimension is shifted to the 1st dimension or not
     subroutine ReadPlot3DGridPlaneGen(fn,idim,jdim,kdim,x,y,z,ish)
       character(*), intent(in) :: fn
       real(8), intent(out), dimension(:,:,:), allocatable :: x, y, z
       integer, intent(in), optional :: ish
       integer, intent(out) :: idim, jdim, kdim
       integer :: nblks, ishift
       integer :: i, j, k
       if (present(ish)) then
         ishift = ish
       else
         ishift = 1
       end if
       open(unit=11,file=fn,form='unformatted',status='old')
       read(11)nblks
       if (nblks.ne.1) then
         print*,'number of blocks must be 1 in file: ', trim(fn)
         stop
       end if
       read(11)idim,jdim,kdim
       print*,'idim,jdim,kdim:',idim,jdim,kdim
       if (allocated(x)) deallocate(x)
       if (allocated(y)) deallocate(y)
       if (allocated(z)) deallocate(z)
       if (ishift.eq.1) then
         allocate(x(kdim,idim,jdim),y(kdim,idim,jdim),z(kdim,idim,jdim))
       else
         allocate(x(idim,jdim,kdim),y(idim,jdim,kdim),z(idim,jdim,kdim))
       end if
       
       do k=1,kdim
         if (ishift.eq.1) then
           read(11)((x(k,i,j),i=1,idim),j=1,jdim)&
                 , ((y(k,i,j),i=1,idim),j=1,jdim)&
                 , ((z(k,i,j),i=1,idim),j=1,jdim)
         else
           read(11)((x(i,j,k),i=1,idim),j=1,jdim)&
                 , ((y(i,j,k),i=1,idim),j=1,jdim)&
                 , ((z(i,j,k),i=1,idim),j=1,jdim)
         end if
       end do
       close(11)
     end subroutine ReadPlot3DGridPlaneGen

     !Write Plot3D plane format grid file
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !kdim (input): size of the 3rd dimension
     !x,y,z(input): 3 dimensional arrays containing grid
     !               dimensions are idimxjdimxkdim
     subroutine WritePlot3DGridPlane(fn,idim,jdim,kdim,x,y,z)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim
       real(8), intent(in), dimension(idim,jdim,kdim) :: x, y, z
       integer :: i, j, k
       open(unit=11,file=fn,form='unformatted',status='unknown')
       write(11)1 !1 block
       write(11)idim,jdim,kdim
       do k=1,kdim
         write(11)((x(i,j,k),i=1,idim),j=1,jdim)&
                , ((y(i,j,k),i=1,idim),j=1,jdim)&
                , ((z(i,j,k),i=1,idim),j=1,jdim)
       end do
       close(11)
     end subroutine WritePlot3DGridPlane

     !Write Plot3D plane format grid file
     !The 3rd dimension of x,y,z are shifted to the 1st dimension
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !kdim (input): size of the 3rd dimension
     !x,y,z(input): 3 dimensional arrays containing grid
     !               dimensions are kdimxidimxjdim
     subroutine WritePlot3DGridPlaneShift(fn,idim,jdim,kdim,x,y,z)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim
       real(8), intent(in), dimension(kdim,idim,jdim) :: x, y, z
       integer :: i, j, k
       open(unit=11,file=fn,form='unformatted',status='unknown')
       write(11)1 !1 block
       write(11)idim,jdim,kdim
       do k=1,kdim
         write(11)((x(k,i,j),i=1,idim),j=1,jdim)&
                , ((y(k,i,j),i=1,idim),j=1,jdim)&
                , ((z(k,i,j),i=1,idim),j=1,jdim)
       end do
       close(11)
     end subroutine WritePlot3DGridPlaneShift

     !Read Plot3D whole format function file
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !kdim (input): size of the 3rd dimension
     !nvar (input): number of variables
     !vars (output): 4 dimensional arrays to store the data
     !               dimensions are idimxjdimxkdimxnvar
     !time (output, optional): if present, an extra double precision number
     !               will be read, which is the time parameter for DNS data files
     subroutine ReadPlot3DSolWhole(fn,idim,jdim,kdim,nvar,vars,time)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim,nvar
       real(8), intent(out), dimension(idim,jdim,kdim,nvar) :: vars
       real(8), intent(out), optional :: time
       real(8) :: rh
       integer :: id,jd,kd, nv,nb
       integer :: i, j, k, n
       open(unit=11,file=fn,form='unformatted',status='old')
       read(11)nb
       if (nb.ne.1) then
         print*,'number of blocks must be 1 in file: ', trim(fn)
         stop
       end if
       read(11)id,jd,kd,nv
       if (id.ne.idim.or.jd.ne.jdim.or.kd.ne.kdim.or.nv.ne.nvar) then
         print*,'size inconsistent with input in file: ', trim(fn)
         print*,'input idim jdim kdim nvar', idim, jdim, kdim, nvar
         print*,'file idim jdim kdim nvar', id, jd, kd, nv
         print*,'can not proceed due to smaller file size'
         stop
       end if
       read(11)((((vars(i,j,k,n),i=1,idim),j=1,jdim),k=1,kdim),n=1,nvar)
       if (present(time)) then
         Read(11)time
       end if
       close(11)
     end subroutine ReadPlot3DSolWhole

     !Write Plot3D whole format function file
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !kdim (input): size of the 3rd dimension
     !nvar (input): number of variables
     !vars (input): 4 dimensional arrays to store the data
     !               dimensions are idimxjdimxkdimxnvar
     !time (input, optional): if present, an extra double precision number
     !               will be written, which is the time parameter for DNS data files
     !apd2 (input, optional): if present, two extra double precision numbers
     !               will be written, including time
     subroutine WritePlot3DSolWhole(fn,idim,jdim,kdim,nvar,vars,time,apd2)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim,nvar
       real(8), intent(in), dimension(idim,jdim,kdim,nvar) :: vars
       real(8), intent(in), optional :: time, apd2
       integer :: i, j, k, n
       open(unit=11,file=fn,form='unformatted',status='unknown')
       write(11)1 !1 block
       write(11)idim,jdim,kdim,nvar
       write(11)((((vars(i,j,k,n),i=1,idim),j=1,jdim),k=1,kdim),n=1,nvar)
       if (present(time)) then
         if (present(apd2)) then
           write(11)time, apd2
         else
           write(11)time
         end if
       end if
       close(11)
     end subroutine WritePlot3DSolWhole

     !Read Plot3D whole format grid file
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !kdim (input): size of the 3rd dimension
     !x,y,z(input): 3 dimensional arrays containing grid
     !               dimensions are idimxjdimxkdim
     subroutine ReadPlot3DGridWhole(fn,idim,jdim,kdim,x,y,z)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim
       real(8), intent(out), dimension(idim,jdim,kdim) :: x, y, z
       real(8) :: rh
       integer :: id,jd,kd,nb
       integer :: i, j, k
       open(unit=11,file=fn,form='unformatted',status='old')
       read(11)nb
       if (nb.ne.1) then
         print*,'number of blocks must be 1 in file: ', trim(fn)
         stop
       end if
       read(11)id,jd,kd
       if (id.ne.idim.or.jd.ne.jdim.or.kd.ne.kdim) then
         print*,'size inconsistent with input in file: ', trim(fn)
         print*,'input idim jdim kdim', idim, jdim, kdim
         print*,'file idim jdim kdim', id, jd, kd
         print*,'can not proceed due to smaller file size'
         stop
       end if
       read(11)(((x(i,j,k),i=1,idim),j=1,jdim),k=1,kdim)&
             , (((y(i,j,k),i=1,idim),j=1,jdim),k=1,kdim)&
             , (((z(i,j,k),i=1,idim),j=1,jdim),k=1,kdim)
       close(11)
     end subroutine ReadPlot3DGridWhole

     !Write Plot3D whole format grid file
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !kdim (input): size of the 3rd dimension
     !x,y,z(input): 3 dimensional arrays containing grid
     !               dimensions are idimxjdimxkdim
     subroutine WritePlot3DGridWhole(fn,idim,jdim,kdim,x,y,z)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim
       real(8), intent(in), dimension(idim,jdim,kdim) :: x, y, z
       integer :: i, j, k
       open(unit=11,file=fn,form='unformatted',status='unknown')
       write(11)1 !1 block
       write(11)idim,jdim,kdim
       write(11)(((x(i,j,k),i=1,idim),j=1,jdim),k=1,kdim)&
              , (((y(i,j,k),i=1,idim),j=1,jdim),k=1,kdim)&
              , (((z(i,j,k),i=1,idim),j=1,jdim),k=1,kdim)
       close(11)
     end subroutine WritePlot3DGridWhole

     !Read Plot2D format function file
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !nvar (input): number of variables
     !vars (output): 3 dimensional arrays to store the data
     !               dimensions are idimxjdimxnvar
     !apd1 (output, optional): if present, an extra double precision number
     !               will be read
     !apd2 (output, optional): if present, two extra double precision numbers
     !               will be read, including apd1
     subroutine ReadPlot2DSolution(fn,idim,jdim,nvar,vars, apd1, apd2)
       implicit none
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, nvar
       real(8), intent(out), optional :: apd1, apd2
       real(8), intent(out), dimension(idim,jdim,nvar) :: vars
       integer :: i, j, n, nblks, id,jd,nd
       open(unit=11,file=fn,form='unformatted',status='unknown')
       read(11)nblks
       if (nblks.ne.1) then
         print*,'number of blocks must be 1 in file: ', trim(fn)
         stop
       end if
       read(11)id,jd,nd
       if (id.ne.idim.or.jd.ne.jdim.or.nd.ne.nvar) then
         print*,'size inconsistent with input in file: ', trim(fn)
         print*,'input idim jdim nvar', idim, jdim, nvar
         print*,'file idim jdim nvar', id, jd, nd
         print*,'can not proceed due to smaller file size'
         stop
       end if
       read(11)(((vars(i,j,n),i=1,idim),j=1,jdim),n=1,nvar)
       if (present(apd1)) then
         if (present(apd2)) then
           read(11)apd1, apd2
         else
           read(11)apd1
         end if
       end if
       close(11)
     end subroutine ReadPlot2DSolution

     !Write Plot2D format function file
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !nvar (input): number of variables
     !vars (output): 3 dimensional arrays to store the data
     !               dimensions are idimxjdimxnvar
     !apd1 (output, optional): if present, an extra double precision number
     !               will be written
     !apd2 (output, optional): if present, two extra double precision numbers
     !               will be written, including apd1
     subroutine WritePlot2DSolution(fn,idim,jdim,nvar,vars, apd1, apd2)
       implicit none
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, nvar
       real(8), intent(in), optional :: apd1, apd2
       real(8), intent(in), dimension(idim,jdim,nvar) :: vars
       integer :: i, j, n
       open(unit=11,file=fn,form='unformatted',status='unknown')
       write(11)1 !1 block
       write(11)idim,jdim,nvar
       write(11)(((vars(i,j,n),i=1,idim),j=1,jdim),n=1,nvar)
       if (present(apd1)) then
         if (present(apd2)) then
           write(11)apd1, apd2
         else
           write(11)apd1
         end if
       end if
       close(11)
     end subroutine WritePlot2DSolution

     !Write Plot2D format function file
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !nvar (input): number of variables
     !vars (output): 3 dimensional arrays to store the data
     !               dimensions are jdim x idim x nvar
     !apd1 (output, optional): if present, an extra double precision number
     !               will be written
     !apd2 (output, optional): if present, two extra double precision numbers
     !               will be written, including apd1
     subroutine WritePlot2DSolutionShift(fn,idim,jdim,nvar,vars, apd1, apd2)
       implicit none
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, nvar
       real(8), intent(in), optional :: apd1, apd2
       real(8), intent(in), dimension(jdim,idim,nvar) :: vars
       integer :: i, j, n
       open(unit=11,file=fn,form='unformatted',status='unknown')
       write(11)1 !1 block
       write(11)idim,jdim,nvar
       write(11)(((vars(j,i,n),i=1,idim),j=1,jdim),n=1,nvar)
       if (present(apd1)) then
         if (present(apd2)) then
           write(11)apd1, apd2
         else
           write(11)apd1
         end if
       end if
       close(11)
     end subroutine WritePlot2DSolutionShift

     !Read Plot2D grid file
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !x, y (output): 2 dimensional arrays to store the grid
     !              dimensions are idimxjdim
     subroutine ReadPlot2DGrid(fn,idim,jdim,x,y)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim
       real(8), intent(out), dimension(idim,jdim) :: x, y
       integer :: i, j,id,jd,nblks
       open(unit=11,file=fn,form='unformatted',status='old')
       read(11)nblks
       if (nblks.ne.1) then
         print*,'number of blocks must be 1 in file:', trim(fn)
         stop
       end if
       read(11)id,jd
       if (id.ne.idim.or.jd.ne.jdim) then
         print*,'dimension inconsistent in file: ', trim(fn)
         print*,'in file idim, jdim: ', id, jd
         print*,'input idim, jdim: ', idim, jdim
       end if
       read(11)((x(i,j),i=1,idim),j=1,jdim)&
             , ((y(i,j),i=1,idim),j=1,jdim)
       close(11)
     end subroutine ReadPlot2DGrid

     !Write Plot2D grid file
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !x, y (input): 2 dimensional arrays containing the grid
     !              dimensions are idimxjdim
     subroutine WritePlot2DGrid(fn,idim,jdim,x,y)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim
       real(8), intent(in), dimension(idim,jdim) :: x, y
       integer :: i, j
       open(unit=11,file=fn,form='unformatted',status='unknown')
       write(11)1 !1 block
       write(11)idim,jdim
       write(11)((x(i,j),i=1,idim),j=1,jdim)&
              , ((y(i,j),i=1,idim),j=1,jdim)
       close(11)
     end subroutine WritePlot2DGrid

     !Write Plot2D grid file
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !x, y (input): 2 dimensional arrays containing the grid
     !              dimensions are jdim x idim
     subroutine WritePlot2DGridShift(fn,idim,jdim,x,y)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim
       real(8), intent(in), dimension(jdim,idim) :: x, y
       integer :: i, j
       open(unit=11,file=fn,form='unformatted',status='unknown')
       write(11)1 !1 block
       write(11)idim,jdim
       write(11)((x(j,i),i=1,idim),j=1,jdim)&
              , ((y(j,i),i=1,idim),j=1,jdim)
       close(11)
     end subroutine WritePlot2DGridShift

     !Read cartesian grid file used in the DNS code
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !kdim (input): size of the 3rd dimension
     !x,y,z(output): 1 dimensional arrays to store the grid
     subroutine ReadCartGrid(fn,idim,jdim,kdim,x,y,z)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim
       real(8), intent(out) :: x(idim), y(jdim), z(kdim)
       integer :: it, id, jd, kd
       open(21,file=fn,form='unformatted',status='old')
       read(21)it
       if (it.ne.1) then
         print*,'type is not 1 in cartesian grid file: ', trim(fn)
         stop
       end if
       read(21)id,jd,kd
       if (id.ne.idim.or.jd.ne.jdim.or.kd.ne.kdim) then
         print*,'size inconsistent with input in file: ', trim(fn)
         print*,'input idim jdim kdim', idim, jdim, kdim
         print*,'file idim jdim kdim', id, jd, kd
         stop
       end if
       read(21)x,y,z
       close(21)
     end subroutine ReadCartGrid

     !Write cartesian grid file used in the DNS code
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !kdim (input): size of the 3rd dimension
     !x,y,z(input): 1 dimensional arrays containing the grid
     subroutine WriteCartGrid(fn,idim,jdim,kdim,x,y,z)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim
       real(8), intent(in) :: x(idim), y(jdim), z(kdim)
       open(unit=11,file=fn,form='unformatted',status='unknown')
       write(11)1 !type1
       write(11)idim,jdim,kdim
       write(11)x,y,z
       close(11)
     end subroutine WriteCartGrid

     !Write cartesian grid file with grid derivatives used in the DNS code
     !fn (input): file name
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !kdim (input): size of the 3rd dimension
     !x,y,z(input): 1 dimensional arrays containing the grid
     !dzdk (input): 1 dimensional arrays containing the 1st order grid derivative dz/dk
     !dzdk2(input): 1 dimensional arrays containing the 2nd order grid eerivative d2z/dk2
     subroutine WriteCartGridEx(fn,idim,jdim,kdim,x,y,z,dzdk,dzdk2)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim
       real(8), intent(in) :: x(idim), y(jdim), z(kdim)
       real(8), intent(in) :: dzdk(kdim), dzdk2(kdim)
       open(unit=11,file=fn,form='unformatted',status='unknown')
       write(11)2 !type2
       write(11)idim,jdim,kdim
       write(11)x,y,z
       write(11)dzdk,dzdk2
       close(11)
     end subroutine WriteCartGridEx

     !Read old_plane.xxx DNS grid
     !fn (input): file path
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !kdim (input): size of the 3rd dimension
     !x,y,z(output): 1 dimensional arrays containing the grid
     subroutine ReadOldPlaneGrid(fn,idim,jdim,kdim,x,y,z)
     character(*), intent(in) :: fn
     integer, intent(in) :: idim, jdim, kdim
     real(8), intent(out) :: x(idim), y(jdim), z(kdim)
     character(300) :: fname
     integer :: i, j, k, n
     character(3) :: fnext
     real(8) :: tempreal
     integer :: tempint
     integer :: id, jd, kd

     write(unit=fnext,fmt='(i3)') 101
     fname=trim(fn)//'REST/old_plane.'//fnext
     open(11,file=fname,form='unformatted',status='old')
     rewind(11) 
     read(11) tempreal, tempint, id, jd, kd
     read(11)
     if (id.ne.idim.or.jd.ne.jdim.or.kd.ne.kdim) then
         print*,'size inconsistent with input in file: ', trim(fn)
         print*,'input idim jdim kdim', idim, jdim, kdim
         print*,'file idim jdim kdim', id, jd, kd
         stop
     end if
     read(11) x, y, z
     close(11)
     end subroutine ReadOldPlaneGrid

     !Read old_plane.xxx DNS flow fields
     !fn (input): file path
     !idim (input): size of the 1st dimension
     !jdim (input): size of the 2nd dimension
     !kdim (input): size of the 3rd dimension
     !nvars (input): number of flow variables (5+ns)
     !time(output): time of the flow volume
     !nstart(output): iteration indexes of the flow volume
     !vars (output): 3d flow field (u,v,w,p,t,rhon)
     subroutine ReadOldPlaneSol(fn,idim,jdim,kdim,nvars,time,nstart,vars)
     character(*), intent(in) :: fn
     integer, intent(in) :: idim, jdim, kdim, nvars
     real(8), intent(out) :: time
     integer, intent(out) :: nstart
     real(8), intent(out) :: vars(idim,jdim,kdim,nvars)
     character(3) :: fnext
     character(300) :: fname
     integer :: i, j, k, n
     real(8) :: temp
     integer :: id, jd, kd

     do k = 1, kdim
        write(unit=fnext,fmt='(i3)') 100+k
        fname=trim(fn)//'REST/old_plane.'//fnext
        open(11,file=fname,form='unformatted',status='old')
        read(11) time, nstart, id, jd, kd
        read(11)
        read(11) 
        if (id.ne.idim.or.jd.ne.jdim.or.kd.ne.kdim) then
           print*,'size inconsistent with input in file: ', trim(fn)
           print*,'input idim jdim kdim', idim, jdim, kdim
           print*,'file idim jdim kdim', id, jd, kd
            stop
        end if
        read(11) (((vars(i,j,k,n),i=1,idim),j=1,jdim),n=1,nvars)
        close(11)
     enddo
     end subroutine ReadOldPlaneSol
end module MFileIO
