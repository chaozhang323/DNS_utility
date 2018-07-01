!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! All subroutines contained here do NOT necessarily need
!! to be in a module. Putting them in a module removes the
!! necessity for declaring interface when using them (otherwise
!! optional argument may not work for some reason).
!!
!! This file works for ifort version 11
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module MFileIO
  implicit none
  contains
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
         if (id.ge.idim.and.jd.ge.jdim.and.kd.ge.kdim.and.nv.ge.nvar) then
           print*,'!!!!!!!!!!!WARNING WARNING WARNING!!!!!!!!!!!!!'
           print*,'less number of points/variables will be read in'
         else
           print*,'can not proceed due to smaller file size'
           stop
         end if
         do k=1,kdim
           read(11) ( ( ( ( (vars(i,j,k,n),i=1,idim)&
                    , (rh,i=idim+1,id) ), j=1,jdim )&
                    , ((rh, i=1,id),j=jdim+1,jd)), n=1,nvar)&
                    , (((rh,i=1,id),j=1,jd),n=nvar+1,nv)
         end do
       else
         do k=1,kdim
           read(11)(((vars(i,j,k,n),i=1,idim),j=1,jdim),n=1,nvar)
         end do
       end if
       if (present(time)) then
         Read(11)time
       end if
       close(11)
     end subroutine ReadPlot3DSolPlane

     subroutine ReadPlot3DSolPlaneShift(fn,idim,jdim,kdim,nvar,vars,time)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim,nvar
       real(8), intent(out), dimension(kdim,idim,jdim,nvar) :: vars
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
         if (id.ge.idim.and.jd.ge.jdim.and.kd.ge.kdim.and.nv.ge.nvar) then
           print*,'!!!!!!!!!!!WARNING WARNING WARNING!!!!!!!!!!!!!'
           print*,'less number of points/variables will be read in'
         else
           print*,'can not proceed due to smaller file size'
           stop
         end if
         do k=1,kdim
           read(11)( ( ( ( (vars(k,i,j,n),i=1,idim)&
                    , (rh,i=idim+1,id) ), j=1,jdim )&
                    , ((rh, i=1,id),j=jdim+1,jd)), n=1,nvar)&
                    , (((rh,i=1,id),j=1,jd),n=nvar+1,nv)
         end do
       else
         do k=1,kdim
           read(11)(((vars(k,i,j,n),i=1,idim),j=1,jdim),n=1,nvar)
         end do
       end if
       if (present(time)) then
         Read(11)time
       end if
       close(11)
     end subroutine ReadPlot3DSolPlaneShift

     subroutine WritePlot3DSolPlane(fn,idim,jdim,kdim,nvar,vars,time)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim,nvar
       real(8), intent(in), dimension(idim,jdim,kdim,nvar) :: vars
       real(8), intent(in), optional :: time
       integer :: i, j, k, n
       open(unit=11,file=fn,form='unformatted',status='unknown')
       write(11)1 !1 block
       write(11)idim,jdim,kdim,nvar
       do k=1,kdim
         write(11)(((vars(i,j,k,n),i=1,idim),j=1,jdim),n=1,nvar)
       end do
       if (present(time)) then
         write(11)time
       end if
       close(11)
     end subroutine WritePlot3DSolPlane

     subroutine WritePlot3DSolPlaneShift(fn,idim,jdim,kdim,nvar,vars,time)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim,nvar
       real(8), intent(in), dimension(kdim,idim,jdim,nvar) :: vars
       real(8), intent(in), optional :: time
       integer :: i, j, k, n
       open(unit=11,file=fn,form='unformatted',status='unknown')
       write(11)1 !1 block
       write(11)idim,jdim,kdim,nvar
       do k=1,kdim
         write(11)(((vars(k,i,j,n),i=1,idim),j=1,jdim),n=1,nvar)
       end do
       if (present(time)) then
         write(11)time
       end if
       close(11)
     end subroutine WritePlot3DSolPlaneShift

     subroutine ReadPlot3DGridPlane(fn,idim,jdim,kdim,x,y,z)
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
         if (id.ge.idim.and.jd.ge.jdim.and.kd.ge.kdim) then
           print*,'!!!!!!!!!!!WARNING WARNING WARNING!!!!!!!!!!!!!'
           print*,'less number of points will be read in'
         else
           print*,'can not proceed due to smaller file size'
           stop
         end if
         do k=1,kdim
           read(11)(((x(i,j,k),i=1,idim),(rh,i=idim+1,id)),j=1,jdim)&
                  ,((rh,i=1,id),j=jdim+1,jd)&
               , (((y(i,j,k),i=1,idim),(rh,i=idim+1,id)),j=1,jdim)&
                  ,((rh,i=1,id),j=jdim+1,jd)&
               , (((z(i,j,k),i=1,idim),(rh,i=idim+1,id)),j=1,jdim)&
                  ,((rh,i=1,id),j=jdim+1,jd)
         end do
       else
         do k=1,kdim
           read(11)((x(i,j,k),i=1,idim),j=1,jdim)&
                 , ((y(i,j,k),i=1,idim),j=1,jdim)&
                 , ((z(i,j,k),i=1,idim),j=1,jdim)
         end do
       end if
       close(11)
     end subroutine ReadPlot3DGridPlane

     subroutine ReadPlot3DGridPlaneShift(fn,idim,jdim,kdim,x,y,z)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim
       real(8), intent(out), dimension(kdim,idim,jdim) :: x, y, z
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
         if (id.ge.idim.and.jd.ge.jdim.and.kd.ge.kdim) then
           print*,'!!!!!!!!!!!WARNING WARNING WARNING!!!!!!!!!!!!!'
           print*,'less number of points will be read in'
         else
           print*,'can not proceed due to smaller file size'
           stop
         end if
         do k=1,kdim
           read(11)(((x(k,i,j),i=1,idim),(rh,i=idim+1,id)),j=1,jdim)&
                  ,((rh,i=1,id),j=jdim+1,jd)&
               , (((y(k,i,j),i=1,idim),(rh,i=idim+1,id)),j=1,jdim)&
                  ,((rh,i=1,id),j=jdim+1,jd)&
               , (((z(k,i,j),i=1,idim),(rh,i=idim+1,id)),j=1,jdim)&
                  ,((rh,i=1,id),j=jdim+1,jd)
         end do
       else
         do k=1,kdim
           read(11)((x(k,i,j),i=1,idim),j=1,jdim)&
                 , ((y(k,i,j),i=1,idim),j=1,jdim)&
                 , ((z(k,i,j),i=1,idim),j=1,jdim)
         end do
       end if
       close(11)
     end subroutine ReadPlot3DGridPlaneShift

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
         if (id.ge.idim.and.jd.ge.jdim.and.kd.ge.kdim.and.nv.ge.nvar) then
           print*,'!!!!!!!!!!!WARNING WARNING WARNING!!!!!!!!!!!!!'
           print*,'less number of points/variables will be read in'
         else
           print*,'can not proceed due to smaller file size'
           stop
         end if
         read(11)( ( ( ( ( ( (vars(i,j,k,n),i=1,idim)&
                    , (rh,i=idim+1,id) ), j=1,jdim )&
                    , ((rh, i=1,id),j=jdim+1,jd)),k=1,kdim)&
                    , (((rh,i=1,id),j=1,jd),k=kdim+1,kd)), n=1,nvar)&
                    , ((((rh,i=1,id),j=1,jd),k=1,kd),n=nvar+1,nv)
       else
         read(11)((((vars(i,j,k,n),i=1,idim),j=1,jdim),k=1,kdim),n=1,nvar)
       end if
       if (present(time)) then
         Read(11)time
       end if
       close(11)
     end subroutine ReadPlot3DSolWhole

     subroutine WritePlot3DSolWhole(fn,idim,jdim,kdim,nvar,vars,time,apd2)
       character(*), intent(in) :: fn
       integer, intent(in) :: idim, jdim, kdim,nvar
       real(8), intent(in), dimension(idim,jdim,kdim,nvar) :: vars
       real(8), intent(in), optional :: time,apd2
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
         if (id.ge.idim.and.jd.ge.jdim.and.kd.ge.kdim) then
           print*,'!!!!!!!!!!!WARNING WARNING WARNING!!!!!!!!!!!!!'
           print*,'less number of points will be read in'
         else
           print*,'can not proceed due to smaller file size'
           stop
         end if
         read(11)(((((x(i,j,k),i=1,idim),(rh,i=idim+1,id)),j=1,jdim)&
                  ,((rh,i=1,id),j=jdim+1,jd)),k=1,kdim), (((rh,i=1,id),j=1,jd),k=kdim+1,kd)&
               , (((((y(i,j,k),i=1,idim),(rh,i=idim+1,id)),j=1,jdim)&
                  ,((rh,i=1,id),j=jdim+1,jd)),k=1,kdim), (((rh,i=1,id),j=1,jd),k=kdim+1,kd)&
               , (((((z(i,j,k),i=1,idim),(rh,i=idim+1,id)),j=1,jdim)&
                  ,((rh,i=1,id),j=jdim+1,jd)),k=1,kdim), (((rh,i=1,id),j=1,jd),k=kdim+1,kd)
       else
         read(11)(((x(i,j,k),i=1,idim),j=1,jdim),k=1,kdim)&
               , (((y(i,j,k),i=1,idim),j=1,jdim),k=1,kdim)&
               , (((z(i,j,k),i=1,idim),j=1,jdim),k=1,kdim)
       end if
       close(11)
     end subroutine ReadPlot3DGridWhole

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

     subroutine WritePlot2DSolution(fn,idim,jdim,nvar,vars, apd1, apd2)
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
end module MFileIO
