module MAverage
implicit none

contains
   pure function ave1D(idim,var1D,dxdi)
      integer, intent(in) :: idim
      real(8), intent(in) :: var1D(idim)
      real(8), intent(in), optional :: dxdi(idim)
      real(8) :: ave1D
      if(present(dxdi)) then
         ave1D = sum(var1D*dxdi)/sum(dxdi)
      else
         ave1D = sum(var1D)/dble(idim)
      endif
   end function ave1D

   pure function ave2D(idim,jdim,var2D,dxdi,dydj)
      integer, intent(in) :: idim, jdim
      real(8), intent(in) :: var2D(idim,jdim)
      real(8), intent(in), optional :: dxdi(idim)
      real(8), intent(in), optional :: dydj(jdim)
      real(8) :: ave2D
      real(8) :: avetemp(jdim)
      integer :: j

      if(present(dxdi)) then
         forall(j=1:jdim)
            avetemp(j) = ave1D(idim,var2D(:,j),dxdi)
         end forall
      else
         forall(j=1:jdim)
            avetemp(j) = ave1D(idim,var2D(:,j))
         end forall
      endif

      if(present(dydj)) then
         ave2D = ave1D(jdim,avetemp,dydj)
      else
         ave2D = ave1D(jdim,avetemp)
      endif
   end function ave2D

   pure function gausswin(ihw,sigma_inp)
      integer, intent(in) :: ihw
      real(8), intent(in), optional :: sigma_inp
      real(8) :: gausswin(-ihw:ihw)
      real(8) :: sigma
      integer :: ctr

      if(present(sigma_inp)) then
          sigma = sigma_inp
      else
          sigma = dble(ihw)/3.d0
      endif
      gausswin = 0.d0; gausswin(0) = 1.d0
      forall(ctr = -ihw:ihw,ihw.gt.0.and.sigma.gt.1.d-16)
         gausswin(ctr) = exp(-0.5d0* (dble(ctr)/sigma)**2.d0 )
      end forall
      gausswin = gausswin/sum(gausswin) ! Normalize gausswin
      end function gausswin

      subroutine compute_average_time(ihw,var,varave)
      implicit none
      integer, intent(in) ::  ihw
      real(8), dimension(:,-ihw:), intent(in) :: var
      real(8), dimension(:), intent(out) :: varave
      integer :: ctr, k
      real(8) :: gausswin(-ihw:ihw)
      real(8) :: sigma

      if(ihw.eq.0) then
          varave = var(:,0)
      else
          sigma = dble(ihw)/3.d0
          forall(ctr = -ihw:ihw)
             gausswin(ctr) = exp(-0.5d0* (dble(ctr)/sigma)**2.d0 )
          end forall
          gausswin = gausswin/sum(gausswin) ! Normalize gausswin
          forall( k=1:size(var,dim=1) )
             varave(k) = sum(gausswin*var(k,:))
          end forall
      endif
   end subroutine compute_average_time

   subroutine compute_average_streamspan(dxdi,dydj,var,varave)
      implicit none
      real(8), dimension(:), intent(in) :: dxdi, dydj
      real(8), dimension(:,:,:), intent(in) :: var
      real(8), dimension(:), intent(out) :: varave
      integer :: imax, jmax, kmax
      real(8), dimension(:,:,:), allocatable :: varm
      real(8), dimension(:,:), allocatable :: area
      real(8) :: c12i
      integer :: i, j

      c12i = 1.d0/12.d0
      imax = size(var,dim=1)
      jmax = size(var,dim=2)
      kmax = size(var,dim=3)
      allocate( varm(imax,jmax,kmax), area(imax,jmax) )
      forall(i=1:imax,j=1:jmax)
         varm(i,j,:) = var(i,j,:)*dxdi(i)*dydj(j)
         area(i,j) = dxdi(i)*dydj(j)
      end forall
      varave = sum(sum(varm,dim=1),dim=1)/sum(sum(area,dim=1))
      deallocate(varm,area)
   end subroutine compute_average_streamspan

   subroutine compute_average_streamspan_plane(dxdi,dydj,var,varave)
      implicit none
      real(8), dimension(:,:), intent(in) :: var
      real(8), intent(in) :: dxdi(size(var,1)), dydj(size(var,2))
      real(8), intent(out) :: varave
      real(8) :: varm(size(var,1),size(var,2)), area(size(var,1),size(var,2))
      integer :: i, j

      forall(i=1:size(var,1),j=1:size(var,2))
         varm(i,j) = var(i,j)*dxdi(i)*dydj(j)
         area(i,j) = dxdi(i)*dydj(j)
      end forall
      varave = sum(sum(varm,dim=1),dim=1)/sum(sum(area,dim=1))
   end subroutine compute_average_streamspan_plane

end module MAverage

