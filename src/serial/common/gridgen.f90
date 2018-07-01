module cm_gridgen
use omp_lib
use modMetrics
implicit none

contains
!------------------------------------------!
!	get_intersect_point:	judge if two straight lines on a plane have a intersection
!	In:	p00,p01:	two end points of line 1
!		p10,p11:	two end points of line 2
!	Out:	inter:	0 if there's no intersection, 1 if there is
!------------------------------------------!
subroutine get_intersect_point(p00, p01, p10, p11, inter, point)
real(8), dimension(2), intent(in) :: p00, p01, p10, p11
integer, intent(out) :: inter
real(8), dimension(2), intent(out) :: point
real(8), dimension(2) :: vec0, vec1
real(8) :: det, one_over_det, l(2), v,t

point = 0.d0
vec0 = p01 - p00
vec1 = p11 - p10
det = -vec0(1)*vec1(2) + vec0(2)*vec1(1)
if( det==0.d0 ) then
	inter = 0
	return
endif

one_over_det = 1.d0 / det
l = p10 - p00
t = (-vec1(2)*l(1) + vec1(1)*l(2) )*one_over_det
v = (-vec0(2)*l(1) + vec0(1)*l(2) )*one_over_det

if( t>=0.d0 .and. t<=1.d0 .and. v>=0.d0 .and. v<=1.d0) then
	inter = 1
	point = p00 + t*vec0
	return
else
	inter = 0
	return
endif

end subroutine
!-------------------------------------------------------------------!
!	uniform_distribution:	to uniformly distribute an array from 0 to 1
!	In:	nx
!	Out:	x, dx
!-------------------------------------------------------------------!
subroutine uniform_distribution(x,dx,nx)
integer, intent(in) :: nx
real(8), intent(out) :: x(nx), dx
integer n

dx = 1.d0 / real(nx-1, 8)
x(1) = 0.d0
x(nx) = 1.d0
do n = 2, nx-1
	x(n) = (n-1)*dx
enddo

end subroutine 
!-------------------------------------------------------------------!
!	init_TFI:	initialize with transfinite interpolation method
!	In:	u_in, v_in: note that the 4 boundaries are given
!		nx,ny:	dimension of the grid
!	Out:	u,v:	results
!-------------------------------------------------------------------!
subroutine init_TFI(u_in, v_in, u, v, nx, ny)
integer, intent(in) :: nx, ny
real(8), dimension(nx,ny), intent(in) :: u_in, v_in
real(8), dimension(nx,ny), intent(out) :: u, v
real(8) :: x(nx), y(ny), dx, dy
integer :: i, j

u = u_in
v = v_in
dx = 1.d0 / real(nx-1, 8)
dy = 1.d0 / real(ny-1, 8)

do i = 2, nx-1
	x(i) = dx*(i-1)
enddo
do j = 2, ny-1
	y(j) = dy*(j-1)
enddo

do j = 2, ny-1
do i = 2, nx-1
	u(i,j) = (1.d0-x(i))*u(1,j) + x(i)*u(nx,j) &
			+ (1.d0-y(j))*u(i,1) + y(j)*u(i,ny) &
			-( (1.d0-x(i))*(1.d0-y(j))*u(1,1) &
			+ x(i)*y(j)*u(nx,ny) &
			+ x(i)*(1.d0-y(j))*u(nx,1) &
			+ (1.d0-x(i))*y(j)*u(1,ny) )
	v(i,j) = (1.d0-x(i))*v(1,j) + x(i)*v(nx,j) &
			+ (1.d0-y(j))*v(i,1) + y(j)*v(i,ny) &
			-( (1.d0-x(i))*(1.d0-y(j))*v(1,1) &
			+ x(i)*y(j)*v(nx,ny) &
			+ x(i)*(1.d0-y(j))*v(nx,1) &
			+ (1.d0-x(i))*y(j)*v(1,ny) )
enddo
enddo

end subroutine
	


!-------------------------------------------------------------------!
subroutine Laplace(u, v, dx,dy, error, nx, ny)
integer, intent(in) :: nx, ny
real(8), dimension(nx,ny), intent(inout) :: u, v
real(8), intent(in) :: dx, dy, error
real(8) :: err, dudx, dudy, dvdx, dvdy, rhs, alpha, beta, gamma
real(8), dimension(nx,ny) :: u1, v1
integer :: i, j
real(8) :: w, a

w = (dx/dy)**2
a = 1.d0 / ((w+1.d0)*2.d0)

u1 = u
v1 = v
err = 1.
do 
	do j = 2, ny-1
	do i = 2, nx-1
		u(i,j) = ( (u(i-1,j)+u(i+1,j)) + w*(u(i,j-1)+u(i,j+1)) ) * a
		v(i,j) = ( (v(i-1,j)+v(i+1,j)) + w*(v(i,j-1)+v(i,j+1)) ) * a
	enddo
	enddo
	
	err = maxval( sqrt((u-u1)**2+(v-v1)**2) )
	u1 = u
	v1 = v
	!print*, err
	if( err<error ) exit
enddo

end subroutine
!-------------------------------------------------------------------!
!	inverse_laplace:	governing equation for orthogonal grid generation
!	In:	nx,ny:	dimension
!		error:	tolerance for the convergency
!		dx, dy: spacing
!	Inout:	u,v
!	Out: iter:	num of iteration 
!-------------------------------------------------------------------!
subroutine inverse_Laplace(u, v, dx,dy, error, nx, ny, maxiter, iter)
integer, intent(in) :: nx, ny, maxiter
real(8), dimension(nx,ny), intent(inout) :: u, v
real(8), intent(in) :: dx, dy, error
integer, intent(out) :: iter
real(8) :: err, dudx, dudy, dvdx, dvdy, rhs1,rhs2, alpha, beta, gamma, Jacobi_square
real(8), dimension(nx,ny) :: u1, v1
real(8), dimension(nx,ny) :: err_mat, w
integer :: i, j, n
real(8) :: lx, ly, lxx, lyy, deno, xdeno, p, q
integer :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS, omp_get_max_threads, tid

lx = 0.5d0/dx
ly = 0.5d0/dy
lxx = 1.d0/dx**2
lyy = 1.d0/dy**2

!w = 1.d0
!w(nx-10:nx,:) = 1.d9
!w = w / (sum(w)/(nx*ny))
	
iter = 0
err = 1.
do n = 1,maxiter
	iter = iter + 1
	u1 = u
	v1 = v
	err_mat = 0.d0
	!$OMP PARALLEL do private(i,j,dudx,dudy,dvdx,dvdy,alpha,beta,gamma,Jacobi_square,deno,xdeno,rhs1,rhs2,p,q)  &
    !$OMP& shared(u,v,lx,ly,lxx,lyy,err)
!	TID = OMP_GET_THREAD_NUM()
!	PRINT *, 'Hello World from thread = ', TID
	do j = ny-1, 2, -1
	do i = nx-1, 2, -1
		dudx = (u(i+1,j)-u(i-1,j)) *lx
		dudy = (u(i,j+1)-u(i,j-1)) *ly
		dvdx = (v(i+1,j)-v(i-1,j)) *lx
		dvdy = (v(i,j+1)-v(i,j-1)) *ly
		alpha = (dudy*dudy + dvdy*dvdy) *lxx
		gamma = (dudx*dudx + dvdx*dvdx) *lyy
		deno = (alpha + gamma)*2.d0
		if(deno > 0.d0) then
			xdeno = 1.d0/deno
			beta = (dudx*dudy + dvdx*dvdy) *lx*ly
			Jacobi_square = (dudx*dvdy-dudy*dvdx)*(dudx*dvdy-dudy*dvdx)
			
			!p = -10.0*sign(x(i)-1.d0, 1.d0)*exp(-abs(y(j)-1.d0))!sqrt((x(i)-1.d0)**2+(y(j)-1.d0)**2))
			!q = -10.0*sign(y(j)-1.d0, 1.d0)*exp(-abs(y(j)-1.d0))!sqrt((x(i)-1.d0)**2+(y(j)-1.d0)**2))
			
			rhs1 = alpha*(u(i+1,j)+u(i-1,j))  + gamma*(u(i,j+1)+u(i,j-1)) &
				- 2.*beta*(u(i+1,j+1)+u(i-1,j-1)-u(i+1,j-1)-u(i-1,j+1)) !&
				!+ Jacobi_square*(p*dudx+q*dudy)
			u(i,j) = rhs1 * xdeno
			rhs2 = alpha*(v(i+1,j)+v(i-1,j))  + gamma*(v(i,j+1)+v(i,j-1)) &
				- 2.*beta*(v(i+1,j+1)+v(i-1,j-1)-v(i+1,j-1)-v(i-1,j+1)) !&
				!+ Jacobi_square*(p*dvdx+q*dvdy)
			v(i,j) = rhs2 * xdeno
			!err_mat(i,j) = abs(rhs1-deno*u(i,j)) + abs(rhs2-deno*v(i,j))
		endif
	enddo
	enddo
	!$OMP END PARALLEL do
	err_mat = abs(u1 - u) + abs(v1 - v)
	
	
	err = maxval( err_mat )
	!err = maxval( sqrt((u-u1)**2+(v-v1)**2) )
	if( err<error ) exit
enddo   ! n=1,maxiter

end subroutine
!-------------------------------------------------------------------!
subroutine pde_RyskinLeal(u, v, dist_f, nx, ny, omega, error, maxiter,iforce, iter)
integer, intent(in) :: nx, ny, maxiter, iforce
real(8), dimension(nx,ny), intent(inout) :: u, v
real(8), dimension(4,nx,ny), intent(in) :: dist_f
real(8), intent(in) :: omega, error
integer, intent(out) :: iter
real(8) :: err, dudx, dudy, dvdx, dvdy, dd,ee
real(8), dimension(nx) :: aa,bb,cc, rhs1,rhs2
real(8), dimension(nx,ny) :: u1, v1, force_x, force_y
integer :: i, j, n
real(8) :: lx, ly, lxx, lyy, deno, xdeno, p, q, alpha

! convergence history
!open(33,file='convergence_inner.dat',status='unknown')
!write(33,*) 'Variables = i, residual'
! source term
alpha = 0.01
force_x = 0.0d0; force_y = 0.0d0
if(iforce.ne.0) call get_pseudo_force(u,v,nx,ny,force_x,force_y)
force_x = alpha*force_x; force_y = alpha*force_y

iter = 0
err = 1.
do n=1,maxiter
    iter = iter + 1
    u1 = u
    v1 = v
!    err_mat = 0.d0
    do j = 2, ny-1
        do i = 2, nx-1
            aa(i) = dist_f(1,i,j)
            cc(i) = dist_f(2,i,j)
            dd = 1.d0/dist_f(3,i,j)
            ee = 1.d0/dist_f(4,i,j)
            bb(i) = -(dist_f(1,i,j) + dist_f(2,i,j) + dd + ee)
            rhs1(i) = force_x(i,j)-dd*u(i,j-1)-ee*u1(i,j+1)
            rhs2(i) = force_y(i,j)-dd*v(i,j-1)-ee*v1(i,j+1)
        enddo
        bb(1) = 1.d0; cc(1) = 0.d0; rhs1(1) = u(1,j); rhs2(1) = v(1,j)
        aa(nx) = 0.d0; bb(nx) = 1.d0; rhs1(nx) = u(nx,j); rhs2(nx) = v(nx,j)
        call tridiagonal_solve(nx,aa,bb,cc,rhs1,u(:,j))
        call tridiagonal_solve(nx,aa,bb,cc,rhs2,v(:,j))
        ! SOR
        u(:,j) = omega*(u(:,j)-u1(:,j)) + u1(:,j)
        v(:,j) = omega*(v(:,j)-v1(:,j)) + v1(:,j)
    enddo
!    err_mat = abs(u1 - u) + abs(v1 - v)

    err = max(maxval( abs(u1-u) ), maxval( abs(v1-v) ))
!    write(33,*) n, err
    !err = maxval( sqrt((u-u1)**2+(v-v1)**2) )
    if( err<error) exit
enddo
!close(33)
end subroutine pde_RyskinLeal
!-------------------------------------------------------------------!
subroutine get_distortion_func(u,v,distort)
real(8),dimension(:,:),intent(in) :: u,v
real(8),dimension(:,:,:),intent(out) :: distort
integer :: i,j,nx,ny
real(8) :: err = 0.d0
nx = size(u,dim=1)
ny = size(u,dim=2)
distort = 1.0
do j=2,ny-1
do i=2,nx-1
    distort(1,i,j) = 0.25d0*sqrt( ((u(i-1,j+1)+u(i,j+1)-u(i-1,j-1)-u(i,j-1))**2+(v(i-1,j+1)+v(i,j+1)-v(i-1,j-1)-v(i,j-1))**2)/ &
                             ((u(i,j)-u(i-1,j))**2+(v(i,j)-v(i-1,j))**2+err) )  ! f_{i-1/2,j}
    distort(2,i,j) = 0.25d0*sqrt( ((u(i+1,j+1)+u(i,j+1)-u(i+1,j-1)-u(i,j-1))**2+(v(i+1,j+1)+v(i,j+1)-v(i+1,j-1)-v(i,j-1))**2)/ &
                             ((u(i+1,j)-u(i,j))**2+(v(i+1,j)-v(i,j))**2+err) )  ! f_{i+1/2,j}
    distort(3,i,j) = 4.d0*sqrt( ((u(i,j)-u(i,j-1))**2+(v(i,j)-v(i,j-1))**2)/&   ! f_{i,j-1/2}
                            ((u(i+1,j)+u(i+1,j-1)-u(i-1,j)-u(i-1,j-1))**2+(v(i+1,j)+v(i+1,j-1)-v(i-1,j)-v(i-1,j-1))**2+err) )
    distort(4,i,j) = 4.d0*sqrt( ((u(i,j+1)-u(i,j))**2+(v(i,j+1)-v(i,j))**2)/&   ! f_{i,j+1/2}
                            ((u(i+1,j)+u(i+1,j+1)-u(i-1,j)-u(i-1,j+1))**2+(v(i+1,j)+v(i+1,j+1)-v(i-1,j)-v(i-1,j+1))**2+err) )
enddo
enddo
distort = distort +err
end subroutine get_distortion_func
!-------------------------------------------------------------------!
subroutine get_pseudo_force(u,v,nx,ny,rhs1,rhs2)
integer,intent(in) :: nx,ny
real(8),dimension(nx,ny),intent(in) :: u,v
real(8),dimension(nx,ny),intent(out) :: rhs1,rhs2
real(8),dimension(nx,ny) :: hksi,heta
real(8) :: hksi_ave(ny),heta_ave(nx)
real(8) :: p0,p1,q0,q1
integer :: i,j, im1,ip1,jm1,jp1
hksi = 0.d0; heta = 0.d0
do j=1,ny-1
do i=1,nx-1
    hksi(i,j) = sqrt( (u(i+1,j)-u(i,j))**2+(v(i+1,j)-v(i,j))**2 )
    heta(i,j) = sqrt( (u(i,j+1)-u(i,j))**2+(v(i,j+1)-v(i,j))**2 )
enddo
enddo

hksi_ave = sum( hksi, dim=1 )/dble(nx-1)
heta_ave = sum( heta, dim=2 )/dble(ny-1)

do j=2,ny-1
do i=2,nx-1
    im1 = i-1
    ip1 = i+1
    jm1 = j-1
    jp1 = j+1
    p0 = hksi(im1,j)-hksi_ave(j)**2/hksi(im1,j)
    p1 = hksi(i,j)-hksi_ave(j)**2/hksi(i,j)
    q0 = heta(i,jm1)-heta_ave(i)**2/heta(i,jm1)
    q1 = heta(i,j)-heta_ave(i)**2/heta(i,j)
    rhs1(i,j) = -(p1*(u(ip1,j)-u(i,j))/hksi(i,j) + p0*(u(i,j)-u(im1,j))/hksi(im1,j) +&
                  q1*(u(i,jp1)-u(i,j))/heta(i,j) + q0*(u(i,j)-u(i,jm1))/heta(i,jm1))
    rhs2(i,j) = -(p1*(v(ip1,j)-v(i,j))/hksi(i,j) + p0*(v(i,j)-v(im1,j))/hksi(im1,j) +&
                  q1*(v(i,jp1)-v(i,j))/heta(i,j) + q0*(v(i,j)-v(i,jm1))/heta(i,jm1))
enddo
enddo

end subroutine get_pseudo_force
!-------------------------------------------------------------------!
!	get_closest_point:	find the closest point from an array
!	In:	choice:	available points to choose from
!		num_choice:	the dimension of choice
!		p:	the point
!	Out:	idx_out:	location of the closest point
!-------------------------------------------------------------------!
subroutine get_closest_point(choice, num_choice, p, idx_out)
integer, intent(in) :: num_choice
real(8), dimension(num_choice, 2), intent(in) :: choice
real(8), dimension(2), intent(in) :: p
integer, intent(out) :: idx_out
real(8), dimension(num_choice) :: s
integer :: n, idx(1:1)

!$omp parallel do private(n) shared(choice,s)
do n = 1, num_choice
	s(n) = sqrt( sum((p-choice(n,:))**2) )
enddo
!$omp end parallel do
idx = minloc(s)
idx_out = idx(1)

end subroutine
!-------------------------------------------------------------------!
!	evaluate_left_boundary:	evaluate orhogonality of left boundary
!	In:	u,v:	the grid
!		nx,ny:	dimension of the grid
!	Out:	unit_vector:	the unit vector of the boundary point
!			dev:	the deviation from what it should be
!-------------------------------------------------------------------!
subroutine evaluate_left_boundary(u, v, nx, ny, unit_vector, dev)
integer, intent(in) :: nx, ny
real(8), dimension(nx,ny), intent(in) :: u, v
real(8), dimension(2:ny-1,2), intent(out) :: unit_vector
real(8), dimension(2:ny-1), intent(out) :: dev
real(8), dimension(2) :: p0, p1
integer :: j

do j = 2, ny-1
	p0(1) = u(1,j+1) - u(1,j-1)
	p0(2) = v(1,j+1) - v(1,j-1)
	p0 = p0 / sqrt( sum(p0**2) )
	p1(1) = u(2,j) - u(1,j)
	p1(2) = v(2,j) - v(1,j)
	
	dev(j) = sum(p0*p1)
	unit_vector(j,:) = p0
enddo
	
end subroutine
!-------------------------------------------------------------------!
subroutine evaluate_right_boundary(u, v, nx, ny, unit_vector, dev)
integer, intent(in) :: nx, ny
real(8), dimension(nx,ny), intent(in) :: u, v
real(8), dimension(2:ny-1,2), intent(out) :: unit_vector
real(8), dimension(2:ny-1), intent(out) :: dev
real(8), dimension(2) :: p0, p1
integer :: j

do j = 2, ny-1
	p0(1) = u(nx,j+1) - u(nx,j-1)
	p0(2) = v(nx,j+1) - v(nx,j-1)
	p0 = p0 / sqrt( sum(p0**2) )
	p1(1) = u(nx-1,j) - u(nx,j)
	p1(2) = v(nx-1,j) - v(nx,j)
	
	dev(j) = sum(p0*p1)
	unit_vector(j,:) = p0
enddo
	
end subroutine
!-------------------------------------------------------------------!
subroutine evaluate_bottom_boundary(u, v, nx, ny, unit_vector, dev)
integer, intent(in) :: nx, ny
real(8), dimension(nx,ny), intent(in) :: u, v
real(8), dimension(2:nx-1,2), intent(out) :: unit_vector
real(8), dimension(2:nx-1), intent(out) :: dev
real(8), dimension(2) :: p0, p1
integer :: i

do i = 2, nx-1
	p0(1) = u(i+1,1) - u(i-1,1)
	p0(2) = v(i+1,1) - v(i-1,1)
	p0 = p0 / sqrt( sum(p0**2) )
	p1(1) = u(i,2) - u(i,1)
	p1(2) = v(i,2) - v(i,1)
	
	dev(i) = sum(p0*p1)
	unit_vector(i,:) = p0
enddo
	
end subroutine
!-------------------------------------------------------------------!
subroutine evaluate_top_boundary(u, v, nx, ny, unit_vector, dev)
integer, intent(in) :: nx, ny
real(8), dimension(nx,ny), intent(in) :: u, v
real(8), dimension(2:nx-1,2), intent(out) :: unit_vector
real(8), dimension(2:nx-1), intent(out) :: dev
real(8), dimension(2) :: p0, p1
integer :: i

do i = 2, nx-1
	p0(1) = u(i+1,ny) - u(i-1,ny)
	p0(2) = v(i+1,ny) - v(i-1,ny)
	p0 = p0 / sqrt( sum(p0**2) )
	p1(1) = u(i,ny-1) - u(i,ny)
	p1(2) = v(i,ny-1) - v(i,ny)
	
	dev(i) = sum(p0*p1)
	unit_vector(i,:) = p0
enddo
	
end subroutine
!-------------------------------------------------------------------!
!	choose_boundary:	choose a better boundary point based on evaluation
!	In:	choice:	available points to choose from
!		num_choice:	the dimension of choice
!		unit_vector:	the unit vector of the boundary point
!		dev:	the deviation from what it should be
!		idx_in:	original index for boundary
!		ny:	dimension for the boundary
!	Out:	idx_out:	the new index
!			u,v:	the redistributed boundary
!-------------------------------------------------------------------!
subroutine choose_boundary(choice, num_choice, unit_vector, dev, idx_in, idx_out, u, v, ny )
integer, intent(in) :: num_choice, ny
real(8), dimension(num_choice,2), intent(in) :: choice
real(8), dimension(2:ny-1,2), intent(in) :: unit_vector
real(8), dimension(2:ny-1), intent(in) :: dev
integer, dimension(ny), intent(in) :: idx_in
integer, dimension(ny), intent(out) :: idx_out
real(8), dimension(2:ny-1), intent(out) :: u,v
real(8), dimension(2) :: p
integer j, idx!, c_index

idx_out = idx_in

do j = 2, ny-1
	p = choice(idx_out(j),:) + dev(j)*unit_vector(j,:)
	call get_closest_point(choice, num_choice, p, idx)
	idx_out(j) = idx
	if( idx_out(j)<=idx_out(j-1) ) idx_out(j) = idx_out(j-1) + 1
	!if( idx_out(j)>=idx_out(j+1) ) idx_out(j) = idx_out(j+1) - 1
	u(j) = choice(idx_out(j),1)
	v(j) = choice(idx_out(j),2)
enddo

if( idx_out(ny-1)>=idx_out(ny) ) then
	idx_out(ny-1) = idx_out(ny) - 1
	u(ny-1) = choice(idx_out(ny-1),1)
	v(ny-1) = choice(idx_out(ny-1),2)
endif

end subroutine
!--------------------------------------------------------------------!
!	init_index:	initialize index
!	In:	choice:	available points to choose from
!		num_choice:	the dimension of choice
!		
!		nx:	dimension for the boundary
!		itype:	type of initilization
!	Out:	idx_out:	the new index
!			u,v:	the redistributed boundary
!--------------------------------------------------------------------!
subroutine init_index(choice, num_choice, idx_out, u,v, nx, itype)
integer, intent(in) :: nx, num_choice
real(8), dimension(num_choice,2), intent(in) :: choice
integer, dimension(nx), intent(out) :: idx_out
real(8), dimension(nx), intent(out) :: u,v
integer, intent(in) :: itype
integer i, nseg

if(num_choice==nx) then
	do i = 1, nx
		idx_out(i) = i
		u(i) = choice(i,1)
		v(i) = choice(i,2)
	enddo
	print*, 'Num_choice equals actual num_points...'
	return
endif

idx_out(1) = 1
idx_out(nx) = num_choice

selectcase(itype)
case(0)
	do i = 2, nx-1
		idx_out(i) = i
	enddo
case(1)
	nseg = num_choice / (nx-1)
	do i = 2, nx-1
		idx_out(i) = nseg + idx_out(i-1)
	enddo
endselect

do i = 1, nx
	u(i) = choice(idx_out(i),1)
	v(i) = choice(idx_out(i),2)
enddo

end subroutine
!--------------------------------------------------------------------!
!	compute_grid:	the routine assembling the whole process
!
!--------------------------------------------------------------------!
subroutine compute_grid(u_in, v_in, u, v, idx_left_in, idx_right_in, &
						idx_bottom_in, idx_top_in, idx_left , idx_right, idx_bottom, idx_top, nx, ny, &
						choice_left, num_choice_left, choice_right, num_choice_right, &
						choice_bottom, num_choice_bottom, choice_top, num_choice_top, &
						error, max_change, maxiter_global,maxiter_inner, nsave, boundary_fixed, boundary_symmetric)!, iconverge)
integer, intent(in) :: nx, ny, num_choice_left, num_choice_right, num_choice_bottom, num_choice_top, nsave, &
                       maxiter_global, maxiter_inner, max_change
real(8), dimension(nx,ny), intent(in) :: u_in, v_in
real(8), dimension(nx,ny), intent(out) :: u, v
integer, dimension(nx), intent(in) :: idx_bottom_in, idx_top_in
integer, dimension(ny), intent(in) :: idx_left_in, idx_right_in
integer, dimension(nx), intent(out) :: idx_bottom, idx_top
integer, dimension(ny), intent(out) :: idx_left, idx_right
integer, dimension(4), intent(in) :: boundary_fixed, boundary_symmetric
real(8), dimension(num_choice_left,2), intent(in) :: choice_left
real(8), dimension(num_choice_right,2), intent(in) :: choice_right
real(8), dimension(num_choice_bottom,2), intent(in) :: choice_bottom
real(8), dimension(num_choice_top,2), intent(in) :: choice_top
real(8), intent(in) :: error
!integer, intent(out) :: iconverge
real(8) :: dx, dy, unit_vector_nx(2:nx-1,2), unit_vector_ny(2:ny-1,2), dev_nx(2:nx-1), dev_ny(2:ny-1)
integer, dimension(nx) :: idx_bottom1, idx_top1
integer, dimension(ny) :: idx_left1, idx_right1
real(8) :: x(nx), y(ny)
real(8) :: lchange(4), lchange_total
integer :: ncount, iter, i
integer :: cr_t,cm_t,c1_t,c2_t
real(8) :: rate

idx_left = idx_left_in
idx_right = idx_right_in
idx_bottom = idx_bottom_in
idx_top = idx_top_in
u = u_in
v = v_in

call uniform_distribution(x,dx,nx)
call uniform_distribution(y,dy,ny)

!iconverge = 0
ncount = 0
! initialize the system clock
call system_clock(count_rate=cr_t)
call system_clock(count_max=cm_t)
rate = dble(cr_t)
call system_clock(c1_t)
do i=1,maxiter_global
	!call Laplace(u, v, dx,dy, error, nx, ny)
	call inverse_Laplace(u, v, dx,dy, error, nx, ny, maxiter_inner, iter)
	
	lchange = 0.d0
	
	if( boundary_fixed(1)==0 ) then
		call evaluate_left_boundary(u, v, nx, ny, unit_vector_ny, dev_ny)
		call choose_boundary(choice_left, num_choice_left, unit_vector_ny, dev_ny, idx_left, idx_left1, &
								u(1,2:ny-1), v(1,2:ny-1), ny)
		lchange(1) = maxval(abs(idx_left-idx_left1))
		idx_left = idx_left1
	endif
	if( boundary_fixed(2)==0 ) then
		call evaluate_right_boundary(u, v, nx, ny, unit_vector_ny, dev_ny)
		call choose_boundary(choice_right, num_choice_right, unit_vector_ny, dev_ny, idx_right, idx_right1, &
								u(nx,2:ny-1), v(nx,2:ny-1), ny)
		lchange(2) = maxval(abs(idx_right-idx_right1))
		idx_right = idx_right1
	else
		!call bend_right_boundary(u, v, nx, ny)
	endif
	
	if( boundary_fixed(3)==0 ) then
		call evaluate_bottom_boundary(u, v, nx, ny, unit_vector_nx, dev_nx)
		call choose_boundary(choice_bottom, num_choice_bottom, unit_vector_nx, dev_nx, idx_bottom, idx_bottom1, &
								u(2:nx-1,1), v(2:nx-1,1), nx)
		lchange(3) = maxval(abs(idx_bottom-idx_bottom1))
		idx_bottom = idx_bottom1
	endif
	if( boundary_fixed(4)==0 ) then
		call evaluate_top_boundary(u, v, nx, ny, unit_vector_nx, dev_nx)
		call choose_boundary(choice_top, num_choice_top, unit_vector_nx, dev_nx, idx_top, idx_top1, &
								u(2:nx-1,ny), v(2:nx-1,ny), nx)
		lchange(4) = maxval(abs(idx_top-idx_top1))
		idx_top = idx_top1
	endif
	lchange_total = maxval(lchange)
	
	if( boundary_symmetric(1)==1 ) then
		call symmetric_left(v, nx,ny)
	endif
	if( boundary_symmetric(2)==1 ) then
		call symmetric_right(v, nx,ny)
	endif
	if( boundary_symmetric(3)==1 ) then
		call symmetric_bottom(u, nx,ny)
	endif
	if( boundary_symmetric(4)==1 ) then
		call symmetric_top(u, nx,ny)
	endif
	
	
	ncount = ncount + 1
!	print*, 'Counting...No.', ncount,',with change of', lchange_total, 'iteration = ', iter
	
	if(mod(ncount,nsave)==0) then
        print*, 'Counting...No.', ncount,',with change of', lchange_total, 'iteration = ', iter
!        exit
	endif
	if( lchange_total<=max_change ) then
!		iconverge = 1
		print*,'Converged! '
		exit
	endif
	
enddo
call system_clock(c2_t)
print*,'compute_grid time [sec]: ',dble(c2_t-c1_t)/rate

end subroutine
!--------------------------------------------------------------------!
! Orthogonal mesh generator with distortion function (weak constraint method):
!
! Luis Eca, 2D Orthogonal Grid Generation with Boundary Point Distribution Control, JCP, 1996
!--------------------------------------------------------------------!
subroutine compute_grid_new(u, v, idx_left , idx_right, idx_bottom, idx_top, nx, ny, &
                        choice_left, num_choice_left, choice_right, num_choice_right, &
                        choice_bottom, num_choice_bottom, choice_top, num_choice_top, &
                        boundary_fixed, boundary_symmetric, error_inner, error_global, maxiter_inner, maxiter_global,max_change,nsave)
integer, intent(in) :: nx, ny, num_choice_left, num_choice_right, num_choice_bottom, num_choice_top,max_change,nsave
real(8), dimension(nx,ny), intent(inout) :: u, v
integer, dimension(nx), intent(inout) :: idx_bottom, idx_top
integer, dimension(ny), intent(inout) :: idx_left, idx_right
integer, dimension(4), intent(in) :: boundary_fixed, boundary_symmetric
real(8), dimension(num_choice_left,2), intent(in) :: choice_left
real(8), dimension(num_choice_right,2), intent(in) :: choice_right
real(8), dimension(num_choice_bottom,2), intent(in) :: choice_bottom
real(8), dimension(num_choice_top,2), intent(in) :: choice_top
real(8), intent(in) :: error_inner, error_global
integer, intent(in) :: maxiter_inner, maxiter_global
real(8) :: dx, dy, unit_vector_nx(2:nx-1,2), unit_vector_ny(2:ny-1,2), dev_nx(2:nx-1), dev_ny(2:ny-1)
integer, dimension(nx) :: idx_bottom1, idx_top1
integer, dimension(ny) :: idx_left1, idx_right1
real(8) :: x(nx), y(ny)
real(8) :: lchange(4), lchange_total, err_distortion
integer :: ncount, iter, i
integer :: cr_t,cm_t,c1_t,c2_t, iforce
real(8) :: rate
real(8), dimension(4,nx,ny) :: distort0,distort1

! convergence history
open(33,file='convergence.dat',status='unknown')
write(33,*) 'Variables = i, residual'

iforce = 0
ncount = 0
! initialize the system clock
call system_clock(count_rate=cr_t)
call system_clock(count_max=cm_t)
rate = dble(cr_t)
call system_clock(c1_t)
do i=1,maxiter_global

    ! initial distortion
    call get_distortion_func(u,v,distort0)
!    print*,'Before: ',distort0(:,2,ny/2)
    call pde_RyskinLeal(u, v, distort0, nx, ny, 1.1, error_inner, maxiter_inner, iforce, iter)

    ! check distortion function
    call get_distortion_func(u,v,distort1)
!    print*,'After: ',distort1(:,2,ny/2)

    lchange = 0.d0

    if( boundary_fixed(1)==0 ) then
        call evaluate_left_boundary(u, v, nx, ny, unit_vector_ny, dev_ny)
        call choose_boundary(choice_left, num_choice_left, unit_vector_ny, dev_ny, idx_left, idx_left1, &
                                u(1,2:ny-1), v(1,2:ny-1), ny)
        lchange(1) = maxval(abs(idx_left-idx_left1))
        idx_left = idx_left1
    endif
    if( boundary_fixed(2)==0 ) then
        call evaluate_right_boundary(u, v, nx, ny, unit_vector_ny, dev_ny)
        call choose_boundary(choice_right, num_choice_right, unit_vector_ny, dev_ny, idx_right, idx_right1, &
                                u(nx,2:ny-1), v(nx,2:ny-1), ny)
        lchange(2) = maxval(abs(idx_right-idx_right1))
        idx_right = idx_right1
    else
        !call bend_right_boundary(u, v, nx, ny)
    endif

    if( boundary_fixed(3)==0 ) then
        call evaluate_bottom_boundary(u, v, nx, ny, unit_vector_nx, dev_nx)
        call choose_boundary(choice_bottom, num_choice_bottom, unit_vector_nx, dev_nx, idx_bottom, idx_bottom1, &
                                u(2:nx-1,1), v(2:nx-1,1), nx)
        lchange(3) = maxval(abs(idx_bottom-idx_bottom1))
        idx_bottom = idx_bottom1
    endif
    if( boundary_fixed(4)==0 ) then
        call evaluate_top_boundary(u, v, nx, ny, unit_vector_nx, dev_nx)
        call choose_boundary(choice_top, num_choice_top, unit_vector_nx, dev_nx, idx_top, idx_top1, &
                                u(2:nx-1,ny), v(2:nx-1,ny), nx)
        lchange(4) = maxval(abs(idx_top-idx_top1))
        idx_top = idx_top1
    endif
    lchange_total = maxval(lchange)

    if( boundary_symmetric(1)==1 ) then
        call symmetric_left(v, nx,ny)
    endif
    if( boundary_symmetric(2)==1 ) then
        call symmetric_right(v, nx,ny)
    endif
    if( boundary_symmetric(3)==1 ) then
        call symmetric_bottom(u, nx,ny)
    endif
    if( boundary_symmetric(4)==1 ) then
        call symmetric_top(u, nx,ny)
    endif

    err_distortion = maxval( abs((distort1-distort0)/(distort1+1.d-20)) )
    write(33,*) i, err_distortion

!    if(mod(i,maxiter_global/2-1).eq.0) then
!        iforce = 1
!    else
!        iforce = 0
!    endif
!    iforce = 1
    ncount = ncount + 1
!   print*, 'Counting...No.', ncount,',with change of', lchange_total, 'iteration = ', iter

    if(mod(ncount,nsave)==0) then
!        print*, 'Counting...No.', ncount,',with change of', lchange_total, ', err_distortion: ',err_distortion
        print*, 'Counting...No.', ncount,',with change of', lchange_total, 'iteration = ', iter
    endif
    if( int(lchange_total)<=max_change .or. err_distortion<=error_global ) then
        print*,'Converged! '
        exit
    endif

enddo
call system_clock(c2_t)
print*,'compute_grid time [sec]: ',dble(c2_t-c1_t)/rate

close(33)
end subroutine
!-------------------------------------------------------!
!	refine_boundary:	refine by adding one point in every 2 points
!	In:	boundary(nx):	original
!	Out:	boundary_new
!---------------------------------------------------!
subroutine refine_boundary(boundary, nx, boundary_new)
integer, intent(in) :: nx
real(8), dimension(nx,2), intent(in) :: boundary
real(8), dimension(nx*2-1,2), intent(out) :: boundary_new
integer i

boundary_new(1,:) = boundary(1,:)
do i = 2, nx
	boundary_new(2*i-2,:) = 5.d-1*( boundary(i,:) + boundary(i-1,:) )
	boundary_new(2*i-1,:) = boundary(i,:)
enddo

end subroutine
!-----------------------------------------------------!
!	symmetric_bottom:	symmetric at bottom boundary
!	In:	nx,ny
!	Inout:	u:	the whole grid-x
!---------------------------------------------------!
subroutine symmetric_bottom(u, nx,ny)
integer, intent(in) :: nx, ny
real(8), dimension(nx,ny), intent(inout) :: u

u(:,1) = u(:,2)

end subroutine
!-----------------------------------------------------!
subroutine symmetric_top(u, nx,ny)
integer, intent(in) :: nx, ny
real(8), dimension(nx,ny), intent(inout) :: u

u(:,ny) = u(:,ny-1)

end subroutine
!-----------------------------------------------------!
subroutine symmetric_left(v, nx,ny)
integer, intent(in) :: nx, ny
real(8), dimension(nx,ny), intent(inout) :: v

v(1,:) = v(2,:)

end subroutine
!-----------------------------------------------------!
subroutine symmetric_right(v, nx,ny)
integer, intent(in) :: nx, ny
real(8), dimension(nx,ny), intent(inout) :: v

v(nx,:) = v(nx-1,:)

end subroutine
!-----------------------------------------------------!
subroutine bend_right_boundary(u, v, nx, ny)
integer, intent(in) :: nx, ny
real(8), dimension(nx,ny), intent(inout) :: u, v
real(8), dimension(2) :: p0, p1
real(8) :: dev
integer :: j

do j = 2, ny-1
	p0(1) = u(nx,j+1) - u(nx,j-1)
	p0(2) = v(nx,j+1) - v(nx,j-1)
	p0 = p0 / sqrt( sum(p0**2) )
	p1(1) = u(nx-1,j) - u(nx,j)
	p1(2) = v(nx-1,j) - v(nx,j)
	
	dev = sum(p0*p1)
	u(nx-1,j) = u(nx-1,j) - dev*p0(1)
	v(nx-1,j) = v(nx-1,j) - dev*p0(2)
enddo
	
end subroutine
!---------------------------------------------------!
!	gridgen_checker:	checker method
!	In:	u_in,v_in(nx,ny):	the grid
!	Out:	u,v
!---------------------------------------------------!
subroutine gridgen_checker(u_in,v_in, u,v,nx,ny)
integer, intent(in) :: nx,ny
real(8), dimension(nx,ny), intent(in) :: u_in,v_in
real(8), dimension(nx,ny), intent(out) :: u,v
real(8), dimension(nx,ny) :: u1, v1
real(8), dimension(2) :: p00,p01,p10,p11, point
integer :: i, j, m, n, inter
integer :: cr_t,cm_t,c1_t,c2_t
real(8) :: rate

u = u_in
v = v_in
call init_TFI(u_in, v_in, u1, v1, nx, ny)

! initialize the system clock
call system_clock(count_rate=cr_t)
call system_clock(count_max=cm_t)
rate = dble(cr_t)
call system_clock(c1_t)
print*,'gridgen_checker ...'
!$OMP PARALLEL &
!$OMP private(i,j,m,n,p00,p01,p10,p11,inter, point)
!$OMP DO
do j = 2, ny-1
do i = 2, nx-1
	do n = 1, ny-1
	do m = 1, nx-1
		p00 = (/ u(i,n),v(i,n) /)
		p01 = (/ u(i,n+1),v(i,n+1) /)
		p10 = (/ u1(m,j),v1(m,j) /)
		p11 = (/ u1(m+1,j),v1(m+1,j) /)
		call get_intersect_point(p00, p01, p10, p11, inter, point)
		if( inter==1 ) then
			u(i,j) = point(1)
			v(i,j) = point(2)
			exit
		endif
	enddo
	if( inter==1 ) exit
	enddo
enddo
enddo
!$OMP END DO
!$OMP END PARALLEL
print*,' '
call system_clock(c2_t)
print*,'gridgen_checker time [sec]: ',dble(c2_t-c1_t)/rate
print*,' '

end subroutine
	
!-------------------------------------------------------!
! check the angle
subroutine check_angle_cell(u,v,nx,ny,theta)
integer,intent(in) :: nx,ny
real(8),dimension(nx,ny),intent(in) :: u,v
real(8),dimension(nx,ny),intent(out) :: theta
integer :: i,j
real(8),dimension(2,2) :: mm
real(8) :: g12,hksi,heta

do j=1,ny
do i=1,nx
    mm = Calmm_2D(i,j,u,v)
    g12 = mm(1,1)*mm(1,2) + mm(2,1)*mm(2,2)
    hksi = sqrt((mm(1,1))**2+(mm(2,1))**2)
    heta = sqrt((mm(1,2))**2+(mm(2,2))**2)
    theta(i,j) = acos(g12/(hksi*heta+1.d-20))
enddo
enddo
end subroutine
!-------------------------------------------------------!
! check the aspect ratio
subroutine check_aspect_ratio(u,v,nx,ny,metric)
integer,intent(in) :: nx,ny
real(8),dimension(nx,ny),intent(in) :: u,v
real(8),dimension(nx,ny),intent(out) :: metric
integer :: i,j
real(8),dimension(2,2) :: mm
real(8) :: hksi,heta

do j=1,ny
do i=1,nx
    mm = Calmm_2D(i,j,u,v)
    hksi = sqrt((mm(1,1))**2+(mm(2,1))**2)
    heta = sqrt((mm(1,2))**2+(mm(2,2))**2)
    metric(i,j) = max(hksi,heta)/(min(hksi,heta)+1.d-20)
enddo
enddo
end subroutine
!-------------------------------------------------------!
! check the angle between vectors
function get_angle_vector(p00, p01, p10, p11) result(angle)
real(8), dimension(2), intent(in) :: p00, p01, p10, p11
real(8) :: angle
real(8) :: det, norm0, norm1
real(8), dimension(2) :: vec0, vec1

vec0 = p01 - p00
vec1 = p11 - p10
det = vec0(1)*vec1(1) + vec0(2)*vec1(2)
norm0 = sqrt(vec0(1)**2+vec0(2)**2)
norm1 = sqrt(vec1(1)**2+vec1(2)**2)
angle = acos(det/(norm0*norm1+1.d-20))

end function
!-------------------------------------------------------!
! check the smoothness of line (p0,p1) and (p1,p2)
function get_smoothness(p0, p1, p2) result(smooth)
real(8), dimension(2), intent(in) :: p0, p1, p2
real(8) :: smooth
real(8) :: det, norm0, norm1
real(8), dimension(2) :: vec0, vec1

vec0 = p0-p1
vec1 = p2-p1
det = vec0(1)*vec1(1) + vec0(2)*vec1(2)
norm0 = sqrt(vec0(1)**2+vec0(2)**2)
norm1 = sqrt(vec1(1)**2+vec1(2)**2)
smooth = 0.5d0*(1.d0-det/(norm0*norm1+1.d-20))

end function

end module cm_gridgen
