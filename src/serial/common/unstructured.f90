module uns2str
use omp_lib
contains
subroutine get_element_center(nodemap, num_elements, np_per_element, coor_point, num_points, coor_element)
implicit none
integer, intent(in) :: num_elements, np_per_element, num_points
integer, dimension(num_elements, np_per_element), intent(in) :: nodemap
real(8), dimension(num_points,2), intent(in) :: coor_point
real(8), dimension(num_elements, 2), intent(out) :: coor_element
integer :: i, n, point(np_per_element)
real(8) :: weight

weight = 1.d0 / real(np_per_element, 8)
coor_element = 0.
do i = 1, num_elements
!	point = nodemap(i, :) + 1
    point = nodemap(i, :)
!	coor_element(i, :) = 2.5d-1* (coor_point(point(1), :)+coor_point(point(2), :)+coor_point(point(3), :)+coor_point(point(4), :))
	do n = 1, np_per_element
		coor_element(i, :) = coor_element(i, :) + coor_point(point(n), :)
	enddo
!	coor_element(i, :) = coor_element(i, :) * 0.25d0
enddo
coor_element= coor_element * weight
end subroutine

!------------------------------------------------------------------!
subroutine get_position(p0, p, pstn)
implicit none
real(8), dimension(2), intent(in) :: p0, p
integer, intent(out) :: pstn
real(8), dimension(2) :: s

s = abs(p-p0)
!print*,'s',s,'p',p,p0
if( s(1)>s(2) ) then
	if( p(1)<p0(1) ) then
		pstn = 1
	else
		pstn = 2
	endif
else
	if( p(2)<p0(2) ) then
		pstn = 3
	else
		pstn = 4
	endif
endif

end subroutine

!----------------------------------------------------------------!
subroutine get_neighbor_element(nodemap, num_elements, no, neighbors, num_neighbors)
implicit none
integer, intent(in) :: num_elements, no
integer, dimension(num_elements,4), intent(in) :: nodemap
integer, intent(out) :: neighbors(4), num_neighbors
integer :: i,j, points_ele(4), n1, n2, p1, p2, sw1, sw2

points_ele = nodemap(no, :)
num_neighbors = 0
neighbors = 0
do n1 = 1, 3
	do n2 = n1+1, 4
		p1 = points_ele(n1)
		p2 = points_ele(n2)
		
		do i = 1, num_elements
			if( i/=no ) then
				sw1 = 0
				sw2 = 0
				do j = 1, 4
					if( nodemap(i,j)==p1 ) sw1 = 1
					if( nodemap(i,j)==p2 ) sw2 = 1
				enddo
				if(sw1==1 .and. sw2==1) then
					num_neighbors = num_neighbors + 1
					neighbors(num_neighbors) = i
					exit
				endif
			endif
		enddo
	enddo
enddo

end subroutine
!----------------------------------------------------------------!
subroutine get_structured_shape(structured_map, num_elements, start_point, shape)
integer, intent(in) :: num_elements
integer, dimension(num_elements, 4), intent(in) :: structured_map
integer, intent(out) :: shape(2), start_point
integer :: n

n = 1
do 
	if( structured_map(n,1)==-1 ) then
		exit
	else
		n = structured_map(n,1)
	endif
enddo
do 
	if( structured_map(n,3)==-1 ) then
		exit
	else
		n = structured_map(n,3)
	endif
enddo
start_point = n
!---shape---!
shape = 1
n = start_point
do 
	if( structured_map(n,2)==-1 ) then
		exit
	else
		n = structured_map(n,2)
		shape(1) = shape(1) + 1
	endif
enddo
n = start_point
do 
	if( structured_map(n,4)==-1 ) then
		exit
	else
		n = structured_map(n,4)
		shape(2) = shape(2) + 1
	endif
enddo


end subroutine
!----------------------------------------------------------------!

!----------------------------------------------------------------!
! build structure map to nodal
subroutine rebuild_structured_nodal(structured_map, nodemap, num_elements, data_uns, num_nodes, start_point, nx,ny,grid_nodal)
integer,intent(in) :: nx,ny,num_elements,num_nodes,start_point
integer,dimension(num_elements,4),intent(in) :: structured_map,nodemap
real(8),dimension(num_nodes),intent(in) :: data_uns
integer,dimension(nx+1,ny+1),intent(out) :: grid_nodal
integer,dimension(nx,ny) :: grid
integer,dimension(num_elements,4) :: nodemap_order

imax = nx+1
jmax = ny+1
! build structure map to cell center
call get_structured_grid(structured_map, num_elements, grid, start_point, nx, ny)
! get nodemap order
call get_nodemap_order(structured_map, nodemap, num_elements, data_uns, num_nodes, grid, nx, ny, nodemap_order)
! build structure map to nodal
call get_structured_grid_nodal(nodemap, nodemap_order, num_elements, grid, imax, jmax, grid_nodal)
end subroutine
!----------------------------------------------------------------!
subroutine rebuild_structured(nodemap, num_elements, coor_element, structured_map, start_point, shape)
implicit none
integer, intent(in) :: num_elements
integer, dimension(num_elements, 4), intent(in) :: nodemap
real(8), dimension(num_elements, 2), intent(in) :: coor_element
integer, dimension(num_elements, 4), intent(out) :: structured_map
integer, intent(out) :: shape(2), start_point
integer :: neighbors(4), num_neighbors
integer :: i, n, pstn
real(8) :: time_begin,time_end

!call get_element_center(nodemap, num_elements, 4, coor_point, num_points, coor_element)
structured_map = -1
time_begin = omp_get_wtime()
!$OMP PARALLEL &
!$OMP private(neighbors,num_neighbors,pstn,i,n)
!$OMP DO
do i = 1, num_elements
	call get_neighbor_element(nodemap, num_elements, i, neighbors, num_neighbors)
	do n = 1, num_neighbors
		call get_position(coor_element(i,:), coor_element(neighbors(n),:), pstn)
		structured_map(i, pstn) = neighbors(n)
	enddo
enddo
!$OMP END DO
!$OMP END PARALLEL
print*,' '
print*,'rebuild_structured OMP time [sec]: ',omp_get_wtime()-time_begin
print*,' '

call get_structured_shape(structured_map, num_elements, start_point, shape)
!print*,shape, 'Start at', start_point

end subroutine
!----------------------------------------------------------------!
! get the order of nodes around an element
! 4      3
! +------+
! - elem -
! -      -
! +------+
! 1      2
subroutine get_nodemap_order(structured_map, nodemap, num_elements, grid_uns, num_nodes, grid, nx, ny, nodemap_order)
integer,intent(in) :: num_elements, num_nodes
integer,dimension(num_elements,4),intent(in) :: structured_map, nodemap
real(8),dimension(num_nodes),intent(in) :: grid_uns
integer,dimension(nx,ny),intent(in) :: grid
integer,dimension(num_elements,4),intent(out) :: nodemap_order
integer :: i,j,ip,m,n,n1,n2,sw1,sw2,p1,p2

nodemap_order = -1
do j = 1,ny
do i = 2,nx
    n = grid(i,j)
    ! left neighbor
    m = structured_map(n,1)
    do n1 = 1, 3
    do n2 = n1+1,4
        p1 = nodemap(m,n1)
        p2 = nodemap(m,n2)
        sw1 = 0; sw2 = 0
        do ip = 1,4
            if(nodemap(n,ip) == p1) sw1 = ip
            if(nodemap(n,ip) == p2) sw2 = ip
        enddo
        if(sw1 > 0 .and. sw2 > 0) then
            if(grid_uns(nodemap(n,sw1)) > grid_uns(nodemap(n,sw2))) then
                nodemap_order(n,sw1) = 4
                nodemap_order(n,sw2) = 1
            else
                nodemap_order(n,sw1) = 1
                nodemap_order(n,sw2) = 4
            endif
        endif
    enddo
    enddo
enddo
enddo
do j = 1,ny
do i = 1,nx-1
    n = grid(i,j)
    ! right neighbor
    m = structured_map(n,2)
    do n1 = 1, 3
    do n2 = n1+1,4
        p1 = nodemap(m,n1)
        p2 = nodemap(m,n2)
        sw1 = 0; sw2 = 0
        do ip = 1,4
            if(nodemap(n,ip) == p1) sw1 = ip
            if(nodemap(n,ip) == p2) sw2 = ip
        enddo
        if(sw1 > 0 .and. sw2 > 0) then
            if(grid_uns(nodemap(n,sw1)) > grid_uns(nodemap(n,sw2))) then
                nodemap_order(n,sw1) = 3
                nodemap_order(n,sw2) = 2
            else
                nodemap_order(n,sw1) = 2
                nodemap_order(n,sw2) = 3
            endif
        endif
    enddo
    enddo
enddo
enddo
do j=1,ny
    ! left boundary
    n = grid(1,j)
    do i=1,4
        if(nodemap_order(n,i) == -1) then
            sw1 = i
            exit
        endif
    enddo
    do i=4,1,-1
        if(nodemap_order(n,i) == -1) then
            sw2 = i
            exit
        endif
    enddo
    if(grid_uns(nodemap(n,sw1)) > grid_uns(nodemap(n,sw2))) then
        nodemap_order(n,sw1) = 4
        nodemap_order(n,sw2) = 1
    else
        nodemap_order(n,sw1) = 1
        nodemap_order(n,sw2) = 4
    endif
    ! right boundary
    n = grid(nx,j)
    do i=1,4
        if(nodemap_order(n,i) == -1) then
            sw1 = i
            exit
        endif
    enddo
    do i=4,1,-1
        if(nodemap_order(n,i) == -1) then
            sw2 = i
            exit
        endif
    enddo
    if(grid_uns(nodemap(n,sw1)) > grid_uns(nodemap(n,sw2))) then
        nodemap_order(n,sw1) = 3
        nodemap_order(n,sw2) = 2
    else
        nodemap_order(n,sw1) = 2
        nodemap_order(n,sw2) = 3
    endif
enddo
end subroutine
!----------------------------------------------------------------!
subroutine get_structured_grid(structured_map, num_elements, grid, start_point, nx, ny)
integer, intent(in) :: num_elements, start_point
integer, intent(in) :: nx, ny
integer, dimension(num_elements,4), intent(in) :: structured_map
integer, dimension(nx, ny), intent(out) :: grid
integer :: i, j, nj, n

nj = start_point
do j = 1, ny
	n = nj
	do i = 1, nx
		grid(i,j) = n
		n = structured_map(n, 2)
	enddo
	nj = structured_map(nj, 4)
enddo

end subroutine
!----------------------------------------------------------------!
! map unstructured quadrilateral mesh to structured mesh by nodal point
subroutine get_structured_grid_nodal(nodemap, nodemap_order, num_elements, grid, nx, ny, grid_nodal)
integer,intent(in) :: num_elements,nx,ny
integer,dimension(num_elements,4),intent(in) :: nodemap, nodemap_order
integer,dimension(nx-1,ny-1),intent(in) :: grid
integer,dimension(nx,ny),intent(out) :: grid_nodal
integer :: i,j,n

do j = 1, ny-1
do i = 1, nx-1
    n = grid(i,j)
    grid_nodal(i,j) = nodemap(n,nodemap_order(n,1)) ! left bottom
enddo
enddo
! right
do j = 1, ny-1
    n = grid(nx-1,j)
    grid_nodal(nx,j) = nodemap(n,nodemap_order(n,2))    ! right bottom
enddo
! top
do i = 1, nx-1
    n = grid(i,ny-1)
    grid_nodal(i,ny) = nodemap(n,nodemap_order(n,4))    ! left top
enddo
! right top
n = grid(nx-1,ny-1)
grid_nodal(nx,ny) = nodemap(n,nodemap_order(n,3))
end subroutine
!----------------------------------------------------------------!
subroutine get_structured_data_nodal(grid_nodal,nx,ny,data_uns,num_nodes,data_str)
integer,intent(in) :: nx,ny,num_nodes
integer,dimension(nx,ny),intent(in) :: grid_nodal
real(8),dimension(num_nodes),intent(in) :: data_uns
real(8),dimension(nx,ny),intent(out) :: data_str
integer :: i,j

do j=1,ny
do i=1,nx
    data_str(i,j) = data_uns(grid_nodal(i,j))
enddo
enddo
end subroutine
!----------------------------------------------------------------!
subroutine get_structured_data(grid, nx, ny, data_uns, data_str, num_elements)
integer, intent(in) :: nx, ny, num_elements
integer, dimension(nx, ny), intent(in) :: grid
real(8), dimension(num_elements), intent(in) :: data_uns
real(8), dimension(nx, ny), intent(out) :: data_str
integer :: i, j

do j = 1, ny
	do i = 1, nx
		data_str(i,j) = data_uns(grid(i,j))
	enddo
enddo

end subroutine
endmodule






