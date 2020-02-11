module uns2str
implicit none
contains
!! --------------------- !!
subroutine wrapper(x_in, y_in, data_in, &
                   x_out, y_out, data_out, &
                   nx,ny, nodemap, num_elements, num_points, num_vars)
integer, intent(in) :: num_elements, num_points, num_vars
real(kind=8), dimension(num_points), intent(in) :: x_in,y_in
real(kind=8), dimension(num_vars, num_points), intent(in) :: data_in
integer, dimension(num_elements,4), intent(in) :: nodemap
integer, intent(out) :: nx, ny
real(kind=8), dimension(num_elements), intent(out) :: x_out, y_out
real(kind=8), dimension(num_vars, num_elements), intent(out) :: data_out

integer, dimension(:,:), allocatable :: grid
real(kind=8), dimension(:,:), allocatable :: x_tmp, y_tmp
real(kind=8), dimension(:,:,:), allocatable :: data_tmp

real(kind=8), dimension(num_points,2) :: coor_point
real(kind=8), dimension(num_elements,2) :: coor_element
integer, dimension(num_elements) :: num_neighbors, elemap
integer num_bdry
logical, dimension(num_elements) :: mask_bdry
integer, dimension(:), allocatable :: index_bdry_ordered, num_neighbors_bdry
integer domain_shape(2)
integer n

call get_element_center(nodemap, num_elements, 4, coor_point, num_points, coor_element)
call get_rough_elemap(nodemap, num_elements, elemap, num_neighbors)
call get_num_boundary_elements(num_neighbors, num_elements, mask_bdry, num_bdry)
call get_boundary_elements(elemap, num_neighbors, num_elements, &
                           mask_bdry, index_bdry_ordered, num_neighbors_bdry, num_bdry)
allocate(index_bdry_ordered(num_bdry), num_neighbors_bdry(num_bdry))
call get_domain_shape(num_neighbors_bdry, num_bdry, domain_shape)
nx = domain_shape(1)
ny = domain_shape(2)
allocate(grid(nx,ny))
allocate(x_tmp(nx,ny), y_tmp(nx,ny))
allocate(data_tmp(num_vars,nx,ny))
call rebuild_structured_grid(elemap, num_neighbors, num_elements, &
                      index_bdry_ordered, num_neighbors_bdry, num_bdry, grid, nx, ny)

call get_structured_data(grid, nx, ny, x_in, x_tmp, num_elements)
call get_structured_data(grid, nx, ny, y_in, y_tmp, num_elements)
do n = 1, num_vars
call get_structured_data(grid, nx, ny, data_in(n,:), data_tmp(n,:,:), num_elements)
enddo

x_out = reshape(x_tmp, (/num_elements/))
y_out = reshape(y_tmp, (/num_elements/))
data_out = reshape(data_tmp, (/num_vars,num_elements/))

end subroutine
!! ------------------------------------------------- !!


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
	point = nodemap(i, :) + 1
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
subroutine rebuild_structured(nodemap, num_elements, coor_element, structured_map, start_point, shape)
implicit none
integer, intent(in) :: num_elements
integer, dimension(num_elements, 4), intent(in) :: nodemap
real(8), dimension(num_elements, 2), intent(in) :: coor_element
integer, dimension(num_elements, 4), intent(out) :: structured_map
integer, intent(out) :: shape(2), start_point
integer :: neighbors(4), num_neighbors
integer :: i, n, pstn

!call get_element_center(nodemap, num_elements, 4, coor_point, num_points, coor_element)
structured_map = -1
do i = 1, num_elements
	call get_neighbor_element(nodemap, num_elements, i, neighbors, num_neighbors)
    do n = 1, num_neighbors
		call get_position(coor_element(i,:), coor_element(neighbors(n),:), pstn)
		structured_map(i, pstn) = neighbors(n)
		
	enddo
enddo

call get_structured_shape(structured_map, num_elements, start_point, shape)

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
!---------------------------------------------------------------!
subroutine get_rough_elemap(nodemap, num_elements, elemap, num_neighbors)
integer, intent(in) :: num_elements
integer, dimension(num_elements,4), intent(in) :: nodemap
integer, dimension(num_elements,4), intent(out) :: elemap
integer, dimension(num_elements), intent(out) :: num_neighbors
integer i

do i = 1, num_elements
    call get_neighbor_element(nodemap, num_elements, i, elemap(i,:), num_neighbors(i))
enddo
end subroutine
!--------------------------------------------------------------!!
!subroutine search_bordered_element(index_ele, not_ele, num_ele, no, ele_found)
!integer, intent(in) :: num_ele, no
!integer, intent(in) :: index_ele(num_ele), not_ele(:)
!integer, intent(out) :: ele_found


subroutine get_num_boundary_elements(num_neighbors, num_elements, mask_bdry, num_bdry)
integer, intent(in) :: num_elements
integer, dimension(num_elements), intent(in) :: num_neighbors
integer, intent(out) :: num_bdry
logical, dimension(num_elements), intent(out) :: mask_bdry

!------get bdry elements-------------!
mask_bdry = (num_neighbors<=3)
num_bdry = count(mask_bdry)
end subroutine
!------------------------------------------------------------!
subroutine get_boundary_elements(elemap, num_neighbors, num_elements, &
                                    mask_bdry, index_bdry_ordered, num_neighbors_bdry, num_bdry)
integer, intent(in) :: num_elements
integer, dimension(num_elements, 4), intent(in) :: elemap
integer, dimension(num_elements), intent(in) :: num_neighbors
integer, intent(in) :: num_bdry
logical, dimension(num_elements), intent(in) :: mask_bdry
integer, dimension(num_bdry), intent(out) :: index_bdry_ordered, num_neighbors_bdry

integer :: num_corners
logical, dimension(num_elements) :: mask_corners
integer :: i, n, n1,n2, m
integer :: index_all(num_elements)
integer :: index_bdry(num_bdry), elemap_bdry(num_bdry,4)
integer, allocatable :: index_corners(:)

index_all = (/ (i, i=1,num_elements) /)
!------get bdry elements-------------!
index_bdry = pack(index_all, mask_bdry)
elemap_bdry = elemap(index_bdry,:)
!--corner_element--!
mask_corners = (num_neighbors==2)
num_corners = count(mask_corners)
allocate(index_corners(num_corners))
index_corners = pack(index_all, mask_corners)
!-----------------------------------------!
index_bdry_ordered(1) = index_corners(1)
index_bdry_ordered(2) = elemap(index_corners(1), 1)
do i = 3, num_bdry
    n1 = index_bdry_ordered(i-1)    !current ele
    n2 = index_bdry_ordered(i-2)    !previous ele
    do m = 1, num_neighbors(n1) !num of neighbors
        n = elemap(n1,m)    !investigated ele
        if( n/=n2 .and.  any(n==index_bdry) ) then  !ele is not the previous one and in boudary
            index_bdry_ordered(i) = n
            exit
        endif
    enddo
enddo
num_neighbors_bdry = num_neighbors(index_bdry_ordered)
end subroutine
!-----------------------------------------------------!
subroutine get_domain_shape(num_neighbors_bdry, num_bdry, domain_shape)
integer, intent(in) :: num_bdry
integer, dimension(num_bdry), intent(in) :: num_neighbors_bdry
integer, intent(out) :: domain_shape(2)

integer i, rem

do i = 2, num_bdry
    if(num_neighbors_bdry(i)==2) then
        domain_shape(1) = i
        exit
    endif
enddo
rem = mod(num_bdry+4, 2)
if(rem/=0) print*, 'Severe warning: the struture building is not working well!'
domain_shape(2) = (num_bdry+4)/2 - domain_shape(1)

end subroutine
!-------------------------------------------------!
subroutine fill_boundary_index(index_bdry, num_bdry, grid, nx, ny)
integer, intent(in) :: num_bdry
integer, dimension(num_bdry), intent(in) :: index_bdry
integer, intent(in) :: nx, ny
integer, dimension(nx,ny), intent(inout) :: grid

integer i,j, n

n = 0
do i = 1, nx
    n = n + 1
    grid(i, 1) = index_bdry(n)
enddo
do j = 2, ny
    n = n + 1
    grid(nx, j) = index_bdry(n)
enddo
do i = nx-1, 1, -1
    n = n + 1
    grid(i, ny) = index_bdry(n)
enddo
do j = ny-1, 2, -1
    n = n + 1
    grid(1, j) = index_bdry(n)
enddo

if(n/=num_bdry) then
    print*, 'Severe warning: Boundary points failed to match!'
endif

end subroutine
!--------------------------------------------!
subroutine fill_interior_index(elemap, num_neighbors, num_elements, &
                                grid, nx, ny)
integer, intent(in) :: num_elements
integer, intent(in) :: nx, ny
integer, intent(in) :: elemap(num_elements,4), num_neighbors(num_elements)
integer, intent(inout) :: grid(nx,ny)
integer i, j, ew, es, m, n, ifilled

!! boundary already filled
!! fill interior
do j = 2, ny-1
do i = 2, nx-1
    ew = grid(i-1,j)    !ele west
    es = grid(i,j-1)    !ele south
    ifilled = 0
    do m = 1, num_neighbors(ew)
    do n = 1, num_neighbors(es)
        if( elemap(ew,m)==elemap(es,n) .and. all(elemap(ew,m)/=grid) ) then
            grid(i,j) = elemap(ew,m)
            ifilled = 1
            exit
        endif
    enddo
        if(ifilled==1) exit
    enddo
    if(ifilled==0) then
1100 format( 'Severe warning: Failed to find a neighbored cell at (',I4,I4,')!')
        write(*,1100) i,j
    endif
enddo
enddo
end subroutine
!---------------------------------!
subroutine regain_structured_boundary(nodemap, elemap, num_neighbors, num_elements, &
                                        mask_bdry, num_bdry)
integer, intent(in) :: num_elements
integer, dimension(num_elements,4), intent(in) :: nodemap
integer, dimension(num_elements), intent(out) :: num_neighbors
integer, dimension(num_elements), intent(out) :: elemap
integer, intent(out) :: num_bdry
logical, dimension(num_elements), intent(out) :: mask_bdry

call get_rough_elemap(nodemap, num_elements, elemap, num_neighbors)
call get_num_boundary_elements(num_neighbors, num_elements, mask_bdry, num_bdry)
end subroutine
!------------------------------!
subroutine regain_structured_shape(elemap, num_neighbors, num_elements, &
                                    mask_bdry, index_bdry_ordered, num_neighbors_bdry, num_bdry,&
                                    domain_shape)
integer, intent(in) :: num_elements
integer, dimension(num_elements), intent(in) :: elemap, num_neighbors
integer, intent(in) :: num_bdry
logical, dimension(num_elements), intent(in) :: mask_bdry
integer, dimension(num_bdry), intent(out) :: index_bdry_ordered, num_neighbors_bdry
integer, intent(out) :: domain_shape(2)

call get_boundary_elements(elemap, num_neighbors, num_elements, &
                                    mask_bdry, index_bdry_ordered, num_neighbors_bdry, num_bdry)
call get_domain_shape(num_neighbors_bdry, num_bdry, domain_shape)

end subroutine
!---------------------------------!
subroutine rebuild_structured_grid(elemap, num_neighbors, num_elements, &
                                    index_bdry, num_neighbors_bdry, num_bdry, grid, nx, ny)
integer, intent(in) :: num_elements
integer, dimension(num_elements,4), intent(in) :: elemap
integer, dimension(num_elements), intent(in) :: num_neighbors
integer, intent(in) :: num_bdry
integer, dimension(num_bdry), intent(in) :: index_bdry, num_neighbors_bdry
integer, intent(in) :: nx, ny
integer, dimension(nx,ny), intent(out) :: grid

grid = 0
call fill_boundary_index(index_bdry, num_bdry, grid, nx, ny)
call fill_interior_index(elemap, num_neighbors, num_elements,  &
                                grid, nx, ny)
end subroutine

endmodule






