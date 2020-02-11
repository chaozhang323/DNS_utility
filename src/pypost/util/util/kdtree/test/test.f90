program test
use hdf5
use kdtree
use kriging

implicit none
integer(hid_t) :: f_id, dset_id, dspace_id, dtype_id
integer :: hdferr, n
integer(hsize_t) :: npoints, shape_(2)
integer(size_t) :: dset_count


character(len=256) :: flnm = "/usr/local/home/yl8bc/dev_local/grid_interpolation/data_fluent.h5"
!character(len=256) :: flnm = "../sub.h5"
character(len=256), dimension(4) :: dsetnm = (/'X', 'Y', 'u', 'p'/)
real(kind=8), dimension(:,:), allocatable :: data_, data_f

type(tree_master_record), pointer :: tree_id
integer,parameter :: n_nearest = 4, nx0 = 2, rank=2
integer, dimension(n_nearest, nx0) :: indexes_nearest
real(kind=8), dimension(n_nearest, nx0) :: distance_nearest
real(kind=8) :: x0(rank, nx0), theta(rank,4-rank,nx0), p(rank,4-rank,nx0), weights(n_nearest, 4-rank,nx0), f_interp(4-rank,nx0)
!! ---
call h5open_f(hdferr)
call h5eset_auto_f(1, hdferr)

call h5fopen_f(flnm, h5f_acc_rdonly_f, f_id, hdferr)
call h5dopen_f(f_id, dsetnm(1), dset_id, hdferr)
call h5dget_type_f(dset_id, dtype_id, hdferr)
call h5dget_space_f(dset_id, dspace_id, hdferr)
call h5sget_simple_extent_npoints_f(dspace_id, npoints, hdferr)
print*, 'Num of points = ', npoints
shape_(1) = npoints
shape_(2) = len(dsetnm)
allocate(data_(shape_(1), shape_(2)))
allocate(data_f(shape_(2), shape_(1)))
data_ = 0.d0
call h5dread_f(dset_id, dtype_id, data_(:,1), shape_(1:1), hdferr, file_space_id=dspace_id)
!print*, data_(1:10,1)

do n = 2, 4
    call h5dopen_f(f_id, dsetnm(n), dset_id, hdferr)
    call h5dread_f(dset_id, dtype_id, data_(:,n), shape_(1:1), hdferr, file_space_id=dspace_id)
enddo


call h5fget_obj_count_f(f_id, h5f_obj_dataset_f, dset_count, hdferr)
print*, 'Count dsets:', dset_count


!! --
data_f = transpose(data_)

tree_id => create_tree(data_(:,1:rank))

!! --
x0(:,1) = (/0.39d0, 0.057d0/)
x0(:,2) = (/0.39d0, 0.025d0/)
do n = 1, nx0
    call n_nearest_to(tree_id, x0(:,n), n_nearest, indexes_nearest(:,n), distance_nearest(:,n))
enddo
!print*, indexes_nearest


call init_params_naive(theta, p)
call weights_ok(data_f(1:2,:), indexes_nearest, x0, weights, theta, p)
!call kernel_ok(data_(indexes_nearest,1:2), n_nearest, x0, 2, w, theta, p)
print*, weights

call interp_ok(data_f(3:4,:), indexes_nearest, weights, f_interp)
print*, f_interp

call grdd_kriging(data_f(1:2,:), data_f(3:4,:), indexes_nearest, theta, p, 1.d-6, 1.d-6, 1000)

call weights_ok(data_f(1:2,:), indexes_nearest, x0, weights, theta, p)
print*, weights
call interp_ok(data_f(3:4,:), indexes_nearest, weights, f_interp)
print*, f_interp
end program
