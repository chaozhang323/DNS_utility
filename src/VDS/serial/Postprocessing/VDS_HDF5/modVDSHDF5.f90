 
  module modVDSHDF5
    use modRWHDF5
    implicit none

    integer(HID_T), private :: file_id, group_id, dset_id, dspace_id, memspace_id, prp_id, src_space_id
    integer, private :: hdferr
contains
     ! Create 3D virtual dataset
     subroutine CreateHDF5_VDS_3D(hslab,num_subset,hslab_subset)
       integer, intent(in) :: num_subset
       type(tp_rdwt_hdf5), intent(in) :: hslab, hslab_subset(num_subset)
       integer :: n, nn

       CALL h5fcreate_f(trim(hslab%fname), H5F_ACC_TRUNC_F, file_id, hdferr)
       CALL h5screate_simple_f(hslab%rank, hslab%dimsf, dspace_id, hdferr)
       CALL H5Pcreate_f(H5P_DATASET_CREATE_F,prp_id,hdferr)
       CALL H5Pset_fill_value_f(prp_id, H5T_NATIVE_DOUBLE, 1.e30,hdferr)

       do n=1, hslab%dnum
         !print *, 'linking dataset...', hslab%dname(n)
         do nn=1,num_subset
           !print *, 'linking file ... ', trim(hslab_subset(nn)%fname)
           CALL h5screate_simple_f(hslab_subset(nn)%rank, hslab_subset(nn)%dimsf, src_space_id, hdferr) ! create source space
           CALL h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hslab_subset(nn)%offset, &
                                      hslab_subset(nn)%count, hdferr, hslab_subset(nn)%stride, hslab_subset(nn)%dimsf) !select location to be linked
           CALL h5pset_virtual_f(prp_id,dspace_id,hslab_subset(nn)%fname,hslab_subset(nn)%dname(n), src_space_id, hdferr) !set link
           CALL h5sclose_f(src_space_id, hdferr)
         enddo!end nn loop
         CALL h5dcreate_f(file_id, hslab%dname(n), H5T_NATIVE_DOUBLE, dspace_id, dset_id, hdferr,prp_id ) !create data set
         CALL h5dclose_f(dset_id, hdferr)
       enddo !end n loop

       CALL h5sclose_f(dspace_id, hdferr)
       Call h5pclose_f(prp_id, hdferr)
       CALL h5fclose_f(file_id, hdferr)

     end subroutine CreateHDF5_VDS_3D


  end module modVDSHDF5
