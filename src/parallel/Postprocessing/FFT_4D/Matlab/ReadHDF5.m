function rdata = ReadHDF5(hslab)


  file = H5F.open(hslab.fname,'H5F_ACC_RDONLY', 'H5P_DEFAULT');
  group = H5G.open(file,hslab.gname);
  hslab.block = hslab.dimsm;
  dset = H5D.open(group,hslab.dname); 
  dspace = H5D.get_space(dset);
  mem_space_id = H5S.create_simple(hslab.rank,fliplr(hslab.dimsm),[]);
  H5S.select_hyperslab(dspace,'H5S_SELECT_SET', fliplr(hslab.offset), ...
                       fliplr(hslab.stride), fliplr(hslab.count), fliplr(hslab.block) );  

  rdata= H5D.read(dset,'H5T_NATIVE_DOUBLE', mem_space_id, dspace,'H5P_DEFAULT');

  H5S.close(mem_space_id);
  H5D.close(dset);
  H5S.close(dspace);

  H5G.close(group);
  H5F.close(file);

end 
