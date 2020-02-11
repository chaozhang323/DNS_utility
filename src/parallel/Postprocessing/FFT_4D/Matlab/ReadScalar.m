

function ttime = ReadScalar(hslab)
  file = H5F.open(hslab.fname, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
  dset = H5D.open(file,hslab.sname);
  space = H5D.get_space(dset);
 
  ttime = H5D.read (dset,'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT');
  H5D.close (dset);
  H5S.close (space);
  H5F.close (file);
end
