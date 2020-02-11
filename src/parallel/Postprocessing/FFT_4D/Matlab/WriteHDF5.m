function  WriteHDF5(hslab, buffer,file_exist,group_exist)
if nargin < 4
    group_exist=0;
end


if file_exist
    file=H5F.open(hslab.fname,'H5F_ACC_RDWR','H5P_DEFAULT');
else
    file=H5F.create(hslab.fname,'H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');    
end

if hslab.gname ~= '/'
    if group_exist
        group=H5G.open(file,hslab.gname);
    else
        group=H5G.create(file,hslab.gname,'H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
    end
end
    
    


space_id=H5S.create_simple(hslab.rank,fliplr(hslab.dimsm),[]);
if hslab.gname ~= '/'
    dset=H5D.create(group,hslab.dname,'H5T_NATIVE_DOUBLE',space_id,'H5P_DEFAULT');
else
    dset=H5D.create(file,hslab.dname,'H5T_NATIVE_DOUBLE',space_id,'H5P_DEFAULT');
end

H5D.write(dset,'H5T_NATIVE_DOUBLE','H5S_ALL','H5S_ALL','H5P_DEFAULT',buffer);

H5D.close(dset)
H5S.close(space_id)

if hslab.gname ~= '/'
    H5G.close(group);
end
H5F.close(file)


