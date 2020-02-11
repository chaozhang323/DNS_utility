function [] = WriteScalar(hslab,scalar,file_exist,group_exist)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

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

space_id=H5S.create_simple(1,1,[]);
if hslab.gname ~= '/'
    dset=H5D.create(group,hslab.sname,'H5T_NATIVE_DOUBLE',space_id,'H5P_DEFAULT');
else
    dset=H5D.create(file,hslab.sname,'H5T_NATIVE_DOUBLE',space_id,'H5P_DEFAULT');
end

H5D.write(dset,'H5T_NATIVE_DOUBLE','H5S_ALL','H5S_ALL','H5P_DEFAULT',scalar);

H5D.close(dset)
H5S.close(space_id)

if hslab.gname ~= '/'
    H5G.close(group);
end
H5F.close(file)



end

