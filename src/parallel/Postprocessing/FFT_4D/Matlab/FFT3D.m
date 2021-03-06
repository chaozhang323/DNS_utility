clc; clear; format long;

%%%%%%%%%%%%%%%% Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This program takes multiple flow data files and 1 grd file in order to
%produce a plt file with a zone for each flow file

%%%%%%%%%%%%%%%% Explanation of Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%File_Dim = how many dimensions that a file has
%FlowFileLoc = location of your flowdata files
%GridFileLoc = location of your grid file
%OutputFileName = What to name output file
%Dnum = how many varibles to read
%Dname = the name of the variables
%max = file dimensions
%st = start of data selection
%end = end of data selection
%fstart = file number you want to start reading at
%fstride = Skip these many files
%fend = stop reading files at this file number


%%%%%%%%%%%%%%%% START OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
File_Dim = 3;

FlowFileLoc='/usr/local/home/czb58/duanforge/duanl/Acoustics/M6_Purdue_Convergence/AIII_2880x320x786_Res_zod0p925_timeseries_OldRHS_NewPHDF5/TIMESERIES/';
GridFileLoc='/usr/local/home/czb58/duanforge/czb58/M6_AIII_2880x320x786/TIMESERIES/';
 
Dnum=5; Dname={'p' 'T' 'u' 'v' 'w' };

imax=1381;jmax=320; kmax=630;
ist =801; jst =  1; kst =301;
iend=830; jend= 20; kend=310;

file_be= 340050; file_end=340050; file_skip=25;

dt = 6.25e-7;

%%%%%%%%%%%%%%%% END OF INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Intialize Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fsol=InitFlowHDF5(File_Dim);

idimsm=iend-ist+1;
jdimsm=jend-jst+1;
kdimsm=kend-kst+1;

fsol.dimsm=[kdimsm, idimsm, jdimsm];
fsol.offset=[kst-1, ist-1, jst-1];

grdname={'x', 'y', 'z'};

nfiles = fix((file_end-file_be)/file_skip+1);
fnum=file_be;

buffer=zeros(kdimsm,idimsm,jdimsm,Dnum,nfiles);
grd=zeros(kdimsm,idimsm,jdimsm,3);

if (idimsm+ist-1 > imax || jdimsm+jst-1 > jmax || kdimsm+kst-1 > kmax)
   error='File selection outside file dimensions'
   return
end

%%%%%%%%%%%%%%% Main Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read grid data
fsol.fname=strcat(GridFileLoc,'timeseriesVol_grid.h5');
for n=1:3
    fsol.dname=grdname{n};
    grd(:,:,:,n)=ReadHDF5(fsol);
end

dx = grd(1,2,1,1) - grd(1,1,1,1);
dy = grd(1,1,2,2) - grd(1,1,1,2);
dz = grd(2,1,1,1) - grd(1,1,1,3);

for i=1:idimsm
    grd(:,i,:,1) = (i-1)*dx; 
end
for j=1:jdimsm
    grd(:,:,j,2) = (j-1)*dy; 
end
for k=1:kdimsm
    grd(k,:,:,3) = (k-1)*dz; 
end

Lx = dx*idimsm;
Ly = dy*jdimsm;
Lz = dz*kdimsm;
TT = dt*nfiles;

kx = zeros(kdimsm,idimsm,jdimsm);
ky = zeros(kdimsm,idimsm,jdimsm);
kz = zeros(kdimsm,idimsm,jdimsm);
omega = zeros(nfiles);

for i=1:idimsm
   kx(:,i,:) = (i-1)*2*pi/Lx; 
end
for j=1:jdimsm
   ky(:,:,j) = (j-1)*2*pi/Ly; 
end
for k=1:kdimsm
   kz(k,:,:) = (k-1)*2*pi/Lz; 
end

for n=1:nfiles
   omega(n) = (n-1)*2*pi/TT;
end

for n=1:nfiles
    t(n) = (n-1)*dt;
end

%read flow data
fsol.gname = '/vol';
for n=1:nfiles
    fsol.fname=sprintf(strcat(FlowFileLoc,'timeseriesVol_','%08u.h5'),fnum);
    disp(['reading file:', fsol.fname]);
       
    for nn=1:1
        fsol.dname=Dname{nn};
        buffer(:,:,:,nn,n)=ReadHDF5(fsol);
    end
    
    fnum=fnum+file_skip;
end

buffer(1,1,1,1,1)


%% part 1
buffer_3D      = zeros(kdimsm,idimsm,jdimsm);
buffer_3D_plot = zeros(kdimsm,idimsm,jdimsm);
for k=1:kdimsm
for i=1:idimsm
    for j=1:jdimsm
        buffer_3D(k,i,j) = buffer(k,i,j,1,1);
    end
end
end
buffer_3D_plot = buffer_3D; % for debug

% mean value
mean = 0;
for k=1:kdimsm
for i=1:idimsm
    for j=1:jdimsm
      mean = mean + buffer_3D(k,i,j);
    end
end
end
mean = mean/(kdimsm*idimsm*jdimsm)
buffer_3D = buffer_3D - mean;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D real-to-complex transform
buffer_comp = fftn(buffer_3D(:,:,:));

buffer_3D_debug = ifftn(buffer_comp);
buffer_3D_debug = buffer_3D_debug + mean;
% normalization
buffer_comp = buffer_comp/(kdimsm*idimsm*jdimsm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
 






 

 






 