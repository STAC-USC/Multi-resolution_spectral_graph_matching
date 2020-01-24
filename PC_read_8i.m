function [ V,C,J ] = PC_read_8i( filename )
% The following funtion reads .ply files and transforms the information
% into three different variables, where V represent the x-y-z coordinates
% This function is used for the following datasets:
% https://jpeg.org/plenodb/pc/8ilabs/
fid=fopen(filename,'r','n', 'UTF-8');
fscanf(fid,'ply\n');
fscanf(fid,'format ascii 1.0\n');
fscanf(fid,'comment Version 2, Copyright 2017, 8i Labs, Inc.\n');
world_scale=fscanf(fid,'comment frame_to_world_scale %g\n',1);
world_translation=fscanf(fid,'comment frame_to_world_translation %g %g %g\n',3);
width=fscanf(fid,'comment width %d\n',1);
N=fscanf(fid,'element vertex %d\n',1);

fscanf(fid,'property float x\n');
fscanf(fid,'property float y\n');
fscanf(fid,'property float z\n');
fscanf(fid,'property uchar red\n');
fscanf(fid,'property uchar green\n');
fscanf(fid,'property uchar blue\n');

s=fgetl(fid);
data=fscanf(fid,'%g %g %g %u %u %u\n',[6,N]);
V=data(1:3,:)';
C=data(4:6,:)';
J=log2(width+1);
end

