% nio_coordinate_transform
% Transform coordinates and convert to 1*1*1 micron.
% 
% [ trees ] = nio_coordinate_transform( trees, options, voxel_size_xyz )
% ----------------------------------------------------------------------
% 
% Transform coordinates because the coordinate system of original swc file 
% and Matlab doesn't match. Convert the data into 1*1*1 micron by providing
% voxel size in all coordinates.
% 
% Input
% -----
% - trees:: cell array, the tree structure 
% - options:: string, inform the two transformed coordinates.{DEFAULT: 'xy'}
%        'xy' : transform x coordinate and y coordinate
%        'yz' : transform y coordinate and z coordinate
%        'xz' : transform x coordinate and z coordinate
% - voxel_size_xyz:: 1*3 linear array, the voxel size of original swc file,
%                    its elements represent voxel size of x, y and 
%                    z coordinate, respectively. 
% 
% Output
% ------
% - trees:: cell array, the transformed tree structure
% 
% Example
% -------
% sample_tree = nio_coordinate_transform( sample_tree, 'xz', [2 2 2] );

function [ trees ] = nio_coordinate_transform( trees, options, varargin )
num_tree = length(trees);

if ~isempty(options)&&~ischar(options)
    varargin{1} = options;
    options = [];
end

if(isempty(options))
    options = 'xy';
end

if(isempty(varargin))
    voxel_size_xyz = [1 1 1];
else
    voxel_size_xyz = varargin{1};
end

for ward = 1 : num_tree
    trees{ward}.X = trees{ward}.X * voxel_size_xyz(1);
    trees{ward}.Y = trees{ward}.Y * voxel_size_xyz(2);
    trees{ward}.Z = trees{ward}.Z * voxel_size_xyz(3);
end

switch options,
    case 'xy'
        for ward = 1 : num_tree
            tmp = trees{ward}.X;
            trees{ward}.X = trees{ward}.Y;
            trees{ward}.Y = tmp;
        end
    case 'yz'
        for ward = 1 : num_tree
            tmp = trees{ward}.Y;
            trees{ward}.Y = trees{ward}.Z;
            trees{ward}.Z = tmp;
        end
    case 'xz'
        for ward = 1 : num_tree
            tmp = trees{ward}.X;
            trees{ward}.X = trees{ward}.Z;
            trees{ward}.Z = tmp;
        end
end
end

