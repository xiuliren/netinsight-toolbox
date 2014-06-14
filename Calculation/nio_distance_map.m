% nio_distance_map
% Calculate a distance map from soma to microvasculature within a
% single block at different cortex depth.
% 
% dtm = nio_distance_map( d, bm, zv )
% ----------------------------------
% 
% Get average distance from soma to the nearest microvessel within a
% single block at different cortex depth by sliding window.
% 
% Input
% -----
% - d:: n*4 array, n equals to number of somas.the first three column 
%       represent x-coordinate¡¢y-coordinate and z-coordinate of a single 
%       soma, respectively. The last column represents distance from the .
%       soma to the nearest microvessel. unit: ¦Ìm
% - bm:: a binary image(m*n matrix whose elements are 0 or 1)which 
%        marks location of a single block.
% - zv:: a set of depth(z-coordinate) for calculation
% 
% Output
% ------
% - dtm:: n*2 array, n equals to length of zv.The first column represents
%         average distance, the latter column represents number of somas in
%         particular window.
% 
% Example
% -------
% dtm{1} = nio_distance_map( distance, barrel_mask, zv );
% 
% See also nio_soma2tree

function [ dtm ] = nio_distance_map( d, bm, zv )
winsize = 25;
N = size(d, 1);
area = length(find(bm == 1));
if area == 0 % empty block
    dtm = zeros(length(zv), 2);
    return;
end
dtm = zeros(length(zv), 2);
d_bm = zeros(N, 1);
for i = 1 : N
    d_bm(i) = bm(ceil(d(i, 1)), ceil(d(i, 2)));
end
del_idx = d_bm == 0;
d(del_idx, :) = [];
floor = 0;
for z = zv
    floor = floor + 1;
    n1 = find(d(:, 3) >= z - winsize);
    n2 = find(d(:, 3) < z + winsize);
    n = intersect(n1, n2);
    if ~isempty(n)
        dtm(floor, 1) = mean(d(n, 4));
    else
        dtm(floor, 1) = 0;
    end
    dtm(floor, 2) = length(n);
end
end

