% nio_intersection
% Insert a new node when vasculature hit the boundary of a block in a one-root tree
% structure.
% 
% [ x_is y_is z_is d_is ] = nio_intersection( vessel, bm, n1, n2 )
% ------------------------------------------------------------------
% 
% Insert a new node when vasculature cross the boundary of a block by
% linear interpolation in a one-root tree structure at xy space.
% It can be helpful to collect vasculature segments within a single block
% which defined in 2-dimension.
% 
% Input
% -----
% - vessel:: the vectorized vessels which has a tree structure
% - bm:: a binary 2-dimension matrix, marking a single block
% - n1:: index of one of terminal points of a vessel segment which cross the boundary of a block
% - n2:: index of another terminal point of a vessel segment which cross the boundary of a block
% 
% Output
% ------
% - x_is:: x-coordinate of the new node
% - y_is:: y-coordinate of the new node
% - y_is:: z-coordinate of the new node
% - d_is:: diameter of the new node. unit: ¦Ìm
% 
% Example
% -------
% [ x_is y_is z_is d_is ] = nio_intersection( sample_tree, barrel_mask, 1, 2 )
% 
% See also nio_intersection_stk nio_vessel_single_block

function [ x_is y_is z_is d_is ] = nio_intersection( vessel, bm, n1, n2 )
%QA_INTERSECTION Summary of this function goes here
%   Detailed explanation goes here
x1 = vessel.X(n1);         y1 = vessel.Y(n1);       z1 = vessel.Z(n1);
d1 = vessel.D(n1);
x2 = vessel.X(n2);         y2 = vessel.Y(n2);       z2 = vessel.Z(n2);
d2 = vessel.D(n2);
xx = abs(x1 - x2);
yy = abs(y1 - y2);
l = cos(atan(yy/xx)); % step size for x coordinate
if(x1 <= x2)
    x = x1 : l : x2;
else
    x = x1 : -l : x2;
end
if xx == 0
    x_is = x1;
    y = min(y1, y2) : max(y1, y2); % interpolation directly at y coordinate
    if y(length(y)) ~= max(y1, y2)
        y = [y max(y1, y2)];
    end
    for n = 2 : length(y)
        if bm(ceil(x_is), ceil(y(n-1))) ~= bm(ceil(x_is), ceil(y(n))) % hit the boundary 
            if bm(ceil(x_is), ceil(y(n-1)))
                y_is = y(n-1);
            else
                y_is = y(n);
            end
            break;
        end
    end
else
    y = interp1([x1 x2], [y1 y2], x, 'linear', 'extrap'); % interpolation at y coordinate through step increase at x
    nodes = [[x' y']; [x2 y2]];
    for n = 2 : length(nodes)
        if bm(ceil(nodes(n-1, 1)), ceil(nodes(n-1, 2))) ~= bm(ceil(nodes(n, 1)), ceil(nodes(n, 2))) % hit the boundary 
            if bm(ceil(nodes(n-1, 1)), ceil(nodes(n-1, 2)))
                x_is = nodes(n-1, 1);
                y_is = nodes(n-1, 2);
            else
                x_is = nodes(n, 1);
                y_is = nodes(n, 2);
            end
            break;
        end
    end
end
d_is1 = sqrt((x_is - x1)^2 + (y_is - y1)^2); % distance from node1
d_is2 = sqrt((x_is - x2)^2 + (y_is - y2)^2); % distance from node2
z_is = z1 + (z2 - z1)*(d_is1 / (d_is1 + d_is2)); % detemine the z coordinate according to distance
d_is = d1 + (d2 - d1)*(d_is1 / (d_is1 + d_is2)); % detemine the diameter according to distance
end

