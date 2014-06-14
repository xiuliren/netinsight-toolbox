% nio_intersection_stk
% Insert a new node when vasculature hit the boundary of a block in a one-root tree
% structure.
% 
% [ x_is y_is z_is d_is ] = nio_intersection_stk( vessel, stack, n1, n2 )
% ------------------------------------------------------------------
% 
% Insert a new node when vasculature cross the boundary of a block by unit vector
% in a one-root tree structure at 3-dimension stack space.
% It can be helpful to collect vasculature segments within a single block
% which defined in 3-dimension.
% 
% Input
% -----
% - vessel:: the vectorized vessels which has a tree structure
% - stack:: a binary 3-dimension array, marking a single block
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
% [ x_is y_is z_is d_is ] = nio_intersection_stk( sample_tree, stack, 1, 2 )
% 
% See also nio_intersection

function [ x_is y_is z_is d_is ] = nio_intersection_stk( vessel, stack, n1, n2 )
x1 = vessel.X(n1);         y1 = vessel.Y(n1);       z1 = vessel.Z(n1);
D1 = vessel.D(n1);
x2 = vessel.X(n2);         y2 = vessel.Y(n2);       z2 = vessel.Z(n2);
D2 = vessel.D(n2);
dis = sqrt( (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2) );
% unit vector
uv = [ (x2-x1) (y2-y1) (z2-z1) ] ./ dis;
x0 = [];
y0 = [];
z0 = [];
R = [];
for k = 0 : dis
    % center
    x0 = [x0 x1 + k * uv(1)];
    y0 = [y0 y1 + k * uv(2)];
    z0 = [z0 z1 + k * uv(3)];
    % radius
    R = [R ( D1 + (D2-D1) * k / dis ) / 2];
end
if dis ~= fix(dis)
    x0 = [x0 x2];
    y0 = [y0 y2];
    z0 = [z0 z2];
    R = [R D2];
end
for n = 2 : length(R)
    if stack(ceil(x0(n-1)), ceil(y0(n-1)), ceil(z0(n-1))) ~= stack(ceil(x0(n)), ceil(y0(n)), ceil(z0(n)))
        if stack(ceil(x0(n-1)), ceil(y0(n-1)), ceil(z0(n-1)))
            x_is = x0(n-1);
            y_is = y0(n-1);
            z_is = z0(n-1);
            d_is = R(n-1);
        else
            x_is = x0(n);
            y_is = y0(n);
            z_is = z0(n);
            d_is = R(n);
        end
        break;
    end
end
end

