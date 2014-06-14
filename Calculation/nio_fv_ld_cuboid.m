% nio_fv_ld_cube
% Get the fractional vascular volume and normalized length density in a
% cuboid.
% 
% [ fvv, nld ] = nio_fv_ld_cuboid ( vessel, idx_b, m1, m2, n1, n2, k1, k2 )
% -----------------------------------------------------------------------
% 
% Get the fractional vascular volume(fvv) and normalized length density(nld)
% in a cuboid. Split the tree structure by add new parent or child segments when 
% hits the cuboid boundary.
% 
% Input
% -----
% - vessel:: the vessel tree, the Amira coordinate system
% - idx_b:: linear array, the nodes index inside the cuboid
% - m1,m2,n1,n2,k1,k2:: define the cuboid boundary, 1<2, the matlab coordinate system
% 
% Output
% ------
% - fvv:: the fractional vascular volume
% - nld:: the normalized length density, unit: m/mm^3
% 
% Attention
% ---------
% 1. the voxel size of cuboid must be 1X1X1 micron for correct density
% 2. the coordinate system: Y--N, X--M, Z--K
% 
% Example
% -------
% [ fvv, nld ] = nio_fv_ld_cuboid ( vessel, 1:50, 50, 100, 250, 300, 75, 175
% );
% 
% See also nio_intersection_stk nio_density_map
% Uses idpar_tree

function [ fvv, nld ] = nio_fv_ld_cuboid ( vessel, idx_b, m1, m2, n1, n2, k1, k2 )
%% initialize the parameters
% fractional vascular volume
fvv = 0;
% normalized length density
nld = 0;

idpar = idpar_tree( vessel, '-0' );
%% the density of all vessels
for ward = 1 : length( idx_b )
    % the node index and parent, children
    idx = idx_b( ward );
    x = vessel.X( idx );   y = vessel.Y( idx );   z = vessel.Z( idx );   D1 = vessel.D(idx);
%    idx_p = find( vessels.dA(idx, :) );
    idx_p = idpar( idx );
    idx_c_list = find( vessel.dA(:,idx) );
    
    %% add the parent connection volume and length
    if idx_p ~= 0
        % there is a parent node
        px = vessel.X( idx_p );    py = vessel.Y( idx_p );    pz = vessel.Z( idx_p );    D2 = vessel.D( idx_p );
        if ( py<m1 ) || ( py>m2 ) || ( px<n1 ) || ( px>n2 ) || ( pz<k1 ) || ( pz>k2 )
            % out of boundary
            len0 = sqrt( (px-x)*(px-x) + (py-y)*(py-y) + (pz-z)*(pz-z) );
            [ px py pz ] = adjust_coordinate(  x, y, z, px, py, pz, m1, m2, n1, n2, k1, k2 );
            % distance
            len = sqrt( (px-x)*(px-x) + (py-y)*(py-y) + (pz-z)*(pz-z) );
            
            % the interprated diameter
%             D2 = D1 + len0/len*( vessels.D( idx_p ) - D1 );
            D2 = ( len*D2 + D1*len0 ) / ( len + len0 ) ;
        else
            % inside the cube
            len = sqrt( (px-x)*(px-x) + (py-y)*(py-y) + (pz-z)*(pz-z) );
        end
        
        % frustum-based volume
        vol = (pi*(D1*D1 + D1*D2 + D2*D2 )*len) / 12;
        fvv = fvv + vol;
        nld = nld + len;
    end
    
    %% add children segments which are outside cube
    % skip the inside cube children to avoid repeated addition
    for nc = 1 : length( idx_c_list )
        idx_c = idx_c_list( nc );
        cx = vessel.X( idx_c );    cy = vessel.Y( idx_c );    cz = vessel.Z( idx_c );    
        if ( cy<m1 ) || ( cy>m2 ) || ( cx<n1 ) || ( cx>n2 ) || ( cz<k1 ) || ( cz>k2 )
            % out of boundary
            D2 = vessel.D( idx_c );
            len0 = sqrt( (cx-x)*(cx-x) + (cy-y)*(cy-y) + (cz-z)*(cz-z) );
            [ cx cy cz ] = adjust_coordinate(  x, y, z, cx, cy, cz, m1, m2, n1, n2, k1, k2 );
            % distance
            len = sqrt( (cx-x)*(cx-x) + (cy-y)*(cy-y) + (cz-z)*(cz-z) );
            
            % the interprated diameter
%             D2 = D1 + len0/len*( D2 - D1 );
            D2 = ( len*D2 + D1*len0 ) / ( len + len0 ) ;
            % frustum-based volume
            vol = (pi*(D1*D1 + D1*D2 + D2*D2 )*len) / 12;
            fvv = fvv + vol;
            nld = nld + len;
        end
    end
end

% divid the cube volume to get density
V_cube = ( m2 - m1 ) * ( n2 - n1 ) * ( k2 - k1 );
fvv = fvv / V_cube;
nld = nld * 1000 / V_cube;
return;

%% interprete the nodes coordinate to locate inside the cube
function [ x2 y2 z2 ] = adjust_coordinate( x, y, z, x2, y2, z2, m1, m2, n1, n2, k1, k2 )
% interprate node by unit vector
L = sqrt( (x2-x)*(x2-x) + (y2-y)*(y2-y) + (z2-z)*(z2-z) );
xuv = ( x2 - x ) / L;   yuv = ( y2 - y ) / L;   zuv = ( z2 - z ) / L;

% initialze the walking point
xt = x; yt = y; zt = z;
for len = 1 : L
    xt = xt + xuv;  yt = yt + yuv;  zt = zt + zuv;
    if ( xt < n1 ) || ( xt > n2 ) || ( yt < m1 ) || ( yt > m2 ) || ( zt < k1 ) || ( zt > k2 )
        x2 = xt - xuv;  y2 = yt - yuv;  z2 = zt - zuv;
        return;
    end
end
return;