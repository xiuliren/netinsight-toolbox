% nio_density_map_V2
% Drawing somas and vessels density map of barrel cortex around L4.
% 
% [ dm_s, dm_fav, dm_fmv, dm_lav, dm_lmv ] = nio_density_map_V2( vessel, T_D, soma, WS, d_L4, t_L4, Ms, Ns )
% ----------------------------------------------------------------------------------------------------------
% 
% Calculate density by nio_fv_ld_cuboid.
% 
% Input
% -----
% - vessel:: the vectorized vessels which has tree structure, one-root
% - T_D:: the threshold of diameter, the vessel segments whose radius is
%        bigger than T_D/2 will be removed, micron
% - soma:: the somas which has tree structure but no real connections
% - WS:: half of the window size in XY plane, micron
% - d_L4, t_L4:: depth and thickness of Layer 4, micron
% - Ms, Ns:: the size of density map
% 
% Output
% ------
% - dm_s:: density map for cell number density, unit: 10^5/mm^3
% - dm_fav:: density map for fractional vascular volume
% - dm_fmv:: density map for fractional microvascular volume
% - dm_lav:: density map for normalized vascular length density, unit:
%            m/mm^3
% - dm_lmv:: density map for normalized microvascular length density, unit:
%            m/mm^3
% 
% Attention
% ---------
% The coordinate of Matlab and Amira are both used.
% X Y is not exchanged, and M--Y, N--X, K--Z.
% It can be applied to other layers.
% 
% Example
% -------
% [ dm_s, dm_fav, dm_fmv, dm_lav, dm_lmv ] = nio_density_map_V2( vessels, 4, soma, 25, 400, 
% 100, 560, 600 );
% 
% Uses nio_fv_ld_cuboid

function [ dm_s, dm_fav, dm_fmv, dm_lav, dm_lmv ] = nio_density_map_V2( vessel, T_D, soma, WS, d_L4, t_L4, Ms, Ns )
%% parameters
vessel_c = nio_extract_microvessels(vessel, T_D);
%% the soma density map
disp('------ getting the soma density map ...');
dm_s = zeros( Ms, Ns, 'double' );

% get the soma density map
k1 = d_L4 - t_L4;    k2 = d_L4 + t_L4;
for m = ( WS + 1 ) : ( Ms - WS )
    m1 = m - WS;    m2 = m + WS;
    for n = ( WS + 1 ) : ( Ns - WS )
        % the block boundary
        n1 = n - WS;    n2 = n + WS;
        
        idx = find ( ( soma.X > n1 ) & ( soma.X <= n2 ) & ...
            ( soma.Y > m1 ) & ( soma.Y <= m2 ) & ...
            ( soma.Z > k1 ) & ( soma.Z <= k2 ) );
        
        % count the somas
        dm_s(m, n) = length( idx );
    end
end

% get the density, 10e5/mm^3
dm_s = dm_s .* 10000 / (m2 - m1 ) / (n2 -n1 ) / ( k2 - k1 );

%% density map of vessels
disp('------ getting density map of vessels... ')
% fractional volume and length density map of all vessels
dm_fav = zeros( Ms, Ns, 'double' );
dm_lav = zeros( Ms, Ns, 'double' );
% fractional volume and length density map of microvessels
dm_fmv = zeros( Ms, Ns, 'double' );
dm_lmv = zeros( Ms, Ns, 'double' );

% get the density map
% the Z range
k1 = d_L4 - t_L4;    k2 = d_L4 + t_L4;
for m = ( WS + 1 ) : ( Ms - WS )
    % the block boundary
    m1 = m - WS;    m2 = m + WS;   
    for n = ( WS + 1 ) : ( Ns - WS )
        % the block boundary
        n1 = n - WS;    n2 = n + WS;
        
        % the nodes inside current block
        idx_v = find ( ( vessel.X > n1 ) & ( vessel.X <= n2 ) & ...
            ( vessel.Y > m1 ) & ( vessel.Y <= m2 ) & ...
            ( vessel.Z > k1 ) & ( vessel.Z <= k2 ) );
        idx_m = find ( ( vessel_c.X > n1 ) & ( vessel_c.X <= n2 ) & ...
            ( vessel_c.Y > m1 ) & ( vessel_c.Y <= m2 ) & ...
            ( vessel_c.Z > k1 ) & ( vessel_c.Z <= k2 ) );
        
        % the density of all vessels
        [ dm_fav( m, n ), dm_lav( m, n ) ] = nio_fv_ld_cuboid ( vessel, idx_v, m1, m2, n1, n2, k1, k2  );
        % the density of microvessels
        [ dm_fmv( m, n ), dm_lmv( m, n ) ] = nio_fv_ld_cuboid ( vessel_c, idx_m, m1, m2, n1, n2, k1, k2  );
    end
end

%% show density map
disp('------ plot the density map...')

dm_s_p = dm_s( (WS+1):( Ms-WS ), (WS+1):(Ns-WS) );
dm_s_p = imadjust( uint8( dm_s_p/max(dm_s_p(:))*255) );
figure, imshow( dm_s_p ),     colormap jet  colorbar
title('Soma Density Map')

dm_fav_p = dm_fav( (WS+1):( Ms-WS ), (WS+1):(Ns-WS) );
dm_fav_p = imadjust( uint8( dm_fav_p/max(dm_fav_p(:)) *255 ) );
figure, imshow( dm_fav_p ),     colormap jet  colorbar
title('Fractional Volume Map of All Vessels')

dm_lav_p = dm_lav( (WS+1):( Ms-WS ), (WS+1):(Ns-WS) );
dm_lav_p = imadjust( uint8( dm_lav_p/max(dm_lav_p(:))*255) );
figure, imshow( dm_lav_p ),     colormap jet  colorbar
title('Length Density Map of All Vessels')

dm_fmv_p = dm_fmv( (WS+1):( Ms-WS ), (WS+1):(Ns-WS) );
dm_fmv_p = imadjust( uint8( dm_fmv_p/max(dm_fmv_p(:))*255) );
figure, imshow( dm_fmv_p ),     colormap jet  colorbar
title('Fractional Volume Map of Microvessels')

dm_lmv_p = dm_lmv( (WS+1):( Ms-WS ), (WS+1):(Ns-WS) );
dm_lmv_p = imadjust( uint8( dm_lmv_p/max(dm_lmv_p(:))*255) );
figure, imshow( dm_lmv_p ),     colormap jet  colorbar
title('Length Density Map of Microvessels')

disp('------ end ------')