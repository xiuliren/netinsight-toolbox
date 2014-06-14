% nio_density_map
% Drawing somas and vessels density map of barrel cortex around L4.
% 
% [ dm_s, dm_fav, dm_fmv ] = nio_density_map( vessels, vessels_c, soma, WS, d_L4, t_L4, Ms, Ns, Ks )
% --------------------------------------------------------------------------------------------------
% 
% Calculate density by build vessel-filled or soma-filled binary stacks.
% 
% Input
% -----
% - vessels:: the vectorized vessels which has tree structure, multi-roots
% - vessels_c:: the vectorized microvessels which has tree structure,
%               multi-roots
% - soma:: the somas which has tree structure but no real connections
% - WS:: half of a cube size, micron
% - d_L4, t_L4:: depth and thickness of Layer 4, micron
% - Ms, Ns, Ks:: the size of image stack
% 
% Output
% ------
% - dm_s:: density map for cell number density, unit: 10^5/mm^3
% - dm_fav:: density map for fractional vascular volume
% - dm_fmv:: density map for fractional microvascular volume
% 
% Attention
% ---------
% It can be applied to other layers.
% 
% Example
% -------
% [ dm_s, dm_fav, dm_fmv ] = nio_density_map( vessels, vessels_c, soma, 25, 400, 
% 100, 560, 600, 880 );
% 
% Uses nio_swc2stk

function [ dm_s, dm_fav, dm_fmv ] = nio_density_map( vessels, vessels_c, soma, WS, d_L4, t_L4, Ms, Ns, Ks )
%% the soma density map
dm_s = zeros( Ms, Ns, 'double' );


% construct soma stack
soma_stk = zeros( Ms*Ns, Ks, 'uint8' );

for n = 1 : length( soma.X )
    soma_stk( ceil( soma.X(n) ), ceil( soma.Y(n) ), ceil( soma.Z(n) ) ) = 1;
end
k1 = d_L4 - t_L4;    k2 = d_L4 + t_L4;
% get the soma density map
for m = ( WS + 1 ) : ( Ms - WS )
    for n = ( WS + 1 ) : ( Ns - WS )
            block = soma_stk( (m - WS): (m + WS), ( n - WS ) : ( n + WS ), k1 : k2 );
            dm_s( m, n ) = sum( block(:) );
    end
end

%% density map of vessels
dm_fav = zeros( Ms, Ns, 'double' );
dm_fmv = zeros( Ms, Ns, 'double' );

% construct stk 
v_stk = nio_swc2stk( vessels , Ms, Ns, Ks);
microv_stk = nio_swc2stk( vessels_c, Ms, Ns, Ks );
k1 = d_L4 - t_L4;    k2 = d_L4 + t_L4;
% transform to binary
v_stk = logical(v_stk);
microv_stk = logical( microv_stk );

% get the vessel density map
for m = ( WS + 1 ) : ( Ms - WS )
    for n = ( WS + 1 ) : ( Ns - WS )
        block = v_stk( (m - WS): (m + WS), ( n - WS ) : ( n + WS ), k1 : k2 );
        dm_fav( m, n ) = sum( block(:) );
            
        block = microv_stk( (m - WS): (m + WS), ( n - WS ) : ( n + WS ), k1 : k2 );
        dm_fmv( m, n ) = sum( block(:) );
    end
end

%% show density map
dm_s_p = dm_s( (WS+1):( Ms-WS ), (WS+1):(Ns-WS) );
dm_s_p = imadjust( uint8( dm_s_p/max(dm_s_p(:))*255) );
figure, imshow( dm_s_p ),     colormap jet  colorbar
title('Soma Density Map')

dm_fav_p = dm_fav( (WS+1):( Ms-WS ), (WS+1):(Ns-WS) );
dm_fav_p = imadjust( uint8( dm_fav_p/max(dm_fav_p(:))*255) );
figure, imshow( dm_fav_p ),     colormap jet  colorbar
title('Fractional Volume Map of All Vessels')

dm_fmv_p = dm_fmv( (WS+1):( Ms-WS ), (WS+1):(Ns-WS) );
dm_fmv_p = imadjust( uint8( dm_fmv_p/max(dm_fmv_p(:))*255) );
figure, imshow( dm_fmv_p ),     colormap jet  colorbar
title('Fractional Volume Map of Microvessels')
