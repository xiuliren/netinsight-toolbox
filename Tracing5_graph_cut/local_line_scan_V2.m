%% local tracing using a stick matching
% by jpwu, 2013/10/17
function [node_list mk_local_stk] = local_line_scan_V2( seed, local_stk, mk_local_stk, T )
% the computed radius should be > RMIN
RMIN = 1;
% the noise level amplifier, default = 255
noiseLevelAmplifier = 20;

% eliminate the traced voxels
% local_stk = local_stk .* uint8( mk_local_stk );

%% get the circles representing the sphere
R = (size(local_stk,1)-1)/2;
circles = get_sphere_circles( seed, R );

%% get the rays% build the n link graph
% build the n link graph
[~, map_id, PAIRWISE]= build_n_link_graph_V4_sphere...
    ( circles, local_stk, mk_local_stk );
if isempty(map_id)
    node_list = [];
    return;
end
[~, rays]= build_t_link_graph( local_stk, seed, map_id, T );

%% find the maximum match of rays
vcorr = zeros( length(rays),1 );
for ridx = 1 : length( rays )
    vray = rays{ ridx };
    if ~isempty( vray )
        vcorr( ridx ) = mean( vray(uint16(end/2):end,4) );
    end
end

% choose the maximum ray
[C, Imax] = max( vcorr );

% estimate the local noise level
[M, N, K] = size( local_stk );
noiseLevel = estimate_noise( reshape(local_stk, M, N*K) ) * noiseLevelAmplifier;
% noiseLevel = median( NoiseLevel( double(local_stk) ));
% noiseLevel = std( local_stk(:) );
disp(['noise level: ' num2str( noiseLevel )])

if C < mean( local_stk(:) ) + noiseLevel
    node_list = [];
    return;
end

% compute the new node coordinate
vray = rays{ Imax };
[~, idxmax] =  max( vray(uint16(end/2):end,4) );
center = vray( uint16(idxmax)+end/2-1, 1:3 );

% compute the radius
radius = get_radius_V3(local_stk, center, T);
radius = max( radius, RMIN );

node_list = [ center, radius ];

mk_local_stk = mark_visited_voxels_cube(  mk_local_stk, seed, node_list );
