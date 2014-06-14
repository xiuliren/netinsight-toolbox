% build the n link graph
% by jpwu, 2013/09/21
function [A] = build_n_link_graph_V3_cube( map_co, map_id, local_stk )
% parameters
delta = 5;

% stk size
[M N K] = size(local_stk);

% initiate the graph model
% initialize the adjacency matrix
Nv = size( map_id, 1 );
A = sparse( Nv, Nv );

% build the lookup table of the boundary energy, weight of n link
wn_tab = exp( -1 * [0:255].*[0:255] / (2*delta*delta) );

% get all the neighboring voxels, and establish the connection
% pad the map coordinate array to avoid array bounds
map_co = padarray( map_co, [M+1, N+1, K+1], 'post' );

for idx0 = 1 : Nv
    m0 = map_id(idx0,1);  n0 = map_id(idx0,2);  
    k0 = map_id(idx0,3);  i0 = double(map_id(idx0,4));
    
    % the neighboring voxels
    idx1 = map_co( m0+1, n0, k0 );
    if idx1 > 0
        i1 = double(local_stk( m0+1, n0, k0 ));
        A( idx0, idx1 ) = wn_tab( abs(i0-i1) + 1 );
        A( idx1, idx0 ) = A( idx0, idx1 );
    end
    
    idx2 = map_co( m0, n0+1, k0 );
    if idx2 > 0
        i2 = double(local_stk( m0, n0+1, k0 ));
        A( idx0, idx2 ) = wn_tab( abs(i0-i2) + 1 );
        A( idx2, idx0 ) = A( idx0, idx2 );
    end
    
    idx3 = map_co( m0, n0, k0+1 );
    if idx3 > 0
        i3 = double(local_stk( m0, n0, k0+1 ));
        A( idx0, idx3 ) = wn_tab( abs(i0-i3) + 1 );
        A( idx3, idx0 ) = A( idx0, idx3 );
    end
end