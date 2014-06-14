%% construct the local graph model
% by jpwu, 2013/08/29
function [node_list mk_local_stk] = local_graph_cut_cube( seed, local_stk, mk_local_stk, T)
% parameters
% stack size
[M N K] = size( local_stk );

% enhance the local stack
% local_stk = reshape( imadjust(local_stk(:)), M,N,K);
% local_stk2 = fastflux3(single(local_stk), 1:double((seed(4)*2+1)), 1);
% local_stk2 = uint8(local_stk*(255/max(local_stk(:))));
% local_stk = local_stk/2 + local_stk2/2;
% eliminate the visited voxels
local_stk = local_stk .* uint8(~mk_local_stk);

% build the index and coordinate map
[map_co map_id] = build_id_co_map( M, N, K, local_stk, mk_local_stk );
if isempty( map_id )
    % all local voxels were visited
    node_list = [];
    return;
end
%% establish the n and t link
% build the n link graph
PAIRWISE = build_n_link_graph_cube_V2( map_co, size(map_id, 1), local_stk );

% build the t link graph
UNARY = build_t_link_graph( local_stk, seed, map_id, T );

%% get the optimized graph cut
[flow,LABELS] = maxflow(PAIRWISE, sparse(double( UNARY' )));

%% connectivity analysis 
[node_list]= get_node_list( LABELS, PAIRWISE, map_id, seed );

%% mark the visited voxels
mk_local_stk = mark_visited_voxels_cube( mk_local_stk, seed, node_list );

%% build the index and coordinate map
function [map_co map_id] = build_id_co_map( M, N, K, local_stk, mk_local_stk )
% number of voxels
Nv = M*N*K - (M-2)*(N-2)*(K-2);
% initialize
map_co = zeros(M, N, K);
map_id = zeros( Nv, 4, 'uint16' );
% 6 surfaces
% the upper and lower surface
sid = 0;
for m = [1 M]
    for n = 1 : N
%         for k = 1 : K
%             idx = sub2ind( [N K], n, k );
%             map_id( idx+sid, : ) = [m, n, k, local_stk(m,n,k) ];
%         end
        idx_list = sub2ind([N K], ones(1,K).*n, 1:K);
        map_id( idx_list + sid, :) = [ones(K,1).*m, ones(K,1).*n, (1:K)', ...
            local_stk( sub2ind(size(local_stk), ones(K,1).*m, ones(K,1).*n,(1:K)') )];
    end
    
    sid = sid + N*K;
end

% the left and right surface, without the upper and lower boundary
for n = [1 N]
    for m = 2 : M-1
%         for k = 1 : K
%             idx = sub2ind( [M-2, K], m-1, k );
%             map_id( idx+sid, : ) = [m, n, k, local_stk(m,n,k) ];
%         end
        idx_list = sub2ind([M-2 K], ones(1,K).*m-1, 1:K);
        map_id( idx_list + sid, :) = [ones(K,1).*m, ones(K,1).*n, (1:K)', ...
            local_stk( sub2ind(size(local_stk), ones(K,1).*m, ones(K,1).*n,(1:K)') )];
    end
    sid = sid + (M-2)*K;
end

% the front and back surface, without the upper and lower,left and right boundary,
for k = [1 K]
    for m = 2 : M-1
%         for n = 2 : N-1
%             idx = sub2ind( [M-2, N-2], m-1, n-1 );
%             map_id( idx+sid, : ) = [m, n, k, local_stk(m,n,k) ];
%         end
        idx_list = sub2ind([M-2 N-2], ones(1,N-2).*m-1, (2:N-1)-1);
        map_id( idx_list + sid, :) = [ones(N-2,1).*m, ones(N-2,1).*n, (1:N-2)', ...
            local_stk( sub2ind(size(local_stk), ones(N-2,1).*m, (2:N-1)',ones(N-2,1).*k) )];
    end
    sid = sid + (M-2)*(N-2);
end

% eliminate the visited voxels
% label the visited index
vil = false( Nv, 1 );
for nv = 1 : Nv
    if mk_local_stk(map_id(nv,1),map_id(nv,2),map_id(nv,3))
        vil(nv) = true;
    end
end
% the new map
map_id = map_id( ~vil,: );
if isempty( map_id )
    % all the voxels were visited, stop tracing
    return;
end
Nv = size(map_id,1);

% the map of coordinate index
for nv = 1 : Nv
    map_co(map_id(nv,1),map_id(nv,2),map_id(nv,3)) = nv;
end
