%% construct the local graph model
% by jpwu, 2013/08/29
function [node_list mk_local_stk] = local_graph_cut_sphere( seed, local_stk, mk_local_stk, T )
% parameters

%% get the circles representing the sphere

circles = get_sphere_circles( seed, (size(local_stk,1)-1)/2 );

%% establish the n and t link
% build the n link graph
[~, map_id, PAIRWISE]= build_n_link_graph_V4_sphere...
    ( circles, local_stk, mk_local_stk );
if isempty(map_id)
    node_list = [];
    return;
end

% build the t link graph
[UNARY rays]= build_t_link_graph( local_stk, seed, map_id, T );

%% get the optimized graph cut
[~,LABELS] = maxflow(PAIRWISE, sparse(double( UNARY' )));

%% connectivity analysis 
[node_list]= get_node_list( LABELS, PAIRWISE, map_id, seed );

%% mark the visited voxels
mk_local_stk = mark_visited_voxels_V3( mk_local_stk, seed, uint16(node_list) );
