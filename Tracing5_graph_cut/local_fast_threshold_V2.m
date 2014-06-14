%% construct the local graph model
% by jpwu, 2013/08/29
function [node_list mk_local_stk] = local_fast_threshold_V2( seed, local_stk, mk_local_stk, T )
% parameters
[M, N, K] = size(local_stk);
%% binarize the image
% threshold the local stack
BW_local_stk = reshape(im2bw(local_stk(:), double(T)/255), [M,N,K]);
BW_local_stk = ~(~BW_local_stk | mk_local_stk);

%% get the circles representing the sphere
circles = get_sphere_circles( seed, (size(local_stk,1)-1)/2 );

%% establish the n and t link
% build the n link graph
[map_co, map_id, PAIRWISE]= build_n_link_graph_V4_sphere...
    ( circles, uint8(BW_local_stk), ~BW_local_stk );
if isempty(map_id)
    node_list = [];
    return;
end

%% connectivity analysis 
[node_list P_list]= get_node_list( PAIRWISE, map_id, seed);

%% mark the visited voxels
mk_local_stk = mark_visited_voxels_V3( mk_local_stk, seed, uint16(node_list) );

% %% build the index and coordinate map
% function circles = get_sphere_circles( seed, local_stk )
% % initialization
% seed = uint16(seed);
% [M,N,K] = size(local_stk);
% % the integer radius of sphere
% R = (M-1)/2;
% circles = cell( M,1 );
% 
% % build the circles in the 
% n0 = seed(2);   k0 = seed(3);
% for h = double(R : -1 : -R)
%     m0 = seed(1) - h; 
% %     I(:,:) = local_stk(m0,:,:);
%     % the radius of the circle
%     r = uint16( sqrt( double(R*R - h*h) ) );
%     % the circular points
%     [nc,kc] = getmidpointcircle(n0,k0,r);
%     % remove the replicate nodes, which is a bug in the function
%     nc2 = nc;   nc2(1:end-1) = nc(2:end);   nc2(end) = nc(1); 
%     kc2 = kc;   kc2(1:end-1) = kc(2:end);   kc2(end) = kc(1);
%     label = (nc == nc2) & (kc==kc2);
%     nc(label) = []; kc(label) = [];
%     
%     circles{R-h+1} = [ repmat(m0,length(nc),1), nc, kc];
% end

%% get node list from the label
function [node_list P_list] = get_node_list( PAIRWISE, map_id, seed )
% minimum number of voxels of each component
Cmin = 3;   % Cmin >= 3, make sure that the radius index is >=1

% connectivity analysis
[Nb, C] = graphconncomp(PAIRWISE, 'Directed', false );

% initialization
node_list = [];
P_list = [];
% compute the center and radius of each node
for nb = 1 : Nb
    % component index
    cidx = find( C==nb );
    
%     disp(['number of voxels: ' num2str( length(cidx) ) ])
    
    if length(cidx) <= Cmin
        % less than Cmin voxels
        continue;
    end
    
    % vector of the coordinates
    vcoord = map_id(cidx, 1:3);
    center = mean( vcoord, 1 );
    
    % compute the distance from center to every component voxel
    vcoord = double(vcoord);
    dis2 = (vcoord(:,1) - center(1)).*(vcoord(:,1) - center(1)) + ...
        (vcoord(:,2) - center(2)).*(vcoord(:,2) - center(2)) + ...
        (vcoord(:,3) - center(3)).*(vcoord(:,3) - center(3));
    % the radius was estimated using the largest distance
    dis2 = sort(dis2, 'descend');
    Ndis = length( dis2 );
    % the position was modeled by a circute, Ndis = pi*r*r, 
    %   Nc = pi*r*r - pi*(r-1)*(r-1)
    ndis = ceil( ( 2*sqrt(Ndis/pi) - pi )/2) + 1;
    r = sqrt( dis2( ndis ) );

    % add the point 
    P_list = [P_list; center r];
    
    % adjust the coordinate according to the algorithm of voxel scooping
    ex = min( r, seed(4) ) / max( r, seed(4) );
    center = double(seed(1:3)) + (0.5^double(ex))*( double(center(1:3))...
        - double(seed(1:3)) );
    
    % add the element
    node_list = [node_list; center r];
end
