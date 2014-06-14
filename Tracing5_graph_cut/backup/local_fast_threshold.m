%% construct the local graph model
% by jpwu, 2013/08/29
function [node_list mk_local_stk] = local_fast_threshold( seed, local_stk, mk_local_stk, T, h )
% parameters
% stack size
[M N K] = size( local_stk );

% enhance the local stack
% local_stk = reshape( imadjust(local_stk(:)), M,N,K);
% eliminate the visited voxels
local_stk = local_stk .* uint8(~mk_local_stk);

% threshold the local stack
BW_local_stk = reshape(im2bw(local_stk(:), T/255), [M,N,K]);
BW_local_stk = ~(~BW_local_stk | mk_local_stk);

% build the index and coordinate map
[map_co map_id] = build_id_co_map( M, N, K, local_stk, ~BW_local_stk );
if isempty( map_id )
    % all local voxels were visited
    node_list = [];
    return;
end

%% establish the n link
% build the n link graph
PAIRWISE = build_n_link_graph( map_co, size(map_id, 1), local_stk );

%% connectivity analysis 
[node_list P_list]= get_node_list( PAIRWISE, map_id, seed);

%% mark the visited voxels
mk_local_stk = mark_visited_voxels( M, N, K, seed, P_list );

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
        for k = 1 : K
            idx = sub2ind( [N K], n, k );
            map_id( idx+sid, : ) = [m, n, k, local_stk(m,n,k) ];
        end
    end
    sid = sid + N*K;
end

% the left and right surface, without the upper and lower boundary
for n = [1 N]
    for m = 2 : M-1
        for k = 1 : K
            idx = sub2ind( [M-2, K], m-1, k );
            map_id( idx+sid, : ) = [m, n, k, local_stk(m,n,k) ];
        end
    end
    sid = sid + (M-2)*K;
end

% the front and back surface, without the upper and lower,left and right boundary,
for k = [1 K]
    for m = 2 : M-1
        for n = 2 : N-1
            idx = sub2ind( [M-2, N-2], m-1, n-1 );
            map_id( idx+sid, : ) = [m, n, k, local_stk(m,n,k) ];
        end
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

%% get node list from the label
function [node_list P_list] = get_node_list( PAIRWISE, map_id, seed )
% minimum number of voxels of each component
Cmin = 3;   % Cmin >= 3, make sure that the radius index is >=1

[Nb, C] = graphconncomp(PAIRWISE, 'Directed', false );

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
%     % eliminate the out of boundary case
%     center = max([ center; 1 1 1 ]);
%     center = min([ center; M N K ]);
    
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
    ndis = ceil( ( 2*sqrt(Ndis/pi) - pi )) + 1;
    r = sqrt( dis2( ndis ) );

    % add the point 
    P_list = [P_list; center r];
    
    % adjust the coordinate according to the algorithm of voxel scooping
    ex = min( r, seed(4) ) / max( r, seed(4) );
    center = seed(1:3) + (0.5^ex)*( center(1:3) - seed(1:3) );
    
    % add the element
    node_list = [node_list; center r];
end

%% mark the visited voxels
function mk_local_stk = mark_visited_voxels( M, N, K, seed, P_list )
% the expansion rate
ER = 1.2;

% initialize the mark
mk_local_stk = false(M,N,K);
% evaluate the node list
if isempty( P_list )
    return;
end

for nn = 1 : size(P_list, 1)
    r1 = uint16( ER * seed(4) );   r2 = uint16( ER * P_list(nn, 4) );
    % vector of interprated 3D points using bresenham's algorithm
    [vm vn vk] = bresenham_line3d( seed(1:3), P_list(nn,1:3) );
    Ni = length(vm);
    if Ni==0
        return;
    end
    if r1 > r2
        vr = [ r1 : -1 : r2 ];
    else
        vr = [ r1 : 1 : r2 ];
    end
    Rv = imresize( vr, [1, Ni], 'nearest' );
    for ni = 1 : 2 : Ni
        m0 = vm(ni); n0 = vn(ni); k0 = vk(ni); r0 = Rv(ni);
        m1 = max(1, m0-r0);   m2 = min(M, m0+r0);
        n1 = max(1, n0-r0);   n2 = min(N, n0+r0);
        k1 = max(1, k0-r0);   k2 = min(K, k0+r0);
        mk_local_stk( m1:m2, n1:n2, k1:k2 ) = true( m2-m1+1, n2-n1+1, k2-k1+1 );
    end
end
return;