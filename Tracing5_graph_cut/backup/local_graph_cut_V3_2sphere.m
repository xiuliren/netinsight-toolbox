%% construct the local graph model
% by jpwu, 2013/08/29
function [node_list mk_local_stk] = local_graph_cut_V3_2sphere( seed, local_stk, mk_local_stk )
% parameters
[M N K] = size( local_stk );
%% get the circles representing the sphere
R = (M-1)/2;
circles1 = get_sphere_circles( seed, R );
% a smaller sphere
circles2 = get_sphere_circles( seed, R-1 );

%% establish the n and t link
% build the n link graph
[map_co, map_id, PAIRWISE]= build_n_link_graph_V3_2sphere...
    ( circles1, circles2, local_stk, mk_local_stk );
if isempty(map_id)
    node_list = [];
    return;
end

% build the t link graph
% binarize the local voxel cube
[T, h]= kmeans_binarize( local_stk );
[UNARY rays]= build_t_link_graph( local_stk, seed, map_id, T );

% inital label cunarylass according to threshold
Nv = size( map_id, 1 );
CLASS = zeros(1,Nv);
% for nv = 1 : Nv
%     CLASS(nv) = (T<local_stk( map_id(nv,1), map_id(nv,2), map_id(nv,3) ));
% end

% the label cost
LABELCOST = single( [0 1; 1 0] );

%% get the optimized graph cut
[LABELS ENERGY ENERGYAFTER] = GCMEX(CLASS, UNARY, PAIRWISE, LABELCOST, 0);
% the 0 represents the forground, have to reverse this representation
LABELS = 1 - LABELS;

%% connectivity analysis 
[node_list]= get_node_list( LABELS, PAIRWISE, map_id, seed );

%% mark the visited voxels
mk_local_stk = mark_visited_voxels_V3( mk_local_stk, seed, uint16(node_list) );

%% line scan for the failure of graph cut
% if isempty( node_list )
%     disp('using line scan!')
%     [node_list mk_local_stk] = local_line_scan(seed, rays, local_stk, mk_local_stk, T);
%     disp(['the seed: ' ]);seed
%     disp('the new node: '); node_list
% end

%% get node list from the label
function [node_list ] = get_node_list(LABELS, PAIRWISE, map_id, seed)
% minimum number of voxels of each component
Cmin = 3;   % Cmin >= 3, make sure that the radius index is >=1

% the forground index
fidx = find( LABELS );
Nf = length( fidx );
% build the foreground graph adjacency matrix
fA = sparse( Nf, Nf );
for nf = 1 : Nf
    fA(nf,:) = PAIRWISE( fidx(nf), fidx );
end
[Nb, C] = graphconncomp(fA, 'Directed', false );

node_list = [];

% mark the component in LABELS
LABELS( fidx ) = C;
% compute the center and radius of each node
for nb = 1 : Nb
    % component index
    cidx = find( LABELS==nb );
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
    % the position was modeled by a circute, Ndis = pi*r*r, 
    %   Nc = pi*r*r - pi*(r-1)*(r-1)
    ndis = ceil( ( 2*sqrt(length( dis2 )/pi) - pi )) + 1;
    r = sqrt( dis2( ndis ) );
    
    % adjust the coordinate according to the algorithm of voxel scooping
    ex = double(min( r, seed(4) ) / max( r, seed(4) ));
    center = uint16( double(seed(1:3)) + ...
        (0.5^ex)*( double(center) - double(seed(1:3)) ));
    
    % add the element
    node_list = [node_list; center r];
end

%% mark the visited voxels, only mark the final point
function mk_local_stk = mark_visited_voxels_V3( mk_local_stk, seed, P_list )
% the expansion rate
ER = 1.2;
% minimum and maximum spaning range
MINSR = 2;  MAXSR = 7;
P_list = uint16(P_list);
[M,N,K] = size(mk_local_stk);

for nn = 1 : size(P_list, 1)   
    r2 = uint16( ER * P_list(nn, 4) );
    r2 = min(r2, P_list(nn,4)+MAXSR);   
    r2 = max(r2, P_list(nn,4)+MINSR);
    
    % fill the sphere of new node
    m0 = P_list(nn,1);  n0 = P_list(nn,2);
    k0 = P_list(nn,3);  r0 = r2;
    mk_local_stk = fill_mark_stack_with_sphere(mk_local_stk, m0,n0,k0,r0);
    
    % fill the sphere of center node 
    m0 = (seed(1) + P_list(nn,1)) / 2;
    n0 = (seed(2) + P_list(nn,2)) / 2;
    k0 = (seed(3) + P_list(nn,3)) / 2;
    r0 = (seed(4) + r2 ) / 2;
    mk_local_stk = fill_mark_stack_with_sphere(mk_local_stk, m0,n0,k0,r0);
end

%% fill the mark stack with sphere
function mk_stk = fill_mark_stack_with_sphere(mk_stk, m0,n0,k0,r0)
% stack size
[M N K] = size(mk_stk);
% create sphere as template
S = hypersphere(double([r0*2+1, r0*2+1, r0*2+1]),'full');

% initialize the coordinate
m1 = m0-r0; m2 = m0+r0;
n1 = n0-r0; n2 = n0+r0;
k1 = k0-r0; k2 = k0+r0;

ms1 = 1;    ms2 = 2*r0+1;
ns1 = 1;    ns2 = 2*r0+1;
ks1 = 1;    ks2 = 2*r0+1;

% the boundary condition
if m0-r0 < 1
    m1 = 1;
    ms1 = r0-m0+2;
end
if m0+r0 > M
    m2 = M;
%         ms2 = r0 + 1 -m0 + M;
    ms2 = m2-m1+ms1;
end

if n0-r0 < 1
    n1 = 1;
    ns1 = r0-n0+2;
end
if n0+r0 > N
    n2 = N;
%         ns2 = r0 + 1 -n0 + N;
    ns2 = n2-n1+ns1;
end

if k0-r0 < 1
    k1 = 1;
    ks1 = r0-k0+2;
end
if k0+r0 > K
    k2 = K;
%         ks2 = r0 + 1 -k0 + K;
    ks2 = k2-k1+ks1;
end

mk_stk( m1:m2, n1:n2, k1:k2 ) = mk_stk( m1:m2, n1:n2, k1:k2 )...
    | S( ms1:ms2, ns1:ns2, ks1:ks2 );
return;