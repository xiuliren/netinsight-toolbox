%% construct the local graph model
% by jpwu, 2013/08/29
function [node_list mk_local_stk] = local_graph_cut_V3_cube( seed, local_stk, mk_local_stk, T, h )
% parameters

% build the index and coordinate map
[map_co map_id] = build_id_co_map_V2_cube(local_stk, mk_local_stk );
if isempty( map_id )
    % all local voxels were visited
    node_list = [];
    return;
end
%% establish the n and t link
% build the n link graph
PAIRWISE = build_n_link_graph_V3_cube( map_co, map_id, local_stk );

% build the t link graph
% binarize the local voxel cube
[T, h]= kmeans_binarize( local_stk );
UNARY = build_t_link_graph( local_stk, seed, map_id, T );

% inital label cunarylass
Nv = size( map_id, 1 );
CLASS = zeros(1,Nv);

% the label cost
LABELCOST = single( [0 1; 1 0] );

%% get the optimized graph cut
[LABELS ENERGY ENERGYAFTER] = GCMEX(CLASS, single(UNARY), PAIRWISE, LABELCOST, 0);
% the 0 represents the forground, have to reverse this representation
LABELS = 1 - LABELS;

%% connectivity analysis 
[node_list P_list]= get_node_list( local_stk, LABELS, PAIRWISE, map_id, seed, T );

%% mark the visited voxels
% mk_local_stk = mark_visited_voxels( M, N, K, seed, P_list );
mk_local_stk = mark_visited_voxels_V3( mk_local_stk, seed, P_list );

%% build the index and coordinate map
function [map_co map_id] = build_id_co_map_V2_cube( local_stk, mk_local_stk )
% number of voxels
[M, N, K] = size(local_stk);

% the index list of unvisited voxels
uvIdx = find( mk_local_stk );
Nv = length( uvIdx );
% initialize
map_co = zeros(M,N,K, 'double');
map_co( uvIdx ) = 1:Nv;
map_id = zeros( Nv, 4, 'uint16' );
for m = 1 : M
    for n = 1 : N
        for k = 1 : K
            if map_co(m,n,k) > 0
                idx = map_co(m,n,k);
                map_id( idx,: ) = [m,n,k, local_stk(m,n,k) ];
            end
        end
    end
end

%% build the t link 
function unary = build_t_link_graph( local_stk, seed, map_id, T )
% the noise level
delta = 5;

% the lower and upper threshold
T1 = T/2;
T2 = (double(T)+double(max(local_stk(:))))/2;
% establish the weight lookup table of forground t link,
% cf = (255+T)/2; cb = T/2;
alpha = (T2+2*double( T1 ))/4;
beta = -T2/4/log(9);   % constant!
lut = 1./(1+exp( ([0:255]-alpha)./beta ) );
lut_wf = 1-lut ;

lut_wb = lut;
% lut_wb = 1- lut_wf;

% initialize the unary matrix, CXN
Nv = size( map_id, 1 );
unary = zeros( 2, Nv );
unary(1,:) = lut_wf( map_id(:,4)+1 )';
unary(2,:) = lut_wb( map_id(:,4)+1 )';

% % the t-link weight computed by law of cosine
% % get the feature vector
% for idx = 1 : size( map_id, 1 )
%     % get the voxel coordinates using bresenham's algorithm
%     [vm vn vk] = bresenham_line3d( double(seed(1 : 3)), double(map_id(idx, 1:3)));
%     if size(vm,1) <= 2
%         continue;
%     end
%     
%     % get the voxel intensity
%     vi = local_stk( sub2ind(size(local_stk), vm, vn, vk ) );
% 
% %     % the inner differece
% %     d = double(max(diff(vi)));   
% %     % the boundary energy
% %     bE = exp( -d.*d./(2*delta*delta) ) ;
% %     unary(1,idx) = (unary(1,idx) + bE )/2;
%     
%     vt = ones( 1,length(vi) ) * mean( seed(1)-1:seed(1)+1, seed(2)-1:seed(2)+1, seed(3)-1:seed(3)+1 );
%     vi = double( vi( ceil(length(vi)/2) : length(vi) ) );
%     
%     %the Jaccard index
%     corr = dot(vt, vi) / ( dot(vt,vt) + dot(vi,vi) - dot(vt,vi));
%          
%     % cosine correlation
% %     if find(vi)   % if vi is not all zeros        
% %         corr = dot(vt, vi) / ( sqrt( dot(vt,vt) ) * sqrt( dot(vi,vi) ) );   
% %     end
% 
%     unary(1,idx) = (unary(1,idx) + corr)/2;
% end

%% get node list from the label
function [node_list P_list] = get_node_list( local_stk, LABELS, PAIRWISE, map_id, seed, T )
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
fA = (fA~=0);
[Nb, C] = graphconncomp(fA, 'Directed', false );

node_list = [];
P_list = [];
% mark the component in LABELS
LABELS( fidx ) = C;
% compute the center and radius of each node
for nb = 1 : Nb
    % component index
    cidx = find( LABELS==nb );
    
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
    % the position was modeled by a cube, Ndis = (2*r+1)^3, 
    %   Nc = 6*(2*)
    ndis = ceil( 3*(Ndis^(2/3))-0.5) + 1;
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
        vr = r1 : -1 : r2;
    else
        vr = r1 : 1 : r2;
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
    r2 = min(r2, P_list(nn,4)+MAXSR);   r2 = max(r2, P_list(nn,4)+MINSR);
    
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