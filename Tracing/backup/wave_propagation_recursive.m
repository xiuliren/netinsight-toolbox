%% wave propagation from a single seed
% by jpwu, 2013/02/28

function network = wave_propagation( seed, network )
global stk;
% global mk_stk;
%% parameters
% the size scale of local cube
P = 2;

%% propagation processing
% section number of current network
sn = length( network.sections );
% add the original seed in a new section
sn = sn + 1;
network.sections{sn} = seed;

% dynamic seeds list
dsl = seed;
section_idx = sn;
% iterative tracing
while ~isempty(dsl)
    % get the last seed and current section number
    seed = dsl( end,: );
    c_sn = section_idx(end);
    % delet from the dynamic list
    dsl(end,:) = [];
    section_idx(end) = [];
    
    % get the local voxel cube
    [local_stk, m1, n1, k1]= get_local_stk(stk, seed, P);
    % binarize the local voxel cube
    [mk_local_stk, T ] = kmeans_binarize( local_stk );
    [Ml, Nl, Kl] = size( mk_local_stk );
    
    %%  wave propagation in the local cube
    local_wave = nio_new_wave(seed - [m1, n1, k1, 0]);
    [wave_front, mk_local_stk] = local_propagation( local_wave, local_stk, mk_local_stk, m1, n1, k1 );
    
    %% connectivity analysis of the wave front
    disp('------ connectivity analysis of the wave front--------')
    node_list = wave_front_analysis( wave_front, mk_local_stk, m1, n1, k1 );
    Nn = size( node_list, 1);
    if Nn == 1
        network.sections{c_sn} = [network.sections{c_sn}; node_list(1,:)];
        dsl = [dsl; node_list(1,:) ];
        section_idx = [section_idx; sn];
    elseif Nn > 1
        % have branches
        dsl = [dsl; node_list];
        sn = length(network.sections);
        section_idx = [section_idx; ((1:Nn)+sn)'];
        for nn = 1 : Nn
            network.sections{sn+nn} = node_list(nn,:);
        end
    else
        % no new nodes
        % nothing to do
    end
end
return;

%% connectivity analysis of wave front
function node_list = wave_front_analysis( wave_front, mk_local_stk, m1, n1, k1 )
[Ml Nl Kl] = size(mk_local_stk);
nb = 0;
node_list = [];
for m = 1 : Ml
    for n = 1 : Nl
        for k = 1 : Kl
            if mk_local_stk(m,n,k)  == 2
                % a new branch
                nb = nb + 1;
                wave = nio_new_wave();
                % mark all the interconnected voxels
                [wave, mk_local_stk] = find_neighbor_front( wave, mk_local_stk, [m,n,k] );
                % compute a new node
                node = wave2node(wave) + [m1, n1, k1];
                node_list = [ node_list; node ]
            end
        end
    end
end
return;

%% compute node by a wave
function node = wave2node( wave )
% compute centroid 
centroid = mean( wave.voxelList );

% estimate the radius by distance from centroid to farthest voxel
dis = sum( abs( node - wave ) ,2);
% find the largest distance
radius = max(dis);
node = [centroid, radius];
return;

%% find connected neighboring wave front voxels
function [wave, mk_local_stk] = find_neighbor_front( wave, mk_local_stk, seed )
% mark the seed
mk_local_stk(seed) = 1;
wave.voxelList = [ wave.voxelList; seed ];

% estimate the searching range
[Ml Nl Kl] = size( mk_local_stk );
m1 = max(1, seed(1) - 1);   m2 = min( Ml, seed(1) + 1);
n1 = max(1, seed(2) - 1);   n2 = min( Nl, seed(2) + 1);
k1 = max(1, seed(3) - 1);   k2 = min( Kl, seed(3) + 1);

% find the neighboring voxels
for m = m1 : m2
    for n = n1 : n2
        for k = k1 : k2
            if mk_local_stk(m,n,k) ==2
                % find a new neighboring voxel
                wave.voxelList = [wave.voxelList; m,n,k ];
                % recursive growing
                [wave, mk_local_stk] = find_neighbor_front( wave, mk_local_stk, [m,n,k] );
            end
        end
    end
end
return;

%% get the local cube
function [local_stk, m1, n1, k1] = get_local_stk( stk, centroid, P )
[M, N, K] = size(stk);
% the coordinate and radius
ms = centroid(1);   ns = centroid(2);
ks = centroid(3);   rs = centroid(4);

m1 = max( 1, ms-P*rs ); m2 = min(M, ms+P*rs);
n1 = max( 1, ns-P*rs ); n2 = min(M, ns+P*rs);
k1 = max( 1, ks-P*rs ); k2 = min(M, ks+P*rs);

local_stk = stk(m1:m2, n1:n2, k1:k2);
return;