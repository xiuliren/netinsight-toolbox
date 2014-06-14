%% construct the local graph model
% by jpwu, 2013/08/29
function [node_list, Nb] = local_graph_cut( seed, local_stk, T, h )
% parameters
% the local standard deviation, reflecting the noise level
delta = 4;

% stack size
M = size( local_stk, 1 );
N = size( local_stk, 2 );
K = size( local_stk, 3 );

% the binary mask 
mk_local_stk = ( local_stk > T );

%% initiate the graph model
% number of voxels
Nv = M*N*K - (M-2)*(N-2)*(K-2);

% get the voxels coordinate
vol = double( local_stk );
vol(2:M-1, 2:N-1, 2:K-1) = zeros(M-2, N-2, K-2) - 1;
indices = find( vol>=0 );
indices_cen = find( vol<0 );

% the coordinate matrix
[mc nc kc] = ind2sub([M N K], indices); 

% initialize the adjacency matrix
A = sparse( M*N*K, M*N*K );

%% establish the n link
% build the lookup table of the boundary energy
b_tab = exp( -1 * [0:255].*[0:255] / (2*delta*delta) );

% search all the nodes
for ind = 1 : Nv
    % get the coordinates
    m = mc(ind);
    n = nc(ind);
    k = kc(ind);
    
    if m > 1
        ind1 = sub2ind([M N K], m-1, n, k);
        A(ind, ind1) = b_tab( abs(local_stk(m,n,k) - local_stk(m-1,n,k)) + 1 );
        A(ind1, ind) = A(ind, ind1);
    end
    
    if m < M
        ind1 = sub2ind([M N K], m+1, n, k);
        A(ind, ind1) = b_tab( abs(local_stk(m,n,k) - local_stk(m+1,n,k)) + 1 );
        A(ind1, ind) = A(ind, ind1);
    end
    
    if n > 1
        ind1 = sub2ind([M N K], m, n-1, k);
        A(ind, ind1) = b_tab( abs(local_stk(m,n,k) - local_stk(m,n-1,k)) + 1 );
        A(ind1, ind) = A(ind, ind1);
    end
    
    if n < N
        ind1 = sub2ind([M N K], m, n+1, k);
        A(ind, ind1) = b_tab( abs(local_stk(m,n,k) - local_stk(m,n+1,k)) + 1 );
        A(ind1, ind) = A(ind, ind1);
    end
    
    if k > 1
        ind1 = sub2ind([M N K], m, n, k-1);
        A(ind, ind1) = b_tab( abs(local_stk(m,n,k) - local_stk(m,n,k-1)) + 1 );
        A(ind1, ind) = A(ind, ind1);
    end
    
    if k < K
        ind1 = sub2ind([M N K], m, n, k+1);
        A(ind, ind1) = b_tab( abs(local_stk(m,n,k) - local_stk(m,n,k+1)) + 1 );
        A(ind1, ind) = A(ind, ind1);
    end
end

% % compress the adjacency matrix, delete the central voxels
% % this step is also time comsuming !!!
% A(indices_cen, :) = [];
% A(:, indices_cen) = [];

%% establish the t link, add the source and sink nodes
% establish the weight lookup table of forground t link,
% cf = (255+T)/2; cb = T/2;
alpha = (255+2*double(T))/4;
beta = -255/4/log(9);   % constant!
w_lut = 1./(1+exp( ([0:255]-alpha)./beta ) );

% find the path vector between the seed and the boundary voxels
for ind = 1 : Nv
    % get the coordinaKenichi Ebinates
    m = mc(ind);
    n = nc(ind);
    k = kc(ind);
    
%     % get the interprated points, consider the voxel vector
%     [Ml,Nl,Kl] = bresenham_line3d(seed(1:3), [m n k]);
%     % get the intensity vector
%     idx = sub2ind([M N K], Ml, Nl, Kl );
%     iv = local_stk( idx );
    
    % weigt of forground t link
    wft = w_lut( local_stk( m,n,k ) + 1 );
    % weight of background t link
    wbt = 1 - wft;
    
    % add to the adjacency matrix
    A(Nv+1, ind) = wft;
    A(ind, Nv+1) = wft;
    A(Nv+2, ind) = wbt;
    A(ind, Nv+2) = wbt;
end

% %% normalized cut
% % Shi and Malik Algorithm: second smallest eigen vector
% disp('Finding Eigen Vector');
% d = sum(A,2); D = diag(d); tic
% B = (D-A); B = (B+B')/2; OPTS.disp = 0;
% [v,d,flag] = eigs(B,D,2,'SA',OPTS); ttime = toc;
% disp(sprintf('Time for finding eigen vector = %f',ttime)); clear OPTS
% y = v(:,2);
% Ncut = reshape(y,M,N,K);
% Ncut = Ncut > 0;

%%
k = 2;
D = diag(sum(A));
L = D-A;

opt = struct('issym', true, 'isreal', true);
[V dummy] = eigs(L, D, k, 'SM', opt);
idx = kmeans(V, k);