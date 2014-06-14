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

% initialize the adjacency matrix
A = sparse( Nv, Nv );

%% establish the n link
% build the lookup table of the boundary energy
b_tab = exp( -1 * [0:255].*[0:255] / (2*delta*delta) );

% traverse the outer surface of the cube
% the inner surface
for m = 2 : M-1
    for n = 2 : N-1
        for k = [1 K]
            ind0 = sub2ind([M N K], m,n,k);
            
            % the surrounding voxels, 4 connection relashionship
            ind1 = sub2ind([M N K], m-1,n,k);
            A(ind0, ind1) = b_tab( abs( local_stk(ind0)-local_stk(ind1) )+1 );
            A(ind1, ind0) = A(ind0, ind1);
            
            ind2 = sub2ind([M N K], m+1,n,k);
            A(ind0, ind2) = b_tab( abs( local_stk(ind0)-local_stk(ind2) )+1 );
            A(ind2, ind0) = A(ind0, ind2);
            
            ind3 = sub2ind([M N K], m,n-1,k);
            A(ind0, ind3) = b_tab( abs( local_stk(ind0)-local_stk(ind3) )+1 );
            A(ind3, ind0) = A(ind0, ind3);
            
            ind4 = sub2ind([M N K], m,n+1,k);
            A(ind0, ind4) = b_tab( abs( local_stk(ind0)-local_stk(ind4) )+1 );
            A(ind4, ind0) = A(ind0, ind4);
        end
    end
end

for m = 2 : M-1
    for n = [1 N]
        for k = 2 : K-1
            ind0 = sub2ind([M N K], m,n,k);
            
            % the surrounding voxels, 4 connection relashionship
            ind1 = sub2ind([M N K], m-1,n,k);
            A(ind0, ind1) = b_tab( abs( local_stk(ind0)-local_stk(ind1) )+1 );
            A(ind1, ind0) = A(ind0, ind1);
            
            ind2 = sub2ind([M N K], m+1,n,k);
            A(ind0, ind2) = b_tab( abs( local_stk(ind0)-local_stk(ind2) )+1 );
            A(ind2, ind0) = A(ind0, ind2);
            
            ind3 = sub2ind([M N K], m,n,k-1);
            A(ind0, ind3) = b_tab( abs( local_stk(ind0)-local_stk(ind3) )+1 );
            A(ind3, ind0) = A(ind0, ind3);
            
            ind4 = sub2ind([M N K], m,n,k+1);
            A(ind0, ind4) = b_tab( abs( local_stk(ind0)-local_stk(ind4) )+1 );
            A(ind4, ind0) = A(ind0, ind4);
        end
    end
end

for m = [1 M]
    for n = 2 : N-1
        for k = 2 : K-1
            ind0 = sub2ind([M N K], m,n,k);
            
            % the surrounding voxels, 4 connection relashionship
            ind1 = sub2ind([M N K], m,n-1,k);
            A(ind0, ind1) = b_tab( abs( local_stk(ind0)-local_stk(ind1) )+1 );
            A(ind1, ind0) = A(ind0, ind1);
            
            ind2 = sub2ind([M N K], m,n+1,k);
            A(ind0, ind2) = b_tab( abs( local_stk(ind0)-local_stk(ind2) )+1 );
            A(ind2, ind0) = A(ind0, ind2);
            
            ind3 = sub2ind([M N K], m,n,k-1);
            A(ind0, ind3) = b_tab( abs( local_stk(ind0)-local_stk(ind3) )+1 );
            A(ind3, ind0) = A(ind0, ind3);
            
            ind4 = sub2ind([M N K], m,n,k+1);
            A(ind0, ind4) = b_tab( abs( local_stk(ind0)-local_stk(ind4) )+1 );
            A(ind4, ind0) = A(ind0, ind4);
        end
    end
end

% the inner boundary
n = 1; k = 1;
for m = 2 : M-1
    ind0 = sub2ind([M N K], m,n,k);
    
    % the surrounding voxels, 4 connection relashionship
    ind1 = sub2ind([M N K], m-1,n,k);
    A(ind0, ind1) = b_tab( abs( local_stk(ind0)-local_stk(ind1) )+1 );
    A(ind1, ind0) = A(ind0, ind1);
            
    ind2 = sub2ind([M N K], m+1,n,k);
    A(ind0, ind2) = b_tab( abs( local_stk(ind0)-local_stk(ind2) )+1 );
    A(ind2, ind0) = A(ind0, ind2);
            
    ind3 = sub2ind([M N K], m,n+1,k);
    A(ind0, ind3) = b_tab( abs( local_stk(ind0)-local_stk(ind3) )+1 );
    A(ind3, ind0) = A(ind0, ind3);
    
    ind4 = sub2ind([M N K], m,n,k+1);
    A(ind0, ind4) = b_tab( abs( local_stk(ind0)-local_stk(ind4) )+1 );
    A(ind4, ind0) = A(ind0, ind4);
end

n = 1; k = K;
for m = 2 : M-1
    ind0 = sub2ind([M N K], m,n,k);
    
    % the surrounding voxels, 4 connection relashionship
    ind1 = sub2ind([M N K], m-1,n,k);
    A(ind0, ind1) = b_tab( abs( local_stk(ind0)-local_stk(ind1) )+1 );
    A(ind1, ind0) = A(ind0, ind1);
            
    ind2 = sub2ind([M N K], m+1,n,k);
    A(ind0, ind2) = b_tab( abs( local_stk(ind0)-local_stk(ind2) )+1 );
    A(ind2, ind0) = A(ind0, ind2);
            
    ind3 = sub2ind([M N K], m,n+1,k);
    A(ind0, ind3) = b_tab( abs( local_stk(ind0)-local_stk(ind3) )+1 );
    A(ind3, ind0) = A(ind0, ind3);
    
    ind4 = sub2ind([M N K], m,n,k-1);
    A(ind0, ind4) = b_tab( abs( local_stk(ind0)-local_stk(ind4) )+1 );
    A(ind4, ind0) = A(ind0, ind4);
end

n = N; k = 1;
for m = 2 : M-1
    ind0 = sub2ind([M N K], m,n,k);
    
    % the surrounding voxels, 4 connection relashionship
    ind1 = sub2ind([M N K], m-1,n,k);
    A(ind0, ind1) = b_tab( abs( local_stk(ind0)-local_stk(ind1) )+1 );
    A(ind1, ind0) = A(ind0, ind1);
            
    ind2 = sub2ind([M N K], m+1,n,k);
    A(ind0, ind2) = b_tab( abs( local_stk(ind0)-local_stk(ind2) )+1 );
    A(ind2, ind0) = A(ind0, ind2);
            
    ind3 = sub2ind([M N K], m,n-1,k);
    A(ind0, ind3) = b_tab( abs( local_stk(ind0)-local_stk(ind3) )+1 );
    A(ind3, ind0) = A(ind0, ind3);
    
    ind4 = sub2ind([M N K], m,n,k+1);
    A(ind0, ind4) = b_tab( abs( local_stk(ind0)-local_stk(ind4) )+1 );
    A(ind4, ind0) = A(ind0, ind4);
end

n = N; k = K;
for m = 2 : M-1
    ind0 = sub2ind([M N K], m,n,k);
    
    % the surrounding voxels, 4 connection relashionship
    ind1 = sub2ind([M N K], m-1,n,k);
    A(ind0, ind1) = b_tab( abs( local_stk(ind0)-local_stk(ind1) )+1 );
    A(ind1, ind0) = A(ind0, ind1);
            
    ind2 = sub2ind([M N K], m+1,n,k);
    A(ind0, ind2) = b_tab( abs( local_stk(ind0)-local_stk(ind2) )+1 );
    A(ind2, ind0) = A(ind0, ind2);
            
    ind3 = sub2ind([M N K], m,n-1,k);
    A(ind0, ind3) = b_tab( abs( local_stk(ind0)-local_stk(ind3) )+1 );
    A(ind3, ind0) = A(ind0, ind3);
    
    ind4 = sub2ind([M N K], m,n,k-1);
    A(ind0, ind4) = b_tab( abs( local_stk(ind0)-local_stk(ind4) )+1 );
    A(ind4, ind0) = A(ind0, ind4);
end

% n
m = 1; k = 1;
for n = 2 : N-1
    ind0 = sub2ind([M N K], m,n,k);
    
    % the surrounding voxels, 4 connection relashionship
    ind1 = sub2ind([M N K], m,n-1,k);
    A(ind0, ind1) = b_tab( abs( local_stk(ind0)-local_stk(ind1) )+1 );
    A(ind1, ind0) = A(ind0, ind1);
            
    ind2 = sub2ind([M N K], m,n+1,k);
    A(ind0, ind2) = b_tab( abs( local_stk(ind0)-local_stk(ind2) )+1 );
    A(ind2, ind0) = A(ind0, ind2);
            
    ind3 = sub2ind([M N K], m+1,n,k);
    A(ind0, ind3) = b_tab( abs( local_stk(ind0)-local_stk(ind3) )+1 );
    A(ind3, ind0) = A(ind0, ind3);
    
    ind4 = sub2ind([M N K], m,n,k+1);
    A(ind0, ind4) = b_tab( abs( local_stk(ind0)-local_stk(ind4) )+1 );
    A(ind4, ind0) = A(ind0, ind4);
end

m = M; k = 1;
for n = 2 : N-1
    ind0 = sub2ind([M N K], m,n,k);
    
    % the surrounding voxels, 4 connection relashionship
    ind1 = sub2ind([M N K], m,n-1,k);
    A(ind0, ind1) = b_tab( abs( local_stk(ind0)-local_stk(ind1) )+1 );
    A(ind1, ind0) = A(ind0, ind1);
            
    ind2 = sub2ind([M N K], m,n+1,k);
    A(ind0, ind2) = b_tab( abs( local_stk(ind0)-local_stk(ind2) )+1 );
    A(ind2, ind0) = A(ind0, ind2);
            
    ind3 = sub2ind([M N K], m-1,n,k);
    A(ind0, ind3) = b_tab( abs( local_stk(ind0)-local_stk(ind3) )+1 );
    A(ind3, ind0) = A(ind0, ind3);
    
    ind4 = sub2ind([M N K], m,n,k+1);
    A(ind0, ind4) = b_tab( abs( local_stk(ind0)-local_stk(ind4) )+1 );
    A(ind4, ind0) = A(ind0, ind4);
end

m = 1; k = K;
for n = 2 : N-1
    ind0 = sub2ind([M N K], m,n,k);
    
    % the surrounding voxels, 4 connection relashionship
    ind1 = sub2ind([M N K], m,n-1,k);
    A(ind0, ind1) = b_tab( abs( local_stk(ind0)-local_stk(ind1) )+1 );
    A(ind1, ind0) = A(ind0, ind1);
            
    ind2 = sub2ind([M N K], m,n+1,k);
    A(ind0, ind2) = b_tab( abs( local_stk(ind0)-local_stk(ind2) )+1 );
    A(ind2, ind0) = A(ind0, ind2);
            
    ind3 = sub2ind([M N K], m+1,n,k);
    A(ind0, ind3) = b_tab( abs( local_stk(ind0)-local_stk(ind3) )+1 );
    A(ind3, ind0) = A(ind0, ind3);
    
    ind4 = sub2ind([M N K], m,n,k-1);
    A(ind0, ind4) = b_tab( abs( local_stk(ind0)-local_stk(ind4) )+1 );
    A(ind4, ind0) = A(ind0, ind4);
end

m = M; k = K;
for n = 2 : N-1
    ind0 = sub2ind([M N K], m,n,k);
    
    % the surrounding voxels, 4 connection relashionship
    ind1 = sub2ind([M N K], m,n-1,k);
    A(ind0, ind1) = b_tab( abs( local_stk(ind0)-local_stk(ind1) )+1 );
    A(ind1, ind0) = A(ind0, ind1);
            
    ind2 = sub2ind([M N K], m,n+1,k);
    A(ind0, ind2) = b_tab( abs( local_stk(ind0)-local_stk(ind2) )+1 );
    A(ind2, ind0) = A(ind0, ind2);
            
    ind3 = sub2ind([M N K], m-1,n,k);
    A(ind0, ind3) = b_tab( abs( local_stk(ind0)-local_stk(ind3) )+1 );
    A(ind3, ind0) = A(ind0, ind3);
    
    ind4 = sub2ind([M N K], m,n,k-1);
    A(ind0, ind4) = b_tab( abs( local_stk(ind0)-local_stk(ind4) )+1 );
    A(ind4, ind0) = A(ind0, ind4);
end

% k
m = 1; n = 1;
for k = 2 : K-1
    ind0 = sub2ind([M N K], m,n,k);
    
    % the surrounding voxels, 4 connection relashionship
    ind1 = sub2ind([M N K], m,n,k-1);
    A(ind0, ind1) = b_tab( abs( local_stk(ind0)-local_stk(ind1) )+1 );
    A(ind1, ind0) = A(ind0, ind1);
            
    ind2 = sub2ind([M N K], m,n,k+1);
    A(ind0, ind2) = b_tab( abs( local_stk(ind0)-local_stk(ind2) )+1 );
    A(ind2, ind0) = A(ind0, ind2);
            
    ind3 = sub2ind([M N K], m+1,n,k);
    A(ind0, ind3) = b_tab( abs( local_stk(ind0)-local_stk(ind3) )+1 );
    A(ind3, ind0) = A(ind0, ind3);
    
    ind4 = sub2ind([M N K], m,n+1,k);
    A(ind0, ind4) = b_tab( abs( local_stk(ind0)-local_stk(ind4) )+1 );
    A(ind4, ind0) = A(ind0, ind4);
end

m = M; n = 1;
for k = 2 : K-1
    ind0 = sub2ind([M N K], m,n,k);
    
    % the surrounding voxels, 4 connection relashionship
    ind1 = sub2ind([M N K], m,n,k-1);
    A(ind0, ind1) = b_tab( abs( local_stk(ind0)-local_stk(ind1) )+1 );
    A(ind1, ind0) = A(ind0, ind1);
            
    ind2 = sub2ind([M N K], m,n,k+1);
    A(ind0, ind2) = b_tab( abs( local_stk(ind0)-local_stk(ind2) )+1 );
    A(ind2, ind0) = A(ind0, ind2);
            
    ind3 = sub2ind([M N K], m-1,n,k);
    A(ind0, ind3) = b_tab( abs( local_stk(ind0)-local_stk(ind3) )+1 );
    A(ind3, ind0) = A(ind0, ind3);
    
    ind4 = sub2ind([M N K], m,n+1,k);
    A(ind0, ind4) = b_tab( abs( local_stk(ind0)-local_stk(ind4) )+1 );
    A(ind4, ind0) = A(ind0, ind4);
end

m = 1; n = N;
for k = 2 : K-1
    ind0 = sub2ind([M N K], m,n,k);
    
    % the surrounding voxels, 4 connection relashionship
    ind1 = sub2ind([M N K], m,n,k-1);
    A(ind0, ind1) = b_tab( abs( local_stk(ind0)-local_stk(ind1) )+1 );
    A(ind1, ind0) = A(ind0, ind1);
            
    ind2 = sub2ind([M N K], m,n,k+1);
    A(ind0, ind2) = b_tab( abs( local_stk(ind0)-local_stk(ind2) )+1 );
    A(ind2, ind0) = A(ind0, ind2);
            
    ind3 = sub2ind([M N K], m+1,n,k);
    A(ind0, ind3) = b_tab( abs( local_stk(ind0)-local_stk(ind3) )+1 );
    A(ind3, ind0) = A(ind0, ind3);
    
    ind4 = sub2ind([M N K], m,n-1,k);
    A(ind0, ind4) = b_tab( abs( local_stk(ind0)-local_stk(ind4) )+1 );
    A(ind4, ind0) = A(ind0, ind4);
end

m = M; n = N;
for k = 2 : K-1
    ind0 = sub2ind([M N K], m,n,k);
    
    % the surrounding voxels, 4 connection relashionship
    ind1 = sub2ind([M N K], m,n,k-1);
    A(ind0, ind1) = b_tab( abs( local_stk(ind0)-local_stk(ind1) )+1 );
    A(ind1, ind0) = A(ind0, ind1);
            
    ind2 = sub2ind([M N K], m,n,k+1);
    A(ind0, ind2) = b_tab( abs( local_stk(ind0)-local_stk(ind2) )+1 );
    A(ind2, ind0) = A(ind0, ind2);
            
    ind3 = sub2ind([M N K], m-1,n,k);
    A(ind0, ind3) = b_tab( abs( local_stk(ind0)-local_stk(ind3) )+1 );
    A(ind3, ind0) = A(ind0, ind3);
    
    ind4 = sub2ind([M N K], m,n-1,k);
    A(ind0, ind4) = b_tab( abs( local_stk(ind0)-local_stk(ind4) )+1 );
    A(ind4, ind0) = A(ind0, ind4);
end

% the corners