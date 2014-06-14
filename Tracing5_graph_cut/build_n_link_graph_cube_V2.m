%% build the n link graph of the cube model
% by jpwu, 2013/09/21
function [A] = build_n_link_graph_cube_V2( map_co, Nv, local_stk )
% parameters
delta = 5;

% stk size
[M, N, K] = size(local_stk);

% initiate the graph model
% initialize the adjacency matrix
% A = sparse( Nv, Nv );
rows = zeros(10 * Nv, 1);
cols = zeros(10 * Nv, 1);
weights = zeros(10 * Nv, 1);
rcNum = 0;

% build the lookup table of the boundary energy, weight of n link
wn_tab = exp( -1 * (0:255).*(0:255) / (2*delta*delta) );

% search all the nodes by searching 6 surfaces, 12 edges, and 8 corners.
% the 6 surfaces
% the upper and lower surface
for m = [1 M]
    for n = 2 : N
        for k = 2 : K
            idx = map_co( m, n, k );
            if idx == 0
                continue;
            end
            
            idx1 = map_co( m, n-1, k);
            if 0 ~= idx1 
                % weight of n link
                rcNum = rcNum + 1;
                weights(rcNum) = wn_tab( abs(double(local_stk(m,n,k))-...
                    double(local_stk(m,n-1,k)))+1 );
                rows(rcNum) = idx;  cols(rcNum) = idx1;
            end

            idx1 = map_co( m, n, k-1);
            if 0 ~= idx1
                % weight of n link
                rcNum = rcNum + 1;
                weights(rcNum) = wn_tab( abs(double(local_stk(m,n,k))-...
                    double(local_stk(m,n,k-1)))+1 );
                rows(rcNum) = idx;  cols(rcNum) = idx1;
            end
        end
    end
end

% the left and right surface
for n = [1 N]
    for m = 2 : M
        for k = 2 : K
            idx = map_co( m, n, k );
            if idx == 0
                continue;
            end
            
            idx1 = map_co( m-1, n, k);
            if 0 ~= idx1
                % weight of n link
                rcNum = rcNum + 1;
                weights(rcNum) = wn_tab( abs(double(local_stk(m,n,k))-...
                    double(local_stk(m-1,n,k)))+1 );
                rows(rcNum) = idx;  cols(rcNum) = idx1;
            end

            idx1 = map_co( m, n, k-1);
            if 0 ~= idx1
                % weight of n link
                rcNum = rcNum + 1;
                weights(rcNum) =  wn_tab( abs(double(local_stk(m,n,k))-...
                    double(local_stk(m,n,k-1)))+1 );
                rows(rcNum) = idx;  cols(rcNum) = idx1;
            end
        end
    end
end

% the front and back surface
for k = [1 K]
    for m = 2 : M
        for n = 2 : N
            idx = map_co( m, n, k );
            if idx == 0
                continue;
            end
            
            idx1 = map_co( m-1, n, k);
            if 0 ~= idx1
                % weight of n link
                rcNum = rcNum + 1;
                weights(rcNum) =  wn_tab( abs(double(local_stk(m,n,k))-...
                    double(local_stk(m-1,n,k)))+1 );
                rows(rcNum) = idx;  cols(rcNum) = idx1;
            end

            idx1 = map_co( m, n-1, k);
            if 0 ~= idx1
                % weight of n link
                rcNum = rcNum + 1;
                weights(rcNum) =  wn_tab( abs(double(local_stk(m,n,k))-...
                    double(local_stk(m,n-1,k)))+1 );
                rows(rcNum) = idx;  cols(rcNum) = idx1;
            end
        end
    end
end

% the 12 boundaries
for n = [1 N]
    for k = [1 K]
        for m = 2 : M
            idx = map_co( m,n,k );
            
            idx1 = map_co(m-1,n,k);
            if 0 ~= idx1 && 0~=idx
                % weight of n link
                rcNum = rcNum + 1;
                weights(rcNum) =  wn_tab( abs(double(local_stk(m,n,k))-...
                    double(local_stk(m-1,n,k)))+1 );
                rows(rcNum) = idx;  cols(rcNum) = idx1;
            end
        end
    end
end

for m = [1 M]
    for k = [1 K]
        for n = 2 : N
            idx = map_co( m,n,k );

            idx1 = map_co(m,n-1,k);
            if 0 ~= idx1 && 0~=idx
                % weight of n link
                rcNum = rcNum + 1;
                weights(rcNum) =  wn_tab( abs(double(local_stk(m,n,k))-...
                    double(local_stk(m,n-1,k)))+1 );
                rows(rcNum) = idx;  cols(rcNum) = idx1;
            end
        end
    end
end

for m = [1 M]
    for n = [1 N]
        for k = 2 : K
            idx = map_co( m,n,k );
            idx1 = map_co(m,n,k-1);
            if 0 ~= idx1 && 0~=idx
                % weight of n link
                rcNum = rcNum + 1;
                weights(rcNum) =  wn_tab( abs(double(local_stk(m,n,k))-...
                    double(local_stk(m,n,k-1)))+1 );
                rows(rcNum) = idx;  cols(rcNum) = idx1;
            end
        end
    end
end

% cut the tail
rows((rcNum +1) : end) = [];
cols((rcNum+1) : end ) = [];
weights( (rcNum+1) : end ) = [];
A = sparse( double([rows; cols]), double([cols; rows;]), [weights; weights], Nv, Nv );
