% build the n link graph in a sphere model !
% the neighbor relationship was defined in circles!
% by jpwu, 2014/04/26
function [map_co, map_id, uA]= build_n_link_graph_V3_sphere( circles, local_stk, mk_local_stk )
% parameters
delta = 5;  % 4, 20

% stk size
[M N K] = size(local_stk);

% the coordinate list
co_list = [];
for ci = 1 : length(circles)
    co_list = [co_list; circles{ci}];
end
Nv = size(co_list,1);

% initialize the coordinate map
map_co = zeros(M,N,K,'uint16');
for nv = 1 : Nv
    map_co(co_list(nv,1), co_list(nv,2), co_list(nv,3)) = nv;
end

% initiate the graph model
% initialize the adjacency matrix
%A = sparse( Nv, Nv );
A = initialize_connection_matrix(circles, map_co, Nv);

% build the lookup table of the boundary energy, weight of n link
wn_tab = exp( -1 * (0:255).*(0:255) / (2*delta*delta) );

% find the neighbors of each voxel
% the neighbor relationship inside the circle
for ci = 1 : length(circles)
    c = circles{ ci };
    % add the last element to establish the connection between start and
    % end
    c = [c; c(1,:)];
    
    % establish the connection between neighboring voxels
    for pi = 1 : length(c)-1
        idx1 = map_co( c(pi,1),  c(pi,2),  c(pi,3) );
        idx2 = map_co( c(pi+1,1),  c(pi+1,2),  c(pi+1,3) );
        
        weight = wn_tab( local_stk( c(pi,1),  c(pi,2),  c(pi,3) ) - ...
            local_stk( c(pi+1,1),  c(pi+1,2),  c(pi+1,3) ) + 1 );
        A(idx2,idx1) = weight;  A(idx1,idx2) = weight;        
    end
end

% the small circle was interpreated to big circle for mapping
for ci = 1 : length(circles)-1
    % c2 is bigger than c1
    if length(circles{ci}) > length( circles{ci+1} )
        c2 = circles{ ci };
        c1 = circles{ ci+1 };
    else
        c1 = circles{ ci };
        c2 = circles{ ci+1 };
    end
    
    % interpreate c1 for mapping
    %map = imresize( c1, size(c2), 'nearest' );
    %map2 = zeros(size(c2), 'uint16');
    lengthRatio = size(c2,1) / size(c1,1);
    position = floor( (1: size(c2,1)) / lengthRatio ) +1 ;
    position(end) = position(end) - 1;
    map = c1( position, : );
    
    % build the upper and lower connection
    for pi = 1 : size(c2,1)
        idx1 = map_co( map(pi,1), map(pi,2), map(pi,3) );
        idx2 = map_co( c2(pi,1),  c2(pi,2),  c2(pi,3) );
        
        weight = wn_tab( local_stk( map(pi,1), map(pi,2), map(pi,3) ) - ...
            local_stk( c2(pi,1),  c2(pi,2),  c2(pi,3) ) + 1 );
        A(idx2,idx1) = weight;  A(idx1,idx2) = weight;
    end
    
    % build the corner connection
    map = [map; map(1,:)]; c2 = [c2; c2(1,:)];
    for pi = 2 : size(c2,1)-1
        idx1 = map_co( map(pi,1), map(pi,2), map(pi,3) );
        idx3 = map_co( c2(pi-1,1),  c2(pi-1,2),  c2(pi-1,3) );
        
        weight = wn_tab( local_stk( map(pi,1), map(pi,2), map(pi,3) ) - ...
            local_stk( c2(pi-1,1),  c2(pi-1,2),  c2(pi-1,3) ) + 1 );
        A(idx3,idx1) = weight;  A(idx1,idx3) = weight;
        
        idx1 = map_co( map(pi,1), map(pi,2), map(pi,3) );
        idx4 = map_co( c2(pi+1,1),  c2(pi+1,2),  c2(pi+1,3) );
        
        weight = wn_tab( local_stk( map(pi,1), map(pi,2), map(pi,3) ) - ...
            local_stk( c2(pi+1,1),  c2(pi+1,2),  c2(pi+1,3) ) + 1 );
        A(idx4,idx1) = weight;  A(idx1,idx4) = weight;
    end
end

% lable the unvisited voxels
luv = false(Nv,1);
for nv = 1 : Nv
    luv(nv) = ~mk_local_stk( co_list(nv,1), co_list(nv,2), co_list(nv,3) );
end

% eliminate the visited voxels
uvIdx = find(luv);
Nu = length( uvIdx );
uA = sparse( Nu, Nu );
for nu = 1 : Nu
    uA(nu,:) = A( uvIdx(nu), uvIdx );
end

% update the co_list
co_list = co_list(luv,:);
Nv = size(co_list,1);
% build the map of coordinate and index
map_co = zeros(M,N,K,'uint16');
map_id = zeros(Nv,4,'uint16');
for nv = 1 : Nv
    map_co(co_list(nv,1), co_list(nv,2), co_list(nv,3)) = nv;
    map_id(nv,:) = [co_list(nv,:), ...
        uint16( local_stk(co_list(nv,1), co_list(nv,2), co_list(nv,3)) )];
end
return;

%% initialize the sparse matrix for efficiency
function A = initialize_connection_matrix(circles, map_co, Nv)
rows = zeros(10 * Nv, 1);
columns = zeros(10 * Nv, 1);
rcNum = 0;

% find the neighbors of each voxel
% the neighbor relationship inside the circle
for ci = 1 : length(circles)
    c = circles{ ci };
    % add the last element to establish the connection between start and
    % end
    c = [c; c(1,:)];
    
    % establish the connection between neighboring voxels
    for pi = 1 : length(c)-1
        idx1 = map_co( c(pi,1),  c(pi,2),  c(pi,3) );
        idx2 = map_co( c(pi+1,1),  c(pi+1,2),  c(pi+1,3) );
        % add this cooridnate
        rcNum = rcNum + 1;
        rows(rcNum) = idx2;    columns(rcNum) = idx1;      
    end
end

% the small circle was interpreated to big circle for mapping
for ci = 1 : length(circles)-1
    % c2 is bigger than c1
    if length(circles{ci}) > length( circles{ci+1} )
        c2 = circles{ ci };
        c1 = circles{ ci+1 };
    else
        c1 = circles{ ci };
        c2 = circles{ ci+1 };
    end
    
    % interpreate c1 for mapping
    %map = imresize( c1, size(c2), 'nearest' );
    %map2 = zeros(size(c2), 'uint16');
    lengthRatio = size(c2,1) / size(c1,1);
    position = floor( (1: size(c2,1)) / lengthRatio ) +1 ;
    position(end) = position(end) - 1;
    map = c1( position, : );
    
    % build the upper and lower connection
    for pi = 1 : size(c2,1)
        idx1 = map_co( map(pi,1), map(pi,2), map(pi,3) );
        idx2 = map_co( c2(pi,1),  c2(pi,2),  c2(pi,3) );
        % add this cooridnate
        rcNum = rcNum + 1;
        rows(rcNum) = idx2;    columns(rcNum) = idx1;
    end
    
    % build the corner connection
    map = [map; map(1,:)]; c2 = [c2; c2(1,:)];
    for pi = 2 : size(c2,1)-1
        idx1 = map_co( map(pi,1), map(pi,2), map(pi,3) );
        idx3 = map_co( c2(pi-1,1),  c2(pi-1,2),  c2(pi-1,3) );
        
        % add this cooridnate
        rcNum = rcNum + 1;
        rows(rcNum) = idx3;    columns(rcNum) = idx1;  
        
        idx1 = map_co( map(pi,1), map(pi,2), map(pi,3) );
        idx4 = map_co( c2(pi+1,1),  c2(pi+1,2),  c2(pi+1,3) );
        % add this cooridnate
        rcNum = rcNum + 1;
        rows(rcNum) = idx4;    columns(rcNum) = idx1;
    end
end

% cut the tail
rows((rcNum +1) : end) = [];
columns( (rcNum+1) : end ) = [];
A = sparse( double([rows; columns]), double([columns; rows;]), 0, Nv, Nv );
