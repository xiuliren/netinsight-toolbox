% build the n link graph in a sphere model !
% the neighbor relationship was defined in circles!
% by jpwu, 2013/10/04
function [map_co, map_id, uA]= build_n_link_graph_V3_2sphere( circles1, circles2, local_stk, mk_local_stk )
%% parameters
delta = 5;  % 4, 20

% stk size
[M N K] = size(local_stk);

%%  the coordinate list
co_list1 = [];
for ci = 1 : length(circles1)
    co_list1 = [co_list1; circles1{ci}];
end
Nv1 = size(co_list1,1);

co_list2 = [];
for ci = 1 : length(circles2)
    co_list2 = [co_list2; circles2{ci}];
end
Nv2 = size(co_list2,1);

%% initialize the coordinate map
map_co = zeros(M,N,K,'uint16');
for nv = 1 : Nv1
    map_co(co_list1(nv,1), co_list1(nv,2), co_list1(nv,3)) = nv;
end

for nv = 1 : Nv2
    map_co(co_list2(nv,1), co_list2(nv,2), co_list2(nv,3)) = nv + Nv1;
end

%% initiate the graph model
% initialize the adjacency matrix
A = sparse( Nv1+Nv2, Nv1+Nv2 );

% build the lookup table of the boundary energy, weight of n link
wn_tab = exp( -1 * [0:255].*[0:255] / (2*delta*delta) );

%% find the neighbors of each voxel
% the neighbor relationship inside the circle
for ci = 1 : length(circles1)
    c = circles1{ ci };
    % connect the neighbors
    A = connect_neighbors_in_circle(c, A, wn_tab, map_co, local_stk);
end

% the neighbor relationship inside the circle
for ci = 1 : length(circles2)
    c = circles2{ ci };
    % connect the neighbors
    A = connect_neighbors_in_circle(c, A, wn_tab, map_co, local_stk);
end

%% the small circle was interpreated to big circle for mapping
for ci = 1 : length(circles1)-1
    % map the circles
    A = map_circles( circles1{ ci }, circles1{ ci+1 }, A, local_stk, map_co, wn_tab );
end

% the small circle was interpreated to big circle for mapping
for ci = 1 : length(circles2)-1
    % map the circles
    A = map_circles( circles2{ ci }, circles2{ ci+1 }, A, local_stk, map_co, wn_tab );
end

%% connect inside and outside spheres
for ci = 1 : length( circles2 )
    % map the circles
    A = map_circles( circles1{ ci+1 }, circles2{ ci }, A, local_stk, map_co, wn_tab );
end

%% lable the unvisited voxels
luv = false(Nv1+Nv2, 1);
for nv = 1 : Nv1
    luv(nv) = ~mk_local_stk( co_list1(nv,1), co_list1(nv,2), co_list1(nv,3) );
end
for nv = 1 : Nv2
    luv(nv+Nv1) = ~mk_local_stk( co_list2(nv,1), co_list2(nv,2), co_list2(nv,3) );
end
%% eliminate the visited voxels
uvIdx = find(luv);
Nu = length( uvIdx );
uA = sparse( Nu, Nu );
for nu = 1 : Nu
    uA(nu,:) = A( uvIdx(nu), uvIdx );
end

%% update the co_list
co_list = [co_list1; co_list2];
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

%% map the circles and build the connection
function A = map_circles( c1, c2, A, local_stk, map_co, wn_tab )
% let c1 be the small one, and c2 be the big one
if size(c1,1) > size(c2,1)
    tmp = c1;    c1 = c2;   c2 = tmp;
end
% interpreate c1 for mapping
map = imresize( c1, size(c2), 'nearest' );
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

%% connect the neighbors inside circle
function A = connect_neighbors_in_circle(c, A, wn_tab, map_co, local_stk)
% add the last element to establish the connection between start and end
c = [c; c(1,:)];

% establish the connection between neighboring voxels
for pi = 1 : length(c)-1
    idx1 = map_co( c(pi,1),  c(pi,2),  c(pi,3) );
    idx2 = map_co( c(pi+1,1),  c(pi+1,2),  c(pi+1,3) );
    
    weight = wn_tab( local_stk( c(pi,1),  c(pi,2),  c(pi,3) ) - ...
        local_stk( c(pi+1,1),  c(pi+1,2),  c(pi+1,3) ) + 1 );
    A(idx2,idx1) = weight;  A(idx1,idx2) = weight;
end