%% test sigmoid function
x = 0 : 255;
alpha = 128;
beta = -128/2/log(99);
% s = 128;
y = 1./(1+exp( (x-alpha)./beta ) );
plot(x,y)

%%
T = 128;
alpha = (255+2*T)/4;
beta = -255/4/log(99);   % constant!
w_lut = 1./(1+exp( ((0:255)-alpha)./beta ) );
plot( 0:255, w_lut)

%%
% parameters
delta = 50;

% build the lookup table of the boundary energy, weight of n link
wn_tab = exp( -1 * (0:255).*(0:255) / (2*delta*delta) )

plot( wn_tab)

%% test midpoint circle
[nc kc] = getmidpointcircle(5,5,3);
% eliminate the redundent points
% nc2 = nc;   nc2(1:end-1) = nc(2:end);   nc2(end) = nc(1); 
% kc2 = kc;   kc2(1:end-1) = kc(2:end);   kc2(end) = kc(1);
% label = (nc == nc2) & (kc==kc2);
% nc(label) = []; kc(label) = [];

for p = 1 : length(nc)
    plot(nc(p), kc(p), '.g');
    hold on
    text(nc(p), kc(p),['\leftarrow' num2str(p)])
end

%% test LABELS
% [M,N,K,It] = max(map_id);
M = 31;
N = 31;
K = 31;
stk = false(M,N,K);

for id = 1 : length(fidx)
    if C( id ) > 0
        stk( map_id(id,1), map_id(id,2), map_id(id,3) ) = true;
    end
end
figure; h = vol3d('CData',stk);

%% test plot marker size
x = 1:50;
y = rand(1,50);
s = 5; %Marker width in units of X

h = scatter(x,y); % Create a scatter plot and return a handle to the 'hggroup' object

%Obtain the axes size (in axpos) in Points
currentunits = get(gca,'Units');
set(gca, 'Units', 'Points');
axpos = get(gca,'Position');
set(gca, 'Units', currentunits);

markerWidth = s/diff(xlim)*axpos(3); % Calculate Marker width in points

set(h, 'SizeData', markerWidth^2)

%% test function speed
tic
map = imresize( c1, size(c2), 'nearest' );
toc

tic
map = imresize( c1, size(c2), 'nearest' );
toc

%% test mapping
map2 = zeros(size(c2), 'uint16');
lengthRatio = size(c2,1) / size(c1,1);
position = floor( (1: size(c2,1)) / lengthRatio ) +1 ;
position(end) = position(end) - 1;
map2 = c1( position, : );

%% test

for pi = 1 : length(c)-1
    % add this cooridnate
    rcNum = rcNum + 1;
    rows(rcNum) = map_co( c(pi,1),  c(pi,2),  c(pi,3) );
    cols(rcNum) = map_co( c(pi+1,1),  c(pi+1,2),  c(pi+1,3) );
    weights(rcNum) = wn_tab( local_stk( c(pi,1),  c(pi,2),  c(pi,3) ) - ...
        local_stk( c(pi+1,1),  c(pi+1,2),  c(pi+1,3) ) + 1 );
end

%% test2
pi = 1 : length(c)-1;
rows(rcNum + pi) = map_co( sub2ind( size(map_co),...
    c(pi,1),  c(pi,2),  c(pi,3) ) );

cols(rcNum + pi) = map_co( sub2ind( size(map_co),...
    c(pi+1,1),  c(pi+1,2),  c(pi+1,3) ) );
idx_list = local_stk( sub2ind( size(local_stk), ...
    c(pi,1),  c(pi,2),  c(pi,3) ) ) - ...
    local_stk( sub2ind( size(local_stk),...
    c(pi+1,1),  c(pi+1,2),  c(pi+1,3) ) )  + 1;
weights(rcNum + (1 : length(c)-1)) = wn_tab( idx_list );
