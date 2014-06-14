%% wave propagation from a single seed
% by jpwu, 2013/02/28

function network = wave_propagation( seed, network )
global stk;
%% parameters
% the size scale of local cube
P = 1.5;
% the minimun voxel number of a wave front
Ta = 0;

%% visualization system
MIP = max( stk, [], 3 );
[M N K] = size(stk);
imshow(MIP); hold on
pause(0.1); % draw this figure

% propagation processing
% adjust the seed radius and position
seed(4) = get_radius( seed(1:3), seed(4) );
seed(1:3) = recenter( seed(1:3), seed(4), P );
seed(4) = get_radius( seed(1:3), seed(4) );

% section number of current network
sn = length( network.sections );

% dynamic seeds list
dsl = seed;
sec_idx = sn;
bbs = seed(4) * 2;

% iterative tracing
while ~isempty(dsl)
    % get the last seed and current section number
    seed = uint16( dsl( end,: ) );
    % current section index
    csi = sec_idx(end);
    % store the parent bounding box size
    ps = bbs(end);
    
    % delet current seed from the dynamic list
    dsl(end, :) = [];
    sec_idx(end) = [];
    bbs(end) = [];
    
    % get the local voxel cube
    [local_stk, m1, n1, k1]= get_local_stk(stk, seed, 2);
    % get threshold in a big local stack
    T = kmeans_binarize( local_stk );
%     
%     % get the local voxel cube and binarize
%     [local_stk, m1, n1, k1]= get_local_stk(stk, seed, P);
    mk_local_stk = ( local_stk > T );
    S = ca_cube( mk_local_stk, m1,n1,k1, seed, P );
    Nb = length(S);

    Nb
    
    %% compute the centroid and radius of each component
    sn = length(network.sections);
    if Nb == 0
        disp('end of propagation!')
        continue;
    elseif Nb == 1
        % only one component, add new node
        % current size of the bounding box
        cs = sqrt( S(1).BoundingBox(4)^2 + S(1).BoundingBox(5)^2 + S(1).BoundingBox(6)^2 );
        
        % initialize the center and radius
        box_center = S(1).Centroid + double( [m1-1, n1-1, k1-1] );
        center = get_node_position( box_center, seed, ps, cs );
        % estimate the radius by a edge length adjustable cube
%         radius = get_radius(center, cs/2);
        radius = get_radius_V2(center, T);
        
        % estimate the centroid by mean shift
        center = recenter( center, radius, P );
        % get the new radius
%         radius = get_radius(center, radius);
        radius = get_radius_V2(center, T);
        % estimate the centroid by mean shift
        center = recenter( center, radius, P );
        % get the new radius
%         radius = get_radius(center, radius);
        radius = get_radius_V2(center, T);
        
        center
        
        % update the dynamic list
        dsl = [ dsl; [center radius] ];
        sec_idx = [ sec_idx; csi ];
        bbs = [bbs; cs];
        % add this new node to current section
        network.sections{csi} = [network.sections{csi}; [center radius ] ];
        
        % update the visualization system
        plot( center(2), center(1), '.',...
            'color',[rand, rand, rand], 'markersize',3*radius );
        plot( [seed(2) center(2)], [seed(1) center(1)],'-r','LineWidth',2 );
        pause(0.1); % draw this figure
    else
        % have branches
        for nb = 1 : Nb
            if S(nb).Area > Ta
                % current size of the bounding box
                cs = sqrt( S(nb).BoundingBox(4)^2 + S(nb).BoundingBox(5)^2 + S(nb).BoundingBox(6)^2 );
                % initialize the center and radius
                box_center = S(nb).Centroid + double( [m1-1, n1-1, k1-1] );
                center = get_node_position( box_center, seed, ps, cs );
                
                % estimate the radius by a edge length adjustable cube
%                 radius = get_radius(center, cs/2);
                radius = get_radius_V2(center, T);
                % estimate the centroid by mean shift
                center = recenter( center, radius, P );
                % get the new radius
%                 radius = get_radius(center, radius);
                radius = get_radius_V2(center, T);
                % estimate the centroid by mean shift
                center = recenter( center, radius, P );
                % get the new radius
%                 radius = get_radius(center, radius);
                radius = get_radius_V2(center, T);
                
                center
                
                % connect with the branch point                
                % the new branches
                sn = length( network.sections ) + 1;
                % update the dynamic list
                dsl = [ dsl; [center radius] ];
                sec_idx = [sec_idx; sn];
                bbs = [bbs; cs];
                % add seed and this new branch node to a new section
                network.sections{sn} = [seed; center radius ];
                
                % update the visualization system
                plot( center(2), center(1), '.', ...
                    'color',[rand, rand, rand], 'markersize', 3*radius );
                plot( [seed(2) center(2)], [seed(1) center(1)],'-r','LineWidth',2 );
                pause(0.1); % draw this figure
            end
        end
    end
end

network.sn = length(network.sections);
network = nio_build_net_connectivity(network);
return;

%% connectivity analysis in a local cube
function S = ca_cube( mk_local_stk, m1,n1,k1, seed, P )
global mk_stk;

% eliminate the visited voxels
[Ml, Nl, Kl] = size( mk_local_stk );

% % morphological operation
% str = ones(3,3,3);
% mk_local_stk = imclose( mk_local_stk, str );

%%  wave propagation in the local cube
% directly remove the propagation progress and get binary boundary
disp('---- the propagation time: ')
% local connectivity label
[Ll Num] = bwlabeln(mk_local_stk);
% get the connectivity component index
if mk_local_stk( seed(1)-m1+1, seed(2)-n1+1, seed(3)-k1+1 )
    cci = Ll( seed(1)-m1+1, seed(2)-n1+1, seed(3)-k1+1 );
else
    % the seed is in background !
    disp('seed in background !');
    S = [];
    return;
end
% eliminate other component
Ll = (Ll==cci);
disp('--------- the connectivity analysis of wave front: ')

% eliminate central voxels
r = uint16( seed(4) );
ml1 = max( 2, Ml/2 - r );   ml2 = min( Ml-1, Ml/2 + r );
nl1 = max( 2, Nl/2 - r);    nl2 = min( Nl-1, Nl/2 + r );
kl1 = max( 2, Kl/2 - r);    kl2 = min( Kl-1, Kl/2 + r );

% mark the traced voxels
mk_stk ( m1+ml1-1:(m1+ml2-1), n1+nl1-1:(n1+nl2-1), k1+kl1-1:(k1+kl2-1) ) = ...
    mk_stk ( m1+ml1-1:(m1+ml2-1), n1+nl1-1:(n1+nl2-1), k1+kl1-1:(k1+kl2-1) ) ...
    | Ll( ml1:ml2, nl1:nl2, kl1:kl2 ) ;

Ll( ml1:ml2, nl1:nl2, kl1:kl2 ) = 0;
% Ll(:(Ml-1), 2:(Nl-1), 2:(Kl-1)) = 0;

% eliminate the traced voxels
Ll = Ll & ~mk_stk ( m1:(m1+Ml-1), n1:(n1+Nl-1), k1:(k1+Kl-1) );

% connectivity analysis of wave front voxels, compute the branch number
Llf = bwlabeln(Ll);
S = regionprops(Llf, 'area', 'centroid','BoundingBox');
% the centroid is in amira coordinate system, transform
Nb = length(S);
for nb = 1 : Nb
    tmp = S(nb).Centroid(1);
    S(nb).Centroid(1) = S(nb).Centroid(2);
    S(nb).Centroid(2) = tmp;
end
return;

%% estimate the new node position according to voxel scooping algorithm
function np = get_node_position( box_center, seed, ps, cs )
ratio = min( ps/cs, cs/ps );
% parent center
pc = double( seed(1:3) );
np = pc + (box_center - pc ) * (0.5^ratio); 
return;

%% relocate the center by mean shift
function c = recenter(c, r, P)
% P = 1.2;
global stk;
% get the local stack
[local_stk, m1, n1, k1] = get_local_stk( stk, [c r], P );
[Ml, Nl, Kl] = size( local_stk );

% compute the center of local stack
si = sum( local_stk(:) );
% % estimate the weight of each intensity of 0-255
% wi_local_stk = double(local_stk) ./ si;
c = double([0 0 0]);
local_stk = double( local_stk );

for m = 1 : Ml
    for n = 1 : Nl
        for k = 1 : Kl
            c(1) = c(1) + double(m)*local_stk(m,n,k);
            c(2) = c(2) + double(n)*local_stk(m,n,k);
            c(3) = c(3) + double(k)*local_stk(m,n,k);
        end
    end
end

c = c./si + double( [ (m1-1), (n1-1), (k1-1) ] );
return;

%% get the local cube
function [local_stk, m1, n1, k1] = get_local_stk( stk, centroid, P )
[M N K] = size(stk);

% the coordinate and radius
ms = centroid(1);   ns = centroid(2);
ks = centroid(3);   rs = centroid(4);

m1 = uint16( max( 1, ms-P*rs ) );   m2 = uint16( min(M, ms+P*rs) );
n1 = uint16( max( 1, ns-P*rs ) );   n2 = uint16( min(N, ns+P*rs) );
k1 = uint16( max( 1, ks-P*rs ) );   k2 = uint16( min(K, ks+P*rs) );

local_stk = stk(m1:m2, n1:n2, k1:k2);
return;