%% connectivity analysis in a local cube
function [network mk_stk]= voxelscooping( stk, mk_stk, balls, sdtable, seed, network )
%global stk;
% global mk_stk
%% parameters

%% evaluate the seed
if mk_stk( uint16(seed(1)), uint16(seed(2)), uint16(seed(3)) )
    % already visited
    return;
end

%% initialization
% current section number
network.sn = size( network.sections, 1 ) + 1;
network.sections{network.sn} = seed;

% the dynamic list of seeds, scooping distance, 
% and axis aligned bounding box
dsl = zeros( 100000, 7 );
dsl(1,:) = [network.sn, seed, 1.3*seed(4), 1.5*seed(4)];
ndsl = 1;

% visualization system
% global stk;
mip = max(stk, [], 3);
imshow( mip );  hold on;
pause(0.1); % draw this figure

%% voxel scooping
while ndsl > 0
    % current section index
    csi = dsl(ndsl, 1);
    % get the seed, scooping distance, and axis aligned bounding box
    seed = dsl(ndsl, 2:5);
    sd = dsl(ndsl, 6);
    paabb = dsl(ndsl, 7);
    dsl(ndsl,:) = [0 0 0 0 0 0 0];
    ndsl = ndsl - 1;
    
%     tic, disp('time of create voxel cluster: ')
    % create voxel cluster from seed voxel
%     [vc mk_stk]= create_voxel_cluster_V1( stk, mk_stk, seed,  sd);
    [vc mk_stk]= create_voxel_cluster_V2( stk, mk_stk, balls, seed, sd );
    if isempty(vc)
        continue;
    end
%     toc
    
%     tic, disp('time of get node list from scooped voxels: ')
    % get node list from scooped voxel
    [nodeList sdList aabbList]= create_nodes(stk, sdtable, vc, seed, paabb);
%     toc
    
    % update dynamic seeds list
    [network, dsl, ndsl] = update_dsl( seed, network, dsl, ndsl, nodeList, csi, sdList, aabbList  );
    
    % update the visualization system
    if seed(4)~=0
        plot( seed(2)+1, seed(1)+1, '.') %, 'color',[rand, rand, rand], ); %'markersize', 3*seed(4) );
%             plot( [seed(2)+1 center(2)], [seed(1)
%             center(1)],'-r','LineWidth',2 );
        pause(0.005); % draw this figure
    end
%     disp(['number of sections:  ' num2str( length(network.sections) )]);
end

% rebuild network
network.sn = length(network.sections);
network = nio_build_net_connectivity(network);
return;

%% update the network and dynamic seeds list
function [network, dsl, ndsl] = update_dsl_V2( seed, network, dsl, ndsl, nodeList, csi, sdList, aabbList  )
% number of branches
Nb = size(nodeList, 1);

% number of sections
sn = length(network.sections);

if Nb == 0
    return;
elseif Nb == 1
    ndsl = ndsl + 1;
    % add in current section
    network.sections{csi} = [ network.sections{csi}; nodeList ];
    % update dynamic list
    dsl(ndsl,:) = [csi, nodeList, sdList, aabbList];
else
    % have several branches, sort acoording to radius
    nodeList = sortrows(nodeList, 4);
    for nb = 1 : Nb
        ndsl = ndsl + 1;
        % section index
        si = sn + nb;
        % add in a new section
        network.sections{si} = [ seed; nodeList(nb,:) ];
        dsl(ndsl,:) = [si, nodeList(nb,:), sdList(nb), aabbList(nb)];
    end
end
return;

%% update the network and dynamic seeds list
function [network, dsl, ndsl] = update_dsl( seed, network, dsl, ndsl, nodeList, csi, sdList, aabbList  )
% number of branches
Nb = size(nodeList, 1);
% number of sections
sn = length(network.sections);

if Nb == 0
    return;
elseif Nb == 1
    ndsl = ndsl + 1;
    % add in current section
    network.sections{csi} = [ network.sections{csi}; nodeList ];
    % update dynamic list
    dsl(ndsl,:) = [csi, nodeList, sdList, aabbList];
else
    % have several branches
    for nb = 1 : Nb
        ndsl = ndsl + 1;
        % section index
        si = sn + nb;
        % add in a new section
        network.sections{si} = [ seed; nodeList(nb,:) ];
        
        % add new seed according to the node diameter, assending order
        ind = find( dsl(1:ndsl,5) > nodeList(nb,4), 1 ) - 1;
        if isempty( ind )
            % the largest one
%             ind = size( dsl, 1 );
            dsl(ndsl,:) = [si, nodeList(nb,:), sdList(nb), aabbList(nb)];
        else
            % not the largest one
            dsl = insertrows( dsl, [si, nodeList(nb,:), sdList(nb), aabbList(nb)], ind );
%             dsl(ndsl,:) = [si, nodeList(nb,:), sdList(nb), aabbList(nb)];
        end        
    end
end
return;

%% create node from a voxel cluster
function [nodeList sdList aabbList]= create_nodes( stk, sdtable, vc, pnode, paabb )
% area threshold
Ta = 3;

% initialize the node list
nodeList = [];
aabbList = [];
sdList = [];
if isempty( vc )
    return;
end
% connectivity analysis
[Ll NL]= bwlabeln( logical(vc.neighboringvoxels) ); 
if NL==0
%     disp(' no any new connectivity regions');
    return;
end

% tic, disp('the time of node computing:    ')

% have connectivity regions
S = regionprops(Ll, 'area', 'centroid','BoundingBox');
for nb = 1 : NL
    if S(nb).Area > Ta
        % current axis aligned bounding box
        caabb = sqrt( S(nb).BoundingBox(4)^2 + S(nb).BoundingBox(5)^2 + S(nb).BoundingBox(6)^2 );
        ratio = min( caabb/paabb, paabb/caabb );
        center =  S(nb).Centroid;
        tmp = center(1);    center(1) = center(2); center(2)=tmp;
        center  = center + double(vc.origin) -1;
        center = pnode(1:3) + (0.5^ratio)*( center -  pnode(1:3));
        radius = get_radius_V2( stk, center, vc.threshold );
%         center = recenter( center, radius, 1.2 );
%         radius = get_radius_V2( center, vc.threshold );
        
        % update the node list
        nodeList = [ nodeList; center, radius ];
        % axis aligned bounding box
        aabbList = [aabbList; caabb];
        % compute the scooping distance
%         sd = get_scooping_distance( center, Ll, vc.origin,  nb);
        sd = get_scooping_distance_V2( sdtable, center, Ll, vc.origin,  nb);
        sdList = [ sdList; sd ];
    end
end
% toc
return;

%% get the scooping distance
function sd = get_scooping_distance_V2( sdtable, c, Ll, origin, nb )
[Ml Nl Kl] = size( Ll );

% transfer center to local coordinate
c = uint16( double(c) - origin +1 );
% squre of the distance
% maxd = 0;

% find the foreground voxels
idx = find( Ll==nb );
cmat = repmat(c, length(idx), 1 );
[m n k] = ind2sub(size(Ll), idx);
dmat = abs( uint16([m n k]) - cmat) + 1;

% % get the max distance
% ind_vec = sub2ind(size(sdtable), dmat(:,1), dmat(:,2), dmat(:,3));
% sd = max( sdtable(ind_vec) );

dmax = 0;
for ni = 1 : length( idx )
    dmax = max( dmax, sdtable(dmat(ni,1), dmat(ni,2), dmat(ni,3)) );
end
sd = dmax;
return;

%% get the scooping distance
function sd = get_scooping_distance( c, Ll, origin, nb )
[Ml Nl Kl] = size( Ll );

% transfer center to local coordinate
c = double(c) - origin +1;
% squre of the distance
maxd2 = 0;

% for m = 1 : Ml
%     for n = 1 : Nl
%         for k = 1 : Kl
%             if Ll(m,n,k)==nb
%                 % compute the distance
%                 d2 = (c(1)-m)*(c(1)-m) + (c(2)-n)*(c(2)-n) + (c(3)-k)*(c(3)-k);
%                 maxd2 = max( maxd2, d2 );
%             end
%         end
%     end
% end

% matrix computing is faster than for curculation
idx = find( Ll==nb );
subarray = ind2sub(size(Ll), idx);
for ni = 1 : length(idx)
%     [m n k] = ind2sub(size(Ll), idx(ni));
    m = subarray(ni,1);
    n = subarray(ni,2);
    k = subarray(ni,3);
    d2 = (c(1)-m)*(c(1)-m) + (c(2)-n)*(c(2)-n) + (c(3)-k)*(c(3)-k);
    maxd2 = max( maxd2, d2 );
end

sd = sqrt( maxd2 );
return;

%% create voxel cluster from seed voxel in cube
function [vc mk_stk]= create_voxel_cluster_V3( stk, mk_stk, seed, sd )
% adjust the scooping distance according to cube 
sd = sd/1.732;

% adjust the scooping distance
sd = max( uint16(sd), 2 );

seed = uint16(seed);
% get local image stack 
[local_stk, m1, n1, k1]= get_local_stk(seed, 1.5*sd);
[Ml Nl Kl] = size( local_stk );
vc.origin = double( [m1 n1 k1] );

% get threshold in a big local stack
T = kmeans_binarize( local_stk );
vc.threshold = T;
% evaluate the intensity of seed voxel
if stk( seed(1), seed(2), seed(3) ) < T
    vc = [];
    return;
end
% binarized image stack
mkl = (local_stk > T);

% local connectivity label
Ll = bwlabeln(mkl);
% get the connectivity component index
if mkl( seed(1)-m1+1, seed(2)-n1+1, seed(3)-k1+1 )
    cci = Ll( seed(1)-m1+1, seed(2)-n1+1, seed(3)-k1+1 );
else
%     disp('the seed is in background !');
    vc = [];
    return;
end
% eliminate other component
Ll = (Ll==cci);

% eliminate the traced voxels
Ll = Ll & ~mk_stk ( m1:(m1+Ml-1), n1:(n1+Nl-1), k1:(k1+Kl-1) );

% local seed
ls= seed(1:3) - [m1 n1 k1] + 1;
mls = double(ls(1));    nls = double( ls(2) );    kls = double( ls(3) );

% tic, disp('the time of distance evalueation:    ')
% neighboring voxels
nv = zeros( Ml,Nl,Kl, 'uint8' );
nv = logical( nv );

% the searching range of outside boundary
mo1 = uint16(max( 1, mls-sd-2 ));   mo2 = uint16(min( Ml, mls+sd+2 ));
no1 = uint16(max( 1, nls-sd-2 ));   no2 = uint16(min( Nl, nls+sd+2 ));
ko1 = uint16(max( 1, kls-sd-2 ));   ko2 = uint16(min( Kl, kls+sd+2 ));
% the searching range of inner boundary
mw1 = uint16(max( 1, mls-sd));    mw2 = uint16(min( Ml, mls+sd ));
nw1 = uint16(max( 1, nls-sd));    nw2 = uint16(min( Nl, nls+sd ));
kw1 = uint16(max( 1, kls-sd));    kw2 = uint16(min( Kl, kls+sd ));

% mark voxels as visited
mk_stk( mo1+m1-1:mo2+m1-1, no1+n1-1:no2+n1-1, ko1+k1-1:ko2+k1-1 ) = ...
    mk_stk( mo1+m1-1:mo2+m1-1, no1+n1-1:no2+n1-1, ko1+k1-1:ko2+k1-1 ) | Ll(mo1:mo2,no1:no2,ko1:ko2);

% matrix implemention
Ll( mw1:mw2, nw1:nw2, kw1:kw2 ) = 0;
Ll(1:mo1, :,:) = 0; Ll(:, 1:no1,:) = 0; Ll(:,:, 1:ko1) = 0;
L1(mo2:Ml,:,:) = 0; Ll(:,no2:Nl,:) = 0; Ll(:,:,ko2:Kl) = 0;   
vc.neighboringvoxels = Ll;
return;

%% create voxel cluster from seed voxel
function [vc mk_stk]= create_voxel_cluster_V2( stk, mk_stk, balls, seed, sd )
% adjust the scooping distance
sd = max( uint16(sd), 2 );

seed = uint16(seed);
% get local image stack 
cuberadius = uint16( 1.5*sd );
[local_stk, m1, n1, k1]= get_local_stk(stk, seed, cuberadius);
% [Ml Nl Kl] = size( local_stk );
% vc.origin = double( [m1 n1 k1] );

% get threshold in a big local stack
T = kmeans_binarize( local_stk );
vc.threshold = T;
% evaluate the intensity of seed voxel
if stk( seed(1), seed(2), seed(3) ) < T
    vc = [];
    return;
end

% get the small stack
cuberadius = uint16( sd+2 );
[local_stk, m1, n1, k1]= get_local_stk(stk, seed, sd+2);
[Ml Nl Kl] = size( local_stk );
vc.origin = double( [m1 n1 k1] );

% binarized image stack
mkl = (local_stk > T);

% local connectivity label
Ll = bwlabeln(mkl);
% get the connectivity component index
if mkl( seed(1)-m1+1, seed(2)-n1+1, seed(3)-k1+1 )
    cci = Ll( seed(1)-m1+1, seed(2)-n1+1, seed(3)-k1+1 );
else
%     disp('the seed is in background !');
    vc = [];
    return;
end
% eliminate other component
Ll = (Ll==cci);

% eliminate the traced voxels
Ll = Ll & ~mk_stk ( m1:(m1+Ml-1), n1:(n1+Nl-1), k1:(k1+Kl-1) );

% local seed
ls= seed(1:3) - [m1 n1 k1] + 1;
mls = double(ls(1));    nls = double( ls(2) );    kls = double( ls(3) );

% tic, disp('the time of distance evalueation:    ')
% neighboring voxels
nv = zeros( Ml,Nl,Kl, 'uint8' );
nv = logical( nv );

% the searching range of outside boundary
mo1 = uint16(max( 1, mls-sd-2 ));   mo2 = uint16(min( Ml, mls+sd+2 ));
no1 = uint16(max( 1, nls-sd-2 ));   no2 = uint16(min( Nl, nls+sd+2 ));
ko1 = uint16(max( 1, kls-sd-2 ));   ko2 = uint16(min( Kl, kls+sd+2 ));
% the searching range of inner boundary
mi1 = uint16(max( 1, mls-sd/1.414));    mi2 = uint16(min( Ml, mls+sd/1.414 ));
ni1 = uint16(max( 1, nls-sd/1.414));    ni2 = uint16(min( Nl, nls+sd/1.414 ));
ki1 = uint16(max( 1, kls-sd/1.414));    ki2 = uint16(min( Kl, kls+sd/1.414 ));

% matrix implemention
% Ll( mi1:mi2, ni1:ni2, ki1:ki2 ) = 0;
% Ll(1:mo1, :,:) = 0; Ll(:, 1:no1,:) = 0; Ll(:,:, 1:ko1) = 0;
% L1(mo2:Ml,:,:) = 0; Ll(:,no2:Nl,:) = 0; Ll(:,:,ko2:Kl) = 0;   
% idx = find( Ll );
% for ni = 1 : length(idx)
%     [m n k] = ind2sub( [Ml Nl Kl], idx(ni) );
%     d2 = (m-mls)*(m-mls) + (n-nls)*(n-nls) + (k-kls)*(k-kls);
%     if d2 < sd*sd
%         % in scooping distance
%         % mark voxel as visited
%         mk_stk(m+m1-1, n+n1-1, k+k1-1) = 1;
%     elseif d2 < (sd+2)*(sd+2)
%         % in neiboring zoon
%         nv(m,n,k) = 1;
%         % mark voxel as visited
%         mk_stk(m+m1-1, n+n1-1, k+k1-1) = 1;
%     end
% end

% lookup table implementation, holp to be the fastest !
cubesize = 2*cuberadius+1;
if Ml==cubesize && Nl==cubesize && Kl==cubesize
    % inside stack
    % range
    r1 = uint16(sd);
    r2 = uint16(sd+2);
    % template
    tplt1 = balls{r1};
    tplt2 = balls{r2};
    
    % the searching range of inner boundary
    mw1 = uint16(max( 1, mls-sd));    mw2 = uint16(min( Ml, mls+sd ));
    nw1 = uint16(max( 1, nls-sd));    nw2 = uint16(min( Nl, nls+sd ));
    kw1 = uint16(max( 1, kls-sd));    kw2 = uint16(min( Kl, kls+sd ));
    
    % get visited forground voxels
    try
        visited = tplt2 & Ll(mo1:mo2, no1:no2, ko1:ko2);
    catch
        disp('error!');
    end
    % mark as visited
    mk_stk( mo1+m1-1:mo2+m1-1, no1+n1-1:no2+n1-1, ko1+k1-1:ko2+k1-1 ) = ...
        visited | mk_stk( mo1+m1-1:mo2+m1-1, no1+n1-1:no2+n1-1, ko1+k1-1:ko2+k1-1 );
    
    % get outside voxels
    % eliminate inside voxels
    Ll(mw1:mw2, nw1:nw2, kw1:kw2) = ~tplt1 & Ll(mw1:mw2, nw1:nw2, kw1:kw2);
    % eliminate outside voxels
    Ll(mo1:mo2, no1:no2, ko1:ko2) = tplt2 & Ll(mo1:mo2, no1:no2, ko1:ko2);
    % eliminate the voxels outside of big box
    Ll(1:mo1, :,:) = 0; Ll(:, 1:no1,:) = 0; Ll(:,:, 1:ko1) = 0;
    L1(mo2:Ml,:,:) = 0; Ll(:,no2:Nl,:) = 0; Ll(:,:,ko2:Kl) = 0;
    nv = Ll;
else
    % near the boundary
    % remove the far away voxels
    % for implementation is faster than matrix implementation !
    for m = 1:Ml % [mo1 : mo2]
        for n = 1:Nl %[no1 : no2]
            for k = 1 : Kl %[ko1 : ko2]
                if Ll(m,n,k)
                    if m>mi1 && m<mi2 && n>ni1 && n<ni2 && k>ki1 && k<ki2
                        mk_stk(m+m1-1, n+n1-1, k+k1-1) = 1;
                        continue;
                    end
                    d2 = (m-mls)*(m-mls) + (n-nls)*(n-nls) + (k-kls)*(k-kls);
                    if d2 < (sd)*(sd)
                        % in scooping distance
                        % mark voxel as visited
                        mk_stk(m+m1-1, n+n1-1, k+k1-1) = 1;
                    elseif d2 < (sd+2)*(sd+2)
                        % in neiboring zoon
                        nv(m,n,k) = 1;
                        % mark voxel as visited
                        mk_stk(m+m1-1, n+n1-1, k+k1-1) = 1;
                    end
                end
            end
        end
    end
end
vc.neighboringvoxels = nv;
return;

%% create voxel cluster from seed voxel
function [vc mk_stk] = create_voxel_cluster_V1(stk, mk_stk, seed, sd )
% global stk;
% global mk_stk;

% adjust the scooping distance
% sd = min( sd, seed(4) );

seed = uint16(seed);
% get local image stack 
[local_stk, m1, n1, k1]= get_local_stk(stk, seed, 1.5*sd);
[Ml Nl Kl] = size( local_stk );
vc.origin = double( [m1 n1 k1] );

% get threshold in a big local stack
T = kmeans_binarize( local_stk );
vc.threshold = T;
% evaluate the intensity of seed voxel
if stk( seed(1), seed(2), seed(3) ) < T
    vc = [];
    return;
end
% binarized image stack
mkl = (local_stk > T);

% local connectivity label
Ll = bwlabeln(mkl);
% get the connectivity component index
if mkl( seed(1)-m1+1, seed(2)-n1+1, seed(3)-k1+1 )
    cci = Ll( seed(1)-m1+1, seed(2)-n1+1, seed(3)-k1+1 );
else
    disp('the seed is in background !');
    vc = [];
    return;
end
% eliminate other component
Ll = (Ll==cci);

% eliminate the traced voxels
Ll = Ll & ~mk_stk ( m1:(m1+Ml-1), n1:(n1+Nl-1), k1:(k1+Kl-1) );

% local seed
ls= seed(1:3) - [m1 n1 k1] + 1;
mls = double(ls(1));    nls = double( ls(2) );    kls = double( ls(3) );

% neighboring voxels
nv = zeros( Ml,Nl,Kl, 'uint8' );
% remove the far away voxels
for m = 1 : Ml
    for n = 1 : Nl
        for k = 1 : Kl
            if Ll(m,n,k)
                d2 = (m-mls)*(m-mls) + (n-nls)*(n-nls) + (k-kls)*(k-kls);
                if d2 < sd*sd
                    % in scooping distance 
                    % mark voxel as visited
                    mk_stk(m+m1-1, n+n1-1, k+k1-1) = 1;
                elseif d2 < (sd+2)*(sd+2)
                    % in neiboring zoon
                    nv(m,n,k) = 1;
                    % mark voxel as visited
                    mk_stk(m+m1-1, n+n1-1, k+k1-1) = 1;
                end
            end
        end
    end
end
vc.neighboringvoxels = nv;
return;

%% get the local cube
function [local_stk, m1, n1, k1] = get_local_stk(stk, centroid, r )
% global stk;
[M N K] = size(stk);

% the coordinate and radius
ms = centroid(1);   ns = centroid(2);   ks = centroid(3);   

m1 = uint16( max( 1, ms-r ) );   m2 = uint16( min(M, ms+r) );
n1 = uint16( max( 1, ns-r ) );   n2 = uint16( min(N, ns+r) );
k1 = uint16( max( 1, ks-r ) );   k2 = uint16( min(K, ks+r) );

local_stk = stk(m1:m2, n1:n2, k1:k2);
return;

%% relocate the center by mean shift
% function c = recenter(c, r, P)
% % P = 1.2;
% % get the local stack
% [local_stk, m1, n1, k1] = get_local_stk( c, P*r );
% [Ml, Nl, Kl] = size( local_stk );
% 
% % compute the center of local stack
% si = sum( local_stk(:) );
% % % estimate the weight of each intensity of 0-255
% % wi_local_stk = double(local_stk) ./ si;
% c = double([0 0 0]);
% local_stk = double( local_stk );
% 
% for m = 1 : Ml
%     for n = 1 : Nl
%         for k = 1 : Kl
%             c(1) = c(1) + double(m)*local_stk(m,n,k);
%             c(2) = c(2) + double(n)*local_stk(m,n,k);
%             c(3) = c(3) + double(k)*local_stk(m,n,k);
%         end
%     end
% end
% 
% c = c./si + double( [ (m1-1), (n1-1), (k1-1) ] );
% return;