%% wave propagation from a single seed
% by jpwu, 2013/02/28

function network = tracing_propagation( seed, network )
global stk;
global mk_stk;
%% parameters
% the radius threshold for selecting threshold tracing method
RThresholdTracing = 0;
% the radius threshold for sphere model of graph cut
RGruphCutSphereTracing = 0;
% the radius threshold for line scan
RLineScan = 1;

% counter interval for sorting seed list
counterITV = 200;

% the length of the scan line
ScanLineLen = 10;

%% traced seed
% if mk_stk( uint32(seed(1)), uint32(seed(2)), uint32(seed(3)) )
%     return;
% end
%% propagation processing
% section number of current network
sn = length( network.sections );
% add the original seed in a new section
sn = sn + 1;
network.sections{sn} = seed;

% dynamic seeds list, m, n, k, r, section index, flag of line scan
% dsl = zeros(10000000, 6);
global dsl;
dsl(1,:) = [seed sn 1];
pdsl = 1;

% iterative tracing
counter = uint8(0);
while 0 ~= pdsl
    pdsl
    counter = counter + 1;
    % sort the seed list, ensures that the big seed was traced first
    if mod(counter, counterITV) == 0
        counter = 0;
        dsl(1:pdsl, :) = sortrows( dsl(1:pdsl,:), 4 );
    end
    % get the last seed and current section number
    seed = round( dsl( pdsl, 1:4 ) );
    currentSecIdx = dsl(pdsl, 5);
    flagLineScan = dsl(pdsl, 6);
    % delet from the dynamic list
    dsl(pdsl, :) = [0 0 0 0 0 0];
    pdsl = pdsl - 1;
    
    % get the local voxel cube
    [local_stk, local_seed, gm1, gn1, gk1, gm2, gn2, gk2, ...
        lm1, ln1,lk1, lm2, ln2, lk2]= get_local_stk(stk, seed, 0);
    
    % binarize the local voxel cube
    [T, h]= kmeans_binarize( local_stk );
    
    % get the local marker stack, the pad value was set to be true
    mk_local_stk = true( size( local_stk ) );
    mk_local_stk(lm1:lm2,ln1:ln2, lk1:lk2) = mk_stk(gm1:gm2,gn1:gn2,gk1:gk2);
    
    % select tracing method according to the radius of the seed 
    node_list = [];
    if seed(4) > RThresholdTracing
%         disp('using threshold tracing!')
        [node_list mk_local_stk] = local_fast_threshold_V2( local_seed, local_stk, mk_local_stk, T );
    elseif seed(4) > RGruphCutSphereTracing
        %disp('using graph cut tracing of sphere model!');
         % get the new node via graph cut method, key step
         [node_list mk_local_stk] = local_graph_cut_sphere( local_seed, local_stk, mk_local_stk, T );
    elseif seed(4) > RLineScan
        %disp('using graph cut tracing of cube model!');
        [node_list mk_local_stk] = local_graph_cut_cube( local_seed, local_stk, mk_local_stk, T );
    end
    
    % when the above methods is not effective, use the line scan method
    if isempty( node_list ) && flagLineScan %|| size(node_list, 1) > 3
        % mark that line scan has been used in this section tracing process
        flagLineScan = 0;
        % retracing using line scan
        %disp('using line scan V2 !')
        seed2 = seed;   seed2(4) = ScanLineLen;
        % get the local voxel cube
        [local_stk, local_seed, gm1, gn1, gk1, gm2, gn2, gk2, ...
            lm1, ln1,lk1, lm2, ln2, lk2]= get_local_stk(stk, seed2, 0);

        % get the local marker stack, the pad value was set to be true
        mk_local_stk = true( size( local_stk ) );
        mk_local_stk(lm1:lm2,ln1:ln2, lk1:lk2) = mk_stk(gm1:gm2,gn1:gn2,gk1:gk2);
        
        % binarize the local voxel cube
        [T, h]= kmeans_binarize( local_stk );
        [node_list mk_local_stk] = local_line_scan_V2...
            ( local_seed, local_stk, mk_local_stk, T );
    end   
    
    if isempty( node_list )
        continue;
    end
    % mark the visited voxels
    mk_stk(gm1:gm2, gn1:gn2, gk1:gk2) = mk_stk(gm1:gm2, gn1:gn2, gk1:gk2) ...
        | mk_local_stk(lm1:lm2,ln1:ln2,lk1:lk2);
    % transform the coordinate to global system
    node_list = double(node_list) + ...
        repmat(double([gm1-1 gn1-1 gk1-1 0]), size(node_list,1), 1);

%     % the radius of rayburst
%     radiusRayburst = get_radius_V3( stk, node_list(1:3), T );
%     node_list(4) = (node_list(4) + radiusRayburst) / 2;
    
%     node_list
    
    % sort the rows in ascending order according to the radius.
    % this operation insures that the big branches were traced first
%     node_list = sortrows( node_list, -4 );
       
    %% add the node list to the network
    Nb = size( node_list, 1 );
    if Nb < 1
        disp('end of propagation!')
    elseif Nb == 1
        network.sections{currentSecIdx} = ...
            [network.sections{currentSecIdx}; node_list ];
        % update the dynamic list
        pdsl = pdsl + 1;
        dsl(pdsl,:) =  [ node_list currentSecIdx flagLineScan];
    else
        % have branches
        % sort the new seeds
        node_list = sortrows(node_list, 4);
        for nb = 1 : Nb 
            % set the line scan flag
            flagLineScan = 1;
            % the new nodes
            newSecIdx = length( network.sections ) + 1;
            network.sections{newSecIdx} = [seed; node_list(nb,:)];
            % update the dynamic list
            pdsl = pdsl + 1;
            dsl(pdsl,:) =  [node_list(nb,:) newSecIdx flagLineScan];
        end
    end
    
    % update the visualization system
    for nb = 1 : size(node_list,1)
        plot( node_list(nb,2)-1, node_list(nb,1)-1,...
            '.', 'color',[rand, rand, rand], 'markersize', 3*node_list(nb,4) );
            plot( [seed(2)-1 node_list(nb,2)-1], [seed(1)-1 node_list(nb,1)-1],'-r','LineWidth',2 );
        pause(0.01); % draw this figure
    end
end
return;

%% get the local cube
function [local_stk, local_seed, gm1, gn1, gk1, gm2, gn2, gk2, lm1, ln1,lk1, lm2, ln2, lk2] = get_local_stk( stk, seed, padValue )
% the size scale of local cube
P = 1.5;   % Pl > 1
% minimum and maximum spaning range
MINSR = 3;  MAXSR = 10;

[M, N, K] = size(stk);
% the coordinate and radius
ms = seed(1);   ns = seed(2);
ks = seed(3);   rs = seed(4);

% the new redius of local stack
rl = uint16(P*rs);
rl = max( rl, rs+MINSR);  rl = min( rl, rs+MAXSR );

% the coordinate range in the global stack
gm1 = uint16(max( 1, ms-rl )); gm2 = uint16(min(M, ms+rl));
gn1 = uint16(max( 1, ns-rl )); gn2 = uint16(min(N, ns+rl));
gk1 = uint16(max( 1, ks-rl )); gk2 = uint16(min(K, ks+rl));

% get the original local stack
local_stk_tmp = stk(gm1:gm2, gn1:gn2, gk1:gk2);
[Mlt,Nlt,Klt] = size(local_stk_tmp);

% the distance from the seed to the boundary
% the starting point
slm = ms - gm1 + 1;  sln = ns - gn1 + 1; slk = ks - gk1 + 1;

% the ideal stack size
Ml = uint16(2*rl + 1);    Nl = Ml;    Kl = Ml;

% the pad value was set according to the parameter of padValue
if padValue == 0
    local_stk = zeros(Ml,Nl,Kl,'uint8');
elseif padValue == true
    local_stk = true(Ml,Nl,Kl);
else
    disp('any other situation? ')
end
% the ideal local center of the seed
mcl = uint16(rl + 1);  ncl = mcl; kcl = mcl;

% put the stack inside the ideal one
lm1 = mcl-slm+1;    lm2 = mcl-slm+Mlt;
ln1 = ncl-sln+1;    ln2 = ncl-sln+Nlt;
lk1 = kcl-slk+1;    lk2 = kcl-slk+Klt;
local_stk( lm1:lm2, ln1:ln2, lk1:lk2 ) = local_stk_tmp; 

% update the seed
local_seed = [mcl ncl kcl seed(4)];
return;