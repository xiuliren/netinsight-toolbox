%% get node list from the label
function [node_list ] = get_node_list(LABELS, PAIRWISE, map_id, seed)
% minimum number of voxels of each component
Cmin = 3;   % Cmin >= 3, make sure that the radius index is >=1

% the forground index
fidx = find( LABELS );
Nf = length( fidx );
% build the foreground graph adjacency matrix
fA = sparse( Nf, Nf );
for nf = 1 : Nf
    fA(nf,:) = PAIRWISE( fidx(nf), fidx );
end
[Nb, C] = graphconncomp(fA, 'Directed', false );

node_list = [];

% mark the component in LABELS
LABELS( fidx ) = C;
% compute the center and radius of each node
for nb = 1 : Nb
    % component index
    cidx = find( LABELS==nb );
    if length(cidx) <= Cmin
        % less than Cmin voxels
        continue;
    end
    
    % vector of the coordinates
    vcoord = map_id(cidx, 1:3);
    center = mean( vcoord, 1 );
    
    % compute the distance from center to every component voxel
    vcoord = double(vcoord);
    dis2 = (vcoord(:,1) - center(1)).*(vcoord(:,1) - center(1)) + ...
        (vcoord(:,2) - center(2)).*(vcoord(:,2) - center(2)) + ...
        (vcoord(:,3) - center(3)).*(vcoord(:,3) - center(3));
    % the radius was estimated using the largest distance
    dis2 = sort(dis2, 'descend');
    % the position was modeled by a circute, Ndis = pi*r*r, 
    %   Nc = pi*r*r - pi*(r-1)*(r-1)
    ndis = ceil( ( 2*sqrt(length( dis2 )/pi) - pi )/2) + 1;
    r = sqrt( dis2( ndis ) );
    
    % adjust the coordinate according to the algorithm of voxel scooping
%     ex = double(min( r, seed(4) ) / max( r, seed(4) ));
%     center = uint16( double(seed(1:3)) + ...
%         (0.5^ex)*( double(center) - double(seed(1:3)) ));
    center = (double(seed(1:3)) + center)/2;
    
    % add the element
    node_list = [node_list; center r];
end