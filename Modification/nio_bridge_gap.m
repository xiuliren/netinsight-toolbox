% nio_bridge_gap
% To bridge gaps caused by vessel tracing in segments of a vessel tree.
% 
% [ vessel ] = nio_bridge_gap( vessel, DstDir )
% --------------------------------------
% 
% Input
% -----
% - vessel:: the vectorized vessels which has a tree structure
%
% Output
% -----
% - vessel:: the vessel. one point may have multiple parent points
% Example
% -------
% vessel = nio_bridge_gap ( vessel );
% 
% Uses idpar_tree swc_tree neuron_tree

function [ vessel ] = nio_bridge_gap( vessel )
%% preparation
% the volume size
Ms = ceil(max(vessel.X));
Ns = ceil(max(vessel.Y));
Ks = ceil(max(vessel.Z));

% revise the  coordinate < 0
vessel.X( vessel.X < 1 ) = 1;
vessel.Y( vessel.Y < 1 ) = 1;

%% build the point index stack
disp('-------- getting binary stack ...')
idx_stk = zeros (Ms, Ns, Ks, 'double');

% points connectivity, child and parent
[chl par] = find(vessel.dA == 1);
for con = 1 : length( chl )
    % the child and parent point
    x1 = vessel.X( chl(con) );
    y1 = vessel.Y( chl(con) );
    z1 = vessel.Z( chl(con) );
    D1 = vessel.D( chl(con) );
    
    x2 = vessel.X( par(con) );
    y2 = vessel.Y( par(con) );
    z2 = vessel.Z( par(con) );
    D2 = vessel.D( par(con) );
    
    % point distance
    dis = sqrt( (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2) );
    if dis < 1
        continue;
    end
    
    % unit vector
    uv = [ (x2-x1) (y2-y1) (z2-z1) ] ./ dis;
    
    % fill ball with point id step by step
    for k = 0 : dis
        % center
        x0 = x1 + k * uv(1);
        y0 = y1 + k * uv(2);
        z0 = z1 + k * uv(3);
        % radius
        R = ( D1 + (D2-D1) * k / dis ) / 2;
        
        % the filling index of point
        if k < dis / 2
            id = chl(con);
        else
            id = par(con);
        end
        
        % fill in the index
        for x = floor(x0 - R) : ceil(x0 + R)
            if (x<1) || (x>Ms)
                continue;
            end
            for y = floor(y0 - R) : ceil(y0 + R)
                if (y<1) || (y>Ns)
                    continue;
                end
                for z = floor(z0 - R) : ceil(z0 + R)
                    if (z<1) || (z>Ks)
                        continue;
                    end
                    if idx_stk(x,y,z) > 0
                        continue;
                    end
                    % check distance
                    d_squr = (x - x0) * (x - x0) +  (y - y0) * (y - y0) + (z - z0) * (z - z0) ;
                    if d_squr <= R * R
                        idx_stk(x,y,z) = id;
                    end
                end
            end
        end
    end
end

%% detect gaps
gaps = [];
for n = 1 : length( vessel.X )
    if (sum( vessel.dA(:,n) ) + sum( vessel.dA(n, :) )) < 2
        gaps = [ gaps n ];
    end
end

%% bridge gaps
idpar = idpar_tree(vessel , '-0');
for k = 1 : length( gaps )
    xg = vessel.X( gaps(k) );
    yg = vessel.Y( gaps(k) );
    zg = vessel.Z( gaps(k) );
    rg = vessel.D( gaps(k) ) / 2;
    pg = idpar( gaps(k) );
    
    % neighber points
    neighbers = [ find( vessel.dA(gaps(k),:) ==1 )  find( vessel.dA(:, gaps(k)) ==1 ) ] ;
    
    % voxel scooping
    vs = [];
    % scooping distance
    sd = 2 * rg;
    for x = floor(xg - sd) : ceil(xg + sd)
        if (x<1) || (x>Ms)
            continue;
        end
        for y = floor(yg - sd) : ceil(yg + sd)
            if (y<1) || (y>Ns)
                continue;
            end
            for z = floor(zg - sd) : ceil(zg + sd)
                if (z<1) || (z>Ks)
                    continue;
                end
                if idx_stk(x,y,z) == 0
                    continue;
                end
                % check distance
                d_squr = (x - xg) * (x - xg) +  (y - yg) * (y - yg) + (z - zg) * (z - zg) ;
                if ( d_squr <= sd * sd) && ( d_squr > rg * rg )
                    % the distance between radius and scooping distance, scoope the voxel
                    vs = [ vs; x y z ];
                end
            end
        end
    end
    
    % find the connecting point by analysis the voxels
    % get point index list
    p_list = [];
    for m = 1 : size( vs, 1)
        idx = idx_stk( vs(m, 1), vs(m, 2), vs(m, 3)  );
        if ( idx ~= gaps(k) ) && ( isempty( find (idx == neighbers, 1 ) ) )
            % add the point index
            p_list = [ p_list idx ];
        end
    end
    
   if ~isempty( p_list )
        % get the most frequent index as the connecting point
        mfi = mode( p_list );
        % connect the points by modify the dA matrix
        vessel.dA( gaps(k), mfi ) = 1;
    end
    
    
end

%% save network as hoc or swc format
% neuron_tree( vessel, [DstDir '.nrn'] );
% swc_tree( vessel, [DstDir '.swc'] );

end

