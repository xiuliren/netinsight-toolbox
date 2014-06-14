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

function network = nio_bridge_gap_net( network, bridge_len )

%% load network to debug
% clc
% clear
% % load( 'network.mat' );
% tree = load_tree( 'sample2.mtr' );
% network = nio_tree2network( tree );
% bridge_len = 10;

%% parameters
% detection range of bridging, the default length
if nargin == 1
    bridge_len = 20;
end

%% bridge the sections by creating a new section
% build an image stack with expanded diameter
stk = nio_network2stk( network, 3 );
[ Ms Ns Ks ] = size( stk );
% the terminal sections
idx_t = find( network.terminalSec );

for idx_sec = 1 : length(idx_t)
    t_sec = network.sections{ idx_sec };
    % the unit vector for step
    len = sqrt( ( t_sec(2,1)-t_sec(1,1))*( t_sec(2,1)-t_sec(1,1)) + ...
        (t_sec(2,2)-t_sec(1,2))*( t_sec(2,2)-t_sec(1,2)) + ...
        (t_sec(2,3)-t_sec(1,3))*( t_sec(2,3)-t_sec(1,3)) );
    if len == 0
        break;
    end
    ux = (t_sec(1,1) - t_sec(2,1))./len;
    uy = (t_sec(1,2) - t_sec(2,2))./len;
    uz = (t_sec(1,3) - t_sec(2,3))./len;
    
    % step for detection
    bl = max( [ t_sec(1,4) bridge_len ] );
    x0 = t_sec(1,1);
    y0 = t_sec(1,2);
    z0 = t_sec(1,3);
    for k = 1 : bl
        % step forward
        x = round( x0 + k*ux );
        y = round( y0 + k*uy );
        z = round( z0 + k*uz );
        
        if (x<1) | (y<1) | (z<1) | (x>Ms) | (y>Ns) | (z>Ks)
            % out of range
            break;
        end
        if stk( x,y,z) ~= 0
            disp('find a gap! congratulations !!')
            % find the connect section !
            con_sec_idx = stk( x,y,z );
            sec = network.sections{ con_sec_idx };
            % find the connect point
            xd = abs( sec(:,1) - x );
            yd = abs( sec(:,2) - y );
            zd = abs( sec(:,3) - z );
            dis = xd + yd + zd;
            p_idx = find( dis == min(dis) );
            p_idx = p_idx(1);
            
            if p_idx < 3
                % the starts points, connect the sections directly
                network.dAs( idx_sec, p_idx ) = 1;
            elseif  p_idx > length(sec) - 3
                % the end points, connect the sections directly
                network.dAe( idx_sec, p_idx ) = 1;
            else
                % in the center of this section
                % break this section to two sections and connect
                network.sn = network.sn + 1;
                network.sections{ end } = sec( p_idx:end );
                network.sections{ con_sec_idx } = sec( 1:p_idx );
                
                % the sections that connect to the end of original section
                con_idx_end = find( network.dAe(:, con_sec_idx) );
                % now, connect to the new section
                network.dAe( :, end ) = network.dAe( :, con_sec_idx );
                network.dAe( :, con_sec_idx ) = 0;
                % connect the new section to the broken section
                network.dAe( end, con_sec_idx ) = 1;
                
                % connect the terminal section to the two new sections
                % this gap was bridged !!!
                network.dAe( idx_sec, con_sec_idx ) = 1;
                network.dAs( idx_sec, end ) = 1;
            end
            % finish searching, break
            break;
        end
    end
end
