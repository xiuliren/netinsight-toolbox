% nio_remove_rb
% Remove redundant(unreal) branches in vessel tree.
% 
% [] = nio_remove_rb ( SrcSWC, LR, DstSWC )
% -----------------------------------------
% 
% Input
% -----
% - SrcSWC:: path of the swc file
% - LR:: length/R ratio, length and radius of a segment
% - DstSWC:: output path of swc file with redundant branches removed
% 
% Example
% -------
% nio_remove_rb ( 'original_data.swc', 2,  'modified_data.swc');
% 
% Uses nio_load_tree idpar_tree B_tree T_tree delete_tree swc_tree

function [ ] = nio_remove_rb ( SrcSWC, LR, DstSWC )
%% read source swc file
[ vessel ~ ] = nio_load_tree( SrcSWC );
vessel_rrb = vessel;
%% detect the terminal branches
for k = 1 : length( vessel )
    % parent node index
    idpar = idpar_tree( vessel{k} );
    
    % get branch and termination points
    B = B_tree( vessel{k} );
%     b_list = find(B);
    
    T = T_tree( vessel{k} );
    t_list = find( T );
    
    for t = 1 : length( t_list )
        % get the termination branch points list
        tb = [];
        len = 0;
        id = t_list(t);
        while( 1 )
            % the next point
            idp = idpar( id );
            tb = [ tb id ]; 
            
            % calculate the length
            x1 = vessel{k}.X( id );
            y1 = vessel{k}.Y( id );
            z1 = vessel{k}.Z( id );
            
            x2 = vessel{k}.X( idp );
            y2 = vessel{k}.Y( idp );
            z2 = vessel{k}.Z( idp );
            
            len = len + sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1) );
            % the next point
            id = idp;
            % encounter the branck point, then break
            if ( B( id ) )
                break;
            end
        end
        if len / ( vessel{k}.D(id) / 2 ) < LR
            % delete the short tree
            delete_tree(vessel_rrb{k}, tb);
        end
    end    
end

%% combine multiple trees to single dA, then write the swc file
% number of points
N = 0;
for k = 1 : length( vessel_rrb )
    N = N + length( vessel_rrb{k}.X );
end

% initiate the destination tree matrix
vessel_dst = {};
vessel_dst.dA = sparse(N,N);
vessel_dst.X = zeros(N,1);
vessel_dst.Y = zeros(N,1);
vessel_dst.Z = zeros(N,1);
vessel_dst.R = zeros(N,1);
vessel_dst.D = zeros(N,1);
vessel_dst.rnames = {'1' 'dendrite'};
vessel_dst.Ri = 100;
vessel_dst.Gm = 0;
vessel_dst.Cm = 0;
vessel_dst.name = 'vessels';

% combine multiple trees
n1 = 1;
n2 = 0;
for k = 1 : length( vessel_rrb )
    % the modification range
    Nk = length( vessel_rrb{k}.X );
    n2 = n2 + Nk;
    
    % modify the matrix
    vessel_dst.X( n1 : n2 ) = vessel_rrb{k}.X(:);
    vessel_dst.Y( n1 : n2 ) = vessel_rrb{k}.Y(:);
    vessel_dst.Z( n1 : n2 ) = vessel_rrb{k}.Z(:);
    vessel_dst.R( n1 : n2 ) = vessel_rrb{k}.R(:);
    vessel_dst.D( n1 : n2 ) = vessel_rrb{k}.D(:);
    
    % modify the adjacency matrix
    [ r, c ] = find( vessel_rrb{k}.dA );
    vessel_dst.dA( r + n1 -1, c + n1 - 1 ) = 1;
end

% write the whole vessel
swc_tree( vessel_dst, DstSWC );
end

