% nio_build_net_connectivity
% build the internal connectivity of this network
% 
% network = nio_build_net_connectivity( network )
% --------------------------------------
% 
% Input
% -----
% - network:: the network that is lack of connectivity
%
% Output
% -----
% - net_unit:: the extracted vascular unit
%
% Example
% -------
% 
% 
% Uses 

function network = nio_build_net_connectivity( network )

% initiate
network.sn = length( network.sections );
sn = network.sn;
network.dAs = sparse( sn, sn );
network.dAe = sparse( sn, sn );

% get the start and end points of every section
[ sps eps ] = nio_get_sps_eps( network );

% connect from start point to start point
for k = 1 : sn
    sp = sps( k, : );
    % find the connectivity relasionship
    con_idx = find( sp(1)==sps(:,1) & sp(2)==sps(:,2) & sp(3)==sps(:,3) );
    % eliminate the self connectivity
    con_idx( con_idx==k ) = [];
    % establish the connectivity
    network.dAs( k, con_idx ) = 1;
    
    ep = eps( k, : );
    % find the connectivity relasionship
    con_idx = find( ep(1)==eps(:,1) & ep(2)==eps(:,2) & ep(3)==eps(:,3) );
    % eliminate the self connectivity
    con_idx( con_idx==k ) = [];
    % establish the connectivity
    network.dAe( k, con_idx ) = 1;
end