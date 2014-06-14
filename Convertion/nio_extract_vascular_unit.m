% nio_bridge_gap
% To bridge gaps caused by vessel tracing in segments of a vessel tree.
% 
% net_unit = nio_extract_unit( network )
% --------------------------------------
% 
% Input
% -----
% - network:: the interconnected network
% - kpt:: the key point in a section that supply blood of this capillary bed
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

function net_unit = nio_extract_vascular_unit( network, kpt )

%% test code for debug only
% clc
% kpt = kpt1;

%% parameters
% the diameter threshold of capillaries
Td = 4;

%% find the section index
ksec = [];
dis = zeros(1, network.sn);
for k = 1 : network.sn
    sec = network.sections{k};
%     if ~isempty( find(sec(:,1)==kpt(2)) ) & ...
%             ~isempty( find(sec(:,2)==kpt(1)) ) & ...
%             ~isempty( find(sec(:,3)==kpt(3)) )
%         % the key section was found
%         ksec = [ ksec k ];
%     end
    d_v = (sec(:,1)-kpt(1)).*(sec(:,1)-kpt(1)) + ...
        (sec(:,2)-kpt(2)).*(sec(:,2)-kpt(2)) + ...    
        (sec(:,3)-kpt(3)).*(sec(:,3)-kpt(3));
    dis(k) = min( d_v );
end
ksec = find( dis == min(dis) );

if isempty( ksec )
    disp( 'the key section was not found ! please check the coordinate of key point!' );
    return;
end

%% compute the average diameter of every section
network.avd = zeros( 1, network.sn, 'double' );
for k = 1 : network.sn
    sec = network.sections{k};
    network.avd(k) = mean( sec(:,4) );
end
if length( ksec ) > 1
    % get the section with minimun diameter
    disp( 'find more than one matched sections !!!' );
    d_v = network.avd( ksec );
    idx = find( d_v == min( d_v ) );
    ksec = ksec( idx );
end

%% start extraction
disp('----- start extracting ......');

% get the start and end points
[ sps eps ] = nio_get_sps_eps( network );

% mark the checked sections to avoid repeat
checked_secs_bin = zeros( network.sn, 1 );
checked_secs_bin( ksec ) = 1;
con_secs_idx = ksec;

% the index of terminated sections
tm_secs_idx = [];
% corrent section diameter
cd = network.avd( ksec );

% find smaller connected sections
[ con_secs_idx tm_secs_idx checked_secs_bin ] = nio_find_small_secs( network, ksec, con_secs_idx, tm_secs_idx, checked_secs_bin, Td );
% find the bigger connected sections
for k = 1 : length( tm_secs_idx )
    key_tm = tm_secs_idx( k );
    [ con_secs_idx tm_secs_idx checked_secs_bin ] = nio_find_big_secs( network, key_tm, con_secs_idx, tm_secs_idx, checked_secs_bin, network.avd(ksec)*1.5 );
end

% add the connected sections to network
net_unit = nio_new_network();
for k = 1 : length( con_secs_idx )
    idx = con_secs_idx(k);
    net_unit.sections = [ net_unit.sections; network.sections( idx ) ];
end

%% estimate the connectivity relationship of this new network
disp('estimate the connectivity relationship of this new network!');
net_unit.sn = length( net_unit.sections );
net_unit = nio_build_net_connectivity( net_unit ); 

