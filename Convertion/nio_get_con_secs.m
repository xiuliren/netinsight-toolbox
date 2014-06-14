% nio_get_con_secs
% get the connected secion index of input section
% 
% network = nio_get_con_sections( network, idx )
% --------------------------------------
% 
% Input
% -----
% - network:: the network 
% - idx:: the index of key section
%
% Output
% -----
% - con_secs_idx:: the index of connected neighbor sections
%
% Example
% -------
% 
% 
% Uses 

function con_secs_idx = nio_get_con_secs( network, idx )

[ sps eps ] = nio_get_sps_eps( network );

% the current section
sec = network.sections{ idx };
sp = sec( 1, : );
ep = sec( end, : );

con_secs_idx = find( ( sp(1)==sps(:,1) & sp(2)==sps(:,2) & sp(3)==sps(:,3) ) ...
            | ( sp(1)==eps(:,1) & sp(2)==eps(:,2) & sp(3)==eps(:,3) ) ...
            | ( ep(1)==sps(:,1) & ep(2)==sps(:,2) & ep(3)==sps(:,3) ) ...
            | ( ep(1)==eps(:,1) & ep(2)==eps(:,2) & ep(3)==eps(:,3) ) );