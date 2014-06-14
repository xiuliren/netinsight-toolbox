% nio_get_sps_eps
% get the start and end points of every section
% 
% [ sps eps ] = nio_get_sps_eps( network )
% --------------------------------------
% 
% Input
% -----
% - network:: the network contains several sections 
%
% Output
% -----
% - sps:: the start points of every section
% - eps:: the end points of every section
%
% Example
% -------
% 
% 
% Uses 

function [ sps eps ] = nio_get_sps_eps( network )
% get the start and end points
sps = [];
eps = [];
for k = 1 : network.sn
    sec = network.sections{k};
    sps = [ sps; sec(1, :) ];
    eps = [ eps; sec(end, :) ];
end