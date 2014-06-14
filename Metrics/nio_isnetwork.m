% distinguish networks
% 
%  bool_network = nio_isnetwork( input )
% --------------------------------------
% 
% Input
% -----
% - network:: the input structure which may not be formated of network
%
% Output
% -----
% - isnetwork:: boolen type, whether the input structure is a network
% 
% Example
% -------
% 
%
% Uses 

function isnetwork = nio_isnetwork( network )

%% load test sample for debug only
% clc
% clear
% tree = load_tree( 'sample2.mtr' );
% network = nio_tree2network( tree );

%% whether the input structure is network formated
isnetwork =  isfield( network, 'dAe' );