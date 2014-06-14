% distinguish tree structure
% 
%  istree = nio_istree( input )
% --------------------------------------
% 
% Input
% -----
% - tree:: the input structure which may not be formated of tree
%
% Output
% -----
% - isntree:: boolen type, whether the input structure is a tree
% 
% Example
% -------
% 
%
% Uses 

function isTree = nio_istree( tree )

%% load test sample for debug only
% clc
% clear
% tree = load_tree( 'sample2.mtr' );
% network = nio_tree2network( tree );

%% whether the input structure is network formated
return isfield( tree, 'dA' );