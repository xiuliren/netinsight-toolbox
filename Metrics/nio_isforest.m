% distinguish forest structure
% 
%  isforest = nio_isforest( input )
% --------------------------------------
% 
% Input
% -----
% - forest:: the input structure which may not be formated of forest
%
% Output
% -----
% - isnforest:: boolen type, whether the input structure is a forest
% 
% Example
% -------
% 
%
% Uses 

function isforest = nio_isforest( forest )

%% load test sample for debug only
% clc
% clear
% forest = load_forest( 'sample2.mtr' );
% network = nio_forest2network( forest );

%% whether the input structure is network formated
if iscell( forest ) & ~isempty( forest )
    return isfield( forest{1}, 'dA' );
else
    return false;
end