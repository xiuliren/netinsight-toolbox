% save trees or network as hoc format
% 
% nio_write_tree_hoc( inputTrees, filename )
% --------------------------------------
% 
% Input
% -----
% - inputTrees:: the vectorized trees which has a tree structure
%               note that the points may have multiple parents
% - filename:: the path and name of output hoc file 
%
% Output
% -----
% save a hoc file which can hold circles
% 
% Example
% -------
% nio_save_hoc( inputTrees, filename )
%
% Uses 
function nio_write_tree_hoc( inputTrees, filename )
% by jpwu@CBMP, 2012/09/15

%% load test data
% clc
% clear
% load matlab.mat
% % transfer to input variables
% inputTrees{1} = vessel1;
% % inputTrees{1} = load_tree('sample2.mtr');
% filename = 'gap_bridged.hoc';


%% transfer inputTree to network
if nio_isnetwork( inputTrees )
    network = inputTrees;
else
    % not network
    if nio_istree( inputTrees ) | nio_isforest( inputTrees )
        % transfer tree or forest to network
        network = nio_tree2network(inputTrees);
    else
        % not a propriate data structure
        disp(' not a propriate data structure, need tree, forest or network ');
        disp(' save failed !! I am very sorry ! ')
        return;
    end
end

nio_write_net_hoc( network, filename );
