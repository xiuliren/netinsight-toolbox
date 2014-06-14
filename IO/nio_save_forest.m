% save multiple trees as a single swc file
% 
% nio_save_forest( forest, filename )
% --------------------------------------
% 
% Input
% -----
% - forest:: a cell containing multiple trees
% - filename:: the path and name of output swc file
%
% Output
% -----
% - a swc file containing multiple root points
% 
% Example
% -------
% nio_save_forest( forest, filename )
%
% Uses 

function nio_save_forest( forest, filename )

%% load sample trees for debug only
% clc 
% clear
% forest{1} = load_tree('sample2.mtr');
% forest{2} = load_tree('sample2.mtr');
% filename = 'forest.swc';

%% open file and write some comments
fid = fopen( filename, 'w' );
fwrite(fid, ['# created by NetInsight toolbox' char(13) char(10)] );
%% save as single swc file
idx = 0;
for k = 1 : length( forest )
    tree = forest{ k };
    if isempty(tree)
        continue;
    end
    start_index = idx;
    for n = 1 : length( tree.X )
        % write every single point 
        idx = idx + 1;
        % the index of parent point, may be multiple parent points !
        idx_p = start_index + find( tree.dA(n,:) );
        if isempty(idx_p)
            % the root point
            idx_p = -1;
        elseif length(idx_p) > 1
            % multiple parent points
            % only write the first parent point, this operation may create gap
            idx_p = idx_p(1);
        end
        
        fwrite(fid, [ num2str(idx), ' ', num2str(tree.R(n)), ' ', num2str(tree.X(n)), ' ',...
             num2str(tree.Y(n)), ' ', num2str(tree.Z(n)), ' ', num2str(tree.D(n)/2), ' ',...
             num2str(idx_p), char(13) char(10) ]);
    end
end

fclose(fid);