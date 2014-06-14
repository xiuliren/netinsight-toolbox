% nio_load_tree   
% Loads tree(s) from the swc format.
% (trees package)
%
% [trees, tree, tname, path] = nio_load_tree (tname, options)
% ----------------------------------------------
%
% Loads the swc file and the corresponding directed adjacency matrix to
% create a tree in the trees structure.
% Modification: it can handle multi-roots swc file, making output to be
% a 1*N cell array, every element is tree structure corresponding to one-root 
% condition. which lengh is equal to number of roots. 
% !Note: if file have multi-roots, then every node's parent node should
% appear before itself appears in file. Othewise the output will make no
% sense.
% 
% Input
% -----
% - tname::string: name of the file to be loaded, including the extension.
%                 {DEFAULT : open gui fileselect, replaces format entry}
% - options::string {DEFAULT : '-r'}
%     '-s' : show
%     '-r' : repair tree, preparing trees for most TREES toolbox functions
%
% Output
% ------
% if no output is declared the tree is added to trees
% - trees:: structured output tree(cell array)
% - tree:: structured output tree(a single struct)
% - tname::string: name of output file; [] no file was selected -> no output
% - path::sting: path of the file, complete string is therefore: [path name]
%
% Example
% -------
% tree = nio_load_tree('1.swc');
%
% See also neuron_tree swc_tree start_trees (neu_tree.hoc)

function varargout = nio_load_tree (tname, options)

% trees : contains the tree structures in the trees package
global trees

if (nargin<1)||isempty(tname),
     [tname path] = uigetfile ({'*swc', 'TREES formats (*.swc)'}, ...
        'Pick a file',...
        'multiselect', 'off');
    if tname == 0,
        varargout {1} = []; varargout {2} = []; varargout {3} = [];
        return
    end
else
    path = '';
end
format = tname  (end - 3 : end); % input format from extension:
% extract a sensible name from the filename string:
nstart = unique ([0 strfind(tname, '/') strfind(tname, '\')]);
name   = tname  (nstart (end) + 1 : end - 4);

if (nargin<2)||isempty(options)
    if strcmp (format, 'swc') || strcmp (format, 'neu'),
        options = '-r';
    else
        options = '';
    end
end

switch format,
    case '.swc' % this is then swc
        if ~exist ([path tname], 'file'),
            error ('no such file...');
        end
        A = textread ([path tname], '%s', 'delimiter', '\n');
        swc_length = 0;
        for ward = 1 : length (A),
            if ~isempty (A {ward}), % allow empty lines in between
                if ~strcmp (A {ward} (1), '#') % allow comments: lines starting with #
                    swc_length = swc_length + 1;
                end
            end
        end
        swc = zeros(swc_length, 7);
        swc_idx = 1;
        for ward = 1 : length (A),
            if ~isempty (A {ward}), % allow empty lines in between
                if ~strcmp (A {ward} (1), '#') % allow comments: lines starting with #
                    swc(swc_idx, :) = str2num(A {ward});
                    swc_idx = swc_idx + 1;
                end
            end
        end
        N       = size (swc, 1);
        if sum (swc (:, 1) ~= (1 : N)'),   % check index in first column
            error ('index needs to be 1 .. n');
        end
        idpar   = swc (:, 7); % vector containing index to direct parent

        dA      = sparse (N, N);
        for ward = 2 : N;
           if ( ward < 1 ) || ( idpar(ward) < 1 )
                continue;
           end
            dA (ward, idpar (ward)) = 1;
        end
        
        roots_trees = find(idpar < 1); % index of root nodes
        num_trees = length(roots_trees);
        roots_trees = [roots_trees; swc_length + 1];
        roots_member = cell(1, num_trees);
        for ward = 2 : num_trees + 1
            roots_member{ward - 1} = [];
            roots_member{ward - 1} = roots_trees(ward-1) : 1 : (roots_trees(ward) - 1); % initiation
        end
        for ward = 2 : num_trees
            idx_tr = 2;
            while idx_tr <= length(roots_member{ward}) % traversing indexes
                if isempty(find(roots_member{ward} == idpar(roots_member{ward}(idx_tr)), 1)) == 1
                    for sss = 1 : ward - 1 
                        if ~isempty(find(roots_member{sss} == idpar(roots_member{ward}(idx_tr)), 1))
                            roots_member{sss} = [roots_member{sss} roots_member{ward}(idx_tr)];
                            break;
                        end
                    end
                    roots_member{ward}(idx_tr) = [];
                else
                    idx_tr = idx_tr + 1;
                end
            end
        end
        
        dA_split = cell(1, num_trees);
        swc_split = cell(1, num_trees);
        for ward = 1 : num_trees
            swc_split{ward} = swc(roots_member{ward}, :);
        end
        for ward = 1 : num_trees
            dA_split{ward} = sparse(length(roots_member{ward}), length(roots_member{ward}));
            for idx_dA_split = 1 : length(roots_member{ward})
                if (idpar(roots_member{ward}(idx_dA_split)) < 1)
                    continue;
                else
                    idx_1 = find(roots_member{ward} == idpar(roots_member{ward}(idx_dA_split)), 1);
                    dA_split{ward}(idx_dA_split, idx_1) = 1;
                end
            end
        end
                    
        trees_split = cell(1, num_trees); % output tree as cell array
        for ward = 1 : num_trees
            trees_split{ward}.dA = dA_split{ward};
            trees_split{ward}.X = swc_split{ward}( :, 3);
            trees_split{ward}.Y = swc_split{ward}( :, 4);
            trees_split{ward}.Z = swc_split{ward}( :, 5);
            trees_split{ward}.D = swc_split{ward}( :, 6) * 2;
            [i1 , i2, i3] = unique(swc_split{ward}( :, 2));
            trees_split{ward}.R  = i3;
            trees_split{ward}.rnames = cellstr (num2str (i1))';
            trees_split{ward}.name = name;
        end
        
        tree.dA = dA;
        tree.X  =            swc (:, 3);     % X-locations of nodes on tree
        tree.Y  =            swc (:, 4);     % Y-locations of nodes on tree
        tree.Z  =            swc (:, 5);     % Z-locations of nodes on tree
        tree.D  =            swc (:, 6) * 2; % local diameter values of nodes on tree
        [i1 , i2, i3] = unique (swc (:, 2));
        tree.R  = i3;
        tree.rnames = cellstr (num2str (i1))';%{};
        tree.name   = name;
    otherwise
        warning ('TREES:IO', 'format unknown'); varargout {1} = [];
        varargout {2} = tname; varargout {3} = path; return
end

if strfind (options, '-r'),
    if iscell (tree),
        for ward = 1 : length (tree),
            if iscell (tree {ward}),
                for te = 1 : length (tree {ward}),
                    tree{ward}{te} = repair_tree (tree{ward}{te});
                end
            else
                tree {ward} = repair_tree (tree{ward});
            end
        end
    else
        tree = repair_tree (tree);
    end
end

if strfind (options, '-s'),
    clf; hold on; title ('loaded trees');
    if iscell (tree),
        for ward = 1 : length (tree),
            if iscell (tree {ward}),
                for te = 1 : length (tree {ward}),
                    plot_tree (tree{ward}{te});
                end
            else
                plot_tree (tree{ward});
            end
        end
    else
        plot_tree (tree);
    end
    xlabel ('x [\mum]'); ylabel ('y [\mum]'); zlabel ('z [\mum]');
    view (3); grid on; axis image;
end

if (nargout > 0)
    varargout {1} = trees_split; % if output is defined then it becomes the tree
    varargout {2} = tree; varargout {3} = tname; varargout {4} = path;
else
    trees {length (trees) + 1} = tree; % otherwise add to end of trees cell array
end
end
