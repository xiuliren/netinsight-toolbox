% nio_read
% Loads tree(s) from the swc format.
% (trees package)
%
% [network, trees] = nio_read_swc (tname, options)
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

function varargout = nio_read (tname, options)

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
        [trees, tree, tname, path] = nio_load_tree( tname );
        network = nio_tree2network( trees );
    case '.hoc'
        if ~exist ([path tname], 'file'),
            error ('no such file...');
        end
        network = nio_read_hoc( tname );
        trees = [];
    otherwise
        warning ('TREES:IO', 'format unknown'); varargout {1} = [];
        varargout {2} = tname; varargout {3} = path; return
end

%% options
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
