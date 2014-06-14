% nio_diameter
% Get diameter distribution in one-root tree structure on particular
% segments.
% 
% dia = nio_diameter( vessel, nodes, sample_step )
% ------------------------------
% 
% Get diameter distribution in one-root tree structure on particular
% segments.
% 
% Input
% -----
% - vessel:: the vectorized vessels which has a tree structure
% - nodes:: indexes of the statistic nodes, to calculate particular
%           segments
% - sample_step:: sampling step along the tree path. {Defaut: 1 ¦Ìm}
% 
% Output
% ------
% - dia:: diameter distribution. unit: ¦Ìm
% 
% Example
% -------
% dia = nio_diameter( sample_tree, 1:1000 )
% 
% See also nio_vessel_single_block
% Uses idpar_tree

function [ dia ] = nio_diameter( vessel, nodes, varargin )
if(isempty(varargin))
    sample_step = 1;
else
    sample_step = varargin{1};
end
idpar = idpar_tree(vessel, '-0'); % get the parent nodes
[~, ~, id_idpar] = intersect(nodes, idpar(nodes)); % delete the nodes which don't have parent
nodes = nodes(id_idpar);
d = zeros(length(nodes), 2);
flag = 1;
for n = nodes
    x1 = vessel.X(n); x2 = vessel.X( idpar(n) );
    y1 = vessel.Y(n); y2 = vessel.Y( idpar(n) );
    z1 = vessel.Z(n); z2 = vessel.Z( idpar(n) );
    D1 = vessel.D(n); D2 = vessel.D( idpar(n) );
    distance = sqrt((x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2);
    d(flag, 1) = (D1 + D2) / 2;
    if round(distance/sample_step)
        d(flag, 2) = round(distance/sample_step); % the weight for a diameter value in one vessel segment
    else
        d(flag, 2) = 1;
    end
    flag = flag + 1;
end
dia = d(:, 1);
for n = 1 : length(nodes)
    if d(n, 2) > 1
        tmp = d(n, 1) * ones(d(n, 2) - 1, 1); % weight greater than 1
        dia = [dia; tmp];
    end
end
end