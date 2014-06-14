% nio_vessel_single_block 
% Calculate parameters of vasculature in a single block.
%
% [fv len_den dia] = nio_vessel_single_block(vessels, bm, z_st, z_en, stpSize)
% -------------------------------------------------------------------------------------
% 
% Calculate fractional vascular volume(fv), vascular length density(len_den)
% and diameter distribution(dia) in a single block which is segmented from
% cortex surface. Split the tree structure by insert new nodes when hits
% the block boundary.
% 
% Input
% -----
% - vessels:: the vectorized vessels which has a tree structure
% - bm:: a binary image(m*n matrix whose elements are 0 or 1)which 
%                 marks location of a single block.
% - z_st:: representing the starting depth(z-coordinates) of real cortex in a single block
%          unit: ��m
% - z_en:: representing the ending depth(z-coordinates) of real cortex in a single block
%          unit: ��m
% - stpSize:: the statistic step length along depth, unit : ��m {DEFAULT: 5 ��m}
% 
% Output
% ------
% - fv:: linear array, fractional vascular volume in different depth
% - len_den:: linear array, vascular length density in different depth
%             unit: m/mm^3
% - dia:: linear array, diameter distribution from vasculature in the block, sample
%         interval along vessel centerline: 1 ��m, unit: ��m
% 
% Example
% -------
% vessels = sample_tree; bm = BM{1};
% [fv len_den dia] = nio_vessel_single_block(vessels, bm, 66, 606);
% 
% Uses idpar_tree nio_intersection insert_tree nio_vessel nio_diameter

function [fv len_den dia] = nio_vessel_single_block(vessels, bm, z_st, z_en, varargin)
disp('analyzing single block ... ...');
if(isempty(varargin))
    stpSize = 5;
else
    stpSize = varargin{1};
end
area = length(find( bm ));
zv = z_st : stpSize : z_en;
if area == 0
    fv = zeros(1, length(zv));
    len_den = zeros(1, length(zv));
    dia = [];
    return;
end
num_tree_vessels = length(vessels);
fv = [];
len_den = [];
dia = [];
for ward = 1 : num_tree_vessels
    %% delete the nodes outside the block
    
    N = length(vessels{ward}.X);
    nodes_idx = zeros(N, 4);
    
    nodes_idx(:, 2) = vessels{ward}.X;
    nodes_idx(:, 3) = vessels{ward}.Y;
    for n = 1 : N
        nodes_idx(n, 1) = n;
        nodes_idx(n, 4) = bm(ceil(max(abs(nodes_idx(n, 2)), 1)), ceil(max(abs(nodes_idx(n, 3)), 1)));
    end
    
    del_idx = nodes_idx(nodes_idx(:, 4) == 0, 1);
    d_N = length(del_idx);
%     disp('delete the nodes outside the block ... ...');
    idpar = idpar_tree(vessels{ward}, '-0');
    for n = 1 : d_N
        idpar_n = idpar(del_idx(n));
        if idpar_n && nodes_idx(idpar_n, 4) % need delete child node, generate a new child node
            [x_is y_is z_is d_is] = nio_intersection(vessels{ward}, bm, del_idx(n), idpar_n);   
            vessels{ward}.dA(del_idx(n), idpar_n) = 0;
            vessels{ward} = insert_tree(vessels{ward}, [1 6 x_is y_is z_is d_is idpar_n]);
        end
        
        child_n = find(vessels{ward}.dA(:, del_idx(n)) == 1);
        if ~isempty(child_n)
            for m = 1 : length(child_n)
                if nodes_idx(child_n(m), 4) % need delete parent node, generate a child node for the child node
                    [x_is y_is z_is d_is] = nio_intersection(vessels{ward}, bm, child_n(m), del_idx(n));
                    vessels{ward} = insert_tree(vessels{ward}, [1 6 x_is y_is z_is d_is child_n(m)]);
                end
            end
        end
    end
    
    N = length(vessels{ward}.X);
    nodes_idx = zeros(N, 4);
    nodes_idx(:, 2) = vessels{ward}.X;
    nodes_idx(:, 3) = vessels{ward}.Y;
    for n = 1 : N
        nodes_idx(n, 1) = n;
        nodes_idx(n, 4) = bm(ceil(max(abs(nodes_idx(n, 2)), 1)), ceil(max(abs(nodes_idx(n, 3)), 1)));
    end
    
    nodes = nodes_idx(nodes_idx(:, 4) == 1, 1);
    %% calculate
    
    [fv_tmp len_den_tmp] = nio_vessel( vessels{ward}, nodes', area, z_st, z_en, stpSize, '');
    fv = [fv; fv_tmp];
    len_den = [len_den; len_den_tmp];
    dia_tmp = nio_diameter( vessels{ward}, nodes');
    dia = [dia; dia_tmp];
end
fv = sum(fv, 1);
len_den = sum(len_den, 1);
end

