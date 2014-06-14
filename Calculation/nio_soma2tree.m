% nio_soma2tree
% Calculate average distance from soma to the nearest microvessel.
% 
% [ m_dtm e_dtm DTM d_cortex dis_field] = 
% nio_soma2tree( vessels, vessels_c, soma, Ms, Ns, Ks, Block_Size, Z_ST, Z_EN, options )
% --------------------------------------------------------------------------------------
% 
% Get a distance map from soma to the nearest microvessel.
% By using vasculature-filled 3-dimension image stack, generate the 
% distance field, then collecting distance information accroding to the 
% location of every soma which inside the cortex.
% The image stack has same size to cortex area.
% 
% Input
% -----
% - vessels:: the vectorized vessels which has tree structure, multi-roots
% - vessels_c:: the vectorized microvessels which has tree structure,
%               multi-roots
% - soma:: the somas which has tree structure but no real connections
% - Ms:: the size of image stack, x-coordinate
% - Ns:: the size of image stack, y-coordinate
% - Ks:: the size of image stack, z-coordinate
% - Block_Size:: the size of a block assigned from cortex outside surface,
%                unit : ��m
% - Z_ST:: linear array, representing the starting depth(z-coordinates) of 
%          real cortex in all blocks. unit: ��m
% - Z_EN:: linear array, representing the ending depth(z-coordinates) of 
%          real cortex in all blocks. unit: ��m
% - options:: string, it can be empty by default.
%        '-a' : the average distance from soma to nearest vessel(all vasculature)
%        '-c' : for calculate correlation between distance and fractional vascular volume
%        '-f' : output the distance field which has same size to image
%        stack.
%        '-s' : plot
% 
% Output
% ------
% - m_dtm_m:: average distance at different depth percentage of cortex
%             unit: ��m
% - e_dtm_m:: std of blocks for distance map, unit: ��m
% - DTM:: cell array, original data from every block for distance map
%         unit: ��m
% - d_cortex:: linear array, incudes distances of all soma within cortex area
% - dis_field:: 3-dimension array, representing distance field -- distance
%               distribution in 3-dimension space
% 
% Example
% -------
% [ m_dtm_a e_dtm_a dtm_a ~ ~] = nio_soma2tree( sample_tree, sample_tree2, 
% soma, 1000, 600, 800, 100, 122, 590, '-c-s');
% 
% Uses nio_generate_block nio_swc2stk nio_distance_map

function [ m_dtm e_dtm DTM d_cortex varargout] = nio_soma2tree( vessels, vessels_c, soma, Ms, Ns, Ks, Block_Size, Z_ST, Z_EN, options )
%% build the point index stack
idx_stk = nio_swc2stk(vessels, Ms, Ns, Ks, '');
idx_stkc = nio_swc2stk(vessels_c, Ms, Ns, Ks, '');
%% compute distance map
if ~isempty(options)&&strfind (options, '-a')
    D = bwdist(idx_stk);
else
    D = bwdist(idx_stkc);
end
%%
num_part = length(Z_ST);
zv = cell(1, num_part);
if ~isempty(options)&&strfind (options, '-c')
    for ward = 1 : num_part
        zv{ward} = Z_ST(ward) : Block_Size : Z_EN(ward);
    end
else
    for ward = 1 : num_part
        zv{ward} = Z_ST(ward) : ((Z_EN(ward) - Z_ST(ward))/100) : Z_EN(ward);
    end
end
%% generate binary images for every block
BM = nio_generate_block(Ms, Ns, Block_Size, '');
% the number of somas
N = length( soma.X );
d = zeros(N,4);
%% get the soma distance, delete the somas within vessel voxel
for n = 1 : N
    d(n, 1) = soma.X(n);
    d(n, 2) = soma.Y(n);
    d(n, 3) = soma.Z(n);
    if idx_stk( ceil( soma.X(n) ), ceil( soma.Y(n) ), ceil( soma.Z(n) ) ) ~= 0
        continue;
    else
        d(n, 4) = D( ceil( soma.X(n) ), ceil( soma.Y(n) ), ceil( soma.Z(n) ) );
    end
end
idx_zero = d(:, 4) == 0;
d(idx_zero, :) = [];
%% distance field
if ~isempty(options)&&strfind (options, '-f')
    dis_field = zeros(Ms, Ns, Ks);
    for n = 1 : N
        dis_field(ceil( soma.X(n) ), ceil( soma.Y(n) ), ceil( soma.Z(n) )) = D( ceil( soma.X(n) ), ceil( soma.Y(n) ), ceil( soma.Z(n) ) );
    end
    varargout{1} = dis_field;
else
    varargout{1} = [];
end
%% delete the somas outside the cortex
n_i = floor(Ms/Block_Size);
n_j = floor(Ns/Block_Size);
d_cortex = d(:, 4);
del_idx = [];
for n = 1 : size(d, 1)
    if ceil(d(n, 1)) <= n_i * Block_Size && ceil(d(n, 2)) <= n_j * Block_Size
        idx_x = ceil(ceil(d(n, 1)) / Block_Size);
        idx_y = ceil(ceil(d(n, 2)) / Block_Size);
        idx_bm = n_j * (idx_x - 1) + idx_y;
        if d(n, 3) > Z_EN(idx_bm) || d(n, 3) < Z_ST(idx_bm)d(n, 3) < Z_ST(idx_bm)
            del_idx = [del_idx n];
        end
    else
        del_idx = [del_idx n];
    end
end
d_cortex(del_idx) = [];
%% get the dtm: distance map
DTM = cell(1, num_part);
for ward = 1 : num_part
    [DTM{ward}] = nio_distance_map(d, BM{ward}, zv{ward});
end
%% normalization
m_dtm = zeros(1, 101);
e_dtm = zeros(1, 101);
for i = 1 : 101
    amp = 0; anmsoma = 0; adis = [];
    for ward = 1 : num_part
        amp = amp + DTM{ward}(i, 1) * DTM{ward}(i, 2); % weight: cell number
        anmsoma = anmsoma + DTM{ward}(i, 2);
        adis = [adis DTM{ward}(i, 1)];
    end
    m_dtm(i) = amp / anmsoma;
    e_dtm(i) = std(adis(adis~=0));
end
%% plot
if ~isempty(options)&&strfind (options, '-s')
    plot(0:100, m_dtm_m, 'b', 'Linewidth', 2);hold on
    box off
    patch([0:100 100:(-1):0], [m_dtm_m - e_dtm_m fliplr(m_dtm_m + e_dtm_m)], 'b', 'EdgeColor', 'none');
    alpha(0.2)
    ylabel('Average distance from soma to the nearest microvessel (\mum)')
    xlabel('Cortical depth (%)');
end
end

