%% quantitative analysis of every single barrel in barrel cortex field 
%% undefined parameters
B_vessel_path = ''; % the path of barrel cortex vessels data(.swc)
B_soma_path = ''; % the path of barrel cortex somas data(.swc)
Block_Size = []; % the size of a block assigned from cortex outside surface, unit : ��m
SegImg = ''; % the path of a cortex template
stpSize = []; % the statistic step length, unit : ��m {DEFAULT: 5 ��m}
Ms = []; % barrel cortex field size(x-coordinate) , unit : ��m
Ns = []; % barrel cortex field size(y-coordinate), unit : ��m
T_B = []; % threshold for microvasculature , unit : ��m
bm_path = ''; % the path of image for barrels
%% load vessels and somas
[vessels vessel1] = nio_load_tree(B_vessel_path);
% soma = load_tree(B_soma_path);
num_tree_vessels = length(vessels); % number of root nodes
%% transform the coordinates for adapting to Matlab
vessels = nio_coordinate_transform( vessels, [0.35 0.35 0.35] );
vessel1 = nio_coordinate_transform( vessel1, [0.35 0.35 0.35] );
% soma = nio_coordinate_transform( soma );
%% generate binary images for every block and determine cortex area
barrel_mask = imread(bm_path);
barrel_mask = im2bw(barrel_mask, 0.7) ; % binaryzation
barrel_mask = bwlabel(barrel_mask, 8); % label every barrel
BM = nio_generate_block( Ms, Ns, Block_Size, '-b', barrel_mask);
num_part = length(BM);
[z_st z_en] = nio_cortex_area( BM, SegImg);
%% microvasculature
vessels_c = cell(1, num_tree_vessels);
for ward = 1 : num_tree_vessels
    vessels_c{ward} = nio_extract_microvessels (vessels{ward}, T_B); % threshold
end
%% fractional vascular volume(fv), vascular length density(ld) and diameter distribution(dia)
% for every block, _a;for barrel field in a block, _b;septa field in a
% block, _s
fv = cell(1, num_part); ld = cell(1, num_part);
fv_c = cell(1, num_part); ld_c = cell(1, num_part);
for ward = 1 : num_part
    ward
    [fv{ward} ld{ward}, ~] = nio_vessel_single_block(vessels, BM_a{ward}, z_st(ward), z_en(ward));
    [fv_c{ward} ld_c{ward} ~] = nio_vessel_single_block(vessels_c, BM_a{ward}, z_st(ward), z_en(ward));
end
%% soma density(ds)
ds = cell(1, num_part);
for ward = 1 : num_part
    ward
    ds{ward} = nio_soma(soma, BM_a{ward}, z_st(ward), z_en(ward));
end