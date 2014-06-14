%% quantitative analysis of visual cortex field 
%% undefined parameters
V_vessel_path = ''; % the path of barrel cortex vessels data(.swc)
V_soma_path = ''; % the path of barrel cortex somas data(.swc)
Block_Size = 0; % the size of a block assigned from cortex outside surface, unit : ��m
SegImg = ''; % the path of a cortex template
stpSize = 5; % the statistic step length, unit : ��m {DEFAULT: 5 ��m}
Ms = 0; % visual cortex field size(x-coordinate) , unit : ��m
Ns = 0; % visual cortex field size(y-coordinate), unit : ��m
T_V = []; % threshold for microvasculature , unit : ��m
%% load vessels and somas
% [vessels vessel1] = nio_load_tree(V_vessel_path);
[~,soma] = nio_load_tree(V_soma_path);
% num_tree_vessels = length(vessels); % number of root nodes
%% transform the coordinates for adapting to Matlab
% vessel = cell(1,1);
% vessel{1} = vessel1;
% vessels = nio_coordinate_transform( vessels, [1 1 1] );
% vessel = nio_coordinate_transform( vessel, [1 1 1] );
% vessel1 = vessel{1};
somas = cell(1,1);
somas{1} = soma;
somas = nio_coordinate_transform( somas, [1 1 1] );
soma = somas{1};
%% generate binary images for every block and determine cortex area
BM = nio_generate_block( Ms, Ns, Block_Size, '');
num_part = length(BM);
%%
[z_st z_en] = nio_cortex_area( BM, SegImg, 1);
z_st = z_st+25;
z_en = z_en-25;
%% microvasculature
% vessels_c = cell(1, num_tree_vessels);
% for ward = 1 : num_tree_vessels
%     vessels_c{ward} = nio_extract_microvessels (vessels{ward}, T_V); % threshold
% end
%% fractional vascular volume(fv), vascular length density(ld) and diameter distribution(dia)
% for every block, _a;for barrel field in a block, _b;septa field in a
% block, _s
% fv_a = cell(1, num_part); ld_a = cell(1, num_part);
% fv_a_c = cell(1, num_part); ld_a_c = cell(1, num_part);
% dia_v = cell(1, num_part);
% for ward = 1 : num_part
%     ward
%     [fv_a{ward} ld_a{ward} dia_v{ward}] = nio_vessel_single_block(vessels, BM{ward}, z_st(ward), z_en(ward));
%     [fv_a_c{ward} ld_a_c{ward}, ~] = nio_vessel_single_block(vessels_c, BM{ward}, z_st(ward), z_en(ward));
% end
%% soma density(ds)
ds_a = cell(1, num_part);
for ward = 1 : num_part
    ds_a{ward} = nio_soma(soma, BM{ward}, z_st(ward), z_en(ward), '');
end
%% distance between soma to microvessel(dtm)
% m_dtm: average distance at different depth percentage of cortex
% e_dtm: std of blocks
% dtm: original data
% Ks = max(ceil(vessel1.Z), ceil(soma.Z));
% [ m_dtm_v e_dtm_v dtm_v d_cortex, ~] = nio_soma2tree( vessels, vessels_c, soma, Ms, Ns, Ks, Block_Size, z_st, z_en);