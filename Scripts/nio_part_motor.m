%% quantitative analysis of motor cortex field 
%% undefined parameters
M_vessel_path = ''; % the path of barrel cortex vessels data(.swc)
M_soma_path = ''; % the path of barrel cortex somas data(.swc)
Block_Size = 0; % the size of a block assigned from cortex outside surface, unit : ��m
SegImg = ''; % the path of a cortex template
stpSize = 5; % the statistic step length, unit : ��m {DEFAULT: 5 ��m}
Ms = 0; % motor cortex field size(x-coordinate) , unit : ��m
Ns = 0; % motor cortex field size(y-coordinate), unit : ��m
T_M = []; % threshold for microvasculature , unit : ��m
%% load vessels and somas
[vessels vessel1] = nio_load_tree(M_vessel_path);
% soma = load_tree(M_soma_path);
num_tree_vessels = length(vessels); % number of root nodes
%% transform the coordinates for adapting to Matlab
vessel = cell(1,1);
vessel{1} = vessel1;
vessels = nio_coordinate_transform( vessels, [1 1 1] );
vessel = nio_coordinate_transform( vessel, [1 1 1] );
vessel1 = vessel{1};
% soma = nio_coordinate_transform( soma );
%% generate binary images for every block and determine cortex area
BM = nio_generate_block( Ms, Ns, Block_Size, '');
num_part = length(BM);
[z_st z_en] = nio_cortex_area( BM, SegImg, 1 );
z_st = z_st+25;
z_en = z_en-25;
%% microvasculature
% vessels_c = cell(1, num_tree_vessels);
% for ward = 1 : num_tree_vessels
%     vessels_c{ward} = nio_extract_microvessels (vessels{ward}, T_M); % threshold
% end
%% fractional vascular volume(fv), vascular length density(ld) and diameter distribution(dia)
% for every block, _a;for barrel field in a block, _b;septa field in a
% block, _s
fv_a = cell(1, num_part); ld_a = cell(1, num_part);
% fv_a_c = cell(1, num_part); ld_a_c = cell(1, num_part);
dia_m = cell(1, num_part);
for ward = 1 : num_part
    ward
    [fv_a{ward} ld_a{ward} dia_m{ward}] = nio_vessel_single_block(vessels, BM{ward}, z_st(ward), z_en(ward));
%     [fv_a_c{ward} ld_a_c{ward}, ~] = nio_vessel_single_block(vessels_c, BM{ward}, z_st(ward), z_en(ward));
end
%% soma density(ds)
% ds_a = cell(1, num_part);
% for ward = 1 : num_part
%     ds_a{ward} = nio_soma(soma, BM{ward}, z_st(ward), z_en(ward));
% end
%% distance between soma to microvessel(dtm)
% m_dtm: average distance at different depth percentage of cortex
% e_dtm: std of blocks
% dtm: original data
% Ks = max(ceil(vessel1.Z), ceil(soma.Z));
% [ m_dtm_m e_dtm_m dtm_m d_cortex, ~] = nio_soma2tree( vessels, vessels_c, soma, Ms, Ns, Ks, Block_Size, z_st, z_en);