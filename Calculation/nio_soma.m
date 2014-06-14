% nio_soma
% Calculate soma number density by sliding window.
% 
% 
% ds = nio_soma( soma, bm, z_st, z_en, options, stpSize)
% -----------------------------------------------------------------------------
%
% Calculate soma number density by sliding window in a single block.
% 
% Input
% -----
% - soma:: the somas which has a  tree structure but no real connections
% - bm:: a binary image(m*n matrix whose elements are 0 or 1)which 
%                 marks location of a single block.
% - z_st:: representing the starting depth(z-coordinates) of real cortex in a single block
%          unit: ��m
% - z_en:: representing the ending depth(z-coordinates) of real cortex in a single block
%          unit: ��m
% - options:: string, it can be empty by default.
%        '-s' : plot
% - stpSize:: the statistic step length along depth, unit : ��m {DEFAULT: 5 ��m}
% 
% Output
% ------
% - ds:: linear array, soma number density in different depth
%        unit: 10^5/mm^3
% 
% Example
% -------
% ds{1} = nio_soma( sample_tree, BM{1}, 66, 606, '-s', 10)

function ds = nio_soma( soma, bm, z_st, z_en, options, varargin)
%% parameters
% the size of slide window and step 
winSize = 25;
if ~isempty(options)&&~ischar(options)
    varargin{1} = options;
    options = [];
end
if(isempty(varargin))
    stpSize = 5;
else
    stpSize = varargin{1};
end
zv = z_st : stpSize : z_en;
area = length(find(bm == 1));
if area == 0
    ds = zeros(1, length(zv));
    return;
end
%% process
disp('------- start analysing somas -----------');
N = size(soma.X,1);
idx = zeros(N, 2);
for ward = 1 : N
    idx(ward, 1) = ward;
    idx(ward, 2) = bm(min(size(bm, 1), ceil(soma.X(ward))), min(size(bm, 2), ceil(soma.Y(ward))));
end
x = find(idx(:, 2) == 1);
N = length(x);
ds = [];
% get the soma number of every slide window
for z = zv
    z
    ns1 = length( find(soma.Z(x) >= z - winSize) );
    ns2 = length( find(soma.Z(x) < z + winSize) );
    ns = ns1 + ns2 - N;
    ds = [ds ns];
end

% transform the unit to 10e5/mm3
ds = ds .* ( 10000 / (area*winSize*2) );
disp('------------ somas analysis end ------------');
%% plot results
if ~isempty(options)&&strfind (options, '-s')
    plot (zv-z_st, ds, 'b','LineWidth',2);
    xlabel('Depth below cortical surface(\mum)');
    ylabel('Cell number density\n(10^5mm^-^3)');
end
end
