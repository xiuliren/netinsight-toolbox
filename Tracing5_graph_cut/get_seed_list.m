%% get seeds list
% by jpwu, 2013/03/22
function seedList = get_seed_list( stk, itv )
% parameters
rmin = 2;

% initialization
[~, ~, K] = size( stk );
seedList = [];

% processing
for k = [1 : itv : K, K]
    % binaryize and connectivity analysis
    I(:,:) = stk(:,:,k);
    level = graythresh(I);
    BW = im2bw( I, level );
    S = regionprops(logical(BW), 'MinorAxisLength', 'centroid');
    
    % add the seeds
    for si = 1 : length(S)
        tmp = S(si).Centroid(1);
        S(si).Centroid(1) = S(si).Centroid(2);
        S(si).Centroid(2) = tmp;
        radius = get_radius_V3( stk, [S(si).Centroid k], 255*level );
        seedList = [ seedList; S(si).Centroid, k, radius ];
    end
end

% eliminate the small seeds
ind =  seedList(:,4) < rmin ;
seedList(ind,:) = [];
% sort the seedList according to diameter
seedList = sortrows( seedList, -4 );
