%% get seeds list
% by jpwu, 2013/03/22
function seedList = get_seed_list( stk, itv )
% parameters
rmin = 4;

% initialization
[M, N, K] = size( stk );
seedList = [];

% processing
for k = [1 : itv : K, K]
    % binaryize and connectivity analysis
    I(:,:) = stk(:,:,k);
    level = graythresh(I);
    BW = im2bw( I, level );
    S = regionprops(bwlabel(BW), 'MinorAxisLength', 'centroid');
    
    % add the seeds
    for si = 1 : length(S)
        tmp = S(si).Centroid(1);
        S(si).Centroid(1) = S(si).Centroid(2);
        S(si).Centroid(2) = tmp;
        radius = get_radius_V2( stk, [S(si).Centroid k], 255*level );
        seedList = [ seedList; S(si).Centroid, k, radius ];
    end
end

% sort the seedList according to diameter
seedList = sortrows( seedList, -4 );
% eliminate the small seeds
ind = find( seedList(:,4) < rmin );
seedList(ind,:) = [];