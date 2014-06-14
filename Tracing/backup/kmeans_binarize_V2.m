%% binarize the image stack by kmeans clustering
% by jpwu, 2013/03/08
function bw_stk = kmeans_binarize(stk)

%% parameters
% % the minimum and maximum threshold
Tmin = 1;
% Tmax = 200;

%% test by global threshold
% bw_stk = ( stk > 20 );
% return;

%% processing
[IDX C] = kmeans(double(stk(:)), 2, 'emptyaction', 'drop');

if isnan(C(1))
    bw_stk = [];
    return;
end
stk(:) = IDX;

% the component
if C(1) > C(2)
    ct = 1;
else
    ct = 2;
end

if ct < Tmin
    bw_stk = [];
    return;
end

bw_stk = (stk == ct);
return;