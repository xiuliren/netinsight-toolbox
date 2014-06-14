%% binarize the image stack by kmeans clustering
% by jpwu, 2013/02/28
function [bw_stk T ] = kmeans_binarize(stk)
%% parameters
% the minimum and maximum threshold
Tmin = 50;
Tmax = 200;
%% processing
h = histc( reshape(stk, [],1), 0:255 );
% compute the mean intensity by histogram vector
ml = round( hist_statistics(h) );
% threshold computed by kmeans clustering
IDX = kmeans(h,2,'emptyaction', 'singleton');
T = find(h, 1,'first');
T
T = max(T, 30); T = min(T, Tmax);
% the estimated threshold must bigger than the mean intensity of local cube
% T = ( uint8(ml) + T )/2;

bw_stk = (stk >= T);
return;

%% compute some statistic parameters from frequency distribution
function ml = hist_statistics( h )
s = double(0);
for n = 1 : length(h)
    s = s + h(n)*(n-1);
end
ml = s / length(h);
return;