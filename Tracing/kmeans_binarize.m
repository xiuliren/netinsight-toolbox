%% binarize the image stack by kmeans clustering
% by jpwu, 2013/03/13
function T = kmeans_binarize(stk)

%% parameters
% % the minimum and maximum threshold
Tmin = 2;
Tmax = 50;
Emax = 7;
% background intensity
global BI;

%% test by global threshold
% T = 5;
% return;

E = entropy(stk);

%% processing
% construct the histogram
h = histc( reshape(stk, [],1), 0:255 )';

% initial threshold by mean intensity
T = uint8 ( sum( double(h).*([0:255]) ) / sum( double(h) ) );

while 1
    % comput the centroid of two clusters, low and high intensity
    cl = sum(( h(1:T)/sum(h(1:T)) ) .* double( (0:T-1) ) ) ;
    ch = sum(( h(uint16(T)+1:256)/sum(h(uint16(T)+1:256)) ) .* double( (T:255) ) );
    
    % the new threshold
    newT = uint8( (cl + ch)/2 );
    
    if newT == T || (cl==0 && ch==255) 
        % the final threshold
        T = (T + BI)/2;
        T = T + T*(E/Emax);
        T = min(T, Tmax);
        T = max(T, Tmin);
        return;
    else
        % a new threshold
        T = newT;
    end
end

return;