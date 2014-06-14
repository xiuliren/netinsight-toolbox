%% build the index and coordinate map
function circles = get_sphere_circles( seed, R )
% initialization
seed = uint32(seed);
% [M,N,K] = size(local_stk);
% % the integer radius of sphere
% R = (M-1)/2;
circles = cell( 2*R+1 ,1 );

% build the circles in the 
n0 = seed(2);   k0 = seed(3);
for h = double(R : -1 : -R)
    m0 = seed(1) - h; 
%     I(:,:) = local_stk(m0,:,:);
    % the radius of the circle
    r = uint32( sqrt( double(R*R - h*h) ) );
    % the circular points
    [nc,kc] = getmidpointcircle(n0,k0,r);
    % remove the replicate nodes, which is a bug in the function
    nc2 = nc;   nc2(1:end-1) = nc(2:end);   nc2(end) = nc(1); 
    kc2 = kc;   kc2(1:end-1) = kc(2:end);   kc2(end) = kc(1);
    label = (nc == nc2) & (kc==kc2);
    nc(label) = []; kc(label) = [];
    
    circles{R-h+1} = uint32( [ repmat(m0,length(nc),1), nc, kc] );
end