%% estimate the radius by a edge length adjustable cube
function r = get_radius(c, r)
global stk;
%% parameters
% step size
ss = 0.1;
% foreground voxel ratio threshold
Tfr = 0.3; % the ideal fraction: (4*pi/3)/(2^3) = 0.52
Tc_fr = 0.01;

% the magnitude of radius
P = 2;

%% processing
% get the local stack
[local_stk, m1, n1, k1] = get_local_stk( stk, [c r], P );
% binarization
T = kmeans_binarize( local_stk );
bw_local_stk = ( local_stk > T );
[Ml, Nl, Kl] = size( local_stk );

% check the boundary
if (Ml~=Nl) || (Nl~=Kl) || (Kl~=Ml)
    % hit the boundary
    %         maxr = min([M, N, K]);
    return;
end

% estimate the radius inside the stack
% count foreground inner voxel ratio
intr = uint16(r);
cl = c- double([m1-1, n1-1, k1-1]);
bw_tmp_stk = bw_local_stk( cl(1)-intr+1:cl(1)+intr-1, ...
    cl(2)-intr+1:cl(2)+intr-1, cl(3)-intr+1:cl(3)+intr-1 );

% initial foreground voxel ratio
[Mr, Nr, Kr] = size( bw_tmp_stk );
fr = length(find(bw_tmp_stk(:))) / (Mr*Nr*Kr);
if fr >Tfr
    % need to expand the radius
    while 1
        % adjust the radius
        r = r*(1 + ss);
        
        % count foreground inner voxel ratio
        intr = uint16(r);
        cl = c-double([m1-1, n1-1, k1-1]);
        bw_tmp_stk = bw_local_stk( cl(1)-intr+1:cl(1)+intr-1, ...
            cl(2)-intr+1:cl(2)+intr-1, cl(3)-intr+1:cl(3)+intr-1 );
        
        % the temporal ratio
        fr_tmp = fr;
        % foreground voxel ratio
        fr = length(find(bw_tmp_stk(:))) / (Mr*Nr*Kr);
        
        if fr < Tfr || ( fr_tmp-fr < Tc_fr )    
            return;
        end
    end
else
    % need to shrink the radius
    while 1
        % adjust the radius
        r = r*(1 - ss);
        % count foreground inner voxel ratio
        intr = uint16(r);
        cl = c-double([m1-1, n1-1, k1-1]);
        bw_tmp_stk = bw_local_stk( cl(1)-intr+1:cl(1)+intr-1, ...
            cl(2)-intr+1:cl(2)+intr-1, cl(3)-intr+1:cl(3)+intr-1 );
        
        % the temporal ratio
        fr_tmp = fr;
        % foreground voxel ratio
        fr = length(find(bw_tmp_stk(:))) / (Mr*Nr*Kr);
        
        if (fr > Tfr) || ( fr_tmp-fr < Tc_fr )    
            return;
        end
    end
end
return;

%% get the local cube
function [local_stk, m1, n1, k1] = get_local_stk( stk, centroid, P )
[M N K] = size(stk);
% the coordinate and radius
ms = centroid(1);   ns = centroid(2);
ks = centroid(3);   rs = centroid(4);

m1 = uint16( max( 1, ms-P*rs ) );   m2 = uint16( min(M, ms+P*rs) );
n1 = uint16( max( 1, ns-P*rs ) );   n2 = uint16( min(N, ns+P*rs) );
k1 = uint16( max( 1, ks-P*rs ) );   k2 = uint16( min(K, ks+P*rs) );

local_stk = stk(m1:m2, n1:n2, k1:k2);
return;