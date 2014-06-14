%% mark the visited voxels for line scan
function mk_local_stk = mark_visited_voxels_cube(  mk_local_stk, seed, P_list )
% the expansion rate
ER = 1.2;
% minimum and maximum spaning range
MINSR = 2;  MAXSR = 7;

% the size of the stack
[M,N,K] = size(mk_local_stk);

% evaluate the node list
if isempty( P_list )
    return;
end

for nn = 1 : size(P_list, 1)
    r1 = uint16( ER * seed(4) );        r2 = uint16( ER * P_list(nn, 4) );
    r1 = min(r1, seed(4)+MAXSR);        r2 = min(r2, P_list(nn,4)+MAXSR);
    r1 = max( r1, seed(4)+MINSR );      r2 = max( r2, P_list(nn,4)+MINSR );
    % vector of interprated 3D points using bresenham's algorithm
    [vm vn vk] = bresenham_line3d( double(seed(1:3)), double(P_list(nn,1:3)) );
    Ni = length(vm);
    if Ni==0
        return;
    end
    if r1 > r2
        vr = r1 : -1 : r2;
    else
        vr = r1 : 1 : r2;
    end
    Rv = imresize( vr, [1, Ni], 'nearest' );
    for ni = 1 : 2 : Ni
        m0 = vm(ni); n0 = vn(ni); k0 = vk(ni); r0 = Rv(ni);
        m1 = max(1, m0-r0);   m2 = min(M, m0+r0);
        n1 = max(1, n0-r0);   n2 = min(N, n0+r0);
        k1 = max(1, k0-r0);   k2 = min(K, k0+r0);
        mk_local_stk( m1:m2, n1:n2, k1:k2 ) = true( m2-m1+1, n2-n1+1, k2-k1+1 );
    end
end