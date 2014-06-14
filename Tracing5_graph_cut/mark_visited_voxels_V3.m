%% mark the visited voxels, only mark the final point
function mk_local_stk = mark_visited_voxels_V3( mk_local_stk, seed, P_list )
% the expansion rate
ER = 1.2;
% minimum and maximum spaning range
MINSR = 2;  MAXSR = 7;
P_list = uint16(P_list);

for nn = 1 : size(P_list, 1)   
    r2 = uint16( ER * P_list(nn, 4) );
    r2 = min(r2, P_list(nn,4)+MAXSR);   
    r2 = max(r2, P_list(nn,4)+MINSR);
    
    % fill the sphere of new node
    m0 = P_list(nn,1);  n0 = P_list(nn,2);
    k0 = P_list(nn,3);  r0 = r2;
    mk_local_stk = fill_mark_stack_with_sphere(mk_local_stk, m0,n0,k0,r0);
    
    % fill the sphere of center node 
    m0 = (seed(1) + P_list(nn,1)) / 2;
    n0 = (seed(2) + P_list(nn,2)) / 2;
    k0 = (seed(3) + P_list(nn,3)) / 2;
    r0 = (seed(4) + r2 ) / 2;
    mk_local_stk = fill_mark_stack_with_sphere(mk_local_stk, m0,n0,k0,r0);
end

%% fill the mark stack with sphere
function mk_stk = fill_mark_stack_with_sphere(mk_stk, m0,n0,k0,r0)
% stack size
[M N K] = size(mk_stk);
% create sphere as template
S = hypersphere(double([r0*2+1, r0*2+1, r0*2+1]),'full');

% initialize the coordinate
m1 = m0-r0; m2 = m0+r0;
n1 = n0-r0; n2 = n0+r0;
k1 = k0-r0; k2 = k0+r0;

ms1 = 1;    ms2 = 2*r0+1;
ns1 = 1;    ns2 = 2*r0+1;
ks1 = 1;    ks2 = 2*r0+1;

% the boundary condition
if m0-r0 < 1
    m1 = 1;
    ms1 = r0-m0+2;
end
if m0+r0 > M
    m2 = M;
%         ms2 = r0 + 1 -m0 + M;
    ms2 = m2-m1+ms1;
end

if n0-r0 < 1
    n1 = 1;
    ns1 = r0-n0+2;
end
if n0+r0 > N
    n2 = N;
%         ns2 = r0 + 1 -n0 + N;
    ns2 = n2-n1+ns1;
end

if k0-r0 < 1
    k1 = 1;
    ks1 = r0-k0+2;
end
if k0+r0 > K
    k2 = K;
%         ks2 = r0 + 1 -k0 + K;
    ks2 = k2-k1+ks1;
end

mk_stk( m1:m2, n1:n2, k1:k2 ) = mk_stk( m1:m2, n1:n2, k1:k2 )...
    | S( ms1:ms2, ns1:ns2, ks1:ks2 );
return;