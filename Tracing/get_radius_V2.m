%% estimate the radius by rayburst sampling
function radius = get_radius_V2(stk, c, T)
% global stk;
[M N K] = size(stk);
%% parameters
% P = 3;

%% get the binarized local stack
% % get the local stack
% % [local_stk, m1, n1, k1] = get_local_stk( stk, [c r], P );
% % binarization
% % T = kmeans_binarize( local_stk );
% bw_local_stk = ( local_stk > T );
% [Ml, Nl, Kl] = size( local_stk );
% 
% % the local center
% [mc nc kc] = c(1:3) - [m1 n1 k1] + 1;
%% boundary checking
c = int16(c);
if c(1)<=0 || c(2)<=0 || c(3)<=0 || c(1)>=M || c(2)>=N || c(3)>=K
%     disp('out of boundary!')
    radius = 0;
    return;
end

%% bursting rays in 26 directions
% establish the step and length
sl =   [ 0  0 1  1;         0  0 -1 1;
         0  1 0  1;         0 -1 0  1;
         0  1 1  1.4142;    0  1 -1 1.4142;     0  -1 1 1.4142;     0 -1 -1  1.4142;
         1  0 0  1;         -1 0  0 1;
         1  0 1  1.4142;    1  0 -1 1.4142;     -1  0 1 1.4142;     -1 0 -1  1.4142;
         1  1 0  1.4142;    1  -1 0 1.4142;     -1  1 0 1.4142;     -1 -1 0  1.4142;
         1  1 1  1.7321;    1  1 -1 1.7321;     1  -1 1 1.7321;     1  -1 -1 1.7321;    -1  1  1  1.7321;   -1  1  -1  1.7321;  -1  -1  1  1.7321;  -1  -1  -1 1.7321;
        ];
   
% the length of casted rays
lr = zeros( 26, 1 );

% direction index
for di = 1 : 26
    % current coordinate
    cc = c;
    while stk(cc(1), cc(2), cc(3)) > T
        cc = cc + int16(sl(di, 1:3));
        lr(di) = lr(di) + sl(di, 4);
        if cc(1)==0 || cc(2)==0 || cc(3)==0 || cc(1)==M || cc(2)==N || cc(3)==K
            break;
        end
    end
end

% sort the length vector
lr = sort(lr);
radius = lr(16); % 26*3/4