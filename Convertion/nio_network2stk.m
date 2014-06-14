% transform network to image stack which contains the section index
% by jpwu, 2010/09/19
function stk = nio_network2stk( network, D_scale )

%% load network for debug only, this code cell must be commented 
% clc
% clear
% load( 'network.mat' );

%% parameters
if nargin == 1
    D_scale = 1;
end
%% estimate the stack size
% initialize the stack size
Ms = 0;
Ns = 0;
Ks = 0;
% compute the maximun coordinate
for idx_s = 1 : network.sn
    sec = network.sections{ idx_s };
    Ms = max( Ms,  max( sec(:,1) ) );
    Ns = max( Ns, ( max(sec(:,2) ) ) );
    Ks = max( Ks, ( max(sec(:,3)) ) );
end
Ms = ceil(Ms);
Ns = ceil(Ns);
Ks = ceil(Ks);

%% transformation
stk = zeros( Ms, Ns, Ks, 'double' );
for idx_s = 1 : network.sn
    sec = network.sections{ idx_s };
    % get the local box range
    for idx_p = 1 : size(sec,1)
        r = (sec(idx_p,4)/2)*D_scale ;
        m1 = ceil( max( 1,   sec(idx_p,1)-r) );
        m2 = ceil( min( Ms,  sec(idx_p,1)-r ) );
        n1 = ceil( max( 1,   sec(idx_p,2)-r ) );
        n2 = ceil( min( Ns,  sec(idx_p,2)-r ) );
        k1 = ceil( max( 1,   sec(idx_p,3)-r ) );
        k2 = ceil( min( Ks,  sec(idx_p,3)-r ) );
    end
    stk( m1:m2, n1:n2, k1:k2 ) = idx_s; % ones( m2-m1+1, n2-n1+1, k2-k1+1, 'double' ).*idx_s;
end

