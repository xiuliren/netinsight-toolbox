% calculate the minimun distance between somas and tree
% first create a distance map in a stack
% INPUTS:
%   Xs, Ys, Zs: the vector of somas coordinate
%   tree: the neurites tree
% OUTPUTS:
%   d: the minimun distance between somas and tree
%       the length of d must equal to somas
% EXAMPLE:
%   [Xs Ys Zs] = genrate_random_somas( Mb, Nb, Kb, 10, '-s' );
%   d = soma2tree_distance( Xs, Ys, Zs, sample_tree );
% See also rpoints_tree eucdist
%
% Uses
%
% by jpwu, 2011.7.5, all rights reserved.
function d = soma2tree_distance_V2( somas, tree, options )

%% parameters
% the origional volume size
Ms = 401;
Ns = 434;
Ks = 831;

% the number of somas
N = length( somas.X );

% initiate the shortest distance from soma to vessel
d = ones(N,1).*Inf;

% node number of a vessel trees
K = length(tree.X);
% the parent index of each node
idpar = idpar_tree(tree);

%% create the distance map
% initiate the stack
% distMap = ones( Ms, Ns, Ks ) * Inf;
bin_stk = zeros( Ms, Ns, Ks );
for k = 1 : K
    k
    % the two points
    a = [ tree.X(k) tree.Y(k) tree.Z(k) ];
    b = [ tree.X( idpar(k) ) tree.Y( idpar(k) ) tree.Z( idpar(k) ) ];
    
    if ( idpar(k)== -1 ) | ( a == b )
        bin_stk( round(a(1)), round(a(2)), round(a(3)) ) = 1;
        continue;
    end

    % the unit vector
    ab = b - a;
    % the length
    len = norm( ab );
    % unit vector
    uv = ab / len;
    for m = 0 : 1 : len
        c = a + m * uv;
        try
            bin_stk( round(c(1)), round(c(2)), round(c(3)) ) = 1;
        catch exception
            idx = find( c < 1 );
            c( idx ) = 1;
            bin_stk( round(c(1)), round(c(2)), round(c(3)) ) = 1;
        end
    end
end

% compute distance map
D = bwdist(bin_stk);

% get the soma distance
for n = 1 : N
    d(n) = D( round( somas.X(n) ), round( somas.Y(n) ), round( somas.Z(n) ) );
end
%% plot
if findstr (options, '-s');
    hist(d);
    title  ('histogram of minimum distance between somas and tree');
    xlabel ('distance [\mum]'); ylabel ('frequency');
end

%% calculate the distance from a point to line segment function
% a and b is two vertices, and c is the isolated point
% reference: http://chensavvy.blog.163.com/blog/static/5715718920090292919632/
function d = p2LS(a,b,c)
ab = b - a;
ac = c - a;
f = dot (ab,ac);
if f < 0
    d = sqrt( (a(1) - c(1))^2 + (a(2) - c(2))^2 + (a(2) - c(2))^2 );
    return;
end
dt = dot(ab,ab);
if (f > dt)
    d = sqrt( (b(1) - c(1))^2 + (b(2) - c(2))^2 + (b(2) - c(2))^2 );
    return;
end
f = f./dt;
D = a + f*ab;
d = sqrt( (D(1) - c(1))^2 + (D(2) - c(2))^2 + (D(2) - c(2))^2 );
return;

%by jpwu, 2011.7.5
% calculate euclidian distance between two 3D points
% function dist = euc_distance( X1, X2, Y1, Y2, Z1, Z2)
% dist = euc_distance(Xs(n), resampled_tree.X, Ys(n), resampled_tree.Y,
% Zs(n), resampled_tree.Z);
% dist = sqrt ( (repmat (X1, 1, length (X2)) - repmat (X2', length (X1), 1)).^2 + ...
%     ( repmat(Y1, 1, length (Y2) ) - repmat(Y2', length (Y1), 1) ).^2  + ...
%     ( repmat(Z1, 1, length (Z2) ) - repmat(Z2', length (Z1), 1) ).^2 );