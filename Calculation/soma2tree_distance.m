% calculate the minimun distance between somas and tree
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
function d = soma2tree_distance( Xs, Ys, Zs, tree, options )

%% parameters
N = length( Xs );
d = ones(N,1).*Inf;
K = length(tree.X);
idpar = idpar_tree(tree);
%% calculate distance
% resample the tree
% resampled_tree = resample_tree (tree, Resample_Length, '-d' );
% get the minimum distance
N
for n = 1 : N
    n
    c = [Xs(n) Ys(n) Zs(n)];
    for k = 1 : K
        a = [ tree.X(k) tree.Y(k) tree.Z(k) ];
        b = [ tree.X(idpar(k)), tree.Y( idpar(k) ), tree.Z( idpar(k) )];
        dtmp = p2LS(a,b,c);
        if d(n) > dtmp
            d(n) = dtmp;
        end
    end
end

%% plot
if findstr (options, '-s');
    hist(d);
    %     legend (HP, 'distance');
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