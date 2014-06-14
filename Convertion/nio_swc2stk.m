% nio_swc2stk
% Convert tree strcture to a image stack.
% 
% stk = nio_swc2stk(vessels, Ms, Ns, Ks, options, sample_step )
% -------------------------------------------------------
% 
% Convert tree strcture derived from swc file to a vessel-filled image 
% stack. The image stack can be filled with node id or 1, which depend on 
% options.
% Vasculature filling mode: spheric filling, center of the sphere is 
% determined by unit vector.
% 
% Input
% -----
% - vessels:: the vectorized vessels which has tree structure, multi-roots
% - Ms:: the size of image stack, x-coordinate
% - Ns:: the size of image stack, y-coordinate
% - Ks:: the size of image stack, z-coordinate          
% - options:: string, it can be empty by default.
%        '-i' : fill the stack by node id
%        {DEFAULT: empty, fill the stack by 1, i.e. a binary stack}
% - sample_step:: sampling step along the tree path. {Defaut: 1 ¦Ìm}
% 
% Output
% ------
% - stk:: a vessel-filled image stack
% 
% Example
% -------
% idx_stk = nio_swc2stk( sample_tree, 1000, 600, 800);
% 
% See also nio_soma2tree

function stk = nio_swc2stk(vessels, Ms, Ns, Ks, options, varargin )
% parameters
if ~isempty(options)&&~ischar(options)
    varargin{1} = options;
    options = [];
end
num_trees = length(vessels);
% the resampled node distance
if(isempty(varargin))
    sample_step = 1;
else
    sample_step = varargin{1};
end

%% build the point index stack
disp('-------- getting stack ...')
stk = zeros (Ms, Ns, Ks);
for ward = 1 : num_trees
    % points connectivity, child and parent
    [chl par] = find(vessels{ward}.dA == 1);
    for con = 1 : length( chl )
        % the child and parent point
        x1 = vessels{ward}.X( chl(con) );
        y1 = vessels{ward}.Y( chl(con) );
        z1 = vessels{ward}.Z( chl(con) );
        D1 = vessels{ward}.D( chl(con) );

        x2 = vessels{ward}.X( par(con) );
        y2 = vessels{ward}.Y( par(con) );
        z2 = vessels{ward}.Z( par(con) );
        D2 = vessels{ward}.D( par(con) );

        % point distance
        dis = sqrt( (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2) );
        if dis < sample_step
            continue;
        end

        % unit vector
        uv = [ (x2-x1) (y2-y1) (z2-z1) ] ./ dis;

        % fill ball with point id step by step
        for k = 0 : sample_step : dis
            % center
            x0 = x1 + k * uv(1);
            y0 = y1 + k * uv(2);
            z0 = z1 + k * uv(3);
            % radius
            R = ( D1 + (D2-D1) * k / dis ) / 2;

            % the filling index of point
            if k < dis / 2
                id = chl(con);
            else
                id = par(con);
            end

            % fill in the index
            for x = floor(x0 - R) : ceil(x0 + R)
                if (x<1) | (x>Ms)
                    continue;
                end
                for y = floor(y0 - R) : ceil(y0 + R)
                    if (y<1) | (y>Ns)
                        continue;
                    end
                    for z = floor(z0 - R) : ceil(z0 + R)
                        if (z<1) | (z>Ks)
                            continue;
                        end
                        if stk(x,y,z) > 0
                            continue;
                        end
                        % check distance
                        d_squr = (x - x0) * (x - x0) +  (y - y0) * (y - y0) + (z - z0) * (z - z0) ;
                        if d_squr <= R * R
                            if ~isempty(options)&&strfind (options, '-i')
                                stk(x,y,z) = id;
                            else
                                stk(x,y,z) = 1;
                            end
                        end
                    end
                end
            end
        end
    end
end
end

