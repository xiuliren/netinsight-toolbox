%% converte swc file to stack, mark the node id
% by jpwu, 2011/12/13

function idx_stk = nio_swc2stkcube(vessel1 , Ms, Ns, Ks  )
% parameters

% the resampled node distance
Resample_Length = 1;

%% build the point index stack
disp('-------- getting binary stack ...')
idx_stk = zeros (Ms, Ns, Ks, 'double');

% points connectivity, child and parent
[chl par] = find(vessel1.dA == 1);

for con = 1 : length( chl )
%     con
    % the child and parent point
    x1 = vessel1.X( chl(con) );
    y1 = vessel1.Y( chl(con) );
    z1 = vessel1.Z( chl(con) );
    D1 = vessel1.D( chl(con) );
    
    x2 = vessel1.X( par(con) );
    y2 = vessel1.Y( par(con) );
    z2 = vessel1.Z( par(con) );
    D2 = vessel1.D( par(con) );
    
    % point distance
    dis = sqrt( (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2) );
    if dis < 1
        continue;
    end
    
    % unit vector
    uv = [ (x2-x1) (y2-y1) (z2-z1) ] ./ dis;
    
    % fill ball with point id step by step
    for k = 0 : dis
        % center
        x0 = x1 + k * uv(1);
        y0 = y1 + k * uv(2);
        z0 = z1 + k * uv(3);
        % radius
        R = ( D1 + (D2-D1) * k / dis ) / 2*1.2;
        if R>50
            R = 50;
        end
        % the filling index of point
        if k < dis / 2
            id = chl(con);
        else
            id = par(con);
        end
        
        % fill in the index
        
        xstart = floor((x0-R));        xend = floor((x0+R));
        ystart = floor((y0-R));        yend = floor((y0+R));
        zstart = floor((z0-R));        zend = floor((z0+R));
        
        if xstart<1    
            xstart=1; 
        end
        if xend>Ms
            xend = Ms;
        end
        if ystart<1    
            ystart=1; 
        end
        if yend>Ns
            yend = Ns;
        end 
        if zstart<1    
            zstart=1; 
        end
        if zend>Ks
            zend = Ks;
        end
        
        idx_stk(xstart:xend,ystart:yend,zstart:zend) = id;
%         for x = floor(x0 - R) : ceil(x0 + R)
%             if (x<1) | (x>Ms)
%                 continue;
%             end
%             for y = floor(y0 - R) : ceil(y0 + R)
%                 if (y<1) | (y>Ns)
%                     continue;
%                 end
%                 for z = floor(z0 - R) : ceil(z0 + R)
%                     if (z<1) | (z>Ks)
%                         continue;
%                     end
%                     if idx_stk(x,y,z) > 0
%                         continue;
%                     end
%                     % check distance
%                     d_squr = (x - x0) * (x - x0) +  (y - y0) * (y - y0) + (z - z0) * (z - z0) ;
%                     if d_squr <= R * R
%                         idx_stk(x,y,z) = id;
%                     end
%                 end
%             end
%         end
    end
end

