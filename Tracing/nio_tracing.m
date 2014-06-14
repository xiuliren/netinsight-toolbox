%% automatic iterative tracing from seeds
% by jpwu, 2013/02/27

function network = nio_tracing( stk, seedList )
%% parameters
% RMIN = 2;
% RMAX = 50;

%% processing
% the marker stack 
% global mk_stk;
mk_stk = zeros( size(stk), 'uint8');
mk_stk = logical(mk_stk);

% generate a new network
network = nio_new_network();

% % create balls as a lookup table
% global balls;
% % balls = create_lookup_balls( RMAX );
load('balls.mat');

% create scooping distance lookup table
sdtable = zeros(200, 200, 200, 'double');
for m = 1 : size(sdtable, 1)
    for n = 1 : size( sdtable, 2 )
        for k = 1 : size( sdtable, 3 )
            sdtable(m,n,k) = sqrt( (m-1)*(m-1) + (n-1)*(n-1) + (k-1)*(k-1) );
        end
    end
end

% wave propagation from seeds
disp('---- voxel scooping from seeds ...')
% matlabpool
for e = 1 : size(seedList,1)
    disp([ 'Tracing seed progress: ' num2str(e) ' in ' num2str( size(seedList,1) )]);
%     network = wave_propagation( seedList(e, :), network );
    [network mk_stk]= voxelscooping( stk, mk_stk, balls, sdtable, seedList(e, :), network);
end
% matlabpool('close');

%% create balls in cube as a lookup table to accellerate computing
% function balls = create_lookup_balls( RMAX )
% for r = 1 : RMAX
%     % cube size
%     cs = 2*r+1;
%     % initialize the cube
%     b = zeros(cs,cs,cs, 'uint8');
%     % filling the ball
%     for m = 1 : cs
%         for n = 1 : cs
%             for k = 1 : cs
%                 d2 = (m-r-1)*(m-r-1) + (n-r-1)*(n-r-1) + (k-r-1)*(k-r-1);
%                 if d2 < r*r
%                     b(m,n,k) = 1;
%                 end
%             end
%         end
%     end
%     % store the ball in a cell
%     balls{r} = logical(b);
% end
% return;