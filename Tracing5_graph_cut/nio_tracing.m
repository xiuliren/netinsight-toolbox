%% automatic iterative tracing from seeds
% by jpwu, 2013/02/27

function network = nio_tracing( stk )
%% parameters
secLenThreshold = 5;
% seed interval
SeedSecItv = 5;

%% initialize the data
% global stk
% the marker stack 
global mk_stk;
mk_stk = false(size(stk));

% generate a new network
network = nio_new_network();

%% visualization system
mip = max(stk, [], 3);
imshow( mip );  hold on;
pause(0.1); % draw this figure

%% seeds detection, x, y, z, radius
seedList = get_seed_list( stk, SeedSecItv );

%% wave propagation from seeds
global dsl;
dsl = zeros(10000000, 6);

disp('---- wave propagation from seeds ...')
for e = 1 : size(seedList,1)
    disp([ 'the sees:   ' num2str(e) '  in  ' num2str(length(seedList)) ]);
    network = tracing_propagation( seedList(e, :), network );
end

network.sn = length(network.sections);
network = nio_build_net_connectivity(network);
%% preprocessing of the network
% short sectiong prunning, may modified to prun short leaves
network = nio_short_prunning( network, secLenThreshold );

% global optimal connction and prunning to achieve maximum connectivity, 
% minimum volume redundency, and minimum curvature.
% I have to design a criterial for optimization.