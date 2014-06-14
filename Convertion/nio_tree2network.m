% transfer datastructure of tree to network
% 
% network = nio_tree2network( inputTrees )
% --------------------------------------
% 
% Input
% -----
% - inputTrees:: the vectorized trees which has a tree structure
%               note that the points may have multiple parents
%
% Output
% -----
% - network:: a new data structure containing sections and connections
% 
% Example
% -------
% network = nio_tree2network( inputTrees )
%
% Uses 

function network = nio_tree2network( inputTrees )

%% load test data only for debug
% clc
% clear
% % load matlab.mat
% % transfer to input variables
% % inputTrees{1} = vessel1;
% inputTrees{1} = load_tree('sample2.mtr');

%% parameters
% minimum length, the length was approximately represented by point number
MinLen = 2;
network = nio_new_network();
%% transfer to cells
if ~iscell( inputTrees )
    tmp = cell(1);
    tmp{1} = inputTrees;
    inputTrees = tmp;
end
%% dipart to sections
% estimate the maximun sections by point number
pn = 0;
for n = 1 : length( inputTrees )
    pn = pn + length( inputTrees{n}.R );
end
% initialize the network
% network = structure();
network.sections = cell( floor(pn), 1 );

% transfer to network
% initialize the section number
sn = 0;
% get the sections and connection relationship
for tn = 1 : length( inputTrees )
    % get a separated tree
    tree = inputTrees{ tn };
    % divide tree to sections and connections
    childNum  = sum( tree.dA, 1 );
    parentNum = sum( tree.dA, 2 )';
    branchPoints = ( childNum >= 2 ) | ( parentNum >= 2 );
    rootPoints = ( parentNum == 0 );
    % parent point index.
    idx_par = idpar_tree( tree, '-0' );
    
    % the terminal sections
    % terminal points index
    idx_t = find( childNum == 0 );
    for n = 1 : length( idx_t )
        % every single terminal section
        idx = idx_t( n );
        sn = sn + 1;
        % inserte the terminal point
        network.sections{ sn } = [ network.sections{sn}; tree.Y(idx)+1, tree.X(idx)+1, tree.Z(idx)+1, tree.D(idx)/2 ];
        % insert other points
        while ~ ( branchPoints(idx) | rootPoints(idx) ) 
            % not a branch point or root
            idx = idx_par( idx );
            network.sections{ sn } = [ network.sections{sn}; tree.Y(idx)+1, tree.X(idx)+1, tree.Z(idx)+1, tree.D(idx)/2  ];
        end
    end
%     network.terminalSec = ones( 1, length(idx_t), 'uint8' );
    
    % the internal branch sections
    % branch points index
    idx_b = find( branchPoints );
    for n = 1 : length( idx_b )
        idx = idx_b( n );
        % parent index, may be multiple parent points !
        idx_ps = find( tree.dA(idx,:) );
        for m = 1 : length( idx_ps )
            % a new section
            sn = sn + 1;
%             network.terminalSec = [ network.terminalSec, 0 ];
            % insert the first branch point and the parent
            network.sections{ sn } = [ network.sections{ sn };  tree.Y(idx)+1, tree.X(idx)+1, tree.Z(idx)+1, tree.D(idx)/2  ];
            % set the corrent index to this branck point
            idx = idx_ps(m);
%             if ( idx ~= 0 )
%                 % not root
%                 network.sections{ sn } = [ network.sections{ sn };   tree.Y(idx)+1, tree.X(idx)+1, tree.Z(idx)+1, tree.D(idx)/2  ];
%             else
%                 % if this is root
%                 continue;
%             end
            % insert other points
            while ~ ( branchPoints(idx) | rootPoints(idx) )
                % not a branch point or root
                idx = idx_par( idx );
                network.sections{ sn } = [ network.sections{ sn };  tree.Y(idx)+1, tree.X(idx)+1, tree.Z(idx)+1, tree.D(idx)/2  ];
            end
        end
    end
end
% 
% eliminate the empty cells
network.sections( (sn+1):(length(network.sections)) ) = [];
% eliminate the short sections, lenth < MinLen 
pointNum = zeros(1,sn, 'uint32');
for k = 1 : sn
    pointNum(k) = size( network.sections{k}, 1 );
end
idx_short = find( pointNum < MinLen );
network.sections( idx_short ) = [];
sn = sn - length( idx_short );

%% estimate the connection relationship, represented by directed adjacency matrix
% % initialize the directed adjacency matrix
% network.dAs = sparse( sn, sn );
% network.dAe = sparse( sn, sn );
% % collect the start and end points
% sps = zeros(sn, 4, 'double');
% eps = zeros(sn, 4, 'double');
% for n = 1 : sn
%     sps(n,:) = network.sections{n}(1,:);
%     eps(n,:) = network.sections{n}(end,:);
% end
% 
% % find the connectivity relationship by overlap points
% for n = 1 : sn
%     % connect start points to start points
%     idx_con_s = find( ( sps(n,1) == sps(:,1) ) ...
%         & ( sps(n,2) == sps(:,2) ) & ( sps(n,3) == sps(:,3) ) );
%     % eliminate self connectivity
%     idx_con_s( find( idx_con_s == n ) ) = [];
%     % connect to sections start point
%     network.dAs( n, idx_con_s ) = 1 ;
%     
% %     % connect end points to start points
% %     idx_con_e = find( ( eps(n,1) == sps(:,1) ) ...
% %         & ( eps(n,2) == sps(:,2) ) & ( eps(n,3) == sps(:,3) ) );
% %     % eliminate self connectivity
% %     idx_con_e( find( idx_con_e == n ) ) = [];
% %     % connect to sections start point
% %     network.dAs( idx_con_s, n ) = 1 ;
%     
%     % connect start points to end points
%     idx_con_s = find( ( sps(n,1) == eps(:,1) ) ...
%         & ( sps(n,2) == eps(:,2) ) & ( sps(n,3) == eps(:,3) ) );
%     % eliminate self connectivity
%     idx_con_s( find( idx_con_s == n ) ) = [];
%     % connect to sections end point
%     network.dAe( n, idx_con_s ) = 1 ;
%     
%     
% %     % connect end points to end points
% %     idx_con_e = find( ( eps(n,1) == eps(:,1) ) ...
% %         & ( eps(n,2) == eps(:,2) ) & ( eps(n,3) == eps(:,3) ) );
% %     % eliminate self connectivity
% %     idx_con_e( find( idx_con_e == n ) ) = [];
% %     % connect to sections end point
% %     network.dAe( idx_con_e, n ) = 1 ;
% end

%% estimate the connectivity relationship among sections 
%% for new network data structure
% collect the start and end points
sps = zeros(sn, 4, 'double');
eps = zeros(sn, 4, 'double');
for n = 1 : sn
    sps(n,:) = network.sections{n}(1,:);
    eps(n,:) = network.sections{n}(end,:);
end

% find the connectivity relationship by overlap points
for n = 1 : sn
    % connect start points to start points
    idx_con_s = find( ( sps(n,1) == sps(:,1) ) ...
        & ( sps(n,2) == sps(:,2) ) & ( sps(n,3) == sps(:,3) ) );
    % eliminate self connectivity
    idx_con_s( idx_con_s == n ) = [];
    % connect to sections start point
    network.con( n, idx_con_s ) = 1 ;
    
%     % connect end points to start points
%     idx_con_e = find( ( eps(n,1) == sps(:,1) ) ...
%         & ( eps(n,2) == sps(:,2) ) & ( eps(n,3) == sps(:,3) ) );
%     % eliminate self connectivity
%     idx_con_e( find( idx_con_e == n ) ) = [];
%     % connect to sections start point
%     network.dAs( idx_con_s, n ) = 1 ;
    
    % connect start points to end points
    idx_con_s = find( ( sps(n,1) == eps(:,1) ) ...
        & ( sps(n,2) == eps(:,2) ) & ( sps(n,3) == eps(:,3) ) );
    % eliminate self connectivity
    idx_con_s( idx_con_s == n ) = [];
    % connect to sections end point
    network.con( n, idx_con_s ) = 2 ;
    
    % connect end points to end points
    idx_con_e = find( ( eps(n,1) == eps(:,1) ) ...
        & ( eps(n,2) == eps(:,2) ) & ( eps(n,3) == eps(:,3) ) );
    % eliminate self connectivity
    idx_con_e( idx_con_e == n ) = [];
    % connect to sections end point
    network.con( idx_con_e, n ) = 4 ;
end