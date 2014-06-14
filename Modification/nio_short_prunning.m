%% eliminate short sections
function network = nio_short_prunning( network, Tl )
% by jpwu, 2013/03/18

% %% parameters
% SrcSWC = '../Data/fMOSTstackBWcropEdge.Smt.SptGraph.swc';
% DstHoc = '../Data/trimed_result.hoc';
% 
% % length threshold
% Tl = 30;

%%% read swc file
% neurons = nio_load_tree( SrcSWC );
% % transform to network data structure
% network = nio_tree2network( neurons );
% 
% load fMOST_neurons.mat
% 
% % transform to matlab coordinat
% network = nio_exchange_net_XY( network );

%% detect the short sections
% the length vector
lv = zeros( 1, network.sn );

% compute the length
for si = 1 : network.sn
    sec = network.sections{ si };
    % the length of current section
    sl = 0;
    for ni = 1 : size(sec,1)-1
        len2 =   (sec(ni,1) - sec(ni+1,1))^2 + ...
                (sec(ni,2) - sec(ni+1,2))^2 + ...
                (sec(ni,3) - sec(ni+1,3))^2;
        sl = sl + sqrt( len2 );
    end
    lv(si) = sl;
end

%% eliminate short sections
% mark the section index to be removed
% eliminated section index
esi = find(lv < Tl);

% remove the marked sections
network.sections(esi) = [];
network.sn = network.sn - length(esi);
network = nio_build_net_connectivity(network);
% 
% %% save the result
% nio_save_hoc( network, DstHoc );