% nio_new_network
% generate a new empty network
% 
% newnet = nio_new_network()
% --------------------------------------
% 
% Input
% -----
%
% Output
% -----
% - newnet:: the generated new network
%
% Example
% -------
% 
% 
% Uses 

function newnet = nio_new_network()

newnet.sections = [];
% newnet.sn = 0;
% newnet.dAe = sparse( 1,1 );
% newnet.dAs = sparse( 1,1 );
% newnet.terminalSec = [];

% new connection structure
newnet.con = sparse( 1,1 );
