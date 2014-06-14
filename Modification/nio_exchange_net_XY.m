% nio_exchange_net_XY
% exchange the x and y coordinate of a network
% 
% network = nio_exchange_net_XY( network )
% --------------------------------------
% 
% Input
% -----
% - network:: the interconnected network
% 
% Output
% -----
% - network:: the new network whose X and Y coordinate were exchanged
%
% Example
% -------
% 
% 
% Uses 

function network = nio_exchange_net_XY( network )

for k = 1 : length(network.sections)
    sec = network.sections{k};
    tmp = sec( :,1 );
    sec(:,1) = sec(:,2);
    sec(:,2) = tmp;
    network.sections{k} = sec;
end