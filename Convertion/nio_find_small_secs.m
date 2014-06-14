% nio_find_small_secs
% recursively find interconnected sections with smaller sections
% 
% [ con_secs_idx tm_secs_idx ] = nio_find_small_secs( network, key_idx,
% con_secs_idx, tm_secs_idx, checked_secs_bin, Td )
% --------------------------------------
% 
% Input
% -----
% - network:: the interconnected network
% - key_idx:: the index of key section 
% - con_secs_idx:: the index of connected sections that is not repeated
% - tm_secs_idx:: the index of terminal capillary sections
% - checked_secs_bin:: a binary vector that marks the sections which are
% already checked
% - Td:: a diameter threshold that is greater than the average diameter of  capilary sections
%
% Output
% -----
% - con_secs_idx:: the index of connected sections that is not repeated
% - tm_secs_idx:: the index of terminal capillary sections
% - checked_secs_bin:: a binary vector that marks the sections which are
% already checked
%
% Example
% -------
% 
% 
% Uses 

function [ con_secs_idx tm_secs_idx checked_secs_bin ] = nio_find_small_secs( network, key_idx, con_secs_idx, tm_secs_idx, checked_secs_bin, Td )

% find candidate connected neighbor sections
cndt = nio_get_con_secs( network, key_idx );

% the new sections
new_secs_idx = [];

% flag of terminal section
tm_flag = true;
for k = 1 : length( cndt )
    idx = cndt(k);
    if ~checked_secs_bin( idx )
        % a section not been checked
        checked_secs_bin( idx ) = 1;
        if (network.avd( idx ) < network.avd( key_idx )) | (network.avd(idx) < Td)
            % this key section is not a terminal section
            tm_flag = false;
            % add this section
            cd = network.avd( idx );
            con_secs_idx = [ con_secs_idx idx ];
            new_secs_idx = [ new_secs_idx idx ];
        end
    end
end

% mark the terminal section
if tm_flag
    tm_secs_idx = [ tm_secs_idx key_idx ];
end

% recurent finding
for k = 1 : length( new_secs_idx )
    [ con_secs_idx tm_secs_idx checked_secs_bin] = nio_find_small_secs( network, new_secs_idx(k), con_secs_idx, tm_secs_idx, checked_secs_bin, Td );
end