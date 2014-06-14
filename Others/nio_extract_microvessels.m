% nio_remove_bigvessels
% Remove big vessels, just keep small vessels(microvasculature)
% 
% vessel_c = nio_extract_microvessels( vessel, T_D )
% ------------------------------------------
% 
% Generate microvasculature, remove big vessels in a one-root tree 
% structure by evaluate 0 to dA matrix.
% 
% Input
% -----
% - vessel:: the vectorized vessels which has a tree structure
% - T_D:: the threshold of diameter, the vessel segments whose radius is
%        bigger than T_D/2 will be removed. unit: ¦Ìm
% 
% Output
% ------
% - vessel_c:: the big vessel segments are removed
% 
% Example
% -------
% sample_tree_c = nio_extract_microvessels( sample_tree, 5 );
% 
% Uses idpar_tree

function vessel_c = nio_extract_microvessels( vessel, T_D )
vessel_c = vessel;
N = length( vessel.X );
idpar = idpar_tree (vessel, '-0'); % vector containing index to direct parent
for n = 1 : N
    if (idpar(n) == 0)
        continue;
    end
    if ( (vessel.D(n) + vessel.D(idpar(n))) / 2 > T_D )
        vessel_c.dA(n, idpar(n)) = 0;
    end
end
end
