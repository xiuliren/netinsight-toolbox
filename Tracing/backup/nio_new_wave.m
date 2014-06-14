%% generate a new wave data structure
% by jpwu, 2013/02/27
function wave = nio_new_wave( seed )

% handle the number of input arguments
if nargin < 1
    wave.voxelList = [];
    wave.radius = 3;
else
    wave.voxelList = seed(1:3);
    wave.radius = seed(4);
end

