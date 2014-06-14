%% build a vascular density map
% by jpwu, 2013/04/02
% 
% dm_ld = nio_density_map_net( network, rcs, voxel_size )
% -----------------------------------------------------------------------
% 
% Get the fractional vascular volume(fvv) and normalized length density(nld)
% in a cuboid. 
% 
% Input
% -----
% - network:: the vectorized network
% - rcs:: the radius of the cubic, 25 means the cubic is 50X50X50
% - voxel_size:: voxel size of the density map
% 
% Output
% ------
% - dm_ld:: the normalized length density map, unit: m/mm^3
% 
% Attention
% ---------
% 1. the voxel size of cuboid must be 1X1X1 micron for correct density
% 2. the coordinate system: Y--N, X--M, Z--K
% 
% Example
% -------
% dm_ld = nio_density_map_net( network, 25, 1 )
% 
% See also nio_intersection_stk nio_density_map
% Uses idpar_tree
function dm_ld = nio_density_map_net_V2( network, rcs, voxel_size )
% clc
% clear
% close all
% 
% %% parameters
% % directory
% SrcNet = '../Data/network.mat';
% SrcHoc = '../Data/tracing_result_section.hoc';
% DstMap = '../Data/section_map.mat';
% 
% % roi cubic size, micron
% rcs = 50/2;
% voxel_size = 1;
% 
% %% load data and initialization
% % load( SrcNet );
% disp('------reading hoc file ...')
% % network = nio_read_hoc( SrcHoc );
% load( SrcNet )

sn = length(network.sections);

%% find the boundingbox of each section
disp('------building the bounding box...')
bb = zeros(sn, 6);
for si = 1 : sn
    sec = network.sections{si};
    bb( si, 1:3 ) = min( sec( :, 1:3 ),[],1 );
    bb( si, 4:6 ) = max( sec( :, 1:3 ),[],1 );
end

bbc = ( bb(:,1:3) + bb(:,4:6) ) ./2;
bbw = bb(4:6);
%% get the total range
tbb = zeros( 1, 6 );
tbb(1:3) = min( bb( :, 1:3 ),[],1 );
tbb(4:6) = max( bb( :, 4:6 ),[],1 );

%% compute the density map
dm_ld = zeros( uint16( (tbb(4)+rcs)/voxel_size), uint16((tbb(5)+rcs)/voxel_size), 'double' );

for m = tbb(1)+rcs*2 : voxel_size : tbb(4)-rcs*2 %rcs*2+1 :rcs*2: M-rcs*2
    m
    for n = tbb(2)+rcs*2 :voxel_size: tbb(5)-rcs*2
        % initialize 
        len = 0;
        
        % find the inner and crossing sections     
        overlap = abs(m-bbc(:,1))<bb(:,4)-bb(:,1)+rcs*2+1 & abs(n-bbc(:,2))<bb(:,5)-bb(:,2)+rcs*2+1;
        % the inner section index
        in_sec_idx = find( overlap );
        
        % compute the density
        for sk = 1 : length( in_sec_idx )
            sidx = in_sec_idx(sk);
            sec = network.sections{sidx};
            
            Nn = size(sec,1);
            
            for idx = 1 : Nn-1
                % get the coordinate of the start and end points
                m1 = sec( idx, 1 ); m2 = sec( idx+1, 1);
                n1 = sec( idx, 2 ); n2 = sec( idx+1, 2);
                k1 = sec( idx, 3 ); k2 = sec( idx+1, 3);
                
                % get the total length
                tl = sqrt( (m2-m1)*(m2-m1) + (n2-n1)*(n2-n1) + (k2-k1)*(k2-k1) );
                
                if (m1>m-rcs && m1<m+rcs && n1>n-rcs && n1<n+rcs...
                        && m2>m-rcs && m2<m+rcs && n2>n+rcs && n2<n+rcs) 
                    % two nodes are all in the sliding window
                    len = len + tl;
                elseif (m1>m-rcs && m1<m+rcs && n1>n-rcs && n1<n+rcs)
                    % the first node is in the sliding window, and the
                    % second is not
                    % get the step size
                    dm = (m2-m1)/tl; dn = (n2-n1)/tl;   dk = (k2-k1)/tl;
                    % start from the first node
                    ms = m1;    ns = n1;    ks = k1;
                    % find the cross point
                    for step = 1 : tl
                        ms = ms + dm;   ns = ns + dn;   ks = ks + dk;
                        if ms<m-rcs || ms>m+rcs || ns<n-rcs || ns>n+rcs
                            break;
                        end
                    end
                    len = len + sqrt( (ms-m1)*(ms-m1) + (ns-n1)*(ns-n1) + (ks-k1)*(ks-k1) );
                else
                    % the second node is in the sliding window
                    % and the first node is not
                    dm = (m1-m2)/tl; dn = (n1-n2)/tl;   dk = (k1-k2)/tl;
                    % start from the second node
                    ms = m2;    ns = n2;    ks = k2;
                    % find the cross point
                    for step = 1 : tl
                        ms = ms + dm;   ns = ns + dn;   ks = ks + dk;
                        if ms<m-rcs || ms>m+rcs || ns<n-rcs || ns>n+rcs
                            break;
                        end
                    end
                    len = len + sqrt( (ms-m2)*(ms-m2) + (ns-n2)*(ns-n2) + (ks-k2)*(ks-k2) );
                end
            end
        end
        dm_ld( floor(m/voxel_size), floor(n/voxel_size) ) = len/((2*rcs+1)*(2*rcs+1));
    end
end

%% render the density map
% imagesc( dm_ld );
