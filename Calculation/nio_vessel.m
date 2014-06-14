% nio_vessel
% Calculate parameters of vasculature based on one-root tree structure in 
% particular segments.
% 
% [fv len_den] = nio_vessel( vessel, nodes, area_xy, z_st, z_en, stpSize, options)
% -----------------------------------------------------------------------------
%
% Calculate fractional vascular volume(fv) and vascular length density(ld)
% in a one-root tree structure by sliding window, needs to appoint indexes 
% of the statistic nodes.
% calculate frustum-based volume:
% vol = (pi.*(D.^2 + D.*D(idpar) + D(idpar).^2).*len) / 12;
% 
% Input
% -----
% - vessel:: the vectorized vessels which has a tree structure
% - nodes:: indexes of the statistic nodes, to calculate particular
%           segments
% - area_xy:: the size of a single block(except the depth), unit: ��m^2
% - z_st:: representing the starting depth(z-coordinates) of real cortex in a single block
%          unit: ��m
% - z_en:: representing the ending depth(z-coordinates) of real cortex in a single block
%          unit: ��m
% - stpSize:: the statistic step length along depth, unit : ��m {DEFAULT: 5 ��m}
% - options:: string, it can be empty by default.
%        '-s' : plot
% 
% Output
% ------
% - fv:: linear array, fractional vascular volume in different depth
% - len_den:: linear array, vascular length density in different depth
%             unit: m/mm^3
% 
% Example
% -------
% [fv{1} len_den{1}] = nio_vessel( sample_tree{1}, 1:1000, 400, 65, 606, 5, '-s');
% 
% See also nio_vessel_single_block
% Uses idpar_tree

function [fv len_den] = nio_vessel( vessel, nodes, area_xy, z_st, z_en, stpSize, options)
%% parameters
winSize = 25;
zv = z_st : stpSize : z_en;
num_cal = length(zv);
idpar = idpar_tree(vessel, '-0');
vv = zeros(1, num_cal);
len_den = zeros(1, num_cal);
[~, ~, id_idpar] = intersect(nodes, idpar(nodes));
nodes = nodes(id_idpar);
%% analysis
for z = zv
    vz = 0;
    for n = nodes
        % get the x y z D
        if( vessel.Z(n) <= vessel.Z( idpar(n) ))
            x1 = vessel.X(n); x2 = vessel.X( idpar(n) );
            y1 = vessel.Y(n); y2 = vessel.Y( idpar(n) );
            z1 = vessel.Z(n); z2 = vessel.Z( idpar(n) );
            D1 = vessel.D(n); D2 = vessel.D( idpar(n) );
        else
            x1 = vessel.X( idpar(n) );    x2 = vessel.X(n);
            y1 = vessel.Y( idpar(n) );    y2 = vessel.Y(n);
            z1 = vessel.Z( idpar(n) );    z2 = vessel.Z(n);
            D1 = vessel.D( idpar(n) );    D2 = vessel.D(n);
        end
        
        if ( z1 >= z + winSize ) || ( z2 <= z - winSize )
            % outside the range
            continue;
        end
        len = sqrt( ( x1 - x2 ).^2 + ( y1 - y2 ) .^2 + ( z1 - z2 ).^2 ) ;
        if ( z1 >= z - winSize ) && ( z2 <= z + winSize )
            % totally inside
            % no operation, x y don't need any change
        elseif ( z1 >= z - winSize ) && ( z2 > z + winSize )
            len = len * ( z + winSize - z1) / (z2 - z1);
            % the following two equations are equal
            D2 = D2 + (D1 - D2)* ( z2 - z - winSize ) / (z2 - z1);
%             D2 = ( ( z2 - z -winSize ) * D1 + ( z + winSize - z1 )*D2 ) /
%             ( z2 -z1 );
        elseif ( z1 < z - winSize ) && ( z2 < z + winSize )
            len = len * ( z2 - z + winSize ) / (z2 - z1);
            % the following two equations are equal
            D1 = D1 + (D2 - D1)* ( z - winSize - z1 ) / (z2 - z1);
%             D1 = ( ( z1 - z -winSize ) * D2 + ( z + winSize - z2 )*D1 ) /
%             ( z1 -z2 );
        elseif ( z1 < z - winSize ) && ( z2 > z + winSize )
            len = len * ( 2 * winSize ) / (z2 - z1);
            D2tmp = D2 + (D1 - D2)* ( z2 - z - winSize ) / (z2 - z1);
            D1 = D1 + (D2 - D1)* ( z - winSize - z1 ) / (z2 - z1);
            D2 = D2tmp; 
        end
        vz = vz + (pi.*(D1.^2 + D1.*D2 + D2.^2).*len) / 12;
        len_den((z - zv(1))/stpSize + 1) = len_den((z - zv(1))/stpSize + 1) + len;
    end
    % add one element
    vv((z - zv(1))/stpSize + 1) = vz;
end
% fractional vascular volume (v/v)
fv = vv ./ (area_xy * 2 * winSize) ;
% Normalized vascular length (m/mm^3)
len_den = (len_den .* 1000) ./ (area_xy * 2 * winSize) ;

%% plot
if ~isempty(options)&&strfind (options, '-s')
    figure
    plot (zv-zv(1), fv, 'b','LineWidth',2);
    xlabel('Depth below cortical surface(\mum)');
    ylabel('Fractional vascular volume(v/v)');
end
end