% save trees or network as hoc format
% 
% nio_write_net_hoc( network, filename )
% --------------------------------------
% 
% Input
% -----
% - net:: the network data structure
% - filename:: the path and name of output hoc file 
%
% Output
% -----
% save a hoc file which can hold circles
% 
% Example
% -------
% nio_write_net_hoc( network, filename )
%
% Uses 
function nio_write_net_hoc( network, filename )
% by jpwu@CBMP, 2013/04/11

%% load test data
% clc
% clear
% load matlab.mat
% % transfer to input variables
% inputTrees{1} = vessel1;
% % inputTrees{1} = load_tree('sample2.mtr');
% filename = 'gap_bridged.hoc';


%% write network as the hoc format
fid = fopen( filename, 'w' );

% write the section coordinate and diameter
sn = length( network.sections )
for k = 1 : sn
    fwrite(fid, [ 'create section_' num2str(k-1) char(13), char(10) ], 'char' );
    fwrite(fid, [ 'section_' num2str(k-1) ' {' char(13), char(10) ], 'char');
    fwrite(fid, [ 'pt3dclear()' char(13), char(10)], 'char');
    for m = 1 : size( network.sections{k}, 1 )
        point = network.sections{k}(m,:);
        % note that the m/n was exchanged for amira coordinate system
        fwrite(fid, [ 'pt3dadd(' num2str(point(2)-1) ',' num2str(point(1)-1) ','...
            num2str(point(3)-1) ',' num2str(point(4)) ')' char(13) char(10) ], 'char');
    end
    fwrite(fid, ['}' char(13) char(10)], 'char');
    fwrite(fid, [char(13) char(10)], 'char');
end

% % write the connectivity relationship
% for  k = 1 : sn
%     % connect the start point of this section to the start point parent sections 
%     idx_con_s = find( network.dAs(k,:) );
%     for n = 1 : length( idx_con_s )
%         fwrite(fid, [ 'connect section_' num2str(k-1) '(0), section_' ...
%             num2str(idx_con_s(n)) '(0)' char(13) char(10) ] );        
%     end
%     
%     % connect the start point of this section to the end point parent sections 
%     idx_con_e = find( network.dAe(k,:) );
%     for n = 1 : length( idx_con_e )
%         fwrite(fid, [ 'connect section_' num2str(k-1) '(0), section_' ...
%             num2str(idx_con_e(n)) '(1)' char(13) char(10) ] );        
%     end
% end

%% write the connectivity relationship
[ c1, c2, be ] = find( network.con );
be = uint8(be);
for idx = 1 : length( be )
    fwrite(fid, [ 'connect section_' num2str(c1(idx)-1) '(' num2str((be(idx)-1)/2) '), section_' ...
            num2str( c2(idx)-1 ) '(' num2str(mod((be(idx)-1),2)) ')' char(13) char(10) ] );
end

fclose(fid);
