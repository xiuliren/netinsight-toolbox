% nio_read_hoc  
% Read network from the hoc format.
% 
% network  = nio_read_hoc (tname)
% ----------------------------------------------
%
% Loads the hoc file and build the network data structure.
% 
% Input
% -----
% - tname::string: name of the file to be loaded, including the extension.
%                 {DEFAULT : open gui fileselect, replaces format entry}
%
% Output
% ------
% - network:: the network data structure
%
% Example
% -------
% network = nio_read_hoc('Y.swc');
%
% See also 

function network = nio_read_hoc (tname)

%% initialization for testing
% clc
% clear
% tname = '../Data/Y.hoc';

%% parameter preperation
if (nargin<1)||isempty(tname),
     [tname path] = uigetfile ({'*hoc', 'TREES formats (*.hoc)'}, ...
        'Pick a file',...
        'multiselect', 'off');
    if tname == 0,
        varargout {1} = []; varargout {2} = []; varargout {3} = [];
        return;
    end
else
    path = '';
end
format = tname  (end - 3 : end); % input format from extension:
% extract a sensible name from the filename string:
nstart = unique ([0 strfind(tname, '/') strfind(tname, '\')]);
name   = tname  (nstart (end) + 1 : end - 4);

if (nargin<2)||isempty(options)
    if strcmp (format, 'hoc')
        options = '-r';
    else
        options = '';
    end
end

%% initialization
network = nio_new_network();

%% start reading
% open file
hoc = textread( tname,'%s','delimiter','\n' );

% scan the string line 
sec_idx = 0;
sct = [];
for ln = 1 : length( hoc )
    line = hoc{ln};
    if ~isempty( strfind( line, '{' ) )
        % a new section
%         p1 = strfind(line, 'section_') + 8;
%         p2 = strfind(line, ' {') - 1;
%         sec_idx = str2num( line( p1:p2 ) ) + 1;
        sec_idx = sec_idx + 1;
        sct = [];
    elseif ~isempty( strfind( line, 'pt3dadd(' ) )
        % a new node in section
        ps = strfind( line, '(' ) + 1;
        comma = strfind( line, ',' ) + 1;
        pe = strfind( line, ')' ) + 1;
        % transform the coordinate for matlab
        n = str2num( line( ps:comma(1)-2 ) ) +1;
        m = str2num( line( comma(1):comma(2)-2 )) +1;
        k = str2num( line( comma(2):comma(3)-2 )) +1;
        r = str2num( line( comma(3):pe-2 ));
        sct = [ sct; m,n,k,r, 0 ];
%     elseif ~isempty( strfind( line, '' ) )
      elseif length(line) == 0;
        % end of section
        network.sections{ sec_idx } = sct;
    elseif ~isempty( strfind( line, 'connect' ) )
        % connection
        p1 = strfind( line, 'section_' ) + 8;
        bral = strfind( line, '(') + 1;
        brar = strfind( line, ')') + 1;
        % get the section index and begin/end
        idx1 = str2num( line( p1(1):bral(1)-2 ) ) +1;
        be1 = uint8( str2num( line( bral(1):brar(1)-2 ) ) );
        idx2 = str2num( line( p1(2):bral(2)-2 ) ) +1;
        be2 = uint8( str2num( line( bral(2):brar(2)-2 ) ) );
        % build the connectivity
        network.con(idx1, idx2) = be1*2 + be2 + 1;
    else
        continue;
    end
    
end
