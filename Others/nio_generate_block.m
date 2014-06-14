% nio_generate_block
% Generate binary images for marking every block.
% 
% [ BM, BM_b, BM_s ] = nio_generate_block( Ms, Ns, Block_Size, options, barrel_mask, barrel_mask_ev )
% --------------------------------------------------------------------------------------------------
% 
% Generate binary images for marking every block. One block can be a
% rectangle, a barrel-shape field,a barrel or septa field in a rectangle field. 
% The shape depend on the options. 
% 
% Input
% -----
% - Ms:: the size of binary image, x-coordinate
% - Ns:: the size of binary image, y-coordinate
% - Block_Size:: the size of a block assigned from cortex outside surface,
%                unit : ¦Ìm
% - options:: string, it can be empty by default.{DEFAULT: rectangle block}
%        '-b' : every block represent a barrel-shape field
%        '-s' : one block represent a rectangle, a barrel field or septa
% - barrel_mask:: a labeled image marking different barrel field(options:'-b')
%                 a binary image marking all barrel fields(options:'-s')
% - barrel_mask_ev:: a binary image marking PMBSF(barrels and septa)(options:'-s')
% 
% Output
% ------
% - BM:: cell array, binary images marking rectangle block or one
%        barrel field(options:'-b')
% - BM_b:: cell array, binary images marking barrel field in a rectangle field.(options:'-s')
% - BM_s:: cell array, binary images marking septa in a rectangle field.(options:'-s')
% 
% Example
% -------
% BM{1} = nio_generate_block( 560, 600, 50, '-b', barrel_mask );
% 
% See also nio_soma2tree

function [ BM, varargout ] = nio_generate_block( Ms, Ns, Block_Size, options, varargin )
if ~isempty(options)&&strfind (options, '-b')
    barrel_mask = varargin{1};
    num_barrel = max(max(barrel_mask));
    BM = cell(1, num_barrel); 
    idx_x = cell(1, num_barrel);
    idx_y = cell(1, num_barrel);
    area = zeros(1, num_barrel);
    for ward  = 1 : num_barrel
        [idx_x{ward}, idx_y{ward}] = find(barrel_mask == ward);
        area(ward) = length(idx_x{ward});
        BM{ward} = zeros(Ms, Ns);
        for i = 1 : area(ward)
            BM{ward}(idx_x{ward}(i), idx_y{ward}(i)) = 1;
        end
    end
else
    n_i = fix(Ms / Block_Size);
    n_j = fix(Ns / Block_Size);
    num_part = n_i * n_j; 
    BM = cell(1, num_part); % whole region
    for i = 1 : n_i
        for j = 1 : n_j
            idx = n_j * ( i - 1 ) + j;
            BM{idx} = zeros(Ms, Ns);
            BM{idx}(((i-1)*Block_Size+1):(i*Block_Size), ((j-1)*Block_Size+1):(j*Block_Size)) = 1; % mark every block
        end
    end
    if ~isempty(options)&&strfind (options, '-s')
        barrel_mask = varargin{1};
        barrel_mask_ev = varargin{2};
        BM_b = cell(1, num_part);
        BM_s = cell(1, num_part);
        for idx = 1 : num_part
            BM_b{idx} = BM{idx} & barrel_mask;
            BM_s{idx} = BM{idx} & ~BM_b{idx} & barrel_mask_ev;
        end
        varargout{1} = BM_b;
        varargout{2} = BM_s;
    end
end
end

