%% wave propagation in a local stack
% by jpwu, 2013/02/05

function [wave_front, mk_local_stk] = local_propagation( local_wave, mk_local_stk, m1, n1, k1 )
global mk_stk;

[Ml, Nl, Kl] = size( mk_local_stk );

% the wave front voxels
wave_front = nio_new_wave();

while( ~isempty( local_wave.voxelList ) )
    % initialize a temporal empty wave
    tmp_wave = nio_new_wave();
    tmp_wave.radius = local_wave.radius;

%         disp([ 'the voxel number of local wave: ' num2str(size(local_wave.voxelList,1)) ]);
   
    % region growing in the local cube
    for vn = 1 : size( local_wave.voxelList,1 )
        % mark the traced voxels to avoid repeated evaluation
        cord = local_wave.voxelList(vn,:) + [ m1, n1, k1 ];
        mk_stk( cord(1), cord(2), cord(3) ) = 1;

        % wave propagation through connectivity analysis, 26 connectivity
        cord = local_wave.voxelList(vn,:) + [0,0,1];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if cord(3) == Kl
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end

        cord = local_wave.voxelList(vn,:) + [0,1,0];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;            
            if local_wave.voxelList(vn, 2) + 1 == Nl
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end

        cord = local_wave.voxelList(vn,:) + [0,1,1];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(2) == Nl) || ( cord(3) == Kl)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end

        cord = local_wave.voxelList(vn,:) + [0,1,-1];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk(  cord(1)+m1, cord(2)+n1, cord(3)+k1  )
            mk_stk(  cord(1)+m1, cord(2)+n1, cord(3)+k1  )= 1;
            if (cord(2) == Nl) || ( cord(3) == 1)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end

        cord = local_wave.voxelList(vn,:) + [0,-1,1];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk(  cord(1)+m1, cord(2)+n1, cord(3)+k1  )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(2) == 1) || ( cord(3) == Kl)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end

        cord = local_wave.voxelList(vn,:) + [1,0,0];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(1) == Ml)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end

        cord = local_wave.voxelList(vn,:) + [1,0,1];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(1) == Ml) || ( cord(3) == Kl)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end

        cord = local_wave.voxelList(vn,:) + [1,0,-1];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(1) == Ml) || ( cord(3) == 1)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end

        cord = local_wave.voxelList(vn,:) + [-1,0,1];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(1) == 1) || ( cord(3) == Kl)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end

        cord = local_wave.voxelList(vn,:) + [1,1,0];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(1) == Ml) || ( cord(2) == Nl)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end

        cord = local_wave.voxelList(vn,:) + [1,-1,0];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(1) == Ml) || ( cord(2) == 1)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end

        cord = local_wave.voxelList(vn,:) + [-1,1,0];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(1) == 1) || ( cord(2) == Nl)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end

        cord = local_wave.voxelList(vn,:) + [1,1,1];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(1) == Ml) || ( cord(2) == Nl) || ( cord(3) == Kl)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end

        cord = local_wave.voxelList(vn,:) + [1,1,-1];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(1) == Ml) || ( cord(2) == Nl) || ( cord(3) == 1)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end

        cord = local_wave.voxelList(vn,:) + [1,-1,1];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(1) == Ml) || ( cord(2) == 1) || ( cord(3) == Kl)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end

        cord = local_wave.voxelList(vn,:) + [1,-1,-1];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(1) == Ml) || ( cord(2) == 1) || ( cord(3) == 1)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end

        cord = local_wave.voxelList(vn,:) + [-1,1,1];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(1) == 1) || ( cord(2) == Nl) || ( cord(3) == Kl)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end

        cord = local_wave.voxelList(vn,:) + [-1,1,-1];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(1) == 1) || ( cord(2) == Nl) || ( cord(3) == 1)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end

        cord = local_wave.voxelList(vn,:) + [-1,-1,1];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(1) == 1) || ( cord(2) == 1) || ( cord(3) == Kl)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end

        % subtract direction
        cord = local_wave.voxelList(vn,:) - [0,0,1];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if cord(3) == 1
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end
        
        cord = local_wave.voxelList(vn,:) - [0,1,0];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if cord(2) == 1
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end
        
        cord = local_wave.voxelList(vn,:) - [0,1,1];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(2) == 1) || ( cord(3) == 1)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end
        
        cord = local_wave.voxelList(vn,:) - [1,0,0];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(1) == 1)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end
        
        cord = local_wave.voxelList(vn,:) - [1,0,1];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(1) == 1) || ( cord(3) == 1)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end
        
        cord = local_wave.voxelList(vn,:) - [1,1,0];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(1) == 1) || ( cord(2) == 1)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end
        
        cord = local_wave.voxelList(vn,:) - [1,1,1];
        if mk_local_stk( cord(1), cord(2), cord(3) ) %& ~mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )
            mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 )= 1;
            if ( cord(1) == 1) || ( cord(2) == 1) || ( cord(3) == 1)
                % hit the boundary, mark the voxel as wave front
                mk_local_stk( cord(1), cord(2), cord(3) ) = 2;
                wave_front.voxelList = [ wave_front.voxelList; cord ];
            else
                tmp_wave.voxelList = [tmp_wave.voxelList; cord ];
                mk_stk( cord(1)+m1, cord(2)+n1, cord(3)+k1 ) = 1;
            end
        end
    end
    local_wave = tmp_wave;
end