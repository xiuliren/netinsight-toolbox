%% write image stack
function nio_write_stk(stk, Pre )
[M N K] = size( stk );

for k = 1 : K
    I(:,:) = uint8( stk(:,:,k) );
    I = imadjust(I);
    imwrite(I,[Pre num2str(k, '%05d') '.tif']);
end