n_circles = 20;
color_length = 100;
image_size = 128;
max_radius = 20;

%% plot simultaneously
[X1 Y1] = getmidpointcircle(0, 0, 10);
for np = 1 : 1 : length(X1)
    plot( X1(np),Y1(np), '.r', 'markersize',10);
    text(X1(np),Y1(np),['\leftarrow' num2str(np)])
    hold on
    xlim([-10 10]);ylim([-10 10]);
%     pause(1)
end

%%
I = zeros(image_size, image_size, 3, 'uint8');
colors = hsv(color_length);

for i = 1 : n_circles
    
    x0 = round( image_size * rand);
    y0 = round( image_size * rand);
    radius = round( max_radius * rand );
    
    [x y] = getmidpointcircle(x0, y0, radius);
    
    index = 1 ;
    for j = 1 : numel(x)
        xp = x(j);
        yp = y(j);
        
        if ( xp < 1 || yp < 1 || xp > image_size || yp > image_size )
            continue
        end
        I(xp, yp, :) = round( 255 * colors(index, :) );
        index = index + 1;
        if index > color_length
            index = 1;
        end
    end
    
end

imshow(I, []);