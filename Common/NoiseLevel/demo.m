img = double(imread('lena.png'));

level = [5,10,20,40];
for i=1:size(level,2);
 noise = img + randn(size(img)) * level(i);
 tic;
 nlevel = NoiseLevel(noise);
 t=toc;
 fprintf('True: %5.2f  R:%5.2f G:%5.2f B:%5.2f\n', level(i), nlevel(1), nlevel(2), nlevel(3) );
 fprintf('Calculation time: %5.2f [sec]\n\n', t );
end
