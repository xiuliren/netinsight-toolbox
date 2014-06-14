%% brain vesssel tracing
% by jpwu, 2013/02/27
clc
clear
close all

%% parameters
SrcImg = '../Data/047_030_croped_enhanced-1.tif';
% seedList = [100, 155, 212, 25];  
SrcImg = '../Data/047_030_croped_enhanced.tif';
% seedList = [69 132 72 25];
seedList = [316 95 8 5];
% SrcImg = '';
% stkSize = [3334 2687 334];
% seedList = [335 107 70 8];
% SrcImg = '../Data/barrel_field_data_microvessels.tif';
% seedList = [122 105 17 5];
% SrcImg = '../Data/evaluate_data/vol_sl.mat';
% SrcImg = '../Data/evaluate_data/vol_sl_di.mat';
% SrcImg = '../Data/evaluate_data/sl_di.tif';
% seedList = [20, 11, 20, 10];  % strait line
% SrcMat = '../Data/evaluate_data/vol_Y.mat';
% SrcImg = '../Data/evaluate_data/Y_060.tif';
% seedList = [ 100, 30, 20, 10 ]; % Y
% SrcImg = '../Data/evaluate_data/circular_tube.tif';
% seedList = [ 29, 10, 9, 5 ]; % circle
% SrcImg = '../Data/evaluate_data/strait_tube_varing_radius.tif';
% seedList = [75 90 75 5];
% section interval of detected seeds
SeedSecItv = 5;

%% load the image stack
global stk;
if strfind(SrcImg,'.tif')
    % image stack
    info = imfinfo(SrcImg);
    M = info(1).Height;
    N = info(1).Width;
    K = length(info);
    stk = zeros(M,N,K,'uint8');
    for k = 1 : K
        I(:,:) = imread(SrcImg, k);
        stk(:,:,k) = I;
    end
elseif strfind(SrcImg,'.mat')
    % the volume was stored in mat file
    load( SrcImg );
    stk = uint8(vol);
    [M N K] = size( stk );
elseif strfind(SrcImg,'.raw')
    fid=fopen(SrcImg, 'r');
    stk=fread(fid,'uint8');
%     stk = reshape(stk, );
    fclose(fid);
else
    disp('wrong source image file name!')
end

%% vessel tracing
disp('--------- start tracing ----------')
tic
network = nio_tracing(stk);
toc

%% save the results
if isempty(network.sections)
    disp('network is still empty!! have a check!')
else
    addpath('../')
    nio_write_net_hoc(network, 'tracing_result.hoc');
end

%% display the result
hold off % eliminate previouse plot
mip(:,:) = max(stk, [], 3);
imshow(mip);
hold on; 
% plot the circles
for si = 1 : network.sn
    sec = network.sections{si};
    for ni = 1 : size(sec, 1)
        plot(sec(ni,2), sec(ni,1), '.', ...
            'markersize', sec(ni, 4)*3, 'color', [rand, rand, rand] ); 
    end
end

% plot the center lines
for si = 1 : network.sn
    sec = network.sections{si};
    for ni = 1 : size(sec, 1)-1
        plot( [sec(ni,2) sec(ni+1,2)], ...
            [sec(ni,1) sec(ni+1,1)],'-r','LineWidth',2 );
    end
end

