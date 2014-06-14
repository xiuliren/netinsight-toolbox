%% brain vesssel tracing
% by jpwu, 2013/02/27
clc
clear
close all

%% parameters
% SrcImg = '../Data/047_030_croped_enhanced-1.tif';
% SrcImg = '../Data/zhixian.tif';
SrcImg = '../Data/y.tif';

DstHoc = '../Data/tracing_result.hoc';
itv = 20;

%% load the image stack
info = imfinfo(SrcImg);
M = info(1).Height;
N = info(1).Width;
K = length(info);

% global stk;
stk = zeros(M,N,K,'uint8');
for k = 1 : K
    stk(:,:, k) = imread(SrcImg,k);
end
% load('stk.mat')

% histogram of the whole stack
hw =  histc( reshape(stk, [], 1), 0:255 )';
% background intensity
global BI;
[tmp BI] = max(hw);
BI = BI - 1;    % 1-256 to 0-255
%% seeds detection, x, y, z, radius
% seedList = [175 228 208 26];
% seedList = [ 50, 50, 50, 3 ];
seedList = get_seed_list( stk, itv );

%% vessel tracing
disp('--------- start tracing ----------')
tic
profile on
network = nio_tracing(stk, seedList);
profile viewer
toc
%% save the results
if isempty(network.sections)
    disp('network is still empty!! have a check!')
else
    nio_save_hoc(network, DstHoc);
end