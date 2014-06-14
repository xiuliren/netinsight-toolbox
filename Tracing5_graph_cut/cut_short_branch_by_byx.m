% by byx, 2014/5/30
clc
clear
close all

%% parameters
% directories
Srcmat = '../Data/MBA-IP11186_network_';
DstHoc = '../Data/hoc_with_10_cut/MBA-IP11186_';
secLenThreshold = 10;
ze = 6051; zs = 6200;
for z = ze : 50 : zs
    load([Srcmat num2str(z,'%05d') '-' num2str(z+49,'%05d')]);
    network = nio_short_prunning( network, secLenThreshold );
    nio_write_net_hoc(network, [DstHoc num2str(z,'%05d') '-' num2str(z+49,'%05d') '.hoc']);
end
