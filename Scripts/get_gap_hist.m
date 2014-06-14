%% estimate the gap distribution histogram
% by jpwu, 2013/05/04

clc
clear
close all

%% parameters
SWCDir = '../Data/PMBSF_Blockface_1_T2_eliminated_Add_eliminated_Add_edit.swc';

%% load data
% tree = nio_load_tree( SWCDir );
% net = nio_tree2network( tree );
load( '../Data/barrel.mat' );

%% get the terminal points
Nsec = length( net.sections );

% get the section end points
eps = zeros( 2*Nsec, 3 );
for sn = 1 : Nsec
    sec = net.sections{sn};
    eps( sn, : ) = sec(1,1:3);
    eps( sn+Nsec, : ) = sec(end, 1:3 );
end

% 