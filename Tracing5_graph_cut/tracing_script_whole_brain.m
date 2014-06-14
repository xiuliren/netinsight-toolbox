%% iterative tracing of the whole mouse brain
% by jpwu, 2013/02/27
clc
clear
close all

%% parameters
% directories
SrcImg = '../../jpwu/Indian_ink_vessel/02_remove_outline_010102/MBA-IP11186_';
DstHoc = '../Data/MBA-IP11186_';

% the range of z
zs = 4901;  ze = 4950;
INT = 50;

%% load the image stack
info = imfinfo([SrcImg num2str(zs,'%05d') '.tif' ]);
Ms = info(1).Height;
Ns = info(1).Width;
Ks = INT;
global stk
stk = zeros(Ms, Ns, 2*Ks+10, 'uint8');

for z = zs : INT : ze
    % initialization        
    z1 = z;% - 5;
    z2 = z + INT-1; % + 5;
    
    I1 = imread([ SrcImg num2str(z1,'%05d') '.tif' ]);
    stk(:,:,1) = I1;    stk(:,:,2) = I1;

    % interpolation
    for k = 2 : Ks
        k
        I = imread([ SrcImg num2str(z1+k-1,'%05d') '.tif' ]);
        stk(:,:, 2*k) = I;
        % previous section
        Ip(:,:) = stk(:,:, 2*k-2);
        stk(:,:, 2*k-1) = I/2 + Ip/2;
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
        save( [DstHoc 'network_' num2str(z1,'%05d') '-' num2str(z2,'%05d') '.mat'], 'network' );
        nio_write_net_hoc(network, [DstHoc num2str(z1,'%05d') '-' num2str(z2,'%05d') '.hoc']);
    end
end