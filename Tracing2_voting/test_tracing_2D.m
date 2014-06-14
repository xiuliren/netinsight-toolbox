% test blood pressure tracing algorithm
% by jpwu, 2013/08/09
clc
clear
close all

%% parameters
SrcImg = '../Data/tracing_test_2D.tif';

% the template radius
R = 10;
% the standard deviation of gaussian function
sig = 5;    % sigma

% initial cerebral perfusion pressure, CPP
% CPP is normally between 70 and 90 mmHg in an adult human,
% and cannot go below 70 mmHg.some regard 50-150
icpp = 80;

%% read image
disp('------- reading image -------')
I = double(imread( SrcImg ));
I = imresize( I, 0.2 );
[M N] = size(I);
% imshow(I)
Im = mean( I(:) );

%% create gaussian template 
disp('------ creating gaussian template --------')
% build the gaussian field
gf1 = zeros( 2*R+1, 2*R+1, 'double' );
gf = zeros( 2*R+1, 2*R+1, 2, 'double');
for m = 1 : R+1
    for n = 1 : R+1
        % distance to center
        dc = sqrt( (m-R-1)^2 + (n-R-1)^2 );
        gf1( m, n ) = gaussmf( dc, [ sig, 0 ] );
        
        gf(m,n,1) = gf1(m,n)/dc*(R+1-m);
        gf(m,n,2) = gf1(m,n)/dc*(R+1-n);
    end
end
gf(R+1, R+1, 1) = 0;    gf(R+1, R+1, 2) = 0;
% copy the left up corner to the other three regions
gf1( R+1:2*R+1, 1:R+1 ) = flipud( gf1( 1:R+1, 1:R+1 ) );
gf1( 1:R+1, R+1:2*R+1 ) = fliplr( gf1( 1:R+1, 1:R+1 ) );
gf1( R+1:2*R+1, R+1:2*R+1 ) = fliplr( flipud( gf1( 1:R+1, 1:R+1 ) ) );

gf( R+1:2*R+1, 1:R+1, 1 ) = -flipud( gf( 1:R+1, 1:R+1, 1 ) );
gf( R+1:2*R+1, 1:R+1, 2 ) = flipud( gf( 1:R+1, 1:R+1, 2 ) );
gf( 1:R+1, R+1:2*R+1, 1 ) = fliplr( gf( 1:R+1, 1:R+1, 1 ) );
gf( 1:R+1, R+1:2*R+1, 2 ) = -fliplr( gf( 1:R+1, 1:R+1, 2 ) );
gf( R+1:2*R+1, R+1:2*R+1, 1 ) = -fliplr( flipud( gf( 1:R+1, 1:R+1, 1 ) ) );
gf( R+1:2*R+1, R+1:2*R+1, 2 ) = -fliplr( flipud( gf( 1:R+1, 1:R+1, 2 ) ) );

% visualize the gaussian field
% imagesc(gf1);
% % hold on;
% for m = 1 : 2*R+1
%     for n = 1 : 2*R+1
%         % coordinate in data space
%         md = [ m,(m+gf(m,n,1)) ];
%         nd = [ n, (n+gf(m,n,2)) ];
%         % coordinate in figure space
%         [mf nf] = ds2nfu( nd, 2*R+1-md );
%         annotation('arrow', mf, nf );
%     end
% end

%% create the resistance field
% the black have resistence field,and the white have attraction field
disp('----------- creating the resistance field --------')

% build the initial resistence field
rf = zeros( M, N, 2 );

for m = R+1 : M-R
    for n = R+1 : N-R
        rf(m-R:m+R, n-R:n+R,:) = (I(m,n)-Im)/255*gf(:,:,:);
    end
end

%% visualize the resistance field
% imagesc(I);
% for m = 1 : 10 : M
%     for n = 1 : 10 : N
%         m
%         n
%         % coordinate in data space
%         md = [ (m-0.5),(m-0.5+rf(m,n,1)) ];
%         nd = [ (n-0.5), (n-0.5+rf(m,n,2)) ];
%         % transform to figure coordinate
%         [mf nf] = ds2nfu( nd, M-md );
%         annotation('arrow', mf, nf);
%     end
% end

%% create the pressure field
% get the seed points
% [ms ns] = ginput(1);
ms = 90.8;  ns = 124.8;

% constant cerebral perfusion pressure
cpp = icpp;

% initial propagation

% estimate the diameter

