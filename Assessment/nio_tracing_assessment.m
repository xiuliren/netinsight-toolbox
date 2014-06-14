% tracing_assessment
% Quantitative assessment of vessel tracing accuracy.
% 
% [recall precision] =  tracing_assessment( SWCStk, ManSegStk,
% seg_idx, SrcImgStk, options )
% ----------------------------------------------------------------
%  
% Input
% -----
% - SWCStk:: the tracing result, turn the vectorized swc file to binary
%           image stack by swc2stk
% - ManSegStk:: manual segmentation of blood vessel sections
% - seg_idx:: the manual segmentation image index ( Z coordinate )
% - SrcImgStk:: the origional source images, for check tracing accuracy
%           subjectly
% - options:: if there is  '-s' , than show the origional image sections
%           and boundaries of traced vessel segments and manual labeled regions.
%           {DEFAULT: empty}
% 
% Output
% ------
% - recall:: the detection rate of blood vessels.how many vessel segments are
%           traced in 100 segments.
% - precision:: the precision of vessel tracing . how many real vessel
%           segments in 100 traced vessel segments.
% 
% Example
% -------
% [recall precision] = tracing_assessment( SWCStk, ManSegStk, 126:30:786,
% SrcImgStk, '-s' );
% 
% See also swc2stk

function [recall precision] = tracing_assessment( SWCStk, ManSegStk, seg_idx, SrcImgStk, options )

%% parameters
% the stack size
Ms = size(SWCStk, 1);   Ns = size(SWCStk, 2);   Ks = size(SWCStk, 3);

% segment image index
Ni = length( seg_idx );
% Mmd = zeros(1,Ni, 'double');

%tracing segmentation result
BW_t = zeros(Ms, Ns, 'uint8');
%manual segmentation result
BW_m = zeros(Ms, Ns, 'uint8');
% area ratio
ar = zeros(1,Ni);
% matched area ratio
mar = zeros(1, Ni);
% matched region ration 
mrr_m = zeros(1, Ni);
% diameter excursion ratio
der = zeros(1, Ni);
% area ratio
ar_micro = zeros(1, Ni);
% Precision
P = zeros(1, Ni);
% number of blood vessel 
num_ar = zeros(1, Ni);

%% compare the tracing segmentation results
disp('------- compare the segmentation results ...')
for ni = 1 : length( seg_idx )
    % automatic and manual tracing sections
    idx = seg_idx( ni );
    % tracing segmentation result
    BW_t(:,:) = SWCStk(:,:,idx);
    BW_t = im2bw( BW_t, 0.5/255 );

    % manual segmentation result
    BW_m = ManSegStk(:,:, idx );
    BW_m = im2bw( BW_m, 0.5/255);    
    % get the shared region
    S_m = regionprops( BW_m, 'All' );
    S_t = regionprops( BW_t, 'All' );
    
    % label matrix of manual
    LM_m = zeros( Ms, Ns );
    for m = 1 : length( S_m )
        for k = 1 : length( S_m(m).PixelIdxList )
            LM_m( S_m(m).PixelIdxList(k) ) = m ;
        end
    end
    
    % label matrix of tracing result
    LM_t = zeros( Ms, Ns );
    for t = 1 : length( S_t )
        for k = 1 : length( S_t(t).PixelIdxList )
            LM_t( S_t(t).PixelIdxList(k) ) = t ;
        end
    end
    
    % overlap matrix
    OverLap = zeros( length(S_m), length(S_t) );
%     OverLap = sparse( length(S_m), length(S_t) );
    for m = 1 : length( S_m )
        for k = 1 : length( S_m(m).PixelIdxList )
            regionIdx_t = LM_t( S_m(m).PixelIdxList(k) );
            if  regionIdx_t ~= 0
                OverLap( m, regionIdx_t ) = 1 ;
            end
        end
    end
    
    % compare the matched region area
    % total area of matched region
    MatchNum_m = 0;
    sum_m = 0;
    sum_t = 0;
    for m = 1 : length( S_m )
        regionIdxList_t = find( OverLap(m, :) == 1 );
        
        % add area of manual
        if ~isempty( regionIdxList_t )
            sum_m = sum_m + S_m(m).Area;
            MatchNum_m = MatchNum_m + 1;
        end
        % add area of tracing
        for k = 1 : length( regionIdxList_t )
            sum_t = sum_t + S_t( regionIdxList_t(k) ).Area;
        end
    end
    % matched area ratio
    mar(ni) = sum_t / sum_m; 
    % matched region ratio
    mrr_m(ni) = MatchNum_m / length( S_m );
    % area ratio
    ar(ni) = sum( BW_t(:) ) / sum( BW_m(:) );
    % Precision
    P(ni) =  length( find( max(OverLap) == 1 ) ) / length( S_t );
    % number of blood vessel 
    num_ar(ni) = length(S_m);
   %%  get boundary and plot
   if ~isempty(options)&&strfind(options, '-s')
       %     the original image
       I(:,:) = SrcImgStk(:, :, idx );
       B_t = bwboundaries(BW_t,'noholes');
       B_m = bwboundaries(BW_m,'noholes');
       figure;
       imshow(I)
       hold on
       for k = 1:length(B_m)
           boundary = B_m{k};
           plot(boundary(:,2), boundary(:,1), 'b', 'LineWidth', 2)
       end
       for k = 1:length(B_t)
           boundary = B_t{k};
           plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
       end
       hold off
   end
end

idx = mrr_m>0;    mrr_m = mrr_m(idx);
idx = P>0;        P = P(idx);
idx = ar>0;       ar = ar(idx);
idx = find(mar>0);       mar = mar(idx);
% idx = find(Mmd>0);       Mmd = Mmd(idx);
%% recall and precision
recall = mean( mrr_m );
precision = mean( P );

disp([ 'Recall: ' num2str( mean(mrr_m) ) ' ¡À  '  num2str( std(mrr_m) ) ]);
disp([ 'Precision: ' num2str( mean(P) ) ' ¡À  '  num2str( std(P) ) ]);
disp([ 'area ratio: ' num2str( mean(ar) ) ' ¡À  '  num2str( std(ar) ) ]);
disp([ 'matched area ratio: ' num2str( mean(mar) ) ' ¡À  '  num2str( std(mar) ) ]);
disp([ 'Diameter ratio: ' num2str( mean(sqrt(mar)) ) ' ¡À  '  num2str( std(sqrt(mar)) ) ]);
% disp([ 'Mmd: ' num2str( mean(Mmd) ) ' ¡À '  num2str( std(Mmd) ) ]);

disp('------------- end -----------')