%% multiscale linear correlation of somas and microvessels density
% load the quantitative analysis results of 'mat' files
% get things together and plot
%% parameters
disp('------ setting parameters ... ')

% block size, micron
Scale = 200;

% volume shrinkage
skg = [];
% lenth density revise coefficients
ld2 = []; 

% intersection area ratio
ar_BC = 0.762;
ar_MC = 0.855;
ar_VC = 0.691;

% quantitative analysis mat files
Mat_MC_V = '';
Mat_BC_V = '';
Mat_VC_V = '';

Mat_MC_S = '';
Mat_BC_S = '';
Mat_VC_S = '';

% x y limit ratio
Limit_R = 10;

% number of normalized vector elements
N_nml = 101;

%% load microvessels and somas density
disp('----- loading density ......')
% motor cortex
% vessel density
load( Mat_MC_V );
num_part_MC = num_part;
fv_a_MC     = fv_a;
fv_a_c_MC   = fv_a_c;
ld_a_MC     = ld_a;
ld_a_c_MC   = ld_a_c;
% soma density
load( Mat_MC_S )
ds_MC = ds_a;
z_en_MC = z_en;
z_st_MC = z_st;

% visual cortex, vessels
load( Mat_VC_V );
num_part_VC = num_part;
fv_a_VC     = fv_a;
fv_a_c_VC   = fv_a_c;
ld_a_VC     = ld_a;
ld_a_c_VC   = ld_a_c;
% visual cortex, soma density
load( Mat_VC_S )
ds_VC = ds_a;
z_en_VC = z_en;
z_st_VC = z_st;

% barrel cortex, vessels
load( Mat_BC_V );
num_part_BC = num_part;
fv_a_BC     = fv_a;
fv_a_c_BC   = fv_a_c;
ld_a_BC     = ld_a;
ld_a_c_BC   = ld_a_c;
% barrel cortex, soma density
load( Mat_BC_S )
ds_BC = ds_a;
z_en_BC = z_en;
z_st_BC = z_st;

%% average density along different distance from pia
disp('------ average density ...... ')
% motor cortex
avg_fv_MC = zeros(1, num_part_MC ); % the average fractional vascular volume along different distance from pia
avg_fvc_MC = zeros(1, num_part_MC ); % microvascular
avg_ld_MC = zeros(1, num_part_MC ); % the average length density for vessels
avg_ldc_MC = zeros(1, num_part_MC ); % microvascular
avg_ds_MC = zeros(1, num_part_MC ); % the average cell density
avg_ads_MC = zeros(1, num_part_MC ); % the average cell areal density

for ward = 1 : num_part_MC
    avg_fv_MC(ward)     = mean(fv_a_MC{ward});
    avg_fvc_MC(ward)    = mean(fv_a_c_MC{ward});
    avg_ld_MC(ward)     = mean(ld_a_MC{ward});
    avg_ldc_MC(ward)    = mean(ld_a_c_MC{ward});
    avg_ds_MC(ward)     = mean(ds_MC{ward});
    avg_ads_MC(ward)    = avg_ds_MC(ward) * (z_en_MC(ward) - z_st_MC(ward)) / 1000;
end

% revise parameters
avg_fv_MC   = avg_fv_MC / ar_MC;
avg_fvc_MC  = avg_fvc_MC / ar_MC;
avg_ld_MC   = avg_ld_MC * ld2;
avg_ldc_MC  = avg_ldc_MC * ld2;
avg_ds_MC   = avg_ds_MC * skg;
avg_ads_MC  = avg_ads_MC * ld2;

% visual cortex
avg_fv_VC = zeros(1, num_part_VC ); % the average fractional vascular volume along different distance from pia
avg_fvc_VC = zeros(1, num_part_VC ); % microvascular
avg_ld_VC = zeros(1, num_part_VC ); % the average length density for vessels
avg_ldc_VC = zeros(1, num_part_VC ); % microvascular
avg_ds_VC = zeros(1, num_part_VC ); % the average cell density
avg_ads_VC = zeros(1, num_part_VC ); % the average cell areal density

for ward = 1 : num_part_VC
    avg_fv_VC(ward) = mean(fv_a_VC{ward});
    avg_fvc_VC(ward) = mean(fv_a_c_VC{ward});
    avg_ld_VC(ward) = mean(ld_a_VC{ward});
    avg_ldc_VC(ward) = mean(ld_a_c_VC{ward});
    avg_ds_VC(ward) = mean(ds_VC{ward});
    avg_ads_VC(ward) = avg_ds_VC(ward) * (z_en_VC(ward) - z_st_VC(ward)) / 1000;
end

% revise parameters
avg_fv_VC   = avg_fv_VC / ar_VC;
avg_fvc_VC  = avg_fvc_VC / ar_VC;
avg_ld_VC   = avg_ld_VC * ld2;
avg_ldc_VC  = avg_ldc_VC * ld2;
avg_ds_VC   = avg_ds_VC * skg;
avg_ads_VC  = avg_ads_VC * ld2;

% barrel cortex
avg_fv_BC = zeros(1, num_part_BC ); % the average fractional vascular volume along different distance from pia
avg_fvc_BC = zeros(1, num_part_BC ); % microvascular
avg_ld_BC = zeros(1, num_part_BC ); % the average length density for vessels
avg_ldc_BC = zeros(1, num_part_BC ); % microvascular
avg_ds_BC = zeros(1, num_part_BC ); % the average cell density
avg_ads_BC = zeros(1, num_part_BC ); % the average cell areal density

for ward = 1 : num_part_BC
    avg_fv_BC(ward)     = mean(fv_a_BC{ward});
    avg_fvc_BC(ward)    = mean(fv_a_c_BC{ward});
    avg_ld_BC(ward)     = mean(ld_a_BC{ward});
    avg_ldc_BC(ward)    = mean(ld_a_c_BC{ward});
    avg_ds_BC(ward)     = mean(ds_BC{ward});
    avg_ads_BC(ward)    = avg_ds_BC(ward) * (z_en_BC(ward) - z_st_BC(ward)) / 1000;
end

% revise parameters
avg_fv_BC   = avg_fv_BC / ar_BC;
avg_fvc_BC  = avg_fvc_BC / ar_BC;
avg_ld_BC   = avg_ld_BC * ld2;
avg_ldc_BC  = avg_ldc_BC * ld2;
avg_ds_BC   = avg_ds_BC * skg;
avg_ads_BC  = avg_ads_BC * ld2;

%% combine multiple cortex
avg_fv  = [ avg_fv_MC  avg_fv_VC    avg_fv_BC ];
avg_fvc = [ avg_fvc_MC avg_fvc_VC   avg_fvc_BC ];
avg_ld  = [ avg_ld_MC  avg_ld_VC    avg_ld_BC ];
avg_ldc = [ avg_ldc_MC avg_ldc_VC   avg_ldc_BC ];
avg_ds  = [ avg_ds_MC  avg_ds_VC    avg_ds_BC ];
avg_ads = [ avg_ads_MC avg_ads_VC   avg_ads_BC ];
avg_dtm = [ m_dtm_a m_dtm_m m_dtm_v ];
num_part = num_part_MC + num_part_VC + num_part_BC;
%% normalize the density parameters of BMV cortex
% motor cortex
nml_ds_MC       = zeros( N_nml, num_part_MC);
nml_fv_a_MC     = zeros( N_nml, num_part_MC);
nml_fv_a_c_MC   = zeros( N_nml, num_part_MC);
nml_ld_a_MC     = zeros( N_nml, num_part_MC);
nml_ld_a_c_MC   = zeros( N_nml, num_part_MC);

% normalize the cortex depth
for ward = 1 : num_part_MC
    nml_ds_MC(:,ward)       = imresize( ds_MC{ ward }', [N_nml,1] );
    nml_fv_a_MC(:,ward)     = imresize( fv_a_MC{ ward }', [N_nml,1] );
    nml_fv_a_c_MC(:,ward)   = imresize( fv_a_c_MC{ ward }', [N_nml,1] );
    nml_ld_a_MC(:,ward)     = imresize( ld_a_MC{ ward }', [N_nml,1] );
    nml_ld_a_c_MC(:,ward)   = imresize( ld_a_c_MC{ ward }', [N_nml,1] );
end

% revise
nml_ds_MC     = nml_ds_MC * skg;
nml_fv_a_MC   = nml_fv_a_MC / ar_MC;
nml_fv_a_c_MC = nml_fv_a_c_MC / ar_MC;
nml_ld_a_MC   = nml_ld_a_MC * ld2;
nml_ld_a_c_MC = nml_ld_a_c_MC * ld2;

% mean
m_nml_ds_MC     = mean( nml_ds_MC' );
m_nml_fv_a_MC   = mean( nml_fv_a_MC' );
m_nml_fv_a_c_MC = mean( nml_fv_a_c_MC' );
m_nml_ld_a_MC   = mean( nml_ld_a_MC' );
m_nml_ld_a_c_MC = mean( nml_ld_a_c_MC' );

% standard deviation
std_nml_ds_MC     = std( nml_ds_MC' );
std_nml_fv_a_MC   = std( nml_fv_a_MC' );
std_nml_fv_a_c_MC = std( nml_fv_a_c_MC' );
std_nml_ld_a_MC   = std( nml_ld_a_MC' );
std_nml_ld_a_c_MC = std( nml_ld_a_c_MC' );

% visual cortex
nml_ds_VC       = zeros( N_nml, num_part_VC);
nml_fv_a_VC     = zeros( N_nml, num_part_VC);
nml_fv_a_c_VC   = zeros( N_nml, num_part_VC);
nml_ld_a_VC     = zeros( N_nml, num_part_VC);
nml_ld_a_c_VC   = zeros( N_nml, num_part_VC);

% normalize the cortex depth
for ward = 1 : num_part_VC
    nml_ds_VC(:,ward)       = imresize( ds_VC{ ward }', [N_nml,1] );
    nml_fv_a_VC(:,ward)     = imresize( fv_a_VC{ ward }', [N_nml,1] );
    nml_fv_a_c_VC(:,ward)   = imresize( fv_a_c_VC{ ward }', [N_nml,1] );
    nml_ld_a_VC(:,ward)     = imresize( ld_a_VC{ ward }', [N_nml,1] );
    nml_ld_a_c_VC(:,ward)   = imresize( ld_a_c_VC{ ward }', [N_nml,1] );
end

% revise
nml_ds_VC     = nml_ds_VC * skg;
nml_fv_a_VC   = nml_fv_a_VC / ar_VC;
nml_fv_a_c_VC = nml_fv_a_c_VC / ar_VC;
nml_ld_a_VC   = nml_ld_a_VC * ld2;
nml_ld_a_c_VC = nml_ld_a_c_VC * ld2;

% mean
m_nml_ds_VC     = mean( nml_ds_VC' );
m_nml_fv_a_VC   = mean( nml_fv_a_VC' );
m_nml_fv_a_c_VC = mean( nml_fv_a_c_VC' );
m_nml_ld_a_VC   = mean( nml_ld_a_VC' );
m_nml_ld_a_c_VC = mean( nml_ld_a_c_VC' );

% standard deviation
std_nml_ds_VC     = std( nml_ds_VC' );
std_nml_fv_a_VC   = std( nml_fv_a_VC' );
std_nml_fv_a_c_VC = std( nml_fv_a_c_VC' );
std_nml_ld_a_VC   = std( nml_ld_a_VC' );
std_nml_ld_a_c_VC = std( nml_ld_a_c_VC' );

% barrel cortex
nml_ds_BC       = zeros( N_nml, num_part_BC);
nml_fv_a_BC     = zeros( N_nml, num_part_BC);
nml_fv_a_c_BC   = zeros( N_nml, num_part_BC);
nml_ld_a_BC     = zeros( N_nml, num_part_BC);
nml_ld_a_c_BC   = zeros( N_nml, num_part_BC);

% normalize the cortex depth
for ward = 1 : num_part_BC
    nml_ds_BC(:,ward)       = imresize( ds_BC{ ward }', [N_nml,1] );
    nml_fv_a_BC(:,ward)     = imresize( fv_a_BC{ ward }', [N_nml,1] );
    nml_fv_a_c_BC(:,ward)   = imresize( fv_a_c_BC{ ward }', [N_nml,1] );
    nml_ld_a_BC(:,ward)     = imresize( ld_a_BC{ ward }', [N_nml,1] );
    nml_ld_a_c_BC(:,ward)   = imresize( ld_a_c_BC{ ward }', [N_nml,1] );
end

% revise
nml_ds_BC     = nml_ds_BC * skg;
nml_fv_a_BC   = nml_fv_a_BC / ar_BC;
nml_fv_a_c_BC = nml_fv_a_c_BC / ar_BC;
nml_ld_a_BC   = nml_ld_a_BC * ld2;
nml_ld_a_c_BC = nml_ld_a_c_BC * ld2;

% mean
m_nml_ds_BC     = mean( nml_ds_BC' );
m_nml_fv_a_BC   = mean( nml_fv_a_BC' );
m_nml_fv_a_c_BC = mean( nml_fv_a_c_BC' );
m_nml_ld_a_BC   = mean( nml_ld_a_BC' );
m_nml_ld_a_c_BC = mean( nml_ld_a_c_BC' );

% standard deviation
std_nml_ds_BC     = std( nml_ds_BC' );
std_nml_fv_a_BC   = std( nml_fv_a_BC' );
std_nml_fv_a_c_BC = std( nml_fv_a_c_BC' );
std_nml_ld_a_BC   = std( nml_ld_a_BC' );
std_nml_ld_a_c_BC = std( nml_ld_a_c_BC' );

%% basic quantitative analysis parameters
disp(['basic quantitative analysis parameters with block size: ' num2str( Scale ) ]);
disp('barrel cortex;    motor cortex;   visual cortex;  total');
% fractional volume density of all vessels
disp([ 'fractional volume density of all vessels: ' ])
disp([ num2str( mean(avg_fv_BC) ) '��' num2str( std(avg_fv_BC) ) '   '...
    num2str( mean(avg_fv_MC) ) '��' num2str( std(avg_fv_MC) ) '   '...
    num2str( mean(avg_fv_VC) ) '��' num2str( std(avg_fv_VC) ) '   ' ...
    num2str( mean(avg_fv) ) '��' num2str( std(avg_fv) ) ])

% fractional volume density of all capillaries 
disp([ 'fractional volume density of all capillaries: ' ])
disp([ num2str( mean(avg_fvc_BC) )  '��'     num2str( std(avg_fvc_BC) ) '   '...
    num2str( mean(avg_fvc_MC) )     '��'     num2str( std(avg_fvc_MC) ) '   '...
    num2str( mean(avg_fvc_VC) )     '��'     num2str( std(avg_fvc_VC) ) '   ' ...
    num2str( mean(avg_fvc) )        '��'     num2str( std(avg_fvc) ) ])

% normalized length density of all vessels 
disp([ 'normalized vascular length of all vessels: ' ])
disp([  num2str( mean(avg_ld_BC) ) '��'  num2str( std(avg_ld_BC) ) '   '...
        num2str( mean(avg_ld_MC) ) '��'  num2str( std(avg_ld_MC) ) '   '...
        num2str( mean(avg_ld_VC) ) '��'  num2str( std(avg_ld_VC) ) '   ' ...
        num2str( mean(avg_ld) )    '��'  num2str( std(avg_ld) ) ])
    
% normalized length density of all capillaries 
disp([ 'normalized vascular length of all capillaries: ' ])
disp([  num2str( mean(avg_ldc_BC) ) '��'  num2str( std(avg_ldc_BC) ) '   '...
        num2str( mean(avg_ldc_MC) ) '��'  num2str( std(avg_ldc_MC) ) '   '...
        num2str( mean(avg_ldc_VC) ) '��'  num2str( std(avg_ldc_VC) ) '   ' ...
        num2str( mean(avg_ldc) )    '��'  num2str( std(avg_ldc) ) ])

% cell number density 
disp([ 'cell number density: ' ])
disp([  num2str( mean(avg_ds_BC) ) '��'  num2str( std(avg_ds_BC) ) '   '...
        num2str( mean(avg_ds_MC) ) '��'  num2str( std(avg_ds_MC) ) '   '...
        num2str( mean(avg_ds_VC) ) '��'  num2str( std(avg_ds_VC) ) '   ' ...
        num2str( mean(avg_ds) )    '��'  num2str( std(avg_ds) ) ])    
    
%% linear regression of all the cortex region
disp('----- linear regression ...... ')
%% fractional vascular volume VS cell number density
[b,bint,r,rint,stats] = regress(avg_fv',[ones(1, num_part)' avg_ds']);
figure
plot(avg_ds_MC, avg_fv_MC, '*');hold on;
plot(avg_ds_VC, avg_fv_VC, 'o');hold on;
plot(avg_ds_BC, avg_fv_BC, 's');hold on;
lsline
X_range = [ mean(avg_ds)-Limit_R * std(avg_ds),   mean(avg_ds)+Limit_R * std(avg_ds) ];
Y_range = [ mean(avg_fv)-Limit_R * std(avg_fv),   mean(avg_fv)+Limit_R * std(avg_fv) ];
xx = X_range(1) : 0.01 : X_range(2); 
yy = b(1) + b(2)*xx;
plot(xx, yy, 'LineStyle', '--');
xlim( X_range )
ylim([ Y_range ])
hold off
xlabel('Cell number density (10^5mm^-^3)');
ylabel('Fractional vascular volume (v/v)');
legend('Motor Cortex', 'Visual Cortex', 'Barrel Cortex','Location','NorthWest');
disp(['Fractional vascular volume VS cell number density ' char(13) 'r = ' num2str(sqrt(stats(1))) '   p = ' num2str(stats(3))])

%% fractional microvascular volume VS cell number density
[b,bint,r,rint,stats] = regress(avg_fvc',[ones(1, num_part)' avg_ds']);
figure
plot(avg_ds_MC, avg_fvc_MC, '*');hold on;
plot(avg_ds_VC, avg_fvc_VC, 'o');hold on;
plot(avg_ds_BC, avg_fvc_BC, 's');hold on;
lsline
X_range = [ mean(avg_ds)-Limit_R * std(avg_ds),   mean(avg_ds)+Limit_R * std(avg_ds) ];
Y_range = [ mean(avg_fvc)-Limit_R * std(avg_fvc),   mean(avg_fvc)+Limit_R * std(avg_fvc) ];
xx = X_range(1) : 0.01 : X_range(2); 
yy = b(1) + b(2)*xx;
plot(xx, yy, 'LineStyle', '--');
xlim( X_range )
ylim([ Y_range ])
hold off
xlabel('Cell number density (10^5mm^-^3)');
ylabel('Fractional microvascular volume (v/v)');
legend('Motor Cortex', 'Visual Cortex', 'Barrel Cortex','Location','NorthWest');
disp(['Fractional microvascular volume VS cell number density ' char(13) 'r = ' num2str(sqrt(stats(1))) '   p = ' num2str(stats(3))])
%% fractional microvascular volume VS cell areal density
[b,bint,r,rint,stats] = regress(avg_fvc',[ones(1, num_part)' avg_ads']);
figure
plot(avg_ads_MC, avg_fvc_MC, '*');hold on;
plot(avg_ads_VC, avg_fvc_VC, 'o');hold on;
plot(avg_ads_BC, avg_fvc_BC, 's');hold on;
lsline
X_range = [ mean(avg_ads)-Limit_R * std(avg_ads),   mean(avg_ads)+Limit_R * std(avg_ads) ];
Y_range = [ mean(avg_fvc)-Limit_R * std(avg_fvc),   mean(avg_fvc)+Limit_R * std(avg_fvc) ];
xx = X_range(1) : 0.01 : X_range(2); 
yy = b(1) + b(2)*xx;
plot(xx, yy, 'LineStyle', '--');
xlim( X_range )
ylim([ Y_range ])
hold off
xlabel('Cell areal density (10^5mm^-^2)');
ylabel('Fractional microvascular volume (v/v)');
legend('Motor Cortex', 'Visual Cortex', 'Barrel Cortex','Location','NorthWest');
disp(['Fractional microvascular volume VS cell areal density ' char(13) 'r = ' num2str(sqrt(stats(1))) '   p = ' num2str(stats(3))])
%% fractional vascular volume VS cell areal density
[b,bint,r,rint,stats] = regress(avg_fv',[ones(1, num_part)' avg_ads']);
figure
plot(avg_ads_MC, avg_fv_MC, '*');hold on;
plot(avg_ads_VC, avg_fv_VC, 'o');hold on;
plot(avg_ads_BC, avg_fv_BC, 's');hold on;
lsline
X_range = [ mean(avg_ads)-Limit_R * std(avg_ads),   mean(avg_ads)+Limit_R * std(avg_ads) ];
Y_range = [ mean(avg_fv)-Limit_R * std(avg_fv),   mean(avg_fv)+Limit_R * std(avg_fv) ];
xx = X_range(1) : 0.01 : X_range(2); 
yy = b(1) + b(2)*xx;
plot(xx, yy, 'LineStyle', '--');
xlim( X_range ); ylim( Y_range );
hold off
xlabel('Cell areal density (10^5mm^-^2)');
ylabel('Fractional vascular volume (v/v)');
legend('Motor Cortex', 'Visual Cortex', 'Barrel Cortex','Location','NorthWest');
disp(['Fractional vascular volume VS Cell areal density ' char(13) 'r = ' num2str(sqrt(stats(1))) '   p = ' num2str(stats(3))])

%% Normalized vascular length VS cell areal density
[b,bint,r,rint,stats] = regress(avg_ld',[ones(1, num_part)' avg_ads']);
figure
plot(avg_ads_MC, avg_ld_MC, '*');hold on;
plot(avg_ads_VC, avg_ld_VC, 'o');hold on;
plot(avg_ads_BC, avg_ld_BC, 's');hold on;
lsline
X_range = [ mean(avg_ads)-Limit_R * std(avg_ads),   mean(avg_ads)+Limit_R * std(avg_ads) ];
Y_range = [ mean(avg_ld)-Limit_R * std(avg_ld),   mean(avg_ld)+Limit_R * std(avg_ld) ];
xx = X_range(1) : 0.01 : X_range(2); 
yy = b(1) + b(2)*xx;
plot(xx, yy, 'LineStyle', '--');
xlim( X_range ); ylim( Y_range );
hold off
xlabel('Cell areal density (10^5mm^-^2)');
ylabel('Normalized vascular length (m/mm^3)');
legend('Motor Cortex', 'Visual Cortex', 'Barrel Cortex','Location','NorthWest');
disp(['Normalized vascular length VS Cell areal density ' char(13) 'r = ' num2str(sqrt(stats(1))) '   p = ' num2str(stats(3))])

%% Normalized microvascular length VS cell areal density
[b,bint,r,rint,stats] = regress(avg_ldc',[ones(1, num_part)' avg_ads']);
figure
plot(avg_ads_MC, avg_ldc_MC, '*');hold on;
plot(avg_ads_VC, avg_ldc_VC, 'o');hold on;
plot(avg_ads_BC, avg_ldc_BC, 's');hold on;
lsline
X_range = [ mean(avg_ads)-Limit_R * std(avg_ads),   mean(avg_ads)+Limit_R * std(avg_ads) ];
Y_range = [ mean(avg_ldc)-Limit_R * std(avg_ldc),   mean(avg_ldc)+Limit_R * std(avg_ldc) ];
xx = X_range(1) : 0.01 : X_range(2); 
yy = b(1) + b(2)*xx;
plot(xx, yy, 'LineStyle', '--');
xlim( X_range ); ylim( Y_range );
hold off
xlabel('Cell areal density (10^5mm^-^2)');
ylabel('Normalized microvascular length (m/mm^3)');
legend('Motor Cortex', 'Visual Cortex', 'Barrel Cortex','Location','NorthWest');
disp(['Normalized microvascular length VS Cell areal density ' char(13) 'r = ' num2str(sqrt(stats(1))) '   p = ' num2str(stats(3))])

%% Normalized microvascular length VS cell number density
[b,bint,r,rint,stats] = regress(avg_ldc',[ones(1, num_part)' avg_ds']);
figure
plot(avg_ds_MC, avg_ldc_MC, '*');hold on;
plot(avg_ds_VC, avg_ldc_VC, 'o');hold on;
plot(avg_ds_BC, avg_ldc_BC, 's');hold on;
lsline
X_range = [ mean(avg_ds)-Limit_R * std(avg_ds),   mean(avg_ds)+Limit_R * std(avg_ds) ];
Y_range = [ mean(avg_ldc)-Limit_R * std(avg_ldc),   mean(avg_ldc)+Limit_R * std(avg_ldc) ];
xx = X_range(1) : 0.01 : X_range(2); 
yy = b(1) + b(2)*xx;
plot(xx, yy, 'LineStyle', '--');
xlim( X_range ); ylim( Y_range );
hold off
xlabel('Cell number density (10^5mm^-^3)');
ylabel('Normalized microvascular length (m/mm^3)');
legend('Motor Cortex', 'Visual Cortex', 'Barrel Cortex','Location','NorthWest');
disp(['Normalized microvascular length VS Cell number density ' char(13) 'r = ' num2str(sqrt(stats(1))) '   p = ' num2str(stats(3))])

%% Normalized vascular length VS cell number density
[b,bint,r,rint,stats] = regress(avg_ld',[ones(1, num_part)' avg_ds']);
figure
plot(avg_ds_MC, avg_ld_MC, '*');hold on;
plot(avg_ds_VC, avg_ld_VC, 'o');hold on;
plot(avg_ds_BC, avg_ld_BC, 's');hold on;
lsline
X_range = [ mean(avg_ds)-Limit_R * std(avg_ds),   mean(avg_ds)+Limit_R * std(avg_ds) ];
Y_range = [ mean(avg_ld)-Limit_R * std(avg_ld),   mean(avg_ld)+Limit_R * std(avg_ld) ];
xx = X_range(1) : 0.01 : X_range(2); 
yy = b(1) + b(2)*xx;
plot(xx, yy, 'LineStyle', '--');
xlim( X_range ); ylim( Y_range );
hold off
xlabel('Cell number density (10^5mm^-^3)');
ylabel('Normalized vascular length (m/mm^3)');
legend('Motor Cortex', 'Visual Cortex', 'Barrel Cortex','Location','NorthWest');
disp(['Normalized vascular length VS Cell number density ' char(13) 'r = ' num2str(sqrt(stats(1))) '   p = ' num2str(stats(3))])

%% plot the normalized BMV cortex density
% cells number density
figure,  hold on
plot(1:N_nml, m_nml_ds_MC, 'r', 'LineWidth', 2)
plot(1:N_nml, m_nml_ds_BC, 'Color', [0 0.7 0], 'LineWidth', 2)
plot(1:N_nml, m_nml_ds_VC, 'b', 'LineWidth', 2)
patch([1:N_nml fliplr(1:N_nml)], [m_nml_ds_MC - std_nml_ds_MC fliplr(m_nml_ds_MC + std_nml_ds_MC)], 'r', 'EdgeColor', 'none');alpha(0.2)
patch([1:N_nml fliplr(1:N_nml)], [m_nml_ds_BC - std_nml_ds_BC fliplr(m_nml_ds_BC + std_nml_ds_BC)], [0 0.7 0], 'EdgeColor', 'none');alpha(0.2)
patch([1:N_nml fliplr(1:N_nml)], [m_nml_ds_VC - std_nml_ds_VC fliplr(m_nml_ds_VC + std_nml_ds_VC)], 'b', 'EdgeColor', 'none');alpha(0.2)
xlim([0 N_nml])
hold off
xlabel('Cortical depth (%)');
ylabel('Cell number density (10^5mm^-^3)');
legend('Motor Cortex', 'Barrel Cortex', 'Visual Cortex');

% Fractional vascular volume (v/v)
figure,  hold on
plot(1:N_nml, m_nml_fv_a_MC, 'r', 'LineWidth', 2)
plot(1:N_nml, m_nml_fv_a_BC, 'Color', [0 0.7 0], 'LineWidth', 2)
plot(1:N_nml, m_nml_fv_a_VC, 'b', 'LineWidth', 2)
patch([1:N_nml fliplr(1:N_nml)], [m_nml_fv_a_MC - std_nml_fv_a_MC fliplr(m_nml_fv_a_MC + std_nml_fv_a_MC)], 'r', 'EdgeColor', 'none');alpha(0.2)
patch([1:N_nml fliplr(1:N_nml)], [m_nml_fv_a_BC - std_nml_fv_a_BC fliplr(m_nml_fv_a_BC + std_nml_fv_a_BC)], [0 0.7 0], 'EdgeColor', 'none');alpha(0.2)
patch([1:N_nml fliplr(1:N_nml)], [m_nml_fv_a_VC - std_nml_fv_a_VC fliplr(m_nml_fv_a_VC + std_nml_fv_a_VC)], 'b', 'EdgeColor', 'none');alpha(0.2)
xlim([0 N_nml])
hold off
xlabel('Cortical depth (%)');
ylabel('Fractional vascular volume (v/v)');
legend('Motor Cortex', 'Barrel Cortex', 'Visual Cortex');

% Fractional microvascular volume (v/v)
figure,  hold on
plot(1:N_nml, m_nml_fv_a_c_MC, 'r', 'LineWidth', 2)
plot(1:N_nml, m_nml_fv_a_c_BC, 'Color', [0 0.7 0], 'LineWidth', 2)
plot(1:N_nml, m_nml_fv_a_c_VC, 'b', 'LineWidth', 2)
patch([1:N_nml fliplr(1:N_nml)], [m_nml_fv_a_c_MC - std_nml_fv_a_c_MC fliplr(m_nml_fv_a_c_MC + std_nml_fv_a_c_MC)], 'r', 'EdgeColor', 'none');alpha(0.2)
patch([1:N_nml fliplr(1:N_nml)], [m_nml_fv_a_c_BC - std_nml_fv_a_c_BC fliplr(m_nml_fv_a_c_BC + std_nml_fv_a_c_BC)], [0 0.7 0], 'EdgeColor', 'none');alpha(0.2)
patch([1:N_nml fliplr(1:N_nml)], [m_nml_fv_a_c_VC - std_nml_fv_a_c_VC fliplr(m_nml_fv_a_c_VC + std_nml_fv_a_c_VC)], 'b', 'EdgeColor', 'none');alpha(0.2)
xlim([0 N_nml])
hold off
xlabel('Cortical depth (%)');
ylabel('Fractional microvascular volume (v/v)');
legend('Motor Cortex', 'Barrel Cortex', 'Visual Cortex');

% Normalized vascular length (m/mm^3)
figure,  hold on
plot(1:N_nml, m_nml_ld_a_MC, 'r', 'LineWidth', 2)
plot(1:N_nml, m_nml_ld_a_BC, 'Color', [0 0.7 0], 'LineWidth', 2)
plot(1:N_nml, m_nml_ld_a_VC, 'b', 'LineWidth', 2)
patch([1:N_nml fliplr(1:N_nml)], [m_nml_ld_a_MC - std_nml_ld_a_MC fliplr(m_nml_ld_a_MC + std_nml_ld_a_MC)], 'r', 'EdgeColor', 'none');alpha(0.2)
patch([1:N_nml fliplr(1:N_nml)], [m_nml_ld_a_BC - std_nml_ld_a_BC fliplr(m_nml_ld_a_BC + std_nml_ld_a_BC)], [0 0.7 0], 'EdgeColor', 'none');alpha(0.2)
patch([1:N_nml fliplr(1:N_nml)], [m_nml_ld_a_VC - std_nml_ld_a_VC fliplr(m_nml_ld_a_VC + std_nml_ld_a_VC)], 'b', 'EdgeColor', 'none');alpha(0.2)
xlim([0 N_nml])
hold off
xlabel('Cortical depth (%)');
ylabel('Normalized vascular length (m/mm^3)');
legend('Motor Cortex', 'Barrel Cortex', 'Visual Cortex');

% Normalized microvascular length (m/mm^3)
figure,  hold on
plot(1:N_nml, m_nml_ld_a_c_MC, 'r', 'LineWidth', 2)
plot(1:N_nml, m_nml_ld_a_c_BC, 'Color', [0 0.7 0], 'LineWidth', 2)
plot(1:N_nml, m_nml_ld_a_c_VC, 'b', 'LineWidth', 2)
patch([1:N_nml fliplr(1:N_nml)], [m_nml_ld_a_c_MC - std_nml_ld_a_c_MC fliplr(m_nml_ld_a_c_MC + std_nml_ld_a_c_MC)], 'r', 'EdgeColor', 'none');alpha(0.2)
patch([1:N_nml fliplr(1:N_nml)], [m_nml_ld_a_c_BC - std_nml_ld_a_c_BC fliplr(m_nml_ld_a_c_BC + std_nml_ld_a_c_BC)], [0 0.7 0], 'EdgeColor', 'none');alpha(0.2)
patch([1:N_nml fliplr(1:N_nml)], [m_nml_ld_a_c_VC - std_nml_ld_a_c_VC fliplr(m_nml_ld_a_c_VC + std_nml_ld_a_c_VC)], 'b', 'EdgeColor', 'none');alpha(0.2)
xlim([0 N_nml])
hold off
xlabel('Cortical depth (%)');
ylabel('Normalized microvascular length (m/mm^3)');
legend('Motor Cortex', 'Barrel Cortex', 'Visual Cortex');