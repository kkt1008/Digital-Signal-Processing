%============================================================================
% Digital Signal Processing: Final Project
%  Prob. 4: SAR signal processing
%============================================================================
% Programmed by Prof. Jihoon Choi
% Intelligent Signal Processing Lab. (ISPL), Korea Aerospace University (KAU)
%----------------------------------------------------------------------------
% This program is the property of ISPL at KAU. Redistribution is not allowed
% without permission.
% 이 코드는 한국항공대학교 지능신호처리연구실의 자산으로 허가없이 배포할 수 없음
%============================================================================
clear
close all

tic

%% Common parameters 4-1
SetTarget = 2;          % 0: one point, 1: multi-point, 2: NavalAirStation, 3: SolarTower

D_az = 2;               % Downsampling ratio in the azimuth direction
D_rg = 2;               % Downsampling ratio in the range direction
%%

InterpType = 2;         % Interpolation type - 1: intfilt(), 2: firpm()

AddNoise = 1;           % Add noise - 0: off, 1: on
SNR_dB = 3;             % Signal-to-Noise ratio (dB)

I2_az = 2;              % azimuth interpolation ratio for noise reduction
I2_rg = 2;              % rge interpolation ratio for noise reduction

%==========================================================================
% SAR paramaters
%==========================================================================
% Sidelooking Airborne SAR
N_az        = 256;              % Number of range lines [samples], {256 or 512}
N_rg        = N_az;             % Samples per range line
c           = 3e8;              % Velocity of light [m/s]
theta_r_c   = 0.0 * pi/180;     % Squint angle [rad]

R_0         = 20e3;             % Slant range of scene center [m]
Vr          = 120;              % Effective radar velocity [m/s]
La          = 3.3;              % Antenna length [m]
f_0         = 9.6e9;            % Range center frequency [Hz]
lambda      = c / f_0;          % Wavelength [m]

Fr          = 10.0e6;           % Range sampling rate [Hz]
BW          = 10.0e6;           % Bandwidth [Hz]
Fa          = 150;              % Azimuth sampling rate or PRF [Hz]
Tr          = N_az*1.0e-7;      % Transmitted pulse duration [sec]

Kr = BW/Tr;                     % Range FM rate [Hz/sec] 
Ka = (2*Vr^2)/(lambda*R_0);     % Azimuth FM rate [Hz/sec]

%==========================================================================
% Define the target
%==========================================================================
if (SetTarget==0) || (SetTarget==1) % point targets
    if SetTarget==0 % single point
        target_pos = [0 0];
    elseif SetTarget==1 % multi-point targets
        target_pos = [-40,-40; -40,10; 10,10; 30,30]*2;
    end

    N_of_target = size(target_pos,1);   % number of point targets
    X_org0 = zeros(256,256);     % matrix for scatterer position
    for n=1:N_of_target
        ind_az = 128+1 + target_pos(n,1);
        ind_rg = 128+1 + target_pos(n,2);
        X_org0(ind_az,ind_rg) = 1;
    end

    % Option for plot
    opt = 'default';

elseif (SetTarget==2) || (SetTarget==3) % SAR images
    % Load an image 
    if SetTarget==2 % Naval Air Station
        im0 = imread('s8_NavalAirStation.png');
        im0 = rgb2gray(im0);       % convert RGB to gray
        im1 = imresize(im0, 0.5);  % resize the image
    elseif SetTarget==3 % Solar Tower
        im0 = imread('s10_solarTower.png');
        im0 = rgb2gray(im0);       % convert RGB to gray
        im1 = imresize(im0, 0.45); % resize the image        
    end

    % Define the original scene by scaling
    X_org0 = double(im1(1:256,1:256))/255.0;
    N_of_target = 256^2;    

    % Option for plot
    opt = 'gray';
else
    error('Invalid target number!');
end

% Define the scene with guard area
X_org = zeros(N_az,N_rg);     % matrix for scatterer position
X_org(N_az/2-128+(1:256),N_rg/2-128+(1:256)) = X_org0;

%----- Plot the original target image -----
plot_image(X_org0,'Original target image',[],opt);

%----- Plot the target scene with guard area -----
plot_image(X_org,'Target scene with guard area',[],opt);

%==========================================================================
% SAR geometry and antenna parameters 
%==========================================================================
% Azimuth time
eta_c     = 0;                              % Beam crossing time [sec]
eta       = (eta_c - N_az/2/Fa):(1/Fa):(eta_c + (N_az/2-1)/Fa);
R_eta     = sqrt(R_0^2 + Vr^2 * eta.^2);    % range equation
theta_eta = atan(Vr*(eta-eta_c)/R_0);       % angle at beam crossing time
BeamWidth = 0.886 * lambda / La;            % 3-dB beamwidth
Ant_pttn  = abs(sinc(0.886/BeamWidth*theta_eta));   % antenna pattern

%--------------------------------------------------------------------------
% Generation of SAR raw data
%--------------------------------------------------------------------------
Y0 = gen_SAR_rawData(X_org,Tr,Fr,R_0,R_eta,f_0,Kr,Ant_pttn,N_rg,N_az);

%==========================================================================
% Define Range and Azimuth Matched Filters for RDA
%==========================================================================
% Generation of matched filter for RDA
tau = (-Tr/2):(1/Fr):(+Tr/2 - 1/Fr);    % range time variable

NoS_rg    = fix(Tr * Fr);
pulse_st = fix((N_rg-NoS_rg)/2);        % pulse start
pulse_ed = pulse_st + NoS_rg;           % pulse end

h_rg = zeros(1,N_rg);                   % matched filter - range
h_az = zeros(N_az,1);                   % matched filter - azimuth

h_rg(1     , pulse_st+1:pulse_ed) = exp(-1j*pi*Kr*(tau).^2 );
h_az(1:N_az, 1                  ) = exp( 1j*pi*(Ka)*(eta).^2);

% Compression process in frequency domain - range direction
H_rg = fft(h_rg,N_rg);                  % freq response of range matched filter
H_az = fft(h_az,N_az);                  % freq response of azimuth matched filter
P_rg = ones(N_az,1)*H_rg;               % freq-domain range matched filter
P_az = H_az*ones(1,N_rg);               % freq-domain azimuth matched filter


%==========================================================================
% Original RDA
%==========================================================================
X_RDA = run_RDA(Y0,P_rg,P_az,N_az,N_rg);

%----- Plot the original RDA image -----
plot_image(X_RDA,'Original RDA image',[],opt);

%%
%==========================================================================
% P4-1: RDA with 2D Nulling
%--------------------------------------------------------------------------
%       (Example 1) Observed ratio = 50% x 50% if D_az=2, D_rg=2
%       (Example 2) Observed ratio = 25% x 25% if D_az=4, D_rg=4
%==========================================================================
% Generate the raw data matrix Y_dn with 2D nulling
Y_dn = zeros(size(Y0));  % Y_dn을 0X0 행렬로 초기화
Y_dn(1:D_az:end, 1:D_rg:end) = Y0(1:D_az:end, 1:D_rg:end);  % nulling / 2열과 2행을 제외하고 1열, 1행부터 2칸씩 집어넣음           
% 
% (Example) When D_az=2 and D_rg=2,
%           Y0 = [-7 -6 -5 -4;     ==>  Y_dn = [-7  0 -5  0; 
%                 -3 -2 -1  0;                   0  0  0  0;   
%                  1  2  3  4;                   1  0  3  0;  
%                  5  6  7  8;]                  0  0  0  0]
%
%% 4-1
if exist('Y_dn','var')
    % Generate RDA image using Y_dn
    X_RDA_dn = run_RDA(Y_dn,P_rg,P_az,N_az,N_rg); % RDA 알고리즘 실행

    %----- Plot the RDA image with 2D nulling -----
    plot_image(X_RDA_dn,'RDA image with nulling',[],opt);    
end

if exist('X_RDA_org', 'var')
    get_psnr_contrast(X_RDA_org, X_RDA_dn,'1');
end

%% 4-2, 4-3
%==========================================================================
% Design 2D interpolation filter
%==========================================================================
I_az = D_az;  % azimuth interpolation ratio
I_rg = D_rg;  % range interpolation ratio

if InterpType==1 % 4-2
    %===================================================================
    % P4-2: Design the 2D interpolation filter using intfilt(). 
    %===================================================================
    L = 5;          
    alpha = 0.5;    

    h_az = intfilt(L, I_az, alpha);   
    h_rg = intfilt(L, I_rg, alpha);   

elseif InterpType==2 % 4-3
    %===================================================================
    % P4-3: Design the 2D interpolation filter using firpm(). 
    %===================================================================
    Rp = 0.75;   As = 40;
    h_az = design_FIR_pm(Rp,As,I_az); % azimuth interpolation filter
    h_rg = design_FIR_pm(Rp,As,I_rg); % range interpolation filter
end

% Define 2D interpolation filter, h_M, using h_az and h_rg
h_M  = h_az(:) * h_rg(:)'; % 거리와 방위로 된 1D 데이터를 2D로 만듬 => 행렬로 만듬

% Perform 2D interpolation using filter2()
Y_itp = filter2(h_M,Y_dn,'same');  % 아까 nulling 된 데이터를 interpolation 시켜주는 작업, same은 Y_dn과 같은 크기로 제작

if exist('Y_itp','var')
    % Run RDA using the interpolated raw data
    X_RDA_itp = run_RDA(Y_itp,P_rg,P_az,N_az,N_rg);

    %----- Plot the downsampled and then interpolated RDA image -----
    plot_image(X_RDA_itp,'Down-and-Interpolated RDA image',[],opt);
end

%%
%==========================================================================
% Add Noise to SAR Raw Data
%==========================================================================
if AddNoise==1
    sigma   = sqrt(10^(-SNR_dB/10)*(norm(Y0,'fro')^2/(N_az*N_rg)));
    noise   = sigma/sqrt(2)*(randn(size(Y0)) + 1j*randn(size(Y0)));
    Y_noise = Y0 + noise;
    
    % Run RDA using the noisy raw data Y_noise
    X_RDA_noise = run_RDA(Y_noise,P_rg,P_az,N_az,N_rg);

    %----- Plot the RDA image with noise -----
    plot_image(X_RDA_noise,'RDA image with noise',[],opt);
end

%==========================================================================
% Noise Reduction through 2D Lowpass Filtering
% P4-4: Design the 2D lowpass filter using firpm()
%==========================================================================
% 2D upsampling
Y_up = upsample(transpose(Y_noise),I2_rg);
Y_up = upsample(transpose(Y_up),I2_az);


% Design azimuth and range lowpass filters
Rp = 0.5;   As = 60; % 내가 수정해야하는 값
h_LPF_az = design_FIR_pm_second(Rp, As, I2_az);   
h_LPF_rg = design_FIR_pm_second(Rp, As, I2_rg);   

% Define 2D lowpass filter, h_LPF_2D, using h_LPF_az and h_LPF_rg
h_LPF_2D = h_LPF_az(:) * h_LPF_rg(:)';             % 2D lowpass filter

% Perform 2D lowpass filtering using filter2()
Y_LPF_up = filter2(h_LPF_2D, Y_up, 'same');

% 2D downsampling of Y_LPF_up to get Y_LPF with 256x256 pixels
Y_col = downsample(Y_LPF_up.', I2_rg).'; % 열방향 down
Y_LPF = downsample(Y_col, I2_az); % 행방향 down

%%
if exist('Y_LPF','var')
    % Run RDA using the filtered raw data Y_LPF
    X_RDA_LPF = run_RDA(Y_LPF,P_rg,P_az,N_az,N_rg);

    %----- Plot the noise-filtered RDA image -----
    plot_image(X_RDA_LPF,'Noise-filtered RDA image',[],opt);
end


%==========================================================================
% Print PSNR and Contrast
%==========================================================================
% Load the true RDA image
f_name = sprintf('X_RDA_org%d.mat',SetTarget);
if exist(f_name,'file')
    load(f_name);
else
    X_RDA_org = X_org0;
end

fprintf('PSNR and Contrast:\n');
get_psnr_contrast(X_RDA_org,X_RDA,'RDA  ');

if exist('X_RDA_dn','var')
    get_psnr_contrast(X_RDA_org,X_RDA_dn,'P4-1 ');
end

if exist('X_RDA_itp','var')
    if InterpType==1
        get_psnr_contrast(X_RDA_org,X_RDA_itp,'P4-2 ');
    elseif InterpType==2
        get_psnr_contrast(X_RDA_org,X_RDA_itp,'P4-3 ');
    end
end

%%
if AddNoise==1
    get_psnr_contrast(X_RDA_org,X_RDA_noise,'Noise');
end

if exist('X_RDA_LPF','var')
    get_psnr_contrast(X_RDA_org,X_RDA_LPF,'P4-4 ');
end
fprintf('\n');

%%
%=====================================================
% Local functions
%=====================================================
function s_0 = gen_SAR_rawData(s_org,Tr,Fr,R_0,R_eta,f_0,Kr,Ant_pttn,N_rg,N_az)
    c = 3e8;   % [m/s]     Velocity of light

    % Range time
    NoS_rg    = fix(Tr * Fr);
    tau0 = - Tr/2 + (2*R_0)/c;          % range start time
    tau_vec = linspace(tau0,tau0+(NoS_rg-1)*1/Fr,NoS_rg);
    
    s_rg = zeros(1,N_rg);
    s_az = zeros(N_az,1);
    
    pulse_st = fix((N_rg-NoS_rg)/2);    % pulse start
    s_rg(1,pulse_st+(1:NoS_rg)) = exp(-1j*4*pi*f_0*(R_eta(1))./c).* exp(1j*pi*Kr*(tau_vec - 2*R_eta(1)/c).^2);
    s_az(1:N_az,1) = exp(-1j*4*pi*f_0*(R_eta)./c).* exp(1j*pi*Kr*(tau0 - 2*R_eta/c).^2) .* Ant_pttn;
    
    s_rg_conv = convmtx(s_rg,N_rg);     % range convolution matrix
    s_az_conv = convmtx(s_az,N_az);     % azimuth convoluation matrix
    
    s_0_full = s_az_conv*s_org*s_rg_conv; % SAR raw data
    s_0 = s_0_full(N_az/2+(1:N_az),N_rg/2+(1:N_rg)-1); % SAR raw data for RDA
end

% M(Y) = estimate of X (SAR image)
function X_hat = run_RDA(Y,P_rg,P_az,N_az,N_rg)
    Y_rg = P_rg.*fft(Y,N_rg,2);         % matched filtering in range direction
    y_rc = ifft(Y_rg,N_rg,2);           % range-compressed time-domain data
    Y_ac = P_az.*fft(y_rc,N_az,1);      % matched filtering in azimuth direction
    y_ac = ifft(Y_ac,N_az,1);           % azimuth-compressed time-domain data
    y_ac2 = 1/sqrt(N_az*N_rg)*y_ac;     % scaling

    X_hat = fftshift(fftshift(y_ac2,2),1); % 2D fft shift
    X_hat = X_hat(N_az/2-128+(1:256),N_rg/2-128+(1:256));
end

% Plot 2D image
function plot_image(s,str,no,opt)
    if nargin<4
        opt = 'default';
    end
    
    if isempty(no)
        figure;
    else 
        figure(no);
    end
    s_max = max(max(abs(s)));
    image(abs(s)/s_max*255);
    colormap(opt);
    xlabel('Range time (sample)','FontSize',12);
    ylabel('Azimuth time (sample)','FontSize',12);
    title(str,'FontSize',12);
    xlim([1 size(s,2)]);
    ylim([1 size(s,1)]);
end

%% 4-3

% Design FIR filter
function h = design_FIR_pm(Rp,As,I)
    %===============================================================
    % P4-3: Design an FIR filter using firpm().
    %       The filter length is odd.
    %       Plot the impulse response and log-magnitude response.
    %===============================================================
    % Compute the passband and stopband ripples
    delta1 = (10^(Rp/20)-1)/(10^(Rp/20)+1);
    delta2 = (1+delta1)*(10^(-As/20));

    wp = (1/I-0.1)*pi;           
    ws = (1/I+0.1)*pi; 
    
    weights = [delta2/delta1 1];
    deltaf = (ws-wp)/(2*pi);

    f = [0 wp/pi ws/pi 1];
    m = [1 1 0 0];
    M_num = (-20*log10(sqrt(delta1*delta2))-13);
    M = ceil(M_num/(14.6*deltaf)+1);
    
    h = firpm(M-1,f,m,weights);
    [db,mag,pha,grd,w] = freqz_m(h,1);
    delta_w = 2*pi/1000;
   
    wsi = ws/delta_w+1;
    wpi = wp/delta_w;

    Asd = -max(db(round(wsi):1:501));
    [Hr,omega,P,L] = ampl_res(h);

    fprintf('Asd = %5.2f dB\n',Asd);    

    figure(101);
    
    subplot(2,1,1);
    Hs = stem(0:1:M-1,h);
    grid on;
    title('Impulse Response of LPF');
    xlabel('n'); ylabel('Amplitude');

    % Plot log-magnitude response
    subplot(2,1,2);
    plot(w/pi,db,'linewidth',1);
    grid on;
    ylim([-100 5]);
    xlim([0 1]);
    xlabel('Normalized Frequency (\times\pi rad/sample)');
    ylabel('Magnitude (dB)');
    set(gca,'YTickMode','manual','YTick',[-60 -50 -40, 0]);
    title('Log-Magnitude Response of LPF');    
    %===============================================================
end
%%
% Compute and display PSNR and contrast
function get_psnr_contrast(X_org,X_hat,name)
    % PSNR
    X_org_abs = abs(X_org)/max(abs(X_org(:)));
    X_hat_abs = abs(X_hat)/max(abs(X_hat(:)));
    PSNR = psnr(X_hat_abs,X_org_abs,1);

    % Contrast
    Cntr = std(X_hat_abs)/mean(X_hat_abs);
    fprintf(' %s: PSNR = %.2f(dB), Contrast = %.4f\n',name,PSNR,Cntr);
end

%%
% Design FIR filter
function h = design_FIR_pm_second(Rp,As,I)
    %===============================================================
    % P4-3: Design an FIR filter using firpm().
    %       The filter length is odd.
    %       Plot the impulse response and log-magnitude response.
    %===============================================================
    % Compute the passband and stopband ripples
    delta1 = (10^(Rp/20)-1)/(10^(Rp/20)+1);
    delta2 = (1+delta1)*(10^(-As/20));
    
    wc = 0.3;
    w_delta = 0.15;
    wp = (wc-w_delta)*pi;           
    ws = (wc+w_delta)*pi; 
    
    weights = [delta2/delta1 1];
    deltaf = (ws-wp)/(2*pi);

    f = [0 wp/pi ws/pi 1];
    m = [1 1 0 0];
    M_num = (-20*log10(sqrt(delta1*delta2))-13);
    M = ceil(M_num/(14.6*deltaf)+1);
    
    h = firpm(M-1,f,m,weights);
    [db,mag,pha,grd,w] = freqz_m(h,1);
    delta_w = 2*pi/1000;
   
    wsi = ws/delta_w+1;
    wpi = wp/delta_w;
    
    Asd = -max(db(round(wsi):1:501));
    [Hr,omega,P,L] = ampl_res(h);

    fprintf('Asd = %5.2f dB\n',Asd);    
    
    figure(100);
    subplot(2,1,1);
    Hs = stem(0:1:M-1,h);
    grid on;
    title('Impulse Response of LPF');
    xlabel('n'); ylabel('Amplitude');

    % Plot log-magnitude response
    subplot(2,1,2);
    plot(w/pi,db,'linewidth',1);
    grid on;
    ylim([-100 5]);
    xlim([0 1]);
    xlabel('Normalized Frequency (\times\pi rad/sample)');
    ylabel('Magnitude (dB)');
    set(gca,'YTickMode','manual','YTick',[-60 -50 -40, 0]);
    title('Log-Magnitude Response of LPF');    
    %===============================================================
end