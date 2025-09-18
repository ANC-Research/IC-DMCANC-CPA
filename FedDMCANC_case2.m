%% case 2 ideal network compare centralized, dfxlms,adfxlms,mgdfxlms compressor noise

clc;clear;close all;

set(groot,'defaultAxesTickLabelInterpreter','latex');

%% configuration
Fs = 16000; % sampling frequency
% T  = 60;     % time
% t  = 0:1/Fs:T;
% N  = length(t);

load("path/simulation path/SecondaryPath_6x6.mat");
load("path/simulation path/PrimaryPath_1x6.mat");

PrimaryPath = Primary_path;
SecondaryPath = Secondary_path;

%% system parameters
wLen = 512;  % local control filter length
sLen = 256;  % secondary path length
Numnode = 6; % number of node
cLen = 33;  % compensate filter length
muw = 3e-6; % control filter step size
muc = 1e-5; % compensate filter step size

%% noise generation
load("compressor_16kHz.mat");
t = tr;
N = length(t);
T = 14;
Ref = yr(1:T*Fs);

for i = 1:Numnode
  Dis(i,:) = filter(PrimaryPath(i,:),1,Ref);   % Disturbance         
end

Ref = awgn(Ref,40,'measured');

%% centralized control
CMcANC = McANC_FxLMS_SIMO(wLen,SecondaryPath,sLen,Numnode,Numnode,N,Dis);
[e_CMANC,CMcANC] = McFxLMS_SIMO_166(CMcANC,Ref,muw);

%% diffusion and augmented diffusion (https://github.com/TianyouLi2023)
% [Diffusion_E_n,ADFxLMS_E_n,~]= Main_JASA_Single_Simulation(1,T,Fs,Ref); % mu_MDLMS = 1e-4;  mu_BDLMS_NBC = 5e-4;

%% MGDFxLMS
MGDFxLMS = DMANC_CompensateSP(wLen,SecondaryPath,sLen,Numnode,N,Dis,cLen);
[~,MGDFxLMS] = CompensateSP(MGDFxLMS,muc);
[e_MGDFxLMS,MGDFxLMS] = DMANC_gradient_166(MGDFxLMS,Ref,muw);

%% proposed method
Wcsubopt = zeros(Numnode,(wLen+cLen-1));
FedDMCANC = FedMCANC(wLen,SecondaryPath,sLen,Numnode,N,Dis,Ref,cLen,Wcsubopt);
[err,FedDMCANC] = CompensateSP(FedDMCANC,muc);
% plot convergence of compensated filter
figure;
index = 0;
for i =1:6
    for k = 1:6
        index = index + 1;
        subplot(6,6,index);
        plot(reshape(err(i,k,:),[1,size(err,3)]))
    end
end

alpha = 1000;
% muw = 1e-6;
[e_FedDMCANC,FedDMCANC] = FedMCANC_166_ideal(FedDMCANC,muw,alpha,Wcsubopt);


%% plot figure

figure;
for i = 1:6
    dis22 = smooth((Dis(i,1:T*Fs).^2),2000);
    ecmanc_adaptive = smooth((e_CMANC(i,1:T*Fs).^2),2000);
    eDFxLMS = smooth((Diffusion_E_n(i,1:T*Fs).^2),2000);
    eADFxLMS = smooth((ADFxLMS_E_n(i,1:T*Fs).^2),2000);
    emgdfxlms = smooth((e_MGDFxLMS(i,1:T*Fs).^2),2000);
    efeddmcanc = smooth((e_FedDMCANC(i,1:T*Fs).^2),2000);
    
    mse_adaptive = 10*log10(ecmanc_adaptive./dis22);
    mse_DFxLMS = 10*log10(eDFxLMS./dis22);
    mse_ADFxLMS = 10*log10(eADFxLMS./dis22);
    mse_mgdfxlms = 10*log10(emgdfxlms./dis22);
    mse_feddmcanc = 10*log10(efeddmcanc./dis22);
    
%     figure;
    subplot(3,2,i);
    plot(smooth(mse_adaptive(100:end-1000,1),5000));
    hold on;
    plot(smooth(mse_DFxLMS(100:end-1000,1),5000));
    hold on;
    plot(smooth(mse_ADFxLMS(100:end-1000,1),5000));
    hold on;
    plot(smooth(mse_mgdfxlms(100:end-1000,1),5000));
    hold on;
    plot(smooth(mse_feddmcanc(100:end-1000,1),5000));
    legend('Centralized','DFxLMS','ADFxLMS','MGDFxLMS','FedDMCANC');
    axis([0 inf -inf 10]);
    grid on;
end

nse_adaptive = zeros(Numnode,T*Fs);
nse_DFxLMS = zeros(Numnode,T*Fs);
nse_ADFxLMS = zeros(Numnode,T*Fs);
nse_mgdfxlms = zeros(Numnode,T*Fs);
nse_feddmcanc = zeros(Numnode,T*Fs);

for i = 1:6
    dis22 = smooth((Dis(i,1:T*Fs).^2),2000);
    ecmanc_adaptive = smooth((e_CMANC(i,1:T*Fs).^2),2000);
    eDFxLMS = smooth((Diffusion_E_n(i,1:T*Fs).^2),2000);
    eADFxLMS = smooth((ADFxLMS_E_n(i,1:T*Fs).^2),2000);
    emgdfxlms = smooth((e_MGDFxLMS(i,1:T*Fs).^2),2000);
    efeddmcanc = smooth((e_FedDMCANC(i,1:T*Fs).^2),2000);
    
    nse_adaptive(i,:)    = 10*log10(ecmanc_adaptive./dis22);
    nse_DFxLMS(i,:)      = 10*log10(eDFxLMS./dis22);
    nse_ADFxLMS(i,:)     = 10*log10(eADFxLMS./dis22);
    nse_mgdfxlms(i,:)    = 10*log10(emgdfxlms./dis22);
    nse_feddmcanc(i,:) = 10*log10(efeddmcanc./dis22);    
end

mse_adaptive = mean(nse_adaptive,1);
mse_DFxLMS = mean(nse_DFxLMS,1);
mse_ADFxLMS = mean(nse_ADFxLMS,1);
mse_mgdfxlms = mean(nse_mgdfxlms,1);
mse_feddmcanc = mean(nse_feddmcanc,1);

figure;
plot(smooth(mse_adaptive(100:end-1000),5000));
hold on;
plot(smooth(mse_DFxLMS(100:end-1000),5000));
hold on;
plot(smooth(mse_ADFxLMS(100:end-1000),5000));
hold on;
plot(smooth(mse_mgdfxlms(100:end-1000),5000));
hold on;
plot(smooth(mse_feddmcanc(100:end-1000),5000));
legend('Centralized','DFxLMS','ADFxLMS','MGDFxLMS','FedDMCANC');
axis([0 inf -inf 10]);
grid on;


% save('Result/Case2','CMcANC','e_CMANC','Diffusion_E_n','ADFxLMS_E_n','MGDFxLMS','e_MGDFxLMS',...
%                     'FedDMCANC','e_FedDMCANC','Ref','Dis','muw','alpha')
