%% case 5 each node has individual comm frequency

clc;clear;close all;

set(groot,'defaultAxesTickLabelInterpreter','latex');

%% configuration
Fs = 16000; % sampling frequency
T  = 90;     % time
t  = 0:1/Fs:T;
N  = length(t);

load("simulation path/SecondaryPath_6x6.mat");
load("simulation path/PrimaryPath_1x6.mat");

PrimaryPath = Primary_path;
SecondaryPath = Secondary_path;

%% system parameters
wLen = 512;  % local control filter length
sLen = 256;  % secondary path length
Numnode = 6; % number of node
cLen = 33;  % compensate filter length
muw = 1e-6; % control filter step size
muc = 1e-5; % compensate filter step size

%% noise generation alpha = 1000
noise = randn(N,1);  % random noise
low = 100;
high = 1000;
fil = fir1(63,[2*low/Fs 2*high/Fs]);
Ref = filter(fil,1,noise);   % reference

for i = 1:Numnode
  Dis(i,:) = filter(PrimaryPath(i,:),1,Ref);   % Disturbance         
end

Ref = awgn(Ref,40,'measured');


%% proposed method
Wcsubopt = zeros(Numnode,(wLen+cLen-1));
FedDMCANC = FedMCANC(wLen,SecondaryPath,sLen,Numnode,N,Dis,Ref,cLen,Wcsubopt);
[err,FedDMCANC] = CompensateSP(FedDMCANC,muc);
% plot convergence of compensated filter
% figure;
% index = 0;
% for i =1:6
%     for k = 1:6
%         index = index + 1;
%         subplot(6,6,index);
%         plot(reshape(err(i,k,:),[1,size(err,3)]))
%     end
% end

alpha = 1000;
% muw = 1e-6;

% ideal network
FedDMCANC_ideal = FedDMCANC;
[e_FedDMCANC_ideal,FedDMCANC_ideal] = FedMCANC_166_ideal(FedDMCANC_ideal,muw,alpha,Wcsubopt);

% PCF
Tc = [0.2 0.3 0.4 0.5 0.6 0.7]; % communication interval (seconds) for each node

FedDMCANC_individual = FedDMCANC; % Preserve the original FedDMCANC
[e_FedDMCANC_individual, FedDMCANC_individual] = FedMCANC_166_PCF_individual(FedDMCANC_individual, muw, alpha, Wcsubopt, Tc);



%% plot figure

figure;
for i = 1:6
    dis22 = smooth((Dis(i,1:T*Fs).^2),2000);
    efeddmcanc = smooth((e_FedDMCANC_ideal(i,1:T*Fs).^2),2000);
    efeddmcanc_individual = smooth((e_FedDMCANC_individual(i,1:T*Fs).^2),2000);
    
    mse_feddmcanc = 10*log10(efeddmcanc./dis22);
    mse_feddmcanc_individual = 10*log10(efeddmcanc_individual./dis22);
    
%     figure;
    subplot(3,2,i);
    plot(smooth(mse_feddmcanc(100:end-1000,1),5000));
    hold on;
    plot(smooth(mse_feddmcanc_individual(100:end-1000,1),5000));
    legend('FedDMCANC','FedDMCANC_individual');
    axis([0 inf -inf 10]);
    grid on;
end


% save('Result/Case5','FedDMCANC_ideal','e_FedDMCANC_ideal','e_FedDMCANC_individual','FedDMCANC_individual','Ref','Dis','muw','alpha','Tc');