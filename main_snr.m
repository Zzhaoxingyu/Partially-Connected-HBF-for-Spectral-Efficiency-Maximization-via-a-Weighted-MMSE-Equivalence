clear;
close all;
SNR_dB = -18:3:0;%range of SNR (dB)
snr = 1./db2pow(SNR_dB);
smax = length(snr);
Nt = 64; % the number of Tx antennas
Nr = 32; % the number of Rx antennas
Ns = 2; % the number of data streams
Nt_rf = 4; % the number of RF chains at the transmitter
Nr_rf = 4; % the number of RF chains at the receiver
Nk = 64; % the number of sub_acrriers
realization = 100; % the number of Monte Carlo simulation

%% The initialization mode of WMMSE-EI and WMMSE-MO
% '1' stands for 'MMSEI-ini' 
% '0' stands for 'Random-ini' 
ini_mode = 1; 


%% Initialize parameter loading

% mmWave channel
if Nt == 64 && Nr ==32
    load('64-32.mat')  %  64x32 MIMO
elseif Nt == 144 && Nr ==32
    load('144-32.mat') %  144x32 MIMO
end


% the initialization of analog precoder
if  Nt==64 && Nt_rf == 4
    load('F_RF(64x4)'); % 64x4 analog precoder
elseif Nt ==144 && Nt_rf ==8
    load('F_RF(144x8)'); % 144x8 analog precoder
end
frf_ini = F_RF;

% the initialization of analog combiner
if Nr ==32 && Nr_rf ==2
    load('W_RF(32x2)'); % 32x2 analog combiner
elseif Nr ==32 && Nr_rf ==4
    load('W_RF(32x4)'); % 32x4 analog combiner
elseif Nr ==32 && Nr_rf == 8
    load('W_RF(32x8)'); % 32x8 analog combiner
end
wrf_ini = W_RF;

% the manifold space of the analog precoder
frf_manifold = complexcirclefactory(Nt*Nt_rf);
% the manifold space of the analog combiner
wrf_manifold = complexcirclefactory(Nr*Nr_rf);

H = HH;


FD_BF_rate = zeros(realization,length(snr));
Yuwei_rate = zeros(realization,length(snr));
MMSE_EI_rate = zeros(realization,length(snr));
WMMSE_MO_rate = zeros(realization,length(snr));
WMMSE_EI_rate = zeros(realization,length(snr));



parfor reali = 1:realization
    for i = 1:smax
        [FD_BF_rate(reali,i),V_ropt,W_ropt] = Mrate_wbmethod(Nt,Nr,Ns,H(:,:,:,reali),snr(i),Nk);
        Yuwei_rate(reali,i) = Yuwei_partially_wbmethod(Nt,Nr,Nt_rf,Nr_rf,snr(i),Ns,H(:,:,:,reali),Nk,frf_ini,wrf_ini);
        [MMSE_EI_rate(reali,i),FD,F_RF,WD,W_RF] = MMSE_EI_wbmethod(Nt,Nr,Nt_rf,Nr_rf,snr(i),Ns,H(:,:,:,reali),Nk,frf_ini,wrf_ini);
        WMMSE_MO_rate(reali,i) = WMMSE_MO_wbmethod(Nt,Nr,Nt_rf,Nr_rf,snr(i),Ns,H(:,:,:,reali),Nk,frf_ini,wrf_ini,frf_manifold,wrf_manifold,ini_mode,F_RF,W_RF,WD,FD);
        WMMSE_EI_rate(reali,i) = WMMSE_EI_wbmethod(Nt,Nr,Nt_rf,Nr_rf,snr(i),Ns,H(:,:,:,reali),Nk,frf_ini,wrf_ini,ini_mode,F_RF,W_RF,WD,FD);
    end
end

FD_BF = real(sum(FD_BF_rate,1))/realization;
Yuwei = real(sum(Yuwei_rate,1))/realization;
MMSE_EI = real(sum(MMSE_EI_rate,1))/realization;
WMMSE_MO = real(sum(WMMSE_MO_rate,1))/realization;
WMMSE_EI = real(sum(WMMSE_EI_rate,1))/realization;

plot(SNR_dB,FD_BF);hold on
plot(SNR_dB,WMMSE_EI)
plot(SNR_dB,WMMSE_MO)
plot(SNR_dB,MMSE_EI)
plot(SNR_dB,Yuwei);
legend('FD-BF','WMMSE-EI','WMMSE-MO','MMSE-EI','HBF in [11]')
xlabel('SNR (dB)')
ylabel('Average Spectral Efficiency (bits/s/Hz)')
grid on

