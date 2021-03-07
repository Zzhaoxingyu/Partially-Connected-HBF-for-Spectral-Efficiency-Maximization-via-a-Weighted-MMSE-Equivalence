clear;
close all;
Nt_rf = 16; % the number of RF chains at the transmitter
Nr_rf = 16; % the number of RF chains at the receiver
Nt = 144; % the number of Tx antennas
Nr = 32; % the number of Rx antennas
F_RF = zeros(Nt,Nt_rf);
W_RF = zeros(Nr,Nr_rf);
for i = 1 : Nt_rf %initialize V_RF
    F_RF((i-1) * Nt/Nt_rf + 1:i*Nt/Nt_rf,i) = exp( 1i*unifrnd(0,2*pi,Nt/Nt_rf,1));
end
for i = 1 : Nr_rf %initialize W_RF
    W_RF((i-1) * Nr/Nr_rf + 1:i*Nr/Nr_rf,i) = exp( 1i*unifrnd(0,2*pi,Nr/Nr_rf,1));
end
