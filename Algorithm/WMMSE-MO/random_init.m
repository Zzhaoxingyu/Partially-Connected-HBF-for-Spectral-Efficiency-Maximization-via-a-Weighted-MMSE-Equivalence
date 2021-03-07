function [W_RF,WD,F_RF,FD,T] = random_init(H,Vn,Nr_rf,Nr,Ns,Nk,frf_ini,wrf_ini,fd_ini)

N= Nr/Nr_rf;
FD = fd_ini;
WD = zeros(Nr_rf,Ns,Nk);
T = zeros(Ns,Ns,Nk);
F_RF=frf_ini;
W_RF=wrf_ini;
for k = 1 : Nk
    Fe = F_RF*FD(:,:,k);
    WD(:,:,k) =  ((W_RF')* H(:,:,k)*Fe*Fe'*H(:,:,k)'* W_RF + Vn*N * eye(Nr_rf))^(-1)*W_RF'*H(:,:,k)*Fe;
    We = W_RF*WD(:,:,k);
    E = eye(Ns) -  Fe'*H(:,:,k)'*We-We'*H(:,:,k)*Fe+Vn*We'*We+We'*H(:,:,k)*Fe*Fe'*H(:,:,k)'*We;
    T(:, :, k) = E^(-1);
end
