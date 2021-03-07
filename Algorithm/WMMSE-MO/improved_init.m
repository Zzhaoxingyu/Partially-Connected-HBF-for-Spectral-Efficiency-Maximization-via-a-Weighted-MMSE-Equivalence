function [W_RF,WD,F_RF,FD, T] = improved_init(H,Vn,Ns,Nk,ini_F_RF,ini_W_RF,ini_W_D,ini_F_D)

W_RF = ini_W_RF;
F_RF = ini_F_RF;
FD = ini_F_D;
WD = ini_W_D;
T = zeros(Ns,Ns,Nk);
for k = 1 : Nk
    Fe = F_RF*FD(:,:,k);
    We = W_RF*WD(:,:,k);
    E = eye(Ns) -  Fe'*H(:,:,k)'*We-We'*H(:,:,k)*Fe+Vn*We'*We+We'*H(:,:,k)*Fe*Fe'*H(:,:,k)'*We;
    T(:, :, k) = E^(-1);
%     He1 = H(:,:,k)*Fe;
%     T(:, :, k)= eye(Ns) + 1/Vn * He1'* W_RF * (W_RF'*W_RF)^(-1)*(W_RF') *He1;
end
