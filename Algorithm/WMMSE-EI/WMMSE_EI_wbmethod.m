function rate = WMMSE_EI_wbmethod(Nt,Nr,Nt_rf,Nr_rf,Vn,Ns,H,Nk,frf_ini,wrf_ini,ifVFD,ini_F_RF,ini_W_RF,ini_W_D,ini_F_D)
% =========================================================
% File Name: WMMSE_EI_wbmethod.m
% Function: The WMMSE-EI Algorithm ；
% Author：Xingyu Zhao
% Contact：xingyuzhao19@fudan.edu.cn
% Creation Date：7/1/2019
% Last Modification Date：8/9/2020
% Paper:Partially-Connected Hybrid Beamforming for Spectral Efficiency Maximization via a Weighted MMSE Equivalence
% =========================================================

M = Nt/Nt_rf;
N = Nr/Nr_rf;
Fu = zeros(Nt_rf,Ns,Nk);
v1 = zeros(Nk,1);
v2 = zeros(Nk,1);
xi = zeros(Nk,1);
He1 = zeros(Nr,Ns,Nk);
He2 = zeros(Nt,Ns,Nk);
T= zeros(Ns,Ns,Nk);
wmmse = zeros(Nk,1); 


if ifVFD
    [W_RF,WD,F_RF,FD, T] = improved_init(H,Vn,Ns,Nk,ini_F_RF,ini_W_RF,ini_W_D,ini_F_D); % MMSE-ini
    for k = 1:Nk
        W_equal(:,:,k) = W_RF*WD(:,:,k);
    end
else
    F_RF = frf_ini;
    W_RF = wrf_ini;
   for i = 1:Nk
       Heff = W_RF'*H(:,:,i)*F_RF;
       [~,tmp,V] = svd(Heff);
       ldpow = calc_waterfilling(tmp(1:Ns, 1:Ns) ,1, Vn);
       FD(:,:,i) = V(:,1:Ns) * diag(ldpow);
       FD(:,:,i) = FD(:,:,i)/(norm(F_RF*FD(:,:,i),'fro'));
       Fe = F_RF*FD(:,:,i);
       WD(:,:,i) = ((W_RF')* H(:,:,i)*Fe*Fe'*H(:,:,i)'* W_RF + Vn*N * eye(Nr_rf))^(-1)*W_RF'*H(:,:,i)*Fe;
       We = W_RF*WD(:,:,i);
       W_equal(:,:,i) = We;
       E = eye(Ns) -  Fe'*H(:,:,i)'*We-We'*H(:,:,i)*Fe+Vn*We'*We+We'*H(:,:,i)*Fe*Fe'*H(:,:,i)'*We;
       T(:, :, i) = E^(-1);
   end     
end

k = 1;
deta = 1;
new_wmmse = 10;
while ( k<20 && deta>=1e-4)
    for i = 1:Nk
        v2(i) = Vn * Nt / Nt_rf * trace(T(:,:,i) * W_equal(:,:,i)' * W_equal(:,:,i));
        He2(:,: ,i) = H(:,:,i)'* W_equal(:,:,i);
    end
    [F_RF] = WMMSE_EI_Frf_algorithm(F_RF,v2,T, He2,M,Nt,Nt_rf,Nk);
    for i = 1:Nk
        Fu(:,:,i) = ((F_RF')* He2(:,:,i) * T(:,:,i) * (He2(:,:,i)')* F_RF + v2(i) * eye(Nt_rf))^(-1)*F_RF'*He2(:,:,i)*T(:,:,i);
        xi(i) = ( norm(F_RF*Fu(:,:,i),'fro')^2)^(-0.5);
        v1(i) = Vn *xi(i)^(-2) * Nr/Nr_rf;
        He1(:,:,i) = H(:,:,i) * F_RF * Fu(:,:,i);
    end
    [W_RF] = WMMSE_EI_Wrf_algorithm(W_RF,v1 ,T, He1,N,Nr_rf,Nk);
    for i = 1:Nk
        WD(:,:,i) = ((W_RF')* He1(:,:,i) * (He1(:,:,i)')* W_RF + v1(i) * eye(Nr_rf))^(-1)*W_RF'*He1(:,:,i);
        W_equal(:,:,i) =  W_RF*WD(:,:,i);
        T(:,:,i) = (eye(Ns) + v1(i)^(-1) * He1(:,:,i)'* W_RF * (W_RF') *He1(:,:,i));
        wmmse(i)=trace(eye(Ns))-log(det(T(:,:,i)));
    end
    old_wmmse = new_wmmse;  
    ave_wmmse = mean(wmmse);
    new_wmmse = ave_wmmse;
    deta = old_wmmse - new_wmmse;
    k = k+1;
end

for i = 1:Nk
    FD(:,:,i) = xi(i)*Fu(:,:,i);
end

rate = get_wbrate(FD,F_RF,WD,W_RF,Vn,H,Nk);
