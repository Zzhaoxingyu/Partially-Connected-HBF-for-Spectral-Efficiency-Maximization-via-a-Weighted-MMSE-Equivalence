function rate = WMMSE_MO_Ite_wbmethod(Nt,Nr,Nt_rf,Nr_rf,Vn,Ns,H,Nk,frf_ini,wrf_ini,frf_manifold,wrf_manifold,ifVFD,ini_F_RF,ini_W_RF,ini_W_D,ini_F_D,fd_ini,Iter_max)
% =========================================================
% File Name: WMMSE_MO_wbmethod.m
% Function: The WMMSE_MO Algorithm ；
% Author：Xingyu Zhao
% Contact：xingyuzhao19@fudan.edu.cn
% Creation Date：7/1/2019
% Last Modification Date：8/9/2020
% Paper:Partially-Connected Hybrid Beamforming for Spectral Efficiency Maximization via a Weighted MMSE Equivalence
% =========================================================


M = Nt/Nt_rf;
Fu = zeros(Nt_rf,Ns,Nk);
v1 = zeros(Nk,1);
v2 = zeros(Nk,1);
xi = zeros(Nk,1);
He1 = zeros(Nr,Ns,Nk);
He2 = zeros(Nt,Ns,Nk);
wmmse = zeros(Nk,1); 
rate = zeros(1,Iter_max);
if ifVFD 
    [W_RF,WD,F_RF,FD, T] = improved_init(H,Vn,Ns,Nk,ini_F_RF,ini_W_RF,ini_W_D,ini_F_D); % MMSE-ini
else
    [W_RF,WD,F_RF,FD, T] = random_init(H,Vn,Nr_rf,Nr,Ns,Nk,frf_ini,wrf_ini,fd_ini); % Random-ini
end

k = 0;
deta = 1;
new_wmmse = 10;
% while ( k<10 && deta>=1e-4)
% while ( k<30 && deta>=1e-4)
while ( k<= Iter_max &&  deta>=1e-4)
    for i = 1:Nk
        v2(i) = Vn * Nt * Nr / Nt_rf / Nr_rf * trace(T(:,:,i) * WD(:,:,i)' * WD(:,:,i));
        He2(:,: ,i) = H(:,:,i)'* W_RF*WD(:,:,i);
    end
    [F_RF] = WMMSE_MO_Frf_algorithm(F_RF,v2 ,T, He2,frf_manifold,Nk);
    for i = 1:Nk
        Fu(:,:,i) = ((F_RF')* He2(:,:,i) * T(:,:,i) * (He2(:,:,i)')* F_RF + v2(i) * eye(Nt_rf))^(-1)*F_RF'*He2(:,:,i)*T(:,:,i);
        xi(i) = (M * norm(Fu(:,:,i),'fro')^2)^(-0.5);
        v1(i) = Vn *xi(i)^(-2) * Nr/Nr_rf;
        He1(:,:,i) = H(:,:,i) * F_RF * Fu(:,:,i);
    end
    [W_RF] =WMMSE_MO_Wrf_algorithm(W_RF,v1 ,T, He1,wrf_manifold,Nk);
    for i = 1:Nk
        WD(:,:,i) = ((W_RF')* He1(:,:,i) * (He1(:,:,i)')* W_RF + v1(i) * eye(Nr_rf))^(-1)*W_RF'*He1(:,:,i)+1e-7;
        T(:,:,i) = (eye(Ns) + v1(i)^(-1) * He1(:,:,i)'* W_RF * (W_RF') *He1(:,:,i));
        wmmse(i)=trace(eye(Ns))-log2(det(T(:,:,i)));
        FD(:,:,i) = xi(i)*Fu(:,:,i);
    end
    old_wmmse = new_wmmse;  
    ave_wmmse = mean(wmmse);
    new_wmmse = ave_wmmse;
    deta = old_wmmse - new_wmmse;
    k = k+1;
    rate(1,k)= get_wbrate(FD,F_RF,WD,W_RF,Vn,H,Nk);
end

% rate = get_wbrate(FD,F_RF,WD,W_RF,Vn,H,Nk);

