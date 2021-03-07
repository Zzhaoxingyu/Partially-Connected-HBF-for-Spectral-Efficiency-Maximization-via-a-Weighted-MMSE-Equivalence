function [W_RF,WD,F_RF,T] = random_init_v3(H,Vn,Nt,Nt_rf,Nr_rf,Nr,Ns,Nk,frf_ini,wrf_ini,W_ropt,F_ropt)
WD = zeros(Nr_rf,Ns,Nk);
N= Nr/Nr_rf;
T= zeros(Ns,Ns,Nk);
FD = zeros(Nt_rf,Ns,Nk);
F_RF=frf_ini;
W_RF=wrf_ini;
Q = (F_RF'*F_RF); 
T_eff = Q^(-0.5);
for k = 1:Nk
    L = H(:,:,k)*F_RF*T_eff;
    [~,D,V] = svd(L);
    [~,IX] = sort(diag(D),'descend');
    M = V(:,IX);
    U = M(:,1:Ns);
    FD(:,:,k) = T_eff*U;
    FD(:,:,k)=FD(:,:,k)/norm(F_RF*FD(:,:,k),'fro');
    J = W_RF'*H(:,:,k)* F_RF * FD(:,:,k)* FD(:,:,k)'*F_RF'*H(:,:,k)'*W_RF+Vn*(W_RF')*W_RF;
    WD(:,:,k) = J^(-1)*W_RF'*H(:,:,k)*F_RF*FD(:,:,k);
    He1 = H(:,:,k) * F_RF * FD(:,:,k);
    T(:,:,k) = (eye(Ns) + 1/Vn/N *(W_RF') *He1* He1'* W_RF ) ;
end




