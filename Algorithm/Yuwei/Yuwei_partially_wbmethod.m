 function [rate]  = Yuwei_partially_wbmethod(Nt,Nr,Nt_rf,Nr_rf,Vn,Ns,H,Nk,frf_ini,wrf_ini)

V_RF = yuwei_partiallyA1(Nt,Nt_rf,Vn,H,Nk,frf_ini);
Q = (V_RF'*V_RF); 
T = Q^(-0.5);
V_D = zeros(Nt_rf,Ns,Nk);
W_D = zeros(Nr_rf,Ns,Nk);
for k = 1:Nk
    L = H(:,:,k)*V_RF*T;
    [~,D,V] = svd(L);
    [~,IX] = sort(diag(D),'descend');
    M = V(:,IX);
    U = M(:,1:Ns);
    V_D(:,:,k) = T*U;
    V_D(:,:,k)=V_D(:,:,k)/norm(V_RF*V_D(:,:,k),'fro');
end
W_RF  = yuwei_partiallyA2(V_D,V_RF,Nr,Nr_rf,H,Vn,Nk,wrf_ini);
for k = 1:Nk
    J = W_RF'*H(:,:,k)*V_RF*V_D(:,:,k)*V_D(:,:,k)'*V_RF'*H(:,:,k)'*W_RF+Vn*W_RF'*W_RF;
    W_D(:,:,k) = J^(-1)*W_RF'*H(:,:,k)*V_RF*V_D(:,:,k);
end

rate = get_wbrate(V_D,V_RF,W_D,W_RF,Vn,H,Nk);

    