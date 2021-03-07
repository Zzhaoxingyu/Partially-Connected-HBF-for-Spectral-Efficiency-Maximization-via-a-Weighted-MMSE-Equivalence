function [W_RF,WD,F_RF,T] = random_init_v2(H,Vn,Nt,Nt_rf,Nr_rf,Nr,Ns,Nk,frf_ini,wrf_ini,W_ropt,F_ropt)
M = Nt/Nt_rf;
N= Nr/Nr_rf;
WD = zeros(Nr_rf,Ns,Nk);
T = zeros(Ns,Ns,Nk);
Fu = zeros(Nt_rf,Ns,Nk);
xi = zeros(Nk,1);
F_RF=frf_ini;
W_RF=wrf_ini;
for k = 1 : Nk
    v2 = Vn * M  * trace( W_ropt(:,:,k)' * W_ropt(:,:,k));
    He2 = H(:,:,k)'* W_ropt(:,:,k);
    Fu(:,:,k) = ((F_RF')* He2  * (He2')* F_RF + v2 * eye(Nt_rf))^(-1)*F_RF'*He2;
    xi(k) = (M * norm(Fu(:,:,k),'fro')^2)^(-0.5);
%     v1 = Vn *xi(k)^(-2) * N;
    v1 = Vn * N;
%     He1 = H(:,:,k) * F_RF * Fu(:,:,k);
    He1 = H(:,:,k) * F_ropt(:,:,k);
    WD(:,:,k) = ((W_RF')* He1 * (He1')* W_RF + v1 * eye(Nr_rf))^(-1)*W_RF'*He1;
    Fe = xi(k)*F_RF*Fu(:,:,k);
    We = W_RF*WD(:,:,k);
%     E = eye(Ns) -  Fe'*H(:,:,k)'*We-We'*H(:,:,k)*Fe+Vn*(We')*We+We'*H(:,:,k)*Fe*(Fe')*H(:,:,k)'*We;
%     T(:, :, k) = E^(-1);
    He1 = H(:,:,k) * F_RF * Fu(:,:,k);
    v1 = Vn *xi(k)^(-2) * N;
    T(:, :, k) = (eye(Ns) + v1^(-1) * He1'* W_RF * (W_RF') *He1);
end




