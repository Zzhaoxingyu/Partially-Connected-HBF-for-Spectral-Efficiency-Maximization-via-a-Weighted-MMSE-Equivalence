function [rate,V_B,V_RF,W_B,W_RF] = MMSE_MO_wbmethod(Nt,Nr,Nt_rf,Nr_rf,Vn,Ns,H,Nk,frf_1,wrf_1,fd_ini,frf_manifold,wrf_manifold)

M = Nt/Nt_rf;
N = Nr/Nr_rf;
V_B = zeros (Nt_rf,Ns,Nk);
V_U = zeros (Nt_rf,Ns,Nk);
W_B = zeros (Nr_rf,Ns,Nk);
w = zeros(1,Nk);
t = zeros(1,Nk);
v = zeros(1,Nk);
H1 = zeros(Nt,Ns,Nk);
H2 = zeros(Nr,Ns,Nk);
V_equal = zeros(Nt,Ns,Nk);
H_equal = zeros(Ns,Ns,Nk);

[~,~,W_equal] = MMSE_wbmethod(Nt,Nr,Ns,H,Vn,Nk);
for k = 1:Nk
    w(k) = trace(W_equal(:,:,k)' * W_equal(:,:,k));
end
W_RF = wrf_1;
V_RF = frf_1;

trigger = 1;
m_MSE_new = 100;
iter = 1;

% while (iter<=20 && trigger>1e-4) 
while (iter<=20) 
% while (iter<=10 && trigger>1e-4) 
% while (iter<=100) 
    for k = 1 : Nk
        H1(:,:,k) = H(:,:,k)' * W_equal(:,:,k);
    end
    %precoding
    [V_RF] = MMSE_MO_Frf_algorithm(V_RF,w,H1,Nk,Vn,Ns,frf_manifold);
    for k = 1 : Nk
       V_U(:,:,k) =(V_RF'*H1(:,:,k) * H1(:,:,k)'* V_RF+  Vn * w(k)*M *eye(Nt_rf))^(-1)*V_RF'*H1(:,:,k);
       V_equal(:,:,k) = V_RF * V_U(:,:,k);
       v(k) = trace(V_equal(:,:,k) * V_equal(:,:,k)');
       H2(:,:,k) = H(:,:,k) * V_equal(:,:,k);
    end
    %combining
     [W_RF] = MMSE_MO_Wrf_algorithm(W_RF,v,H2,Nk,Vn,Ns,wrf_manifold);
     for k = 1:Nk   
        W_B(:,:,k) =(W_RF'*H2(:,:,k) * H2(:,:,k)'* W_RF+  Vn * v(k)*N *eye(Nr_rf))^(-1)*W_RF'*H2(:,:,k);
        W_equal(:,:,k) = W_RF * W_B(:,:,k);
        w(k) = trace(W_equal(:,:,k)' * W_equal(:,:,k));
        H_equal(:,:,k) = W_equal(:,:,k)' * H2(:,:,k);
    end
    m_MSE_old = m_MSE_new;
    for k = 1:Nk
        t(k) = trace(H_equal(:,:,k) * H_equal(:,:,k)' - H_equal(:,:,k) - H_equal(:,:,k)') + Vn  * v(k) * w(k);
    end
    m_MSE_new = real(sum(t))/Nk;
    trigger = m_MSE_old - m_MSE_new;
    iter = iter + 1;
end
for k = 1:Nk
    V_B(:,:,k) = V_U(:,:,k) / sqrt(v(k));
end

rate = get_wbrate(V_B,V_RF,W_B,W_RF,Vn,H,Nk);

        
    
    
    
    
    
    
        