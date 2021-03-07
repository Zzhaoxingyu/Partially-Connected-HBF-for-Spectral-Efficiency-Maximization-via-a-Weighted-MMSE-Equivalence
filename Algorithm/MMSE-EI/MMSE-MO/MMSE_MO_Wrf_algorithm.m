function [W_RF, iter] = MMSE_MO_Wrf_algorithm(W_RF,w, H1,Nk,Vn,Ns,wrf_manifold)

[Nr, Nrf] = size(W_RF);
m = Nr/Nrf;
problem.M = wrf_manifold;
fai = zeros(1,Nk);
for i = 1:Nk
    fai(i) = w(i)*Vn*m;
end

problem.cost = @(x)wrf_MMSE_cost(x,fai,H1,Nrf,Nr,Nk,Ns);
problem.egrad = @(x)wrf_MMSE_egrad(x,fai,H1,Nrf,Nr,Nk,Ns);

L = zeros(Nr,Nrf);
for i = 1 : Nrf %initialize V_RF
    L((i-1) * m + 1:i*m,i) = ones(m,1);
end
L = L(:);
[x,iter] = conjugategradient(problem,W_RF(:),L);

W_RF = reshape(x,Nr,Nrf);