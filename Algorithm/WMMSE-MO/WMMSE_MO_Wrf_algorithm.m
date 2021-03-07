function [W_RF, iter] = WMMSE_MO_Wrf_algorithm(W_RF,v1 ,T, He1,wrf_manifold,Nk)

[Nr, Nrf] = size(W_RF);
m = Nr/Nrf;
problem.M = wrf_manifold;

problem.cost = @(x)wrf_cost(x,v1 ,T, He1,Nrf,Nr,Nk);
problem.egrad = @(x)wrf_egrad(x,v1 ,T, He1,Nrf,Nr,Nk);

L = zeros(Nr,Nrf);
for i = 1 : Nrf %initialize V_RF
    L((i-1) * m + 1:i*m,i) = ones(m,1);
end
L = L(:);
[x,iter] = conjugategradient(problem,W_RF(:),L);

W_RF = reshape(x,Nr,Nrf);