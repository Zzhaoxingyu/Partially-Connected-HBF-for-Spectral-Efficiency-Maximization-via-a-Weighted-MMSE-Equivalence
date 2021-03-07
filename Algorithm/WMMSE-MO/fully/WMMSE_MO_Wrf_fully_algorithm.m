function [W_RF, iter] = WMMSE_MO_Wrf_fully_algorithm(W_RF,v1 ,T, He1,wrf_manifold,Nk,Ns)

[Nr, Nrf] = size(W_RF);
problem.M = wrf_manifold;

problem.cost = @(x)wrf_fully_cost(x,v1 ,T, He1,Nrf,Nr,Nk,Ns); 
problem.egrad = @(x)wrf_fully_egrad(x,v1 ,T, He1,Nrf,Nk,Ns);  

L = ones(Nr,Nrf);
L = L(:);
[x,iter] = conjugategradient(problem,W_RF(:),L);

W_RF = reshape(x,Nr,Nrf);