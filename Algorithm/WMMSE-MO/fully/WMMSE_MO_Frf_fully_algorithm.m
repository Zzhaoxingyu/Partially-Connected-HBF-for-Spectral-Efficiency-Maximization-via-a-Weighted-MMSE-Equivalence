function [F_RF, iter] = WMMSE_MO_Frf_fully_algorithm(F_RF,v2 ,T, He2,frf_manifold,Nk)

[Nt, Nrf] = size(F_RF);
problem.M = frf_manifold;

problem.cost = @(x)frf_fully_cost(x,v2 ,T, He2,Nrf,Nt,Nk);
problem.egrad = @(x)frf_fully_egrad(x,v2 ,T, He2,Nrf,Nk );

L = ones(Nt,Nrf);
L = L(:);
[x,~,iter] = conjugategradient(problem,F_RF(:),L);
F_RF = reshape(x,Nt,Nrf);
