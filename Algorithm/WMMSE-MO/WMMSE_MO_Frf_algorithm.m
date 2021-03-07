function [F_RF, iter] = WMMSE_MO_Frf_algorithm(F_RF,v2 ,T, He2,frf_manifold,Nk)

[Nt, Nrf] = size(F_RF);
m = Nt/Nrf;
problem.M = frf_manifold;

problem.cost = @(x)frf_cost(x,v2 ,T, He2,Nrf,Nt,Nk);
problem.egrad = @(x)frf_egrad(x,v2 ,T, He2,Nrf,Nt,Nk );

L = zeros(Nt,Nrf);
for i = 1 : Nrf %initialize V_RF
    L((i-1) * m + 1:i*m,i) = ones(m,1);
end
L = L(:);
[x,~,iter] = conjugategradient(problem,F_RF(:),L);
F_RF = reshape(x,Nt,Nrf);