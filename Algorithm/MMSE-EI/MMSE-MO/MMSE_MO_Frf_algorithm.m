function [F_RF, iter] = MMSE_MO_Frf_algorithm(F_RF, w, H1,Nk,Vn,Ns,frf_manifold)

[Nt, Nrf] = size(F_RF);
m = Nt/Nrf;
fai = zeros(1,Nk);

for i = 1:Nk
    fai(i) = w(i)*Vn*m;
end
problem.M = frf_manifold;

problem.cost = @(x)frf_MMSE_cost(x,fai,H1,Ns,Nrf,Nt,Nk);
problem.egrad = @(x)frf_MMSE_egrad(x,fai,H1,Ns,Nrf,Nt,Nk );

L = zeros(Nt,Nrf);
for i = 1 : Nrf %initialize V_RF
    L((i-1) * m + 1:i*m,i) = ones(m,1);
end
L = L(:);
[x,~,iter] = conjugategradient(problem,F_RF(:),L);
F_RF = reshape(x,Nt,Nrf);

