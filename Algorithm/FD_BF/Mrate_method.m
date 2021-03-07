function [V_ropt, W_ropt] = Mrate_method(Ns,H,Vn,k)
% traditional SVD algorithm for rate maximization

[U,E,V] = svd(H(:,:,k));
ldpow = calc_waterfilling(E(1:Ns, 1:Ns) ,1, Vn);
V_ropt = V(:,1:Ns) * diag(ldpow);
%power constraint
V_ropt = V_ropt / norm(V_ropt,'fro');

W_ropt = U(:,1:Ns);
