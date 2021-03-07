function [rate,V_ropt,W_ropt] = Mrate_wbmethod(Nt,Nr,Ns,H,Vn,Nk)

V_ropt = zeros(Nt,Ns,Nk);
W_ropt = zeros(Nr,Ns,Nk);
rate = 0;
for k = 1:Nk
    [V_ropt(:,:,k),W_ropt(:,:,k)] = Mrate_method(Ns,H,Vn,k);
    V_equal = V_ropt(:,:,k);
    W_equal = W_ropt(:,:,k);
    rate = rate + log2(det(eye(Ns) + 1/Vn * pinv(W_equal) *H(:,:,k) * V_equal * V_equal' * H(:,:,k)'*W_equal ));
end
rate =rate/Nk;
