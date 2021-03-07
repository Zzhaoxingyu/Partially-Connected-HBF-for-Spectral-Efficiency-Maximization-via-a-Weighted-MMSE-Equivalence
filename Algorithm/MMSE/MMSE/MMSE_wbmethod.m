function [rate,V_mopt,W_mopt] = MMSE_wbmethod(Nt,Nr,Ns,H,Vn,Nk )

%the algorithm for fully_digital design for WB case
%regard every sub_carrier problem as a narrowband case
%cite the paper Generalized linear precoder and decoder design for MIMO channels using the
%weighted mmse criterion



V_mopt = zeros(Nt,Ns,Nk);
W_mopt = zeros(Nr,Ns,Nk);
rate = 0;

for k = 1:Nk
    [V_mopt(:,:,k),W_mopt(:,:,k)] = MMSE_method(Ns,Vn,H,k);
    V_equal = V_mopt(:,:,k);
    W_equal = W_mopt(:,:,k);
    rate = rate + log2(det(eye(Ns) + 1/Vn * pinv(W_equal) *H(:,:,k) * V_equal * V_equal' * H(:,:,k)'*W_equal ));
end
rate =rate/Nk;
