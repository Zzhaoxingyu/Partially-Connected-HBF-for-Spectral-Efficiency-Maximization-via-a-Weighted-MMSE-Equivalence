function [W_RF] = MMSE_EI_Wrf_algorithm(W_RF, w, H1,M,Nk,Vn,Nr_rf,Ns )

[nt,~]=size(H1);
A1 = zeros(nt,nt,Nk);
for k = 1:Nk
     fi = w(k)*Vn*M;
     A1(:,:,k) =1/fi* H1(:,:,k)*(1/M*eye(Ns) + 1/fi*H1(:,:,k)'*H1(:,:,k))^(-1)*H1(:,:,k)';
end
G = sum(A1,3);
func_old = trace(W_RF'*G*W_RF);
deta = 1;
ite = 1;
while (deta>1e-6)
% Nloop = 10;
% for loop = 1:Nloop %inner iteration 
    for n = 1 : Nr_rf % columns
        for i = ((n-1)*M+1) : (n*M) % rows
            vrf = W_RF(:,n);
            vrf(i,1) = 0;
            A = vrf'*G(:,i);
            W_RF(i,n) = exp(1j*(-angle(A)));
        end
    end
    func_new = trace(W_RF'*G*W_RF);
    deta = func_new-func_old;
    func_old = func_new;
    ite = ite+1;
end
    
% end

