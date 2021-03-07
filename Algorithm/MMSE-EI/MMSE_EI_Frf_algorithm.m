function [V_RF] = MMSE_EI_Frf_algorithm(V_RF, w, H1,M,Nk,Vn,Nt_rf,Ns)

[nt,~]=size(H1);
A1 = zeros(nt,nt,Nk);
for k = 1:Nk
     fi = w(k)*Vn*M;
     A1(:,:,k) =1/fi* H1(:,:,k)*(1/M*eye(Ns) + 1/fi*H1(:,:,k)'*H1(:,:,k))^(-1)*H1(:,:,k)';
end
G = sum(A1,3);
func_old = trace(V_RF'*G*V_RF);
deta = 1;
ite = 1;
while (deta>1e-6)
% Nloop = 10;
% for loop = 1:Nloop %inner iteration 
    for n = 1 : Nt_rf % columns
        for i = ((n-1)*M+1) : (n*M) % rows
            vrf = V_RF(:,n);
            vrf(i,1) = 0;
            A = vrf'*G(:,i);
            V_RF(i,n) = exp(1j*(-angle(A)));
        end
    end
    func_new = trace(V_RF'*G*V_RF);
    deta = func_new-func_old;
    func_old = func_new;
    ite = ite+1;
end
    
% end
