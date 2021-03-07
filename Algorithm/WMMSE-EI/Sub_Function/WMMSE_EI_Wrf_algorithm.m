 function [W_RF] = WMMSE_EI_Wrf_algorithm(W_RF,v1 ,T, He1 ,N,Nr_rf,Nk)

[s,~] = size(W_RF);
Nloop = 3;
b1 = zeros(1,Nk);
b2 = zeros(1,Nk);
a1 = zeros(1,Nk);
a2 = zeros(1,Nk);
for loop = 1:Nloop
    for  n = 1 : Nr_rf % column index
        for m = (n-1) * N + 1 :n*N % row index
            for i = 1:Nk
                V_m = W_RF;
                V_m(:,n) = [];
                On = T(:,:,i)^(-1) + v1(i)^(-1) * T(:,:,i)^(-1) * He1(:,:,i)' * V_m * (V_m') * He1(:,:,i);
                An = v1(i)^(-1) * He1(:,:,i) * On^(-2) * T(:,:,i)^(-1) * He1(:,:,i)';
                Bn = 1/N * eye(s) + v1(i)^(-1) * He1(:,:,i) * On^(-1)*T(:,:,i)^(-1) * He1(:,:,i)';
                vrf = W_RF(:,n);
                vrf(m,1) = 0;
                b1(i) = vrf'*An(:,m);
                b2(i) = vrf'*Bn(:,m);
                a1(i) = vrf'*An*vrf +W_RF(m,n)'*An(m,m)*W_RF(m,n);
                a2(i) = vrf'*Bn*vrf +W_RF(m,n)'*Bn(m,m)*W_RF(m,n);               
            end
%             a = annealing_algorithm(a1,a2,b1,b2);
             fun = @(x)(sum(-1*(a1 + 2*real(b1*exp(1i*x)))./(a2 + 2*real(b2*exp(1i*x)))/64));
             a = fminbnd(fun,0,2*pi);
%             a = test1(a1,a2,b1,b2);
            W_RF(m,n)=exp(1i*a);
            b1 = zeros(1,Nk);
            b2 = zeros(1,Nk);
            a1 = zeros(1,Nk);
            a2 = zeros(1,Nk);
        end
    end
end

        
        
        
                
        
    