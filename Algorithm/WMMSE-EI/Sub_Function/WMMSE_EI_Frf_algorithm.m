function   [V_RF] = WMMSE_EI_Frf_algorithm(V_RF,v2,T, He2,M,Nt,Nt_rf,Nk)

[s,~] = size(V_RF);
Nloop = 3;
b1 = zeros(1,Nk);
b2 = zeros(1,Nk);
a1 = zeros(1,Nk);
a2 = zeros(1,Nk);
iter = zeros(Nt,Nt_rf);
for loop = 1:Nloop
    for  q = 1 : Nt_rf % column index
        for p = (q-1) * M + 1 :q*M % row index
            for k = 1:Nk
                V_m = V_RF;
                V_m(:,q) = [];
                Oq = T(:,:,k)^(-1) + v2(k)^(-1) * He2(:,:,k)' * V_m * (V_m') * He2(:,:,k);
                Aq = v2(k)^(-1) * He2(:,:,k) * Oq^(-2) * He2(:,:,k)';
                Bq = 1/M * eye(s) + v2(k)^(-1) * He2(:,:,k) * Oq^(-1) * He2(:,:,k)';
                vrf = V_RF(:,q);
                vrf(p,1) = 0;
                b1(k) = vrf'*Aq(:,p);
                b2(k) = vrf'*Bq(:,p);
                a1(k) = vrf'*Aq*vrf +V_RF(p,q)'*Aq(p,p)*V_RF(p,q);
                a2(k) = vrf'*Bq*vrf +V_RF(p,q)'*Bq(p,p)*V_RF(p,q);
            end
            % c = test(a1,b1,a2,b2);
%            a = annealing_algorithm(a1,a2,b1,b2);
            fun = @(x)(sum(-1*(a1 + 2*real(b1*exp(1i*x)))./(a2 + 2*real(b2*exp(1i*x)))/64));
            [a,~,~,output] = fminbnd(fun,0,2*pi);
            iter(p,q)=iter(p,q)+output.iterations;
            V_RF(p,q)=exp(1i*a);
            b1 = zeros(1,Nk);
            b2 = zeros(1,Nk);
            a1 = zeros(1,Nk);
            a2 = zeros(1,Nk);
        end
    end
end


        
        
        
                
        
    