function W_RF  = yuwei_partiallyA2(V_D,V_RF,Nr,Nr_rf,H,Vn,Nk,wrf_ini)

%% partially initialize V_RF
N = Nr/Nr_rf;
% for i = 1 : Nrf
%     W_RF((i-1) * N + 1:i*N,i) = exp( 1i*unifrnd(0,2*pi,N,1));
% end


W_RF = wrf_ini;

for k = 1:Nk
    F(:,:,k) = H(:,:,k)*V_RF*V_D(:,:,k)*V_D(:,:,k)'*V_RF'*H(:,:,k)';
end
F = sum(F,3)/Nk;
g = Nr_rf/Nr;
a = g/Vn;
t = 0;
for Nloop = 1:20
% for Nloop = 1:2
    for j = 1:Nr_rf
        WRF = W_RF;
        WRF(:,j)=[];
        C = eye(Nr_rf-1)+a*WRF'*F*WRF;
        G = a*F-(a^2)*F*WRF*(C^(-1))*WRF'*F;
        for i = ((j-1)*N+1) : (j*N)
            for l = ((j-1)*N+1) : (j*N)
                if i~=l
                    t = t + G(i,l)*W_RF(l,j);
                end
            end
            if t ==0
                W_RF(i,j)=1;
            else
                W_RF(i,j)=t/abs(t);
            end
            t = 0;
        end
    end
end