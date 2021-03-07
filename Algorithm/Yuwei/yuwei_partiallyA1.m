function V_RF  = yuwei_partiallyA1(Nt,Nt_rf,Vn,H,Nk,frf_ini)

    
%% partially initialize 
M = Nt/Nt_rf;
% for i = 1 : Nrf %initialize V_RF
%     V_RF((i-1) * M + 1:i*M,i) = exp( 1i*unifrnd(0,2*pi,M,1));
% end


V_RF = frf_ini;

for k = 1:Nk
    F(:,:,k) = H(:,:,k)'*H(:,:,k);
end
F = sum(F,3)/Nk;
g = 1/Nt;
a = g/Vn;
t = 0;
for Nloop = 1:20
% for Nloop = 1:2
    for j = 1:Nt_rf
        VRF = V_RF;
        VRF(:,j)=[];
        C = eye(Nt_rf-1)+a*VRF'*F*VRF;
        G = a*F-(a^2)*F*VRF*(C^(-1))*VRF'*F;
        for i = ((j-1)*M+1) : (j*M)
            for l = ((j-1)*M+1) : (j*M)
                if i~=l
                    t = t + G(i,l)*V_RF(l,j);
                end
            end
            if t ==0
                V_RF(i,j)=1;
            else
                V_RF(i,j)=t/abs(t);
            end
            t = 0;
        end
    end
end