function egrad = wrf_fully_egrad(x,Vn ,O, H,Nrf,Nk,Ns)

Nr = size(H,1);
egrad = zeros(Nr,Nrf,Nk);
W = reshape(x,Nr,Nrf);
WW = (W'*W)^(-1);

for i = 1:Nk
    A = (1/Vn(i)*H(:,:,i)'*W*WW*W'*H(:,:,i)+eye(Ns))^(-1) * O(:,:,i) *...
        (1/Vn(i)*H(:,:,i)'*W*WW*W'*H(:,:,i)+eye(Ns))^(-1);
    B = H(:,:,i)'*W*WW;
    C = B';
    M = 1/Vn(i)*W*C*A*B;
    N = 1/Vn(i)*H(:,:,i)*A*B;
    egrad(:,:,i) = M-N;
end

egrad = sum(egrad,3);
egrad = egrad(:);