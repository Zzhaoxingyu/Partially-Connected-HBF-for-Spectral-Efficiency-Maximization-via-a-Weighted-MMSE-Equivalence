function egrad = wrf_egrad(x,v1 ,T, He1,Nrf,Nr,Nk)

W = reshape(x,Nr,Nrf);
for i = 1:Nk
    A = T(:,:,i)^(-1)+1/v1(i)*T(:,:,i)^(-1)*He1(:,:,i)'*W*(W')*He1(:,:,i);
    egrad(:,:,i) = -1/v1(i)*He1(:,:,i)*A^(-2)*T(:,:,i)^(-1)*He1(:,:,i)'*W;
end
egrad = sum(egrad,3);
egrad = egrad(:);
