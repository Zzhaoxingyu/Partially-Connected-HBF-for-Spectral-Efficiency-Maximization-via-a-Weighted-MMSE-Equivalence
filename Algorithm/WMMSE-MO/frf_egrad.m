function egrad = frf_egrad(x,v2 ,T, He2,Nrf,Nt,Nk )


W = reshape(x,Nt,Nrf);
for i = 1:Nk
    A = T(:,:,i)^(-1)+1/v2(i)*He2(:,:,i)'*W*(W')*He2(:,:,i);
    egrad(:,:,i) = -1/v2(i)*He2(:,:,i)*A^(-2)*He2(:,:,i)'*W;
end
egrad = sum(egrad,3);
egrad = egrad(:);