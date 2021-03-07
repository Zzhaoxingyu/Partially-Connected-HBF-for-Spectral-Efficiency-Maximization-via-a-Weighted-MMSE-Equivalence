function egrad = frf_MMSE_egrad(x,fai,H1 ,Ns,Nrf,Nt,Nk )


W = reshape(x,Nt,Nrf);
for i = 1:Nk
    A = eye(Ns)+1/fai(i)*H1(:,:,i)'*W*(W')*H1(:,:,i);
    egrad(:,:,i) = -1/fai(i)*H1(:,:,i)*A^(-2)*H1(:,:,i)'*W;
end
egrad = sum(egrad,3);
egrad = egrad(:);