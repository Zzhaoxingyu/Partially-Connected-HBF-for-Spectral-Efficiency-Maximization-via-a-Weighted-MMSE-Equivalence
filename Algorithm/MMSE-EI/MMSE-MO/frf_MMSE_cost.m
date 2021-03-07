function cost = frf_MMSE_cost(x,fai,H1,Ns,Nrf,Nt,Nk)


x = reshape(x,Nt,Nrf);
for i = 1:Nk
    cost(i) = trace((eye(Ns)+1/fai(i)*H1(:,:,i)'*x*x'*H1(:,:,i))^(-1));
end
cost = sum(cost);

%cost = trace((H1'*x*(x'*x)^(-1)*x'*H1/Vn+eye(Ns))^(-1));
end