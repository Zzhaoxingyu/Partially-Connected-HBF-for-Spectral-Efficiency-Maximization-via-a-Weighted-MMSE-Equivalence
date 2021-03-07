function cost = frf_cost(x,v2 ,T, He2,Nrf,Nt,Nk)



x = reshape(x,Nt,Nrf);
for i = 1:Nk
    cost(i) = trace((T(:,:,i)^(-1)+1/v2(i)*He2(:,:,i)'*x*x'*He2(:,:,i))^(-1));
end
cost = sum(cost);

%cost = trace((H1'*x*(x'*x)^(-1)*x'*H1/Vn+eye(Ns))^(-1));
end