function cost = wrf_cost(x,v1 ,T, He1,Nrf,Nr,Nk)


x = reshape(x,Nr,Nrf);
for i = 1:Nk
    cost(i) = trace((T(:,:,i)^(-1)+1/v1(i)*T(:,:,i)^(-1)*He1(:,:,i)'*x*x'*He1(:,:,i))^(-1));
end
cost = sum(cost);



%cost = trace((H1'*x*(x'*x)^(-1)*x'*H1/Vn+eye(Ns))^(-1));
end