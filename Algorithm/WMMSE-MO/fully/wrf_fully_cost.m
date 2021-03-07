function cost = wrf_fully_cost(x,Vn ,O, H,Nrf,Nr,Nk,Ns)



x = reshape(x,Nr,Nrf);
for i = 1:Nk
    cost(i) = trace(O(:,:,i)*(H(:,:,i)'*x*(x'*x)^(-1)*x'*H(:,:,i)/Vn(i)+eye(Ns))^(-1));
end
cost = sum(cost);

%cost = trace((H1'*x*(x'*x)^(-1)*x'*H1/Vn+eye(Ns))^(-1));
end