function y=equation(x,a1,a2,b1,b2)
m = size(x,2);
a1 = repmat(a1,m,1);
a2 = repmat(a2,m,1);
x = x';
y=sum((a1 + 2*real(kron(b1,exp(1i*x)))./(a2 + 2*real(kron(b2,exp(1i*x))))),2);
y = y';
end