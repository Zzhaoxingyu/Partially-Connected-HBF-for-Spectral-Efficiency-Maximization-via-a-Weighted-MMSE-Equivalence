function c = test(a1,b1,a2,b2)

fun = @(x)(sum(-1*(a1 + 2*real(b1*exp(1i*x)))./(a2 + 2*real(b2*exp(1i*x)))));
a = 0:0.001:2;
for i = 1:length(a)
    y(i)=fun(a(i)*pi);
end
plot(a*pi,y)
c=1;