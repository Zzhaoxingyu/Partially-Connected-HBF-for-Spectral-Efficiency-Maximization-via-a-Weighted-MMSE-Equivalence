function ldpow = calc_waterfilling(sv,P,sigma2)
% sv: singular values of the channel
% P: total input power
% sigma2: noise power
% 

sv2 = diag(sv).^2;
K = length(sv2);
for k = K:-1:1
    mu = P/k+sum(sigma2./sv2(1:k))/k;
    if mu > sigma2/sv2(k)
        break
    end
end
ldpow = zeros(K,1);

ldpow(1:k) = mu-sigma2./sv2(1:k);
