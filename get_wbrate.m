% function Rate = get_wbrate(V_D,V_RF,W_D,W_RF,Vn,H,Nk)
% 
% [~,Ns] = size(W_D(:,:,1));
% rate = zeros(1,Nk);
% for k = 1:Nk
%     V_equal = V_RF * V_D(:,:,k);
%     W_equal = W_RF * W_D(:,:,k);
%     if rank(W_D(:,:,k)'*W_D(:,:,k))~=4
%         rate(1,k) = log2(det(eye(Ns) + 1/Vn *V_equal' * H(:,:,k)'* W_equal* (W_equal'*W_equal+1e-7)^(-1)*(W_equal') *H(:,:,k) * V_equal));
%     else
%         rate(1,k) = log2(det(eye(Ns) + 1/Vn *V_equal' * H(:,:,k)'* W_equal* (W_equal'*W_equal)^(-1)*(W_equal') *H(:,:,k) * V_equal));
%     end
%     
% end
% Rate = real(sum(rate))/Nk;



function Rate = get_wbrate(V_D,V_RF,W_D,W_RF,Vn,H,Nk)

[~,Ns] = size(W_D(:,:,1));
rate = zeros(1,Nk);
for k = 1:Nk
    V_equal = V_RF * V_D(:,:,k);
%     W_equal = W_RF * W_D(:,:,k);
    W_equal = W_RF;
%     rate(1,k) = log2(det(eye(Ns) + 1/Vn * pinv(W_equal) *H(:,:,k) * V_equal * V_equal' * H(:,:,k)'*W_equal ));
    rate(1,k) = log2(det(eye(Ns) + 1/Vn *V_equal' * H(:,:,k)'* W_equal* (W_equal'*W_equal)^(-1)*(W_equal') *H(:,:,k) * V_equal));
    
end
Rate = real(sum(rate))/Nk;
    
    