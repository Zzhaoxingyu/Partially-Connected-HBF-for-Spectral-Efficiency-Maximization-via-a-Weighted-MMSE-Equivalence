clear;
clc;
Nloop = 100;
Nt = 144;
Nr = 32;
for i = 1:Nloop
   [H ,Codebook_v, Codebook_w]  = OMPHWB(Nt,Nr);
   HH(:,:,:,i) = H;
end
save('channel64.mat','H1')
