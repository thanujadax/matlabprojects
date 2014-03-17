function image_vq_pds_param(H,N,k,T,P)

% image vector quantization 
% 
% function image_vq_pds_param(H,N,k,T,P); 
% 
% codebook : exhaustive search algorithm 
% quanitzation : partial distance search algorithm 
% 
% k dimension of the quantizer 
% N size of the codebook 
% H number of iterations 
% T training images 
% P image 
% 
% written by Stefan Gachter, EPFL

% trainings set 
% convert into vector

vect_t=block_to_vect(T,k);

% create inital codebook 
C_init=init_codebook(vect_t,N);

save(['init_codebook' num2str(H) '_' num2str(N) '_' num2str(k) '.mat'], 'C_init')

C_old=C_init; 
D_init=zeros(H,1);

for h=1:H

[C_new,D(h)]=codebook_lbg(C_old,vect_t); 
C_old=C_new;
end
C=C_new;

save(['codebook' num2str(H) '_' num2str(N) '_' num2str(k) '.mat'], 'C') 
save(['distortion' num2str(H) '_' num2str(N) '_' num2str(k) '.mat'], 'D')

% quantize

% quantizing set 
% convert into vector 
vect_p=block_to_vect(P,k);

[Y,S_ops,P_ops,C_ops]=quantizer_pds(vect_p,C);

save(['operations' num2str(H) '_' num2str(N) '_' num2str(k) '.mat'], 'S_ops', 'P_ops', 'C_ops')

% convert into block

[Mp,Np]=size(P); 
I=vect_to_block(Y,Mp);

save(['image' num2str(H) '_' num2str(N) '_' num2str(k) '.mat'], 'I') 