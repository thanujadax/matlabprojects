% parameter estimation 
% 
% written by Stefan Gachter, EPFL

new;

tic

N=[64]; 
% N size of the codebook 
k=[64]; 
% k dimension of the quantizer 
H=10; 
% H number of iterations

load IMAGES.MAT lena peppers baboon 
load IMG.MAT barb gold

T=[peppers baboon barb gold]; 
% T training images 
P=lena; 
% P image

for i=1:length(N)

for j=1:length(k)
image_vq_pds_param(H,N(i),k(j),T,P);
end
end
toc 
  