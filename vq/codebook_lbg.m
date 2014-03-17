function [C_new,D]=codebook_lbg(C_old,X)

% [C_new,D]=codebook_lbg(C_old,X) 
% 
% C_old old codebook 
% X trainings set 
% C_new new codebook 
% D mean squared error 
% 
% written by Stefan Gachter, EPFL 
  
 

[N,k]=size(C_old); 
% N size of codebook 
% k dimension of the quantizer

M=size(X,1); 
% M number of training vectors

Vi_sets=zeros(N,M); 
% Vi Voronoi regions

Di=zeros(N,1); 
% Di Distortion

for m=1:M

i=1; 
d2_ref=sq_euk_dist(X(m,:),C_old(i,:)); 
for j=1:N
d2=sq_euk_dist(X(m,:),C_old(j,:)); 
if d2<d2_ref
d2_ref=d2; 
i=j;
end
end 
Vi_sets(i,m)=1; 
Di(i)=Di(i)+d2_ref;
end
C_new=zeros(N,k); 
for i=1:N

Vi=[]; 
for m=1:M
if Vi_sets(i,m)
Vi=[Vi;X(m,:)];
end
end 
g=size(Vi,1); 
if g~=0
C_new(i,:)=(ones(1,g)*Vi)./g;
end
end
D=(ones(1,N)*Di)/N; 
% D mean squared error 
