function [Y,S_ops,P_ops,C_ops]=quantizer_pds(X,C)

% [Y,S_ops,P_ops,C_ops]=quantizer_pds(X,C) 
% 
% X input vectors 
% C codebook 
% Y output vectors 
% S_ops nb of additions per pixel 
% P_ops nb of multiplications per pixel 
% C_ops nb of comparisons per pixel 
% 
% written by Stefan Gachter, EPFL

[N,k]=size(C); 
% N size of codebook 
% k dimesnion of the quantizer

M=size(X,1); 
% M number of input vectors

Y=zeros(M,k); 
% Y output vectors

S_ops=zeros(M,N); 
P_ops=zeros(M,N); 
C_ops=zeros(M,N);

s_ops_k=zeros(1,k); 
p_ops_k=zeros(1,k); 
c_ops_k=zeros(1,k);

% operations

for m=1:M

d2_min=inf; 
i_min=0; 
for i=1:N 
 
s_ops_k=zeros(1,k); 
p_ops_k=zeros(1,k); 
c_ops_k=zeros(1,k);
d2=0; 
j=1; 
while j<=k
a=X(m,j)-C(i,j); 
b=a*a; 
d2=d2+b;
s_ops_k(j)=2; 
p_ops_k(j)=1; 
c_ops_k(j)=1;

if d2>d2_min
j=k;
end 
j=j+1;
end 
if d2<=d2_min
d2_min=d2; 
i_min=i;
end
S_ops(m,i)=mean(s_ops_k); 
P_ops(m,i)=mean(p_ops_k); 
C_ops(m,i)=mean(c_ops_k);

end 
Y(m,:)=C(i_min,:);

end
S_ops=S_ops./M; 
P_ops=P_ops./M; 
C_ops=C_ops./M;