function [B]=vect_to_block(V,M)

% [B]=vect_to_block(V,M) 
% 
% written by Stefan Gachter, EPFL

[Mv,Nv]=size(V);

T=floor(Mv/M); 
B=zeros(M,Nv*T);

for t=1:T

B(1:M,(t-1)*Nv+1:t*Nv)=V((t-1)*M+1:t*M,1:Nv);
end 