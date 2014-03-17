function [V]=block_to_vect(B,k)

% [v]=block_to_vect(B,k) 
% 
% written by Stefan Gachter, EPFL

[Mb,Nb]=size(B);

S=floor(Nb/k); 
V=zeros(Mb*S,k);

for s=1:S

V((s-1)*Mb+1:s*Mb,1:k)=B(1:Mb,(s-1)*k+1:s*k);
end 