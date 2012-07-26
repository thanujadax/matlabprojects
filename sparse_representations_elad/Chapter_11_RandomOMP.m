function [A,S0,x0,y,x_randomp,x_sp]=...
                    Chapter_11_RandomOMP(n,m,k,OMP_exp,sigma)

%================================================
% This function performs a synthetic test of the Random OMP, to demonstrate
% the kind of estimate it produces.
%
% Inputs:
%   m,n - size of the dictionary
%   k - cardinality of the representations
%   OMP_exp - number of Random-OMP runs to average
%
% Output:
%   A - the dictionary
%   S0 - the true support
%   x0 - the true solution
%   y - the given signal
%   x_randomp - the Random-OMP estimate
%   x_sp - the Sparsified Random-OMP estimate
% ================================================

% Setting parameters, the dictionary, and some other things
sigma_x=1;
% A=randn(n,m); W=diag(1./sqrt(diag(A'*A))); A=A*W;
A=dctmtx(m);
A=A(1:2:end,:);
W=diag(1./sqrt(diag(A'*A)));
A=A*W;

% Creation of the signal example
S0=randperm(m);
S0=S0(1:k); % support
x0=zeros(m,1);
x0(S0)=randn(k,1)*sigma_x;
noise=randn(n,1);
y=A*x0+noise*sigma;

% Random-OMP results
C=2*sigma^2*(1+sigma^2/sigma_x^2);
C=1/C;
x_temp=zeros(m,OMP_exp);
for kk=1:1:OMP_exp
    temp=RandOMP(A,y,k,C);
    Sest=find(temp);
    AS=A(:,Sest);
    Qs=AS'*AS/sigma^2+eye(k)/sigma_x^2;
    IQs=inv(Qs);
    x_temp(Sest,kk)=IQs*AS'*y/sigma^2;
end;
x_randomp=mean(x_temp,2);

% Sparsified Random-OMP
[x_sorted,pos]=sort(abs(x_randomp),'descend');
Ssp=pos(1:k);
AS=A(:,Ssp);
Qs=AS'*AS/sigma^2+eye(k)/sigma_x^2;
IQs=inv(Qs);
x_sp=zeros(m,1);
x_sp(Ssp)=IQs*AS'*y/sigma^2;

% Creating a figure to summarize the results

return;

% ================================================
% ================================================
% ================================================

function a=RandOMP(D,x,L,c)

% Orthonormal Matching Pursuit with L non-zeros

[n,K]=size(D);
a=[];
residual=x;
indx=zeros(L,1);
for j=1:1:L,
    proj=D'*residual;
    proj=abs(proj);
    proj=exp(min(c*proj.^2,100));
    proj(indx(1:j-1))=0; % no double choice of atoms
    mm=random_choice(proj/sum(proj));
    indx(j)=mm;
    a=pinv(D(:,indx(1:j)))*x;
    residual=x-D(:,indx(1:j))*a;
end;
temp=zeros(K,1);
temp(indx)=a;
a=sparse(temp);

return;

% ================================================

function m=random_choice(prob)

Ref=cumsum(prob);
x=rand(1);
m=find(x-Ref<0,1);

return;
