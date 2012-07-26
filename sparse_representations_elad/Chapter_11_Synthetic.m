function [Results]=Chapter_11_Synthetic(m,n,k,Exper,OMP_exp,Mask)

%================================================
% This function performs a synthetic test for testing various estimators. We 
% assume that the support cardinality is fixed and known (k). The estimators
% tested are:
% 0. The minimum L2-estomator [1]
% 1. The oracle, error computed empirically [2] and by the formula [3]
% 2. The MAP, error computed empirically [4] and by the formula [5]
% 3. The MMSE, error computed empirically [7] and by the formula [8]
% 4. Approximated MAP by OMP [6]
% 5. Approximated MMSE by Random-OMP [9] and also with a 
%     sparsified version of it [10]. 
% 
% Inputs: 
%   m,n - size of the dictionary
%   k - cardinality of the representations
%   Exper - number of experiment to average
%   OMP_exp - number of Random-OMP runs to average
%   Mask - a binary array stating which estimator to apply 
% 
% Output:
%   Results - an array that contains all the reuslts. 
% ================================================

% Setting parameters, the dictionary, and some other things
if nargin==0
    n=20; m=30; %dictionary size
    k=3; % cardinality of the representations
    OMP_exp=25; % runs of the Rand-OMP
    Exper=1000; % number or experiments to average over
    Mask=ones(1,10);
end;
sigma_x=1; %
A=randn(n,m);
W=diag(1./sqrt(diag(A'*A)));
A=A*W;
PA=pinv(A);
if Mask(4)+Mask(5)+Mask(7)+Mask(8)>=1,
    Omega=GatherSupports([],m,k,[]); %gathering all possible supports
    LL=size(Omega,1);
end;

% Experimenting
Results=zeros(20,Exper,10);
hh=waitbar(0,'Experimenting ... ');
for sigma=0.1:0.1:2,
    for count=1:1:Exper,
        waitbar(((sigma*10-1)*Exper+count)/Exper/20);

        % Creation of the signal example
        S0=randperm(m);
        S0=S0(1:k); % support
        x0=zeros(m,1);
        x0(S0)=randn(k,1)*sigma_x;
        noise=randn(n,1);
        y=A*x0+noise*sigma;
        if Mask(1)==1,
            Results(round(sigma*10),count,1)=(PA*y-x0)'*(PA*y-x0);
        end;

        % Oracle estimation
        if Mask(2)==1,
            AS=A(:,S0);
            Qs=AS'*AS/sigma^2+eye(k)/sigma_x^2;
            IQs=inv(Qs);
            x_oracle=zeros(m,1);
            x_oracle(S0)=IQs*AS'*y/sigma^2;
            Results(round(sigma*10),count,2)=(x_oracle-x0)'*(x_oracle-x0);
        end;
        % Oracle error expression
        if Mask(3)==1,
            AS=A(:,S0);
            Qs=AS'*AS/sigma^2+eye(k)/sigma_x^2;
            IQs=inv(Qs);
            Results(round(sigma*10),count,3)=trace(IQs);
        end;
        
        % Exhaustive MAP 
        if Mask(4)==1
            quality=zeros(LL,1);
            for index=1:1:LL
                S=find(Omega(index,:));
                AS=A(:,S);
                Qs=AS'*AS/sigma^2+eye(k)/sigma_x^2;
                IQs=inv(Qs);
                quality(index)=y'*AS*IQs*AS'*y;
            end;
            pos=find(quality==max(quality),1);
            Smap=find(Omega(pos,:));
            AS=A(:,Smap);
            Qs=AS'*AS/sigma^2+eye(k)/sigma_x^2;
            IQs=inv(Qs);
            x_map=zeros(m,1);
            x_map(Smap)=IQs*AS'*y/sigma^2;
            Results(round(sigma*10),count,4)=(x_map-x0)'*(x_map-x0);
        end;
        
        % Approximated MAP (OMP+adjustment) Result
        if Mask(6)==1
            x_mapest=OMP(A,y,k);
            Smapest=find(x_mapest);
            AS=A(:,Smapest);
            Qs=AS'*AS/sigma^2+eye(k)/sigma_x^2;
            IQs=inv(Qs);
            x_mapest=zeros(m,1);
            x_mapest(Smapest)=IQs*AS'*y/sigma^2;
            Results(round(sigma*10),count,6)=(x_mapest-x0)'*(x_mapest-x0);
        end;
        
        % Exhaustive MMSE 
        if Mask(7)==1
            weight=zeros(LL,1);
            traceIQs=zeros(1,LL);
            x_temp=zeros(m,LL);
            for index=1:1:LL
                S=find(Omega(index,:));
                AS=A(:,S);
                Qs=AS'*AS/sigma^2+eye(k)/sigma_x^2;
                IQs=inv(Qs);
                traceIQs(index)=trace(IQs);
                weight(index)=y'*AS*IQs*AS'*y/2/sigma^4+0.5*log(det(IQs));
                x_temp(S,index)=IQs*AS'*y/sigma^2;
            end;
            weight=weight-max(weight); % because of numerical issues
            weight=exp(weight);
            weight=weight/sum(weight);
            x_mmse=x_temp*weight;
            if isnan((x_mmse-x0)'*(x_mmse-x0))
                pause;
            end;
            Results(round(sigma*10),count,7)=(x_mmse-x0)'*(x_mmse-x0);
        end;
        
        % MMSE error formula
        if Mask(8)==1
            Results(round(sigma*10),count,8)=...
                (sum((x_mmse*ones(1,LL)-x_temp).^2,1)+traceIQs)*weight;
        end;
        
        % MAP error formula
        if Mask(5)==1
            Results(round(sigma*10),count,5)=Results(round(sigma*10),count,7)+...
                (x_map-x_mmse)'*(x_map-x_mmse);
        end;
        
        % Random-OMP results
        if Mask(9)==1
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
            Results(round(sigma*10),count,9)=(x_randomp-x0)'*(x_randomp-x0);
        end;
        
        % Sparsified Random-OMP
        if Mask(10)==1
            [x_sorted,pos]=sort(abs(x_randomp),'descend');
            Ssp=pos(1:k);
            AS=A(:,Ssp);
            Qs=AS'*AS/sigma^2+eye(k)/sigma_x^2;
            IQs=inv(Qs);
            x_sp=zeros(m,1);
            x_sp(Ssp)=IQs*AS'*y/sigma^2;
            Results(round(sigma*10),count,10)=(x_sp-x0)'*(x_sp-x0);
        end;
    end;

    % Creating a figure to summarize the results
    figure(1); clf;
    h=plot(0.1:0.1:2,mean(Results(:,:,2),2)./mean(Results(:,:,1),2),'b');
    % set(h,'LineWidth',2);
    hold on;
    h=plot(0.1:0.1:2,mean(Results(:,:,3),2)./mean(Results(:,:,1),2),'b--');
    % set(h,'LineWidth',2);
    h=plot(0.1:0.1:2,mean(Results(:,:,4),2)./mean(Results(:,:,1),2),'r');
    % set(h,'LineWidth',2);
    h=plot(0.1:0.1:2,mean(Results(:,:,5),2)./mean(Results(:,:,1),2),'r--');
    % set(h,'LineWidth',2);
    h=plot(0.1:0.1:2,mean(Results(:,:,6),2)./mean(Results(:,:,1),2),'rs');
    set(h,'MarkerFaceColor','r');
    h=plot(0.1:0.1:2,mean(Results(:,:,7),2)./mean(Results(:,:,1),2),'g');
    % set(h,'LineWidth',2);
    h=plot(0.1:0.1:2,mean(Results(:,:,8),2)./mean(Results(:,:,1),2),'g--');
    % set(h,'LineWidth',2);
    h=plot(0.1:0.1:2,mean(Results(:,:,9),2)./mean(Results(:,:,1),2),'go');
    set(h,'MarkerFaceColor','g');
    h=plot(0.1:0.1:2,mean(Results(:,:,10),2)./mean(Results(:,:,1),2),'ms');
    set(h,'MarkerFaceColor','m');
    h=xlabel('\sigma');
    set(h,'FontSize',14);
    h=ylabel('Relative Mean-Squared-Error');
    set(h,'FontSize',14);
    set(gca,'FontSize',14);

end;
close(hh);

return;

% ================================================
% ================================================
% ================================================

function [Omega]=GatherSupports(a,m,k,Omega)

% This is a recusrive function that accumulates all the k supports from a
% set of m elements.
if k==0
    ll=size(Omega,1);
    Omega(ll+1,1:length(a))=a;
end;
if m==0 return; end;
if k>0
    Omega=GatherSupports([a,0],m-1,k,Omega);
    Omega=GatherSupports([a,1],m-1,k-1,Omega);
end;

return;

% ================================================

function a=OMP(D,x,L)

% Orthonormal Matching Pursuit with L non-zeros 

[n,K]=size(D);
a=[];
residual=x;
indx=zeros(L,1);
for j=1:1:L,
    proj=D'*residual;
    pos=find(abs(proj)==max(abs(proj)));
    pos=pos(1);
    indx(j)=pos;
    a=pinv(D(:,indx(1:j)))*x;
    residual=x-D(:,indx(1:j))*a;
end;
temp=zeros(K,1); 
temp(indx)=a;
a=sparse(temp);

return;

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
