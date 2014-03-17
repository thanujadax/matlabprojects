% clusterdemo.m - clustering algorithm demonstration program
% (C) copyright 2001 by Yu Hen Hu
% mfiles used: datagen.m
clear all, clf
% generate some data using datagen.m, kmeansf.m
% Nvec: (1xclass) # data in each of the c gaussian distr.
% mean_var: 
%   3 x class: mean (1st 2 rows) and variance of each class.
%   4 x class: mean (1st 2 rows) and variance 1 and variance 2
%   5 x class: mean (1st 2 rows), var 1, var 2, and rotation angle
%      rotation angle is in [0  90)
Nvec=[100 100 100];
mean_var=[...
   0.2   0.2   1.0
   0.2   1.2   0.8
   0.02  0.02  0.03
   0.1   0.1   0.03
   60    60     0];
x=datagen(Nvec,mean_var); % x is 150 x 2
Wtrue=mean_var(1:2,:)';
figure(1),plot(x(:,1),x(:,2),'.r',Wtrue(:,1),Wtrue(:,2),'ob')
legend('data points','True cluster centers')

% apply kmeansf.m to the data
% Input -  W: initial weight vectors  c by N matrix, c: # of clusters
%       -  X: K by N input vectors to be clustered
%       -  er: 0<er<1, fractional error between successive 
%              distortion for convergence. If not specified, default = 0.01
%       -  itmax: maximum iterations before terminate iterations
%              if not specified, default value = c: # of clusters
% Output - W: final weight vectors (code book)
%        - iter: actual number of iterations to converge

c=length(Nvec); [K,N]=size(x); er = 1e-5; itmax=100*c;
done=0;
xmean=mean(x);

while done==0, % while not done yet,
   disp(['current cluster # is set to ' int2str(c)])
   c0=c;
   c=input('Enter a new number or return for no change: ');
   if isempty(c), c=c0; else c0=c; end
   W0=0.1*randn(c,N)+ones(c,1)*xmean; 
   [W,iter]=kmeansf(x,W0,er,itmax);
   if iter<itmax,
      disp(['Kmeans algorithm converges in ' int2str(iter) ' iterations!']);
   end
   
   W1=[W;-5 -5;-5 5;5 -5;5 5];
   [vx,vy]=voronoi(W1(:,1),W1(:,2));
   figure,
   plot(x(:,1),x(:,2),'.g',W(:,1),W(:,2),'*r',Wtrue(:,1),Wtrue(:,2),'ob',vx,vy,'-')
   axis([-1 2 -.5 2]) 
   legend('data samples','converged centers','true centers','cluster boundary')
   done=input('Enter 1 to terminate program, return to continue: ');
   if isempty(done), done=0; end
end

