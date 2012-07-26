function Iout = denoise(I,sigma,Io)
% Inputs:
% I, input noisy image, double, values between 0 and 1
% sigma
% Io, original image (optional) only used for printing the PSNR after each step.
if nargin==2
   Io=I;
end

% output:
% Iout, denoised image

% model parameters
K=512;  % number of dictionary elements
C=80;   
if sigma <= 25/255
   n=9;
elseif sigma <= 50/255
   n=12;
else
   n=16;
end
window_size=32;  % window size for the block matching
thrs=(32*sigma).^2

% load the dictionary
load(sprintf('dicts/dict_n%d.mat',n));

% optimization parameter
threads=4; % number of threads  
J1=20;     % number of passes on the data, first step   
J2=5;      % number of passes on the data, second step
% These two parameters can be greatly reduced for faster speed, for instance, J2=0 or J1=5;

[Iout Dout]=mexDenoise(I,Io,D,sigma,n,C,J1,J2,thrs,window_size,10000,threads,1,C);


