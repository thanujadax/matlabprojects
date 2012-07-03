inputfile1 = '/home/thanuja/matlabprojects/data/curves2/curve1.jpg';
inputfile2 = '/home/thanuja/matlabprojects/data/curves2/curve2.jpg';
inputfile3 = '/home/thanuja/matlabprojects/data/curves2/curve3.jpg';

I = imread(inputfile1,'jpg');
I2 = imread(inputfile2,'jpg');
I3 = imread(inputfile3,'jpg');

I_normalized = double(I)./255.0;
I2_normalized = double(I2)./255.0;
I3_normalized = double(I3)./255.0;

% Inverting black and white
I_normalized = (I_normalized - 255.0) .* (-1.0);
I2_normalized = (I2_normalized - 255.0) .* (-1.0);
I3_normalized = (I3_normalized - 255.0) .* (-1.0);

% figure(1);
% imshow(I);
% 
% figure(2);
% imshow(I2);
% 
% figure(3);
% imshow(I3);

% J1 = conv2(I3_normalized, I_normalized, 'valid');
J1 = normxcorr2(I_normalized, I3_normalized);

% maxJ1 = max(J1);
% 
% J2 = J1 ./ maxJ1(1); % rescale

% rescaling
J2 = J1*1/(max(max(J1)));


figure(4);
imshow(J2);

figure(5);
surf(J2)
%image(J1);