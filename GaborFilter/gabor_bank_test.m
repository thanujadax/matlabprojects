% Thanuja 03.07.2012
% Gabor filter bank
clear;
%% Parameters
inputfile = '/home/thanuja/matlabprojects/data/mitoData/stem1.tiff';
variance = 0.1;
freqAmplitude = 0.01;
numorientations = 32;
% orientation = [0 pi/8 pi/4 3*pi/8 pi/2 5*pi/8 6*pi/8 7*pi/8];  % pi/2
phase = 0;
linearthresh = 50; % thresholding the Gabor output (magnitude) for denoising
%% Filtering
I0=imread(inputfile,'tiff');
I0 = uint8((double(I0)-255) * -1);
%numorientations = size(orientation,2);
orientation = 1:numorientations;
orientation = orientation * pi/numorientations;
whitearea = [];      % variable used to store white spots

% for loop starts here
for i = 1 : numorientations
    [G,GABOUT]=gaborfilter(I0,variance,freqAmplitude,orientation(i),phase);

    R=real(GABOUT);
    I=imag(GABOUT);
    M=abs(GABOUT);
    P=angle(GABOUT);

    %clear GABOUT;

    % post processing
    k=255/max(max(M));
    absgaborout = k*M;  % Gabor magnitude

    % thresholding
    temp_whitearea = find(absgaborout>linearthresh);
    whitearea = union(whitearea,temp_whitearea);

end
% parfor loop ends here

% returns gabor threshold outputs
%whitearea = union_several(temp_whitearea);
threshold_gabor_aggregate = ones(size(I0)) * 255;
threshold_gabor_aggregate(whitearea)=0;

%% Plot results
% original image (input)
figure();
colormap(gray);
imshow(I0);
title('original image (input)');

% Gabor output - magnitude
figure();
colormap(gray);
image(threshold_gabor_aggregate);
title('aggregate of gabor bank (thresh 50)');