% Template matching code based on correlation 
% Threshold should be in [0,1]
% YRP 2005

function corr = TemplateMatching(ima1,ima2,threshold)

s = size(ima1);
% Performs Correlation
corr = StandardCorrelation(ima1,ima2);

% Dislays image
figure(1); imagesc(ima1); axis image;
% Displays correlation
figure(2); imagesc(corr); axis image;
% Display detection results
% Creates an RGB image
ima3 = zeros(s(1),s(2),3);
% Allocates Red plane to scene
ima1 = double(ima1);
ima3(:,:,1) = double(ima1-min(min(ima1)))/(max(max(ima1))-min(min(ima1))); %/512 as data must be in [0,1]
% Allocates Green and Blue to correlation
max(max(corr))
ima3(:,:,2) = (corr > threshold);
ima3(:,:,3) = (corr > threshold);
figure(3); imagesc(ima3); axis image;