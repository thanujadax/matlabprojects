function h = overlayLabelOnImage(backgroundImage,foregroundImage)

% Inputs:
%   backgroundImage - e.g. raw EM image. should be properly normalized
%   since we use imshow to display the images.
%   foregroundImage - e.g. segmentation to visualize overlaid on EM image


% outputPath = '/home/thanuja/projects/data/myelin/ssSEM/s909/overlaid'
% backgroundImageFileName = '/home/thanuja/projects/data/myelin/ssSEM/s909/raw20150428/raw03.png';
% foregroundImageFileName = '/home/thanuja/projects/data/myelin/ssSEM/s909/output20150428/segmentedMyelin3.png';


% backgroundImage = imread(backgroundImage);
% foregroundImage = imread(foregroundImage);

backgroundImage = invertImage(backgroundImage);
backgroundImage = backgroundImage./255;

% figure;imshowpair(backgroundImage,foregroundImage)
% 
% fusedimg = imfuse(backgroundImage,foregroundImage);
% figure;imshow(fusedimg);

figure;
 imshow(backgroundImage, 'InitialMag', 'fit')
 % Make a truecolor all-green image.
 green = cat(3, zeros(size(backgroundImage)),ones(size(backgroundImage)), zeros(size(backgroundImage)));
 hold 
 h = imshow(green); 
 
 hold off
 
  % Use our influence map as the 
 % AlphaData for the solid green image.
 % normalize forground image
 
 foregroundImage = double(foregroundImage);
 foregroundImage = foregroundImage./(max(max(foregroundImage)));
%  foregroundImage(foregroundImage>0) = 0.1;
 set(h, 'AlphaData', foregroundImage);
 set(gca,'position',[0 0 1 1],'units','normalized');
 
 % save figure