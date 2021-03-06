% Hough: detecting oriented bars

%% parameters
savefilepath = '/home/thanuja/Dropbox/RESULTS/hough/houghBars/multijob/';
displayIntermediateFigures=1;
%imagePath = '/home/thanuja/matlabprojects/data/mitoData/stem1_48.png';
imagePath = '/home/thanuja/matlabprojects/data/mitoData/stem1_256by256.png';
% imagePath = 'testImgLines.png';
 
invertImg = 1;      % 1 for membrane images that have to be inverted for Hough transform calculation

grayThresholding = 0;       % 1 if the inverted image should be thresholded
grayThreshold = 0.5;

gaussianFiltering = 0;      % 1 if gaussian filtering should be performed on the input image
sigma = 1;
maskSize = 5;

slidingDist = 1;           % the number of pixels to jump


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% very important tuning parameter
thresholdFraction = 0.5;    % fraction of max(H) to be used as a threshold for peaks
                    % consider the fact that max(H) refers to a global
                    % maximum of H which might overlook smaller line
                    % segments in some patches with less support. 0.5 is
                    % recommended
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

barLength = 13; % should be odd
barWidth = 5; % should be odd

orientations = 0:10:170;    
withBackground = 0;     % plot the detected bars with the original image in the background

%% input preprocessing
imgIn = double(imread(imagePath))/255;
%imgIn = imgIn(1:128,1:128);

if(size(size(imgIn),2)>2)
    img = imgIn(:,:,1);
else
    img = imgIn;
end
if(displayIntermediateFigures)
    figure(1);
    imagesc(img);
    colormap('gray');
    title('original')
end
% invert
if(invertImg)
    imgInv = invertImage(img);
else
    imgInv = img;
end
% rescale 0 - 1
imgInv = imgInv/max(max(imgInv));
if(displayIntermediateFigures)
    figure(2);
    imagesc(imgInv);
    title('inverted input')
    colormap('gray');
end
% thresholding
if(grayThresholding == 1)
 imgInv = simpleThreshold(imgInv,grayThreshold);
 if(displayIntermediateFigures)
     figure(8);
     imagesc(imgInv);
     title('inverted input after thresholding')
     colormap('gray');
 end
end

% gaussian smoothening
if(gaussianFiltering==1)
    imgInv = gaussianFilter(imgInv,sigma,maskSize);
    if(displayIntermediateFigures)
        figure(5);
        imagesc(imgInv);
        title('gaussian smoothening');
        colormap('gray');
    end
end
%% Hough
% perform Hough type processing (voting) for oriented bars
for barWidth = 3:2:5
    for barLength = 5:2:25
        str = sprintf('barWidth = %d ; barLength = %d',barWidth,barLength);
        disp(str);
        display('Computing 3D hough space...');
        t0 = cputime;
        houghSpace3D = houghBars_P(imgInv,barLength,barWidth,orientations,slidingDist);
        t1 = cputime;
        display('3D hough space computed!');
        dt = t1 - t0;
        str = sprintf('Time taken for hough space creation = %0.5f s',dt);
        disp(str);
        % houghSpace3D [row col orientation]

        % peak detection
        peaks3D = houghBarPeaks(houghSpace3D,orientations,thresholdFraction...
                                    ,slidingDist,barLength,barWidth);  

        % draw the detected bars on the image
        display('Reconstructing interpreted outline of the image...');
        t2 = cputime;
        output = reconstructHoughBars_P(peaks3D,orientations,barLength,barWidth);
        t3 = cputime;
        display('Reconstruction completed!');
        titlestr = sprintf('Outline L= %d : W = %d',barLength,barWidth);
        h = figure(101);imagesc(output);title(titlestr);
        dt = t3-t2;
        str = sprintf('Time taken for reconstruction = %0.5f s',dt);
        disp(str);
        % save plots
        savefilename = sprintf('L%d-W%d.png',barLength,barWidth);
        savefilename = strcat(savefilepath,savefilename);
        saveas(h,savefilename);
    end
end
