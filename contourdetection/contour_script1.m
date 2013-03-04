% contour detection
inputfile = '/home/thanuja/Dropbox/data/testImg/stem_256x_t02_V.png'; % BW image
% contains the magnitude of the OFBR
% inputfile = '/home/thanuja/Dropbox/data/testImg/testMem1_V.png';
ofbrImg = '/home/thanuja/Dropbox/data/testImg/stem_256x_t02_RGB.png';

%% get an oversegmentation from the orientation filter bank response (OFBR)
I=imread(inputfile);
%thresh = graythresh(I); % Otsu's method to calculate grey threshold to minimize inter class variance
thresh = 0.3;
L = getWatershedSegmentation(I,thresh);

%% extract the segments that correspond to the potential contour segments
% detected by OFBR

