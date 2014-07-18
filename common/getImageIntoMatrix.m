function IMin = getImageIntoMatrix(imageFileName)

%% Get input image
[IMin_0,pp]=imread(imageFileName);

%% preprocessing (to remove the dark edges at the bottom and to the right)
if (1)
    imSize = size(IMin_0);
    resizeDim = [imSize(1) imSize(2)];
    IMin = IMin_0(1:resizeDim(1), 1:resizeDim(2), 1);
end

%% fixing input image format
IMin=im2double(IMin);
if (length(size(IMin))>2)
    IMin = rgb2gray(IMin);
end
if (max(IMin(:))<2)
    IMin = IMin*255;
end