function rescaledImage = rescaleImg0_1(originalImg)
% rescale image into the range 0 to 1

originalImg = originalImg - min(min(originalImg));
rescaledImage = originalImg./(max(max(originalImg)));