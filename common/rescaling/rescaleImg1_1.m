function rescaledImg = rescaleImg1_1(inputImg)
% rescale image into the range -1 to +1

inputImg = inputImg - min(min(inputImg));
inputImg = inputImg./(max(max(inputImg)));

rescaledImg = inputImg * 2 - 1;