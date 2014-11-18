function deviation = getIntensityDeviationShiftedImg(imageFileName,maxShift)

I = double(imread(imageFileName));

deviation = zeros(1,maxShift);

for g=1:maxShift
    d1I = (I(1+g:size(I,1),:)-I(1:size(I,1)-g,:));
    d2I = (I(:,1+g:size(I,2))-I(:,1:size(I,2)-g));
    deviation(g) = std([d1I(:);d2I(:)]);    
end