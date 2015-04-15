function deviationSigma = getPixIntensityDeviationSigma(I1,I2)
if(size(I1)==sizse(I2))
    d1I = (I2(1:size(I1,1),:)- I1(1:size(I1,1),:));
    d2I = (I2(:,1:size(I1,2))-I1(:,1:size(I,2)));
    deviationSigma = std([d1I(:);d2I(:)]);
else
    error('I1 and I2 should have the same dimensions')
end