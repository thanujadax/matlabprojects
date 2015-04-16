function deviationSigma = getPixIntensityDeviationSigma(I1,I2)
if(size(I1)==size(I2))
    
    dI = I2 - I1;
    deviationSigma = std(dI(:));
    
else
    error('I1 and I2 should have the same dimensions')
end