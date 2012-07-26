function psn = psnr(I1,I2);
eqm=mean((I1(:)-I2(:)).^2);
psn=10*log10(1.0/eqm);
