function [nF1] = f1Measure(nPrec,nRec)

nF1 = 2*nPrec*nRec/(nPrec+nRec);

end

