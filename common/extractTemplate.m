% script to exctract templates from image give x,y start stop coordinates
% originally used to extract templates of vesicles

function outputTemplate = extractTemplate(inputImageMat,rstart,rstop,cstart,cstop)

outputTemplate = inputImageMat(rstart:rstop,cstart:cstop);

