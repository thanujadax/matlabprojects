% rescale z curve using x and y curves

% for each distance point in z (dz), 
% get c.o.c
% find the same distance point in x,y by interpolating the x,y curve (d)
% calculate the multiplication factor where d = k.dz

matFilePath = '/home/thanuja/projects/tests/thickness/similarityCurves/20150512/s108';
fileStr = 'xcorrMat'; % general string that defines the .mat file

distMin = 0;
distMax = 10;

% get x,y cc curve
zDirection = 0; % 0 for x,y direction decay curves
[xyMeanCC,xyCCSD] = makeEnsembleDecayCurveForVolume(matFilePath,fileStr,zDirection);

zDirection = 1; % 0 for x,y direction decay curves
[zMeanCC,zCCSD] = makeEnsembleDecayCurveForVolume(matFilePath,fileStr,zDirection);

method = 'spline';

z1_cc = zMeanCC(3);
z1_map_xy1 = interp1(xyMeanCC,distMin:distMax,z1_cc,method);

k = z1_map_xy1 / z1_cc;

% vq = interp1(x,v,xq,method)