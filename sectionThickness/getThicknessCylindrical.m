function thickness = getThicknessCylindrical(diameterPix,numSections,xyRes)
% Inputs:
% diameterPix - diameter of mitochondria in pixels in xy plane
% numSections - number of sections in which the mitochondria is contained
% xyResolution - pixel width in xy plane (nm): typically 5nm 

% returns the average thickness of a section that contains the parallel
% mitochondria (parallel to the cutting plane) in nm, the same unit as
% xyRes

thickness = diameterPix * xyRes / numSections;