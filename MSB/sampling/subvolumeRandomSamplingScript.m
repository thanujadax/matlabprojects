% random sample subvolumes

% Inputs:
%   spreadsheet - containing msbs
%   size of subvolume: x,y,z
%   how many subvolumes

% Outputs
%   numbers of msbs of 4 equal subvolumes

spreadSheetFileName = '/home/thanuja/Dropbox/PROJECTS/MSB/annotations/xls/s108_msb_mss.xls';

xsize = 600;
ysize = 600;
zsize = 300;
numSubVolumes = 8;

halfSizeX = 800;
halfSizeY = 800;
halfSizeZ = 400;

marginX1 = 100;
marginX2 =100;
marginY1 = 100;
marginY2 =100;
marginZ1 = 50;
marginZ2 =50;

% generate origin for sub volume
x0 = rand(1)*(halfSizeX-marginX1) + marginX2;
y0 = rand(1)*(halfSizeY-marginY1) + marginY2;
z0 = rand(1)*(halfSizeZ-marginZ1) + marginZ2;

xend = x0 + xsize -1;
yend = y0 + ysize -1;
zend = z0 + zsize -1;

% read spreadsheet
% x,y,z of each msb
[num,txt,raw] = xlsread(spreadSheetFileName);

for i=2:size(num,2)
    z = num(i,2);
    x = num(i,3);
    y = num(i,4);
    
    if(inSubVolume)
        % add to the list of point in this subvolume
end