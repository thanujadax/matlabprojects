% sample subvolumes

% Inputs:
%   spreadsheet - containing msbs

% Outputs
%   numbers of msbs of 4 equal subvolumes

spreadSheetFileName = '/home/thanuja/Dropbox/PROJECTS/MSB/annotations/xls/s108_msb_mss.xls';

% read spreadsheet
% x,y,z of each msb
[num,txt,raw] = xlsread(spreadSheetFileName);

% define boundaries of subvolumes
% z
v1234_z_start = 1;
v1234_z_stop = 400;

v5678_z_start = 401;
v5678_z_stop = 800;

%x
v13_x_start = 1;
v13_x_stop = 800;

v24_x_start = 801;
v24_x_stop = 1600;

v57_x_start = 1;
v57_x_stop = 800;

v68_x_start = 801;
v68_x_stop = 1600;

%y
v12_y_start = 1;
v12_y_stop = 800;
v34_y_start = 801;
v34_y_stop = 1600;

v56_y_start = 1;
v56_y_stop = 800;
v78_y_start = 801;
v78_y_stop = 1600;



% sample subvolumes
