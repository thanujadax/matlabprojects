% sample subvolumes

% Inputs:
%   spreadsheet - containing msbs

% Outputs
%   numbers of msbs of 4 equal subvolumes

% _ _        _ _
%|1|2|      |5|6|
% - -   ,    - -
%|3|4|      |7|8|
% - -        - -

spreadSheetFileName = '/home/thanuja/Dropbox/PROJECTS/MSB/annotations/xls/s108_msb_mss.xls';

% read spreadsheet
% x,y,z of each msb
[num,txt,raw] = xlsread(spreadSheetFileName);

%% define key points of subvolumes (in pixel coordinates)
% (local origins)

% z
z1 = 1;
z2 = 401;
z3 = 800;

% x
x1 = 1;
x2 = 1;
x3 = 1;

% y
y1 = 1;
y2 = 1;
y3 = 1;

% general image width 
xRange = 1600;
yRange = 1600;

%% sample subvolumes
for i=2:size(num,1)
    % iterate through the list of coordinates given in the xls file
    z = num(i,2);
    x = num(i,3);
    y = num(i,4);
    
    if(z<z2)
        % should be 1 or 2 or 3 or 4
        xL = (x1+x2)/2;
        xR = xL + xRange;
        xM = (xL+xR)/2; % TODO replace with weighted avg depending on current section
        yM = ((y1+y2)/2 + yRange)/2;
        
        if(x<xM)
            % should be 1 or 3
             
            
        else
            % should be 2 or 4
        end
    else
        % should be 5 or 6 or 7 or 8
        
    end
    
end