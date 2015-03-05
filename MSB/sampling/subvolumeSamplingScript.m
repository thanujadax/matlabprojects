function cells_msbSubVolumes = subvolumeSamplingScript()

% sample subvolumes
% key points of each volume have to be defined before using the script

% Inputs:
%   spreadsheet - containing msbs

% Outputs
%   msbs (IDs) of 4 equal subvolumes

% _ _        _ _
%|1|2|      |5|6|
% - -   ,    - -
%|3|4|      |7|8|
% - -        - -
disp('Sub volume analysis..')
disp('CAUTION: make sure local origins are defined properly for this volume!')

disp('Reading xls file')
spreadSheetFileName = '/home/thanuja/Dropbox/PROJECTS/MSB/annotations/xls/s108_msb_mss.xls';
disp(spreadSheetFileName)

% read spreadsheet
% x,y,z of each msb
[num,txt,raw] = xlsread(spreadSheetFileName);
str1 = sprintf('loaded %d points from the xls file',size(num,1));
disp(str1)
% Output data structures
cells_msbSubVolumes{8} = [];
%% define key points of subvolumes (in pixel coordinates)
% NEED TO BE DEFINED FOR EACH VOLUME
% (local origins: 1,2,3 of the volume in terms of global pixel coordinates)

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
disp('Iterating through the list of points ..')
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
            if(y<yM)
                % 1
                %subVol1 = [subVol1; (i-1)];
                cells_msbSubVolumes{1}(end+1) = (i-1); 
            else
                % 3
                cells_msbSubVolumes{3}(end+1) = (i-1);
            end
        else
            % should be 2 or 4
            if(y<yM)
                % 2
                cells_msbSubVolumes{2}(end+1) = (i-1);
            else
                % 4
                cells_msbSubVolumes{4}(end+1) = (i-1);
            end
        end
    else
        % should be 5 or 6 or 7 or 8
        xL = (x3+x2)/2;
        xR = xL + xRange;
        xM = (xL+xR)/2; % TODO replace with weighted avg depending on current section
        yM = ((y3+y2)/2 + yRange)/2;        

        if(x<xM)
            % should be 5 or 7
            if(y<yM)
                % 5
                cells_msbSubVolumes{5}(end+1) = (i-1);
            else
                % 7
                cells_msbSubVolumes{7}(end+1) = (i-1);
            end
        else
            % should be 6 or 8
            if(y<yM)
                % 6
                cells_msbSubVolumes{6}(end+1) = (i-1);
            else
                % 8
                cells_msbSubVolumes{8}(end+1) = (i-1);
            end
        end        
    end
    
end
disp('.. finished iterating the list')