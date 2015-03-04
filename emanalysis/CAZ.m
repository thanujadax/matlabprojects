%Coordinate Analysis Ziqiang
%%
function CAZ(path, mode, origin)

    % Define default parameters:
    % If not specified otherwise, voxel size is 5nm*5nm*10nm = 250 nm^3
    % !!!Be aware that some stacks have pixel size 3.85*3.85     
    % (e.g.:sample901 sample909)
    voxelSize = [5 5 10]; unit1 = 'nm'; unit2 = 'nm^3'; unit3 = 'µm^3';
    DEFAULT_MODE = 'static';  %default partition mode is using fix-sized grid;
    DEFAULT_ORIGIN = 1;     %default origin is Left-Up-Front corner; 
    DEFAULT_GRID = 3000;    %default grid size is 2500 nm;
    DEFAULT_FACTOR = 2;     %default partition factor is 2;

	% INPUT ARGUMENT CONTROL:    
    % if not specified otherwise, get default parameters:
    if nargin < 2
        mode = DEFAULT_MODE;
        origin = DEFAULT_ORIGIN;
    elseif nargin < 3
        origin = DEFAULT_ORIGIN;       
    end
    
    % get path-specified sample code:
    codestr = inputname(1);
    
    % Load data from local disk	
	Asypath = strcat(path,'\Asy.txt');
    Sympath = strcat(path,'\Sym.txt');
    
    % Read the coordinates from '.txt' file with space separated value   
    [Asy_ID, Asy(:,1), Asy(:,2), Asy(:,3)] = coordinateLoad (Asypath);
	[Sym_ID, Sym(:,1), Sym(:,2), Sym(:,3)] = coordinateLoad (Sympath);
    %read the coordinates from '.txt' file with space separated value

    % Get the range of the coordinates
    minX = min(min(Asy(:,1)), min(Sym(:,1))); maxX = max(max(Asy(:,1)), max(Sym(:,1))); rangeX = maxX - minX;
    minY = min(min(Asy(:,2)), min(Sym(:,2))); maxY = max(max(Asy(:,2)), max(Sym(:,2))); rangeY = maxY - minY;
    minZ = min(min(Asy(:,3)), min(Sym(:,3))); maxZ = max(max(Asy(:,3)), max(Sym(:,3))); rangeZ = maxZ - minZ;
    
    xPhysicalRange = voxelSize(1)*rangeX;
    yPhysicalRange = voxelSize(2)*rangeY;
    zPhysicalRange = voxelSize(3)*rangeZ;
    V = (xPhysicalRange* yPhysicalRange*zPhysicalRange)/(10^9);
     
	fprintf ([' Coordinates loaded for ', num2str(length(Asy_ID)+length(Sym_ID)), ' points',...
             '\n With ', num2str(length(Asy_ID)), ' Asymmetric and ', num2str(length(Sym_ID)), ' Symmetric synapses',...
             '\n Physical ranges of each dimension are:\n',...
             ' X:', num2str(xPhysicalRange), ' ', unit1,...
             ' Y:', num2str(yPhysicalRange), ' ', unit1,...
             ' Z:', num2str(zPhysicalRange), ' ', unit1,...
             '\n with voxel size:', num2str(voxelSize(1)), '*', num2str(voxelSize(2)), '*', num2str(voxelSize(3)), ' ', unit2,...
             '\n The total volume taken up by these points are:', num2str(V), ' ', unit3, '\n']);

         
    %MakeMyVar('Total_Volume',V);
   
         
         
	switch mode
        case 'static';    %partition the volume using fix-sized grid (default value is 2500nm^3);
            
            fprintf (['\n partition mode: STATIC'...
                      '\n default grid size is 3000 nm^3']);
            prompt = '\n Please indicate the grid size (with unit:nm^3): ';
            str = input(prompt);
            if isempty(str)
                grid = DEFAULT_GRID;
                fprintf (['\n Statically partition volume with grid size: ', num2str(DEFAULT_GRID), ' nm^3\n']); 
            elseif isnumeric(str);
                grid = str;
                fprintf (['\n Statically partition volume with grid size: ', num2str(str), ' nm^3\n']);
            end
        
            % calculate sub-cube size of each dimension in pixels     
            cube_X_size = grid/voxelSize(1);
            cube_Y_size = grid/voxelSize(2);
            cube_Z_size = grid/voxelSize(3);
            % calculate sub-cube number of each dimension
            cube_X_number = floor(rangeX / cube_X_size);
            cube_Y_number = floor(rangeY / cube_Y_size);
            cube_Z_number = floor(rangeZ / cube_Z_size);
        
        
        case 'dynamic';   %partition the volume by divide each dimension using the same factor (e.g.:2^3, 3^3);  
            fprintf (['\n partition mode: DYNAMIC'...
                      '\n default partition factor is 2: which means 8 sub-cubes']);
            prompt = '\n Please indicate the partition factor: ';
            str = input(prompt);
            if isempty(str)
                factor = DEFAULT_FACTOR;
                fprintf (['\n Dynamically partition volume with factor: ', num2str(DEFAULT_FACTOR),...
                          '\n This will resulting in 8 sub-cubes\n']);
            elseif isnumeric(str);
                factor = str;
                fprintf (['\n Dynamically partition volume with factor: ', num2str(str),...
                          '\n This will resulting in ', num2str(str^3), ' sub-cubes\n']);
            end
            
            % calculate sub-cube size of each dimension in pixels     
            cube_X_size = rangeX / factor;
            cube_Y_size = rangeY / factor;
            cube_Z_size = rangeZ / factor;
            % calculate sub-cube number of each dimension
            cube_X_number = factor; cube_Y_number = factor; cube_Z_number = factor;
	end
    
    %cube_size = single(zeros(3,1));
    cube_size = horzcat(cube_X_size, cube_Y_size, cube_Z_size);
    MakeMyVar('cube_size',cube_size);
         
    % get the sub-cube value (symmetric synapse ratio), depends on the spatial origin
    
    subNumber_A = single(zeros(cube_X_number,cube_Y_number,cube_Z_number));
    subNumber_S = single(zeros(cube_X_number,cube_Y_number,cube_Z_number));
    
	switch origin
        case 0;     %0: get all 8 corners measurements;   
            % ?still working on this part?
            return;
            
        case 1;     %1: origin is Left-Up-Front corner; 
            for k = 1 : cube_Z_number
                zOnset = minZ + (k-1)*cube_Z_size;
                zOffset = zOnset + cube_Z_size;
                for j = 1 : cube_Y_number
                    yOnset = minY + (j-1)*cube_Y_size;
                    yOffset = yOnset + cube_Y_size;
                    for i = 1 : cube_X_number
                        xOnset = minX + (i-1)*cube_X_size;
                        xOffset = xOnset + cube_X_size;
                        subNumber_A(i,j,k) = getSubNum(Asy(:,1),Asy(:,2),Asy(:,3),xOnset,xOffset,yOnset,yOffset,zOnset,zOffset);
                    end
                end
            end
            for k = 1 : cube_Z_number
                zOnset = minZ + (k-1)*cube_Z_size;
                zOffset = zOnset + cube_Z_size;
                for j = 1 : cube_Y_number
                    yOnset = minY + (j-1)*cube_Y_size;
                    yOffset = yOnset + cube_Y_size;
                    for i = 1 : cube_X_number
                        xOnset = minX + (i-1)*cube_X_size;
                        xOffset = xOnset + cube_X_size;
                        subNumber_S(i,j,k) = getSubNum(Sym(:,1),Sym(:,2),Sym(:,3),xOnset,xOffset,yOnset,yOffset,zOnset,zOffset);
                    end
                end
            end

        case 2;     %2: origin is Right-Up-Front corner; 
            for k = 1 : cube_Z_number
                zOnset = minZ + (k-1)*cube_Z_size;
                zOffset = zOnset + cube_Z_size;
                for j = 1 : cube_Y_number
                    yOnset = minY + (j-1)*cube_Y_size;
                    yOffset = yOnset + cube_Y_size;
                    for i = cube_X_number : -1 : 1
                        xOffset = maxX - (i-1)*cube_X_size;
                        xOnset = xOffset - cube_X_size;
                        subNumber_A(i,j,k) = getSubNum(Asy(:,1),Asy(:,2),Asy(:,3),xOnset,xOffset,yOnset,yOffset,zOnset,zOffset);
                    end
                end
            end
            for k = 1 : cube_Z_number
                zOnset = minZ + (k-1)*cube_Z_size;
                zOffset = zOnset + cube_Z_size;
                for j = 1 : cube_Y_number
                    yOnset = minY + (j-1)*cube_Y_size;
                    yOffset = yOnset + cube_Y_size;
                    for i = cube_X_number : -1 : 1
                        xOffset = maxX - (i-1)*cube_X_size;
                        xOnset = xOffset - cube_X_size;
                        subNumber_S(i,j,k) = getSubNum(Sym(:,1),Sym(:,2),Sym(:,3),xOnset,xOffset,yOnset,yOffset,zOnset,zOffset);
                    end
                end
            end
            
        case 3;     %3: origin is Left-Down-Front corner;   
            for k = 1 : cube_Z_number
                zOnset = minZ + (k-1)*cube_Z_size;
                zOffset = zOnset + cube_Z_size;
                for j = cube_Y_number : -1 : 1
                    yOffset = maxY - (j-1)*cube_Y_size;
                    yOnset = yOffset - cube_Y_size;
                    for i = 1 : cube_X_number
                        xOnset = minX + (i-1)*cube_X_size;
                        xOffset = xOnset + cube_X_size;
                        subNumber_A(i,j,k) = getSubNum(Asy(:,1),Asy(:,2),Asy(:,3),xOnset,xOffset,yOnset,yOffset,zOnset,zOffset);
                    end
                end
            end
            for k = 1 : cube_Z_number
                zOnset = minZ + (k-1)*cube_Z_size;
                zOffset = zOnset + cube_Z_size;
                for j = cube_Y_number : -1 : 1
                    yOffset = maxY - (j-1)*cube_Y_size;
                    yOnset = yOffset - cube_Y_size;
                    for i = 1 : cube_X_number
                        xOnset = minX + (i-1)*cube_X_size;
                        xOffset = xOnset + cube_X_size;
                        subNumber_S(i,j,k) = getSubNum(Sym(:,1),Sym(:,2),Sym(:,3),xOnset,xOffset,yOnset,yOffset,zOnset,zOffset);
                    end
                end
            end
            
        case 4;     %4: origin is Right-Down-Front corner;  
            for k = 1 : cube_Z_number
                zOnset = minZ + (k-1)*cube_Z_size;
                zOffset = zOnset + cube_Z_size;
                for j = cube_Y_number : -1 : 1
                    yOffset = maxY - (j-1)*cube_Y_size;
                    yOnset = yOffset - cube_Y_size;
                    for i = cube_X_number : -1 : 1
                        xOffset = maxX - (i-1)*cube_X_size;
                        xOnset = xOffset - cube_X_size;
                        subNumber_A(i,j,k) = getSubNum(Asy(:,1),Asy(:,2),Asy(:,3),xOnset,xOffset,yOnset,yOffset,zOnset,zOffset);
                    end
                end
            end
            for k = 1 : cube_Z_number
                zOnset = minZ + (k-1)*cube_Z_size;
                zOffset = zOnset + cube_Z_size;
                for j = cube_Y_number : -1 : 1
                    yOffset = maxY - (j-1)*cube_Y_size;
                    yOnset = yOffset - cube_Y_size;
                    for i = cube_X_number : -1 : 1
                        xOffset = maxX - (i-1)*cube_X_size;
                        xOnset = xOffset - cube_X_size;
                        subNumber_S(i,j,k) = getSubNum(Sym(:,1),Sym(:,2),Sym(:,3),xOnset,xOffset,yOnset,yOffset,zOnset,zOffset);
                    end
                end
            end
            
        case 5;     %5: origin is Left-Up-Back corner; 
            for k = cube_Z_number : -1 : 1
                zOffset = maxZ - (k-1)*cube_Z_size;
                zOnset = zOffset - cube_Z_size;
                for j = 1 : cube_Y_number
                    yOnset = minY + (j-1)*cube_Y_size;
                    yOffset = yOnset + cube_Y_size;
                    for i = 1 : cube_X_number
                        xOnset = minX + (i-1)*cube_X_size;
                        xOffset = xOnset + cube_X_size;
                        subNumber_A(i,j,k) = getSubNum(Asy(:,1),Asy(:,2),Asy(:,3),xOnset,xOffset,yOnset,yOffset,zOnset,zOffset);
                    end
                end
            end
            for k = cube_Z_number : -1 : 1
                zOffset = maxZ - (k-1)*cube_Z_size;
                zOnset = zOffset - cube_Z_size;
                for j = 1 : cube_Y_number
                    yOnset = minY + (j-1)*cube_Y_size;
                    yOffset = yOnset + cube_Y_size;
                    for i = 1 : cube_X_number
                        xOnset = minX + (i-1)*cube_X_size;
                        xOffset = xOnset + cube_X_size;
                        subNumber_S(i,j,k) = getSubNum(Sym(:,1),Sym(:,2),Sym(:,3),xOnset,xOffset,yOnset,yOffset,zOnset,zOffset);
                    end
                end
            end
            
        case 6;     %6: origin is Right-Up-Back corner; 
            for k = cube_Z_number : -1 : 1
                zOffset = maxZ - (k-1)*cube_Z_size;
                zOnset = zOffset - cube_Z_size;
                for j = 1 : cube_Y_number
                    yOnset = minY + (j-1)*cube_Y_size;
                    yOffset = yOnset + cube_Y_size;
                    for i = cube_X_number : -1 : 1
                        xOffset = maxX - (i-1)*cube_X_size;
                        xOnset = xOffset - cube_X_size;
                        subNumber_A(i,j,k) = getSubNum(Asy(:,1),Asy(:,2),Asy(:,3),xOnset,xOffset,yOnset,yOffset,zOnset,zOffset);
                    end
                end
            end
            for k = cube_Z_number : -1 : 1
                zOffset = maxZ - (k-1)*cube_Z_size;
                zOnset = zOffset - cube_Z_size;
                for j = 1 : cube_Y_number
                    yOnset = minY + (j-1)*cube_Y_size;
                    yOffset = yOnset + cube_Y_size;
                    for i = cube_X_number : -1 : 1
                        xOffset = maxX - (i-1)*cube_X_size;
                        xOnset = xOffset - cube_X_size;
                        subNumber_S(i,j,k) = getSubNum(Sym(:,1),Sym(:,2),Sym(:,3),xOnset,xOffset,yOnset,yOffset,zOnset,zOffset);
                    end
                end
            end
            
        case 7;     %7: origin is Left-Down-Back corner;                       
            for k = cube_Z_number : -1 : 1
                zOffset = maxZ - (k-1)*cube_Z_size;
                zOnset = zOffset - cube_Z_size;
                for j = cube_Y_number : -1 : 1
                    yOffset = maxY - (j-1)*cube_Y_size;
                    yOnset = yOffset - cube_Y_size;
                    for i = 1 : cube_X_number
                        xOnset = minX + (i-1)*cube_X_size;
                        xOffset = xOnset + cube_X_size;
                        subNumber_A(i,j,k) = getSubNum(Asy(:,1),Asy(:,2),Asy(:,3),xOnset,xOffset,yOnset,yOffset,zOnset,zOffset);
                    end
                end
            end
            for k = cube_Z_number : -1 : 1
                zOffset = maxZ - (k-1)*cube_Z_size;
                zOnset = zOffset - cube_Z_size;
                for j = cube_Y_number : -1 : 1
                    yOffset = maxY - (j-1)*cube_Y_size;
                    yOnset = yOffset - cube_Y_size;
                    for i = 1 : cube_X_number
                        xOnset = minX + (i-1)*cube_X_size;
                        xOffset = xOnset + cube_X_size;
                        subNumber_S(i,j,k) = getSubNum(Sym(:,1),Sym(:,2),Sym(:,3),xOnset,xOffset,yOnset,yOffset,zOnset,zOffset);
                    end
                end
            end 
            
        case 8;     %8: origin is Right-Down-Back corner;                                 
            for k = cube_Z_number : -1 : 1
                zOffset = maxZ - (k-1)*cube_Z_size;
                zOnset = zOffset - cube_Z_size;
                for j = cube_Y_number : -1 : 1
                    yOffset = maxY - (j-1)*cube_Y_size;
                    yOnset = yOffset - cube_Y_size;
                    for i = cube_X_number : -1 : 1
                        xOffset = maxX - (i-1)*cube_X_size;
                        xOnset = xOffset - cube_X_size;
                        subNumber_A(i,j,k) = getSubNum(Asy(:,1),Asy(:,2),Asy(:,3),xOnset,xOffset,yOnset,yOffset,zOnset,zOffset);
                    end
                end
            end            
            for k = cube_Z_number : -1 : 1
                zOffset = maxZ - (k-1)*cube_Z_size;
                zOnset = zOffset - cube_Z_size;
                for j = cube_Y_number : -1 : 1
                    yOffset = maxY - (j-1)*cube_Y_size;
                    yOnset = yOffset - cube_Y_size;
                    for i = cube_X_number : -1 : 1
                        xOffset = maxX - (i-1)*cube_X_size;
                        xOnset = xOffset - cube_X_size;
                        subNumber_S(i,j,k) = getSubNum(Sym(:,1),Sym(:,2),Sym(:,3),xOnset,xOffset,yOnset,yOffset,zOnset,zOffset);
                    end
                end
            end           
	end
    
    subRatio = subNumber_S./(subNumber_A + subNumber_S);               
    returnValue(subRatio,codestr);
% 	MakeMyVar('subNumber_A',subNumber_A);
%	MakeMyVar('subNumber_S',subNumber_S); 
%	MakeMyVar('subRatio',subRatio); 
            
         

         
 

    
%     subRatio = uint8(zeros(cube_X_number, cube_Y_number, cube_Z_number);
%     
%     
%     for i = 1 : length(Asy_ID)+length(Sym_ID)
%         
%     
%     
%     
%     
%     for i = 1 : cube_Z_number
%         for j = 1 : cube_Y_number
%             for k = 1 : cube_X_number
%                 cube
%                 
%                 subRatio(i,j,k) = 
%     
%     xCut = minX + rangeX/G; yCut = minY + rangeY/G; zCut = minZ + rangeZ/G;
%     
%     
%     
%     
%     
%     
%     [Asy_subNum] = getSubNum(Asy_X,Asy_Y,Asy_Z,xCut,yCut,zCut);
%     [Sym_subNum] = getSubNum(Sym_X,Sym_Y,Sym_Z,xCut,yCut,zCut);
%          
%     subRatio =  Sym_subNum ./(Asy_subNum + Sym_subNum);
%     
%     MakeMyVar('subRatio',subRatio);

% 	while (1)
%         prompt = 'Providing sub-volume analysis grid number (e.g.:2, 3, 4, etc): ';
%         str = input(prompt,'s');
%         
%         if isnumeric(str)
%             G = str2num(str);      
%         else
%             error('ErrorTAG:TagName', strcat (' Provide only number for sub-volume analysis.'));
%         end
%        
%         if (G == 2);
%             disp('calculating connected-component from segmentation stack.');
%             getCC;
%             disp('done!');
%             break;
%         else
%             disp('connected-component is not generated. You can manually generate it later with getCC(segStack) command.');
%             break;
%         end
%     end
         
    
%     
%     for i = 1:nargin-1
%         if isnumeric(varargin{i})
%             element(i) = varargin{i};      
%         else
%             error('ErrorTAG:TagName', strcat (' Provide only numeric ID number of merge candidates.'));
%         end       
% 	end
         
         
    %get the range of the coordinates
    
    
    %G=[3 3 3]; % Grid (number of grid lines in X, Y, and Z direction
    
    
    
    
    
%     MakeMyVar('ID',Asy_ID);
%     MakeMyVar('X',Asy_X);
%     MakeMyVar('Y',Asy_Y);
%     MakeMyVar('Z',Asy_Z);
    
    
    
    %MakeMyVar('k',k);
end

function [ID,X,Y,Z] = coordinateLoad (cpath)

    fileID = fopen(cpath);
    Content = textscan(fileID,'%s %f %d8 %f %f %f %f %f');
    ID = Content{3};
    X  = Content{4};
    Y  = Content{5};
    Z  = Content{6};

end

function [subNum] = getSubNum (X,Y,Z,xOnset, xOffset, yOnset, yOffset,zOnset,zOffset)

	
    subNum = 0;
    num = length(X);
    
    for i = 1:num
        if X(i) >= xOnset && X(i) < xOffset
            if Y(i) >= yOnset && Y(i) < yOffset
                if Z(i) >= zOnset && Z(i) < zOffset
                    subNum = subNum + 1;
                end
            end
        end
    end

end

% function [subNum] = getSubNum (X,Y,Z,xCut,yCut,zCut)
% 
% 	
%     subNum = double(zeros(2,2,2));
%     num = length(X);
%     
%     for i = 1:num
%         if X(i) < xCut
%             if Y(i) < yCut
%                 if Z(i) < zCut
%                     subNum (1,1,1) = subNum (1,1,1) + 1;
%                 else
%                     subNum (1,1,2) = subNum (1,1,2) + 1;
%                 end
%             else
%                 if Z(i) < zCut
%                     subNum (1,2,1) = subNum (1,2,1) + 1;
%                 else
%                     subNum (1,2,2) = subNum (1,2,2) + 1;
%                 end
%             end
%         else
%             if Y(i) < yCut
%                 if Z(i) < zCut
%                     subNum (2,1,1) = subNum (2,1,1) + 1;
%                 else
%                     subNum (2,1,2) = subNum (2,1,2) + 1;
%                 end
%             else
%                 if Z(i) < zCut
%                     subNum (2,2,1) = subNum (2,2,1) + 1;
%                 else
%                     subNum (2,2,2) = subNum (2,2,2) + 1;
%                 end
%             end
%         end
%     end
%     
% 
% end


function returnValue(ratio,codestr)

    
    
    sizeOfRatio = size(ratio);
    
    if length(sizeOfRatio) == 3
        NumX = sizeOfRatio(1);
        NumY = sizeOfRatio(2);
        NumZ = sizeOfRatio(3);
    else
        NumX = sizeOfRatio(1);
        NumY = sizeOfRatio(2);
        NumZ = 1;
    end
    
    NumOfSubcube = NumX*NumY*NumZ;
    
    dRatio = reshape(ratio,1,NumOfSubcube);
    
    dRatio(dRatio==0) = [];
    
    %value_mean = mean(dRatio)    %mean
    sem = std(dRatio)/sqrt(length(dRatio)); % standard error of the mean
    
    str_mean = strcat(codestr,'_mean');
    str_sem = strcat(codestr,'_sem');
    
    MakeMyVar(str_mean,mean(dRatio));
    MakeMyVar(str_sem,sem);
   
    
    
    
end

function MakeMyVar(VarName,VarValue)
    assignin('base',VarName,VarValue);
end