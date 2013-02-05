% testing the rainbow with a circle
radius1 = 70;
radius2 = 73;
rshift = 20;
img = zeros(200,200);
% make a thick circle

% take a hough vote for the different orientations

% visualize the votes in HSV and with lines in HSV

for c = 1:200
    r = sqrt(radius2^2 - (c-100)^2);
    r = floor(r);
    if(r>1)
        r = r + 20;
        img(1:r,c) = 1;
    end    
end

for c = 1:200
    r = sqrt(radius1^2 - (c-100)^2);
    r = floor(r);
    if(r>1)
        r = r + 20;
        img(1:r,c) = 0;
    end    
end
imwrite(img,'testCirc.png','png')