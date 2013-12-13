function demoOnImageClick
clc;clear;
imObj = rand(500,500);
figure;
hAxes = axes();
imageHandle = imshow(imObj);
set(imageHandle,'ButtonDownFcn',@ImageClickCallback);

function ImageClickCallback ( objectHandle , eventData )
axesHandle  = get(objectHandle,'Parent');
coordinates = get(axesHandle,'CurrentPoint'); 
coordinates = coordinates(1,1:2);
message     = sprintf('x: %.1f , y: %.1f',coordinates (1) ,coordinates (2));
helpdlg(message);
end

end