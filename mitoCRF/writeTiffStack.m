% Write 3D array stack to outputFileName.tiff
function output = writeTiffStack(stack, outputFileName)

imwrite(stack(:,:,1), outputFileName);
for k = 2:size(stack,3)
    imwrite(stack(:,:,k), outputFileName, 'writemode', 'append');
end

output = 0;