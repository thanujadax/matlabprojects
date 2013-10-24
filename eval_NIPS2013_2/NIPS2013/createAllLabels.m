function [mOutput,vLevels] = createAllLabels(cImages)

% nInt = round(255.0/length(cImages));
vLevels = zeros(1,length(cImages));

for i=1:length(cImages)
    if(i==1)
        mOutput = double(zeros(size(cImages{i})));
    end
    
    mOutput(cImages{i}>0) = i; %(i-1)*nInt + 1;
    vLevels(i) = i; %(i-1)*nInt + 1;
end
end

