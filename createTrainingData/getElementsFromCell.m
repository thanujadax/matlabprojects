function allElementsInList = getElementsFromCell(c_input)

allElementsInList = [];

numCells = numel(c_input);


for i=1:numCells
    elements_i = c_input{i};
    if(~isempty(elements_i))
        [sizeR, sizeC] = size(elements_i);
        if(sizeC==1)
            allElementsInList = [allElementsInList; elements_i];
        elseif(sizeR==1)
            allElementsInList = [allElementsInList; elements_i'];
        else
            disp('*******************************************************')
            disp('ERROR in createTrainingData/getElementsFromCell.m******')
            disp('*******************************************************')
        end
    end
end

allElementsInList = unique(allElementsInList);