function allElementsInList = getElementsFromCell(c_input)

allElementsInList = [];

numCells = numel(c_input);


for i=1:numCells
    elements_i = c_input{i};
    allElementsInList = [allElementsInList; elements_i];
end