function newWS = replaceWSIDs(newWS,ws,oldRegionIDsToChange,newID)

[pixelsToChange_logical,~] = ismember(ws,oldRegionIDsToChange);

newWS(pixelsToChange_logical) = newID;