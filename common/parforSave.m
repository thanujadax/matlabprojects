% wrapper around SAVE function
function parforSave(fname,x,opt1)
    save( fname,'x',opt1);
end