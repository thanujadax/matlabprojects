% MRF learning
% Thanuja

% parameters and variables
hasInitDict = 1;
pathToDict = 'Dictionary_256x_200w_16bb.mat';

% get the dictionary words
if(hasInitDict)
    load(pathToDict);       % loads the structure output (from KSVD)
    Dictionary = output.D;  % copy the dictionary
    clear output;           % output from KSVD is no more required 
else
    % generateDictionary()
end

% build the association matrices 
HorizontalAssociations = getHorizontalAssociations(image,Dictionary);
% (remember to add 1)

% Build MRF
% define potentials
%   1-clique:   Vc(f_i) = (pixel error)
%   pairwise:   Vc(f_i,f_i') = (from associations matrix)

% run energy minimization

% sample MRF to generate sample image

