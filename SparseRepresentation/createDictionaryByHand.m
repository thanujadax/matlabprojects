function Dictionary = createDictionaryByHand(bb)
% creates a synthetic image for EM images based on hand designed words that
% are intended to approximate EM image patches.

% inputs:
% bb - block dimension (width = height)

%% parameters
G0 = 0;    % BLACK
G1 = 50;
G2 = 100;
G3 = 150;
G4 = 200;   
G5 = 250;  % WHITE

wordSize = bb * bb;
numWords = 16;      % hard coded for the time being

Dictionary = zeros(wordSize,numWords);
% wordBlocks = zeros(64,64);
% tmpWordBlock = zeros(bb,bb);

%% plane blocks
% word 1
% plain black
i = 1; % word index
tmpWordBlock = ones(bb,bb).* G1;
Dictionary(:,i) = reshape(wordSize,1,tmpWordBlock); % insert word to Dictionary as col vector

% word 2
% plane white
i = 2;
tmpWordBlock = ones(bb,bb) .* G4;
Dictionary(:,i) = reshape(wordSize,1,tmpWordBlock);

%% half blocks
i = 3; % left vertical half black
tmpWordBlock = ones(bb,bb).*G1;
halfpointHorizontal = bb/2;


