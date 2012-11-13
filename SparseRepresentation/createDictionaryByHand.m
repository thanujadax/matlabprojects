function Dictionary = createDictionaryByHand(bb)
% creates a synthetic dictionary for EM images based on hand designed words that
% are intended to approximate EM image patches.

% inputs:
% bb - block dimension (width = height)

%% parameters
G0 = 0;    % BLACK
G1 = 50;
G2 = 100;
G3 = 150;
G4 = 200;   
G5 = 255;  % WHITE

wordSize = bb * bb;
numWords = 16;      % hard coded for the time being
lineWidth = 4;

Dictionary = zeros(wordSize,numWords);
% wordBlocks = zeros(64,64);
% tmpWordBlock = zeros(bb,bb);

%% plane blocks
% word 1
% plain black
i = 1; % word index
tmpWordBlock = ones(bb,bb).* G1;
Dictionary(:,i) = reshape(tmpWordBlock,wordSize,1); % insert word to Dictionary as col vector

% word 2
% plane white
i = i + 1;
tmpWordBlock = ones(bb,bb) .* G4;
Dictionary(:,i) = reshape(tmpWordBlock,wordSize,1);

% gray word
i = i + 1;
tmpWordBlock = ones(bb,bb) .* G3;
Dictionary(:,i) = reshape(tmpWordBlock,wordSize,1);

%% half blocks - horizontal and vertical
i = i + 1; % left vertical half black
tmpWordBlock = ones(bb,bb).*G1;
halfHorizontal = bb/2;
tmpWordBlock(:,halfHorizontal+1:bb) = ones(bb,halfHorizontal).*G4;
Dictionary(:,i) = reshape(tmpWordBlock,wordSize,1);

i = i + 1; % right vertical half black
tmpWordBlock = ones(bb,bb).*G1;
halfHorizontal = bb/2;
tmpWordBlock(:,1:halfHorizontal) = ones(bb,halfHorizontal).*G4;
Dictionary(:,i) = reshape(tmpWordBlock,wordSize,1);

i = i + 1; % upper horizontal half black
tmpWordBlock = rot90(tmpWordBlock);
Dictionary(:,i) = reshape(tmpWordBlock,wordSize,1);

i = i + 1; % lower horizontal half black
tmpWordBlock = rot90(tmpWordBlock,2);
Dictionary(:,i) = reshape(tmpWordBlock,wordSize,1);

%% half blocks - diagonals
i = i + 1; % upper left triangle black
tmpWordBlock = fliplr(triu(ones(bb,bb),0)).*G1 + fliplr(tril(ones(bb,bb),-1)).*G4;
Dictionary(:,i) = reshape(tmpWordBlock,wordSize,1);

i = i + 1; % lower right triangle black
tmpWordBlock = fliplr(triu(ones(bb,bb),0)).*G4 + fliplr(tril(ones(bb,bb),-1)).*G1;
Dictionary(:,i) = reshape(tmpWordBlock,wordSize,1);

i = i + 1; % lower left triangle black
tmpWordBlock = triu(ones(bb,bb),0).*G4 + tril(ones(bb,bb),-1).*G1;
Dictionary(:,i) = reshape(tmpWordBlock,wordSize,1);

i = i + 1; % upper right triangle black
tmpWordBlock = triu(ones(bb,bb),0).*G1 + tril(ones(bb,bb),-1).*G4;
Dictionary(:,i) = reshape(tmpWordBlock,wordSize,1);

%% Line segments
i = i + 1; % horizontal black line segment
tmpWordBlock = ones(bb,bb).*G4;
horizontalDarkStripe = ones(lineWidth,bb).*G1;
startRow = bb/2 - lineWidth/2 +1;
tmpWordBlock(startRow:(startRow+lineWidth-1),:) = horizontalDarkStripe;
Dictionary(:,i) = reshape(tmpWordBlock,wordSize,1);

i = i + 1; % vertical black line segment
tmpWordBlock = rot90(tmpWordBlock);
Dictionary(:,i) = reshape(tmpWordBlock,wordSize,1);

i = i + 1; % digoanal black strip starting from bottom left
Diagonal_0 = ones(1,bb).*G1;
Diagonal_1 = ones(1,bb-1).*G1;
Diagonal_2 = ones(1,bb-2).*G2;

diagonalStripe1 = diag(Diagonal_1,-1);
diagonalStripe2 = diag(Diagonal_0,0);
diagonalStripe3 = diag(Diagonal_1,1);
diagonalStripe4 = diag(Diagonal_2,2);
diagonalStripe5 = diag(Diagonal_2,-2);

tmpWordBlock = ones(bb,bb).*G4;
tmpWordBlock = tmpWordBlock - diag(diag(tmpWordBlock,0),0)...
                    - diag(diag(tmpWordBlock,1),1) - diag(diag(tmpWordBlock,-1),-1)...
                    - diag(diag(tmpWordBlock,2),2) - diag(diag(tmpWordBlock,-2),-2)...
                    + diagonalStripe1 + diagonalStripe2 + diagonalStripe3...
                    + diagonalStripe4 + diagonalStripe5;
Dictionary(:,i) = reshape(tmpWordBlock,wordSize,1);
diagonalBlock1 = tmpWordBlock;

i = i + 1; % diagonal black strip starting from top left
tmpWordBlock = fliplr(tmpWordBlock);
Dictionary(:,i) = reshape(tmpWordBlock,wordSize,1);

%% Y shapes
i = i + 1; % Y downwards
tmpWordBlock = zeros(bb,bb).*G4;
center = bb/2;
Segment1 = ones(center,lineWidth).*G1;

Diagonal_0 = ones(1,center).*G1;
Diagonal_1 = ones(1,center-1).*G1;
Diagonal_2 = ones(1,center-2).*G2;

diagonalStripe1 = diag(Diagonal_1,-1);
diagonalStripe2 = diag(Diagonal_0,0);
diagonalStripe3 = diag(Diagonal_1,1);
diagonalStripe4 = diag(Diagonal_2,2);
diagonalStripe5 = diag(Diagonal_2,-2);      

Segment2 = ones(center,center).*G4;
Segment2 = Segment2 - diag(diag(Segment2,0),0)...
                    - diag(diag(Segment2,1),1) - diag(diag(Segment2,-1),-1)...
                    - diag(diag(Segment2,2),2) - diag(diag(Segment2,-2),-2)...
                    + diagonalStripe1 + diagonalStripe2 + diagonalStripe3...
                    + diagonalStripe4 + diagonalStripe5;
                
Segment3 = fliplr(Segment2);
startCol = center-lineWidth/2 +1;
tmpWordBlock(1:center,startCol:startCol+lineWidth-1) = Segment1;
tmpWordBlock(center+1:bb,1:center) = Segment3;
tmpWordBlock(center+1:bb,center+1:bb) = Segment2;

zeroInd = find(tmpWordBlock==0);
tmpWordBlock(zeroInd) = G4;
Dictionary(:,i) = reshape(tmpWordBlock,wordSize,1);

i = i + 1; % Y pointing to the right
tmpWordBlock = rot90(tmpWordBlock);
Dictionary(:,i) = reshape(tmpWordBlock,wordSize,1);

i = i + 1; % Y normal way
tmpWordBlock = rot90(tmpWordBlock);
Dictionary(:,i) = reshape(tmpWordBlock,wordSize,1);

i = i + 1; % Y pointing to the left
tmpWordBlock = rot90(tmpWordBlock);
Dictionary(:,i) = reshape(tmpWordBlock,wordSize,1);


%% Return value - normalized dictionary
Dictionary = Dictionary./G5;
