function neighbors = get8Neighbors(pixID,sizeR,sizeC)
[r0,c0] = ind2sub([sizeR sizeC],pixID);
% TODO: sanity checks when entering values into neighbors matrix
numPixInput = numel(pixID);

neighbors = zeros(numPixInput,8);

i = 0;
% N1
r = r0 - 1;
c = c0;

i = i + 1;
neighbors((r>0),i) = sub2ind([sizeR sizeC],r(r>0),c(r>0));


% N2
r = r0 + 1;
c = c0;

i = i + 1;
neighbors((r<=sizeR),i) = sub2ind([sizeR sizeC],r,c);

% N3
r = r0;
c = c0 - 1;

i = i + 1;
neighbors((c>0),i) = sub2ind([sizeR sizeC],r(c>0),c(c>0));


% N4
r = r0;
c = c0 + 1;

i = i + 1;
neighbors((c<=sizeC),i) = sub2ind([sizeR sizeC],r(c<=sizeC),c(c<=sizeC));

% N5
r = r0 - 1;
c = c0 + 1;

i = i + 1;
neighbors((r>0),i) = sub2ind([sizeR sizeC],r(r>0),c(r>0));

% N6
r = r0 + 1;
c = c0 + 1;

    i = i + 1;
    neighbors(:,i) = sub2ind([sizeR sizeC],r,c);


% N7
r = r0 - 1;
c = c0 - 1;

    i = i + 1;
    neighbors(:,i) = sub2ind([sizeR sizeC],r,c);


% N8
r = r0 + 1;
c = c0 - 1;

    i = i + 1;
    neighbors(:,i) = sub2ind([sizeR sizeC],r,c);
