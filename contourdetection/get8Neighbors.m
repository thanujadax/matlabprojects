function neighbors = get8Neighbors(pixID,sizeR,sizeC)

% works with pixID being a vector of pixel indices

% warning: contains zero for non existent neighbors

[r0,c0] = ind2sub([sizeR sizeC],pixID);
% TODO: sanity checks when entering values into neighbors matrix
numPixInput = numel(pixID);

neighbors = zeros(numPixInput,8);

i = 0;
% N1
r = r0 - 1;
c = c0;

i = i + 1;
if(sum(r>0)>0)
neighbors((r>0),i) = sub2ind([sizeR sizeC],r(r>0),c(r>0));
end

% N2
r = r0 + 1;
c = c0;

i = i + 1;
if(sum(r<=sizeR)>0)
neighbors((r<=sizeR),i) = sub2ind([sizeR sizeC],r(r<=sizeR),c(r<=sizeR));
end

% N3
r = r0;
c = c0 - 1;

i = i + 1;
if(sum(c>0)>0)
neighbors((c>0),i) = sub2ind([sizeR sizeC],r(c>0),c(c>0));
end

% N4
r = r0;
c = c0 + 1;

i = i + 1;
if(sum(c<=sizeC)>0)
neighbors((c<=sizeC),i) = sub2ind([sizeR sizeC],r(c<=sizeC),c(c<=sizeC));
end

% N5
r = r0 - 1;
c = c0 + 1;

i = i + 1;
if(sum(r>0 & c<=sizeC)>0)
neighbors((r>0 & c<=sizeC),i) = sub2ind([sizeR sizeC],r(r>0 & c<=sizeC),c(r>0 & c<=sizeC));
end
% N6
r = r0 + 1;
c = c0 + 1;

i = i + 1;
if(sum(r<=sizeR & c<=sizeC)>0)
neighbors((r<=sizeR & c<=sizeC),i) = ...
        sub2ind([sizeR sizeC],r(r<=sizeR & c<=sizeC),c(r<=sizeR & c<=sizeC));
end

% N7
r = r0 - 1;
c = c0 - 1;

i = i + 1;
if(sum(r>0 & c>0)>0)
neighbors((r>0 & c>0),i) = sub2ind([sizeR sizeC],r(r>0 & c>0),c(r>0 & c>0));
end

% N8
r = r0 + 1;
c = c0 - 1;

i = i + 1;
if(sum(r<=sizeR & c>0)>0)
neighbors((r<=sizeR & c>0),i) = ...
        sub2ind([sizeR sizeC],r(r<=sizeR & c>0),c(r<=sizeR & c>0));
end