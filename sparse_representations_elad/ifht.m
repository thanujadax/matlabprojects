function [Out]=ifht(In);

N=length(In);
if N==2,
    Out=[1 -1 ; 1 1]/sqrt(2)*In;
else,
    Out1=ifht(In(1:N/2));
    Out2=ifht(In(N/2+1:N));
    Out=[Out1-Out2; Out1+Out2];
    Out=Out/sqrt(2);
end;
return;
