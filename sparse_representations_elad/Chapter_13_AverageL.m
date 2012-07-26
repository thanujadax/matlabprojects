% Figure - 13.3
% =========================================
% This small program generates the figure that demonstrates the 
% cardinality required for coding in the facial compression chapter 

P=39600; 
m=512; 
n=(8:1:30).^2; 

figure(1); clf; 
for B=500:100:1000
    Bb=B*8; 
    h=plot(sqrt(n),Bb*n/P/(log(m)/log(2)+7),'k');
    set(h,'LineWidth',2);
    if B==500, 
        hold on;
    end;
end;
set(gca,'FontSize',14);
axis([8 30 0 12]);
h=xlabel('Block size n^{0.5}');
set(h,'FontSize',14);
h=ylabel('Average cardinality k_0'); 
set(h,'FontSize',14);

h=text(23,3,'B=500 Bytes');
set(h,'FontSize',14);
h=text(20,9,'B=1000 Bytes');
set(h,'FontSize',14);

% add arrow manually
print -depsc2 Chapter_13_AverageL.eps