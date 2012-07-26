% Figures - 15.27 - 15.31
% =========================================
% This program presents treates the inpainting problem locally as 
% the two previous algorithms, but introduces a leearned dictionary. 
% The results built here are used to fill Table 15.2


% First test - with a text image

[Train,TestL,TestH,Ideal,Result]=...
                    Chapter_15_ScaleUp_Training_and_Testing...
                    ('TextImage3.png','TextImage4.png');

close all; 

figure(1); clf; 
Train(:,1)=0; Train(:,end)=0; 
Train(1,:)=0; Train(end,:)=0;
image(Train); axis image; axis off; truesize; colormap(gray(256)); 
% print -depsc2 Chapter_15_ScaleUp_Text_TrainImage.eps

figure(2); clf; 
TestL(:,1)=0; TestL(:,end)=0; 
TestL(1,:)=0; TestL(end,:)=0;
image(TestL); axis image; axis off; truesize; colormap(gray(256)); 
% print -depsc2 Chapter_15_ScaleUp_Text_TestLImage.eps

figure(3); clf; 
TestH(:,1)=0; TestH(:,end)=0; 
TestH(1,:)=0; TestH(end,:)=0;
image(TestH); axis image; axis off; truesize; colormap(gray(256)); 
% print -depsc2 Chapter_15_ScaleUp_Text_TestHImage.eps

figure(4); clf; 
Ideal(:,1)=0; Ideal(:,end)=0; 
Ideal(1,:)=0; Ideal(end,:)=0;
image(Ideal); axis image; axis off; truesize; colormap(gray(256)); 
% print -depsc2 Chapter_15_ScaleUp_Text_TestIdealImage.eps

figure(5); clf; 
Result(:,1)=0; Result(:,end)=0; 
Result(1,:)=0; Result(end,:)=0;
image(Result); axis image; axis off; truesize; colormap(gray(256)); 
% print -depsc2 Chapter_15_ScaleUp_Text_ResultImage.eps

% Second test - Child

Chapter_15_ScaleUp_Training_and_Testing2('BuildingImage1.png'); 

figure(1); 
% print -depsc2 Chapter_15_ScaleUp_Building.eps

figure(2); 
% print -depsc2 Chapter_15_ScaleUp_Building_Part1.eps

figure(3); 
% print -depsc2 Chapter_15_ScaleUp_Building_Part2.eps

