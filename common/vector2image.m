function wordimg = vector2image(imgAsVector,bb)
wordimg = reshape(imgAsVector,bb,bb);
%imagesc(wordimg(1:bb,1:bb));figure(gcf);
imshow(wordimg);figure(gcf);